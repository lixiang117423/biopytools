"""PSVCP 泛基因组构建编排|PSVCP pangenome construction orchestration

PangenomeConstructor 用 Python 重做 PSVCP 1Genome_construct_Pangenome.py /
Refgenome_update_by_quest.sh 的逐轮循环(步骤序列与命名 1:1 对齐原始),外部工具走
build_conda_command + 命令日志,vendored helper 原样调用,支持断点续传。
|Reimplements the PSVCP per-round loop (step sequence & naming mirror the original),
wraps external tools with build_conda_command + command logging, invokes vendored
helpers verbatim, and supports checkpoint resume.
"""

import logging
import os
from typing import List, Optional, Tuple

from .utils import (
    build_conda_command, check_assemblytics_numpy, run_command,
    run_shell_command,
)


class PangenomeConstructor:
    """泛基因组构建器|Pangenome constructor (per-round MUMmer→Assemblytics→PAV incorporate)"""

    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger
        # vendored helper 脚本目录|vendored helper script dir
        self.scripts_dir = os.path.join(os.path.dirname(__file__), 'scripts', 'construct_pan')

    # ------------------------------------------------------------------ #
    # 路径/命名辅助(照搬原始)|path/naming helpers (mirror original)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _stem(name: str) -> str:
        """基因组名去 .fa 后缀|strip .fa suffix"""
        return name[:-3] if name.endswith('.fa') else name

    @staticmethod
    def _gff_name(name: str) -> str:
        """基因组名 → 对应 .gff 名(原始 re.sub('.fa','.gff'))|name→.gff name"""
        return name[:-3] + '.gff' if name.endswith('.fa') else name + '.gff'

    def _link(self, target: str, linkname: str):
        """强制符号链接(覆盖已存在)|force symlink (overwrite)"""
        link_dir = os.path.dirname(linkname)
        if link_dir:
            os.makedirs(link_dir, exist_ok=True)
        if os.path.lexists(linkname):
            os.remove(linkname)
        os.symlink(target, linkname)

    def _is_nonempty(self, path: str) -> bool:
        """文件存在且非空(原始 [ -s file ])|exists and non-empty"""
        return os.path.isfile(path) and os.path.getsize(path) > 0

    # ------------------------------------------------------------------ #
    # 命令执行封装|command execution wrappers
    # ------------------------------------------------------------------ #
    def _run_tool(self, tool_path: str, args: List[str], desc: str,
                  cwd: Optional[str] = None, capture: bool = True) -> Tuple[str, str]:
        """执行外部工具(conda 包装,列表+shell=False)|run external tool (conda-wrapped)"""
        cmd = build_conda_command(tool_path, args)
        return run_command(cmd, self.logger, desc, cwd=cwd, capture=capture)

    def _run_py(self, script_name: str, args: List[str], desc: str,
                cwd: Optional[str] = None, capture: bool = True) -> Tuple[str, str]:
        """执行 vendored python helper(conda run python3)|run vendored python helper"""
        script = os.path.join(self.scripts_dir, script_name)
        cmd = build_conda_command(self.config.python_path, [script] + args)
        return run_command(cmd, self.logger, desc, cwd=cwd, capture=capture)

    def _run_r(self, script_name: str, args: List[str], desc: str,
               cwd: Optional[str] = None) -> Tuple[str, str]:
        """执行 vendored R helper(conda run Rscript)|run vendored R helper"""
        script = os.path.join(self.scripts_dir, script_name)
        cmd = build_conda_command(self.config.rscript_path, [script] + args)
        return run_command(cmd, self.logger, desc, cwd=cwd)

    def _run_shell(self, cmd_str: str, desc: str, cwd: Optional[str] = None) -> bool:
        """执行 shell 字符串(awk/sort 管道)|run shell string (awk/sort pipeline)"""
        return run_shell_command(cmd_str, self.logger, desc, cwd=cwd)

    def _awk_emitter(self, script_name: str, args: List[str], desc: str, cwd: str) -> bool:
        """awk 发射器:python helper 打印 awk 到 stdout,捕获后执行|awk emitter pattern"""
        awk_str, _ = self._run_py(script_name, args, desc, cwd=cwd, capture=True)
        awk_str = awk_str.strip()
        if not awk_str:
            raise RuntimeError(f"发射器未产出 awk 命令|emitter produced no awk: {script_name}")
        return self._run_shell(awk_str, f"{desc}(awk)", cwd=cwd)

    def _ensure_faidx(self, fa: str):
        """确保 .fai 存在(bedtools getfasta 需要);idempotent|ensure .fai exists"""
        if not os.path.exists(fa + '.fai'):
            self._run_tool(self.config.samtools_path, ['faidx', fa], f"samtools faidx {os.path.basename(fa)}")

    def _grep_to_file(self, pattern: str, in_file: str, out_file: str):
        """Python 实现 grep(行内含 pattern 则写出)|grep substring to file"""
        with open(in_file, 'r', encoding='utf-8') as fin, \
                open(out_file, 'w', encoding='utf-8') as fout:
            for line in fin:
                if pattern in line:
                    fout.write(line)

    # ------------------------------------------------------------------ #
    # 主流程|main pipeline
    # ------------------------------------------------------------------ #
    def run(self) -> bool:
        """运行泛基因组构建|Run pangenome construction"""
        cfg = self.config
        try:
            # 预检:assemblytics 依赖 numpy(BLOCKER 1 早暴露)|pre-check numpy
            if not check_assemblytics_numpy(cfg.python_path, self.logger):
                return False

            os.makedirs(cfg.pan_dir, exist_ok=True)
            os.makedirs(os.path.join(cfg.pan_dir, 'temp'), exist_ok=True)

            # 整体断点续传|overall checkpoint
            pan_fa = os.path.join(cfg.output_dir, 'pan.fa')
            if self._is_nonempty(pan_fa) and not cfg.force:
                self.logger.info(f"跳过已完成流程(pan.fa 已存在)|Skipping: pan.fa already exists: {pan_fa}")
                return True

            names = cfg.read_genome_list()
            ref_name = names[0]
            queries = names[1:]
            self.logger.info(
                f"参考基因组|Reference: {ref_name}; query 数|queries: {len(queries)}; "
                f"总轮数|rounds: {len(queries)}"
            )

            # 初始化 ref0 = 链接到 genome_dir/ref|init ref0
            ref0_fa = os.path.join(cfg.pan_dir, 'ref0.fa')
            ref0_gff = os.path.join(cfg.pan_dir, 'ref0.gff')
            self._link(os.path.join(cfg.genome_dir, ref_name), ref0_fa)
            self._link(os.path.join(cfg.genome_dir, self._gff_name(ref_name)), ref0_gff)

            ref_fa = ref0_fa
            prev_gff = ref0_gff
            prev_pav: Optional[str] = None  # ref0 无 pav|ref0 has no pav
            self._total_rounds = len(queries)

            for n, q_name in enumerate(queries, start=1):
                ref_stem = f"ref{n - 1}"
                q_stem = self._stem(q_name)
                done = self._run_round(n, ref_stem, q_stem, q_name,
                                       ref_fa, prev_gff, prev_pav)
                if not done:
                    return False
                # 推进指针|advance pointers
                ref_fa = os.path.join(cfg.pan_dir, f'ref{n}.fa')
                prev_gff = os.path.join(cfg.pan_dir, f'ref{n}.gff')
                prev_pav = os.path.join(cfg.pan_dir, f'ref{n}.pav.gff')

            self._finalize(len(queries))
            return True

        except Exception as e:
            self.logger.error(f"泛基因组构建失败|Pangenome construction failed: {e}")
            return False

    def _run_round(self, n: int, ref_stem: str, q_stem: str, q_name: str,
                   ref_fa: str, prev_gff: str, prev_pav: Optional[str]) -> bool:
        """单轮:ref{N-1} + query → ref{N}|one round: ref{N-1} + query → ref{N}"""
        cfg = self.config
        chain = f"{ref_stem}{q_stem}"            # 链式名(无下划线,原始约定)|chain name (no underscore)
        prefix = f"{ref_stem}_{q_stem}"          # 工作子目录名(带下划线)|work subdir (underscore)
        work = os.path.join(cfg.pan_dir, prefix)
        # 链式产物放 pan_dir(原始)|chain artifacts in pan_dir (original)
        chain_fa = os.path.join(cfg.pan_dir, f'{chain}.fa')
        chain_gff = os.path.join(cfg.pan_dir, f'{chain}.gff')
        chain_pav = os.path.join(cfg.pan_dir, f'{chain}.pav.gff')

        refN_fa = os.path.join(cfg.pan_dir, f'ref{n}.fa')
        # 逐轮断点续传|per-round checkpoint
        if self._is_nonempty(refN_fa) and not cfg.force:
            self.logger.info(f"跳过已完成轮|Skipping completed round {n}: {prefix}")
            return True

        os.makedirs(work, exist_ok=True)
        os.makedirs(os.path.join(work, 'temp'), exist_ok=True)
        query_fa = os.path.join(cfg.genome_dir, q_name)
        query_gff = os.path.join(cfg.genome_dir, self._gff_name(q_name))
        self.logger.info(
            f"===== 第 {n}/{self._total_rounds} 轮|Round {n}: "
            f"{ref_stem} + {q_stem} ====="
        )

        # 1. samtools faidx(确保 .fai)|ensure .fai
        self._ensure_faidx(ref_fa)
        self._ensure_faidx(query_fa)

        # 2. nucmer|nucmer
        delta = os.path.join(work, f'{prefix}.delta')
        self._run_tool(cfg.nucmer_path,
                       ['-t', str(cfg.threads), '--maxgap', str(cfg.NUCMER_MAXGAP),
                        '--mincluster', str(cfg.NUCMER_MINCLUSTER), '--diagdiff', str(cfg.NUCMER_DIAGDIFF),
                        ref_fa, query_fa, '--prefix', os.path.join(work, prefix)],
                       f"nucmer {ref_stem} vs {q_stem}")

        # 3. assemblytics 2.0.1(新接口 -d -o -l -min -max)|assemblytics 2.0.1
        asm_dir = os.path.join(work, 'assemblytics')
        os.makedirs(asm_dir, exist_ok=True)
        self._run_tool(cfg.assemblytics_path,
                       ['-d', delta, '-o', asm_dir,
                        '-l', str(cfg.ASSEMBLYTICS_UNIQUE_LENGTH),
                        '-min', str(cfg.ASSEMBLYTICS_MIN_SIZE),
                        '-max', str(cfg.ASSEMBLYTICS_MAX_SIZE)],
                       f"assemblytics {ref_stem} vs {q_stem}")
        sv_bed = os.path.join(asm_dir, 'assemblytics_structural_variants.bed')

        # 4. 2variants_to_coords.bed.py(awk 发射器)|awk emitter → coords.bed
        coords_bed = os.path.join(work, f'{prefix}.coords.bed')
        self._awk_emitter('2variants_to_coords.bed.py',
                          [q_stem, sv_bed, coords_bed],
                          f"2variants {ref_stem}+{q_stem}", cwd=work)

        # 5. 3replace(stdout=final.bed;cwd=work 写 temp/)|3replace
        final_bed = os.path.join(work, f'{prefix}.final.bed')
        stdout, _ = self._run_py('3replace_regions_with_nucleotides1.2.py',
                                 [coords_bed, ref_fa, query_fa, chain],
                                 f"3replace {ref_stem}+{q_stem}", cwd=work)
        with open(final_bed, 'w', encoding='utf-8') as f:
            f.write(stdout)

        # 6. grep Insertion/Deletion(Python 实现)|grep
        ins_bed = os.path.join(work, 'ins.bed')
        del_bed = os.path.join(work, 'del.bed')
        self._grep_to_file('Insertion', final_bed, ins_bed)
        self._grep_to_file('Deletion', final_bed, del_bed)

        # 7. 4ins_more_50.py(awk 发射器)|awk emitter → ins.more_50.bed
        ins50 = os.path.join(work, 'ins.more_50.bed')
        if self._is_nonempty(ins_bed):
            self._awk_emitter('4ins_more_50.py', [ins_bed, ins50],
                              f"4ins_more_50 {q_stem}", cwd=work)
        else:
            open(ins50, 'w').close()

        # 8. 分支|branch on insertions
        if self._is_nonempty(ins50):
            self._incorporate_insertions(
                n, work, chain, q_stem, chain_fa, chain_gff, chain_pav,
                ref_fa, query_fa, query_gff, prev_gff, prev_pav, ins50
            )
            with open(ins50, 'r', encoding='utf-8') as _f:
                pav_count = sum(1 for _ in _f)
            self.logger.info(f"本轮并入 PAV 数|PAVs incorporated this round: {pav_count}")
        else:
            self.logger.info(f"本轮无 >{cfg.MIN_INSERTION_SIZE}bp 插入,ref 不变|No insertions >{cfg.MIN_INSERTION_SIZE}bp, ref unchanged")
            # 空分支(照搬原始)|empty branch
            if n == 1:
                open(chain_pav, 'w').close()           # touch
                self._link(prev_gff, chain_gff)
            else:
                self._link(prev_pav, chain_pav)
                self._link(prev_gff, chain_gff)
            self._link(ref_fa, chain_fa)

        # ref{N}.* = 链接到 chain.*|ref{N}.* link to chain.*
        self._link(chain_fa, refN_fa)
        self._link(chain_gff, os.path.join(cfg.pan_dir, f'ref{n}.gff'))
        self._link(chain_pav, os.path.join(cfg.pan_dir, f'ref{n}.pav.gff'))
        return True

    def _incorporate_insertions(self, n: int, work: str, chain: str, q_stem: str,
                                chain_fa: str, chain_gff: str, chain_pav: str,
                                ref_fa: str, query_fa: str, query_gff: str,
                                prev_gff: str, prev_pav: Optional[str], ins50: str):
        """非空分支:5→6→7→7.2→pav 更新→8→9→10→11→12(照搬原始)|non-empty branch"""
        bed2 = ins50 + '2'                 # 5ins 输出 argv[1]+'2' → ins.more_50.bed2
        bed3 = ins50 + '3'                 # 7bed2.R gsub 'bed2$'→'bed3' → ins.more_50.bed3

        # 5. ins.bed_to_bed2 → bed2|bed2
        self._run_py('5ins.bed_to_bed2.py', [ins50], f"5ins_to_bed2 {chain}")
        # 6. update_ref_by_nucmer → chain.fa(拼插入序列进 ref)|splice insertions into ref
        self._run_py('6update_ref_by_nucmer.py', [ref_fa, bed2, chain_fa],
                     f"6update_ref {chain}")
        # 7. bed2_update_bed3.R → bed3|bed3
        self._run_r('7bed2_update_bed3.R', [bed2], f"7bed2_to_bed3 {chain}")
        # 7.2 bed3_to_gff → bed3.pav.gff(放 work)|bed3.pav.gff (in work)
        bed3_pav_gff = os.path.join(work, f'{chain}.bed3.pav.gff')
        self._run_py('7.2bed3_to_gff.py', [bed3, bed3_pav_gff],
                     f"7.2bed3_to_gff {chain}")

        # pav.gff 更新|update pav.gff
        if n == 1:
            self._link(bed3_pav_gff, chain_pav)
        else:
            update1_pav = os.path.join(work, 'update1.pav.gff')
            self._run_r('8.2pav_gff_update_by_bed2info_parLapply.R',
                        [prev_pav, bed2, update1_pav], f"8.2pav_gff_update {chain}")
            self._cat([update1_pav, bed3_pav_gff], chain_pav)

        # 8. gff_update_by_bed2info → update1.gff|update1.gff
        update1_gff = os.path.join(work, 'update1.gff')
        self._run_r('8gff_update_by_bed2info_parLapply.R',
                    [prev_gff, bed2, update1_gff], f"8gff_update {chain}")
        # 9. gff_split_by_bed3 → update2.gff|update2.gff
        update2_gff = os.path.join(work, 'update2.gff')
        self._run_r('9gff_split_by_bed3_6.R',
                    [update1_gff, bed3, update2_gff], f"9gff_split {chain}")
        # 10. gene_in_pv_from_gff → {q_stem}.gff.in_pv|genes in PAV
        # 注:R 脚本仅当有基因与 PAV 重叠(nrow>0)时才写文件,否则不产出 .in_pv
        # |R script writes .in_pv only if genes overlap PAVs; otherwise no file
        in_pv = os.path.join(work, f'{q_stem}.gff.in_pv')
        self._run_r('10gene_in_pv_from_gff_parLapply2.R',
                    [query_gff, bed3, in_pv], f"10gene_in_pv {chain}")

        # 原始三级条件(Refgenome_update_by_quest.sh):无重叠 / 无完全落 PAV 基因 → link update2.gff
        # |original 3-level guard: no overlap / no fully-in-PAV gene → link update2.gff
        if os.path.isfile(in_pv):
            # 11. gene_in_pv_screen → absolutly|absolutely-in-pv genes
            absolutly = os.path.join(work, f'{q_stem}.gff_gene_absolutly.in_pv')
            self._run_py('11gene_in_pv_screen.py', [in_pv, absolutly],
                         f"11gene_in_pv_screen {chain}")
            if self._is_nonempty(absolutly):
                # 12. bed3_update_gff → gene_absolutly_in_pv.gff|update gene coords
                gene_in_pv_gff = os.path.join(work, f'{chain}.gene_absolutly_in_pv.gff')
                self._run_r('12bed3_update_gff.R', [bed3, absolutly, gene_in_pv_gff],
                            f"12bed3_update_gff {chain}")
                if os.path.isfile(gene_in_pv_gff):
                    self._cat([update2_gff, gene_in_pv_gff], chain_gff)
                else:
                    self._link(update2_gff, chain_gff)
            else:
                self._link(update2_gff, chain_gff)
        else:
            self.logger.info(
                f"无基因与 PAV 重叠,link update2.gff|"
                f"No gene overlaps PAV, linking update2.gff: {chain}")
            self._link(update2_gff, chain_gff)

    @staticmethod
    def _cat(files: List[str], out_file: str):
        """拼接文件(cat)|concatenate files"""
        with open(out_file, 'w', encoding='utf-8') as fout:
            for fp in files:
                with open(fp, 'r', encoding='utf-8') as fin:
                    fout.write(fin.read())

    def _finalize(self, n_rounds: int):
        """终化:pan.* 链接 + pav 排序|finalize: pan.* links + sort pav"""
        cfg = self.config
        refN = n_rounds  # 最后一轮的 N|last round N
        pan_fa = os.path.join(cfg.output_dir, 'pan.fa')
        pan_gff = os.path.join(cfg.output_dir, 'pan.gff')
        pan_pav = os.path.join(cfg.output_dir, 'pan.pav.gff')
        self._link(os.path.join(cfg.pan_dir, f'ref{refN}.fa'), pan_fa)
        self._link(os.path.join(cfg.pan_dir, f'ref{refN}.gff'), pan_gff)
        self._link(os.path.join(cfg.pan_dir, f'ref{refN}.pav.gff'), pan_pav)
        # 排序 pav(原始 sort -k1.4,1n -k4,4n)|sort pav
        pan_pav_sorted = os.path.join(cfg.output_dir, 'pan.pav.sorted.gff')
        self._run_shell(
            f"sort -k1.4,1n -k4,4n {pan_pav} > {pan_pav_sorted}",
            "sort pan.pav.gff"
        )
        self.logger.info(f"泛基因组构建完成|Pangenome construction completed")
        self.logger.info(f"  pan.fa         : {pan_fa}")
        self.logger.info(f"  pan.gff        : {pan_gff}")
        self.logger.info(f"  pan.pav.gff    : {pan_pav}")
        self.logger.info(f"  pan.pav.sorted : {pan_pav_sorted}")

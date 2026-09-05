"""
选择性扫荡统计合并模块|Selective Sweep Statistics Merging Module

移植自用户脚本merge_stats.py,修正:RAiSD真实列序/POS分箱、//分隔行解析、
2列与7列报告自适应|Ported from merge_stats.py; fixes: true RAiSD column order,
// separator parsing, 2-col/7-col report auto-detection
"""

from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def bin_pos(pos, win):
    """窗口对齐|Bin position to window grid"""
    return (pos // win) * win


def rank_pct(series, ascending=True):
    """经验百分位排名[0,1],NaN安全|Empirical percentile rank in [0,1], NaN-safe"""
    return series.rank(pct=True, ascending=ascending, na_option='keep')


def load_pi(path, pop, win) -> pd.DataFrame:
    """加载vcftools --window-pi输出|Load vcftools --window-pi output

    输入列|Input cols: CHROM BIN_START [BIN_END] N_VARIANTS PI
    """
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'CHROM': 'CHR', 'BIN_START': 'WIN'})
    df['WIN'] = bin_pos(df['WIN'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[['PI']].mean()
    df = df.rename(columns={'PI': f'PI_{pop}'})
    return df


def load_tajd(path, pop, win) -> pd.DataFrame:
    """加载vcftools --TajimaD输出|Load vcftools --TajimaD output

    输入列|Input cols: CHROM BIN_START N_SNPS TajimaD
    """
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'CHROM': 'CHR', 'BIN_START': 'WIN'})
    df['WIN'] = bin_pos(df['WIN'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[['TajimaD']].mean()
    df = df.rename(columns={'TajimaD': f'TajD_{pop}'})
    return df


def load_fst(path, label, win) -> pd.DataFrame:
    """加载vcftools --weir-fst-pop窗口输出|Load vcftools weir-fst windowed output

    输入列|Input cols: CHROM BIN_START [BIN_END] N_VARIANTS WEIGHTED_FST MEAN_FST
    """
    df = pd.read_csv(path, sep='\t')
    fst_col = 'WEIGHTED_FST' if 'WEIGHTED_FST' in df.columns else 'MEAN_FST'
    df = df.rename(columns={'CHROM': 'CHR', 'BIN_START': 'WIN'})
    df['WIN'] = bin_pos(df['WIN'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[[fst_col]].mean()
    df = df.rename(columns={fst_col: f'Fst_{label}'})
    return df


def load_raisd(files, pop, win) -> pd.DataFrame:
    """加载RAiSD报告文件|Load RAiSD report files

    - 2列(无-R): [POS, MU];7列(-R): [POS, START, END, VAR, SFS, LD, MU]
    - 跳过//注释行;若含"// <名字>"节分隔符,用<名字>作染色体名(优先于文件名后缀)
    - 异常列数的文件抛错(由调用方捕获降级)|Files with unexpected column counts
      raise ValueError (caught and downgraded by caller)

    Args:
        files: 报告文件路径列表|Report file paths
        pop: 群体名|Population name
        win: 窗口大小|Window size

    Returns:
        DataFrame[CHR, WIN, MU_{pop}],每窗口取MU最大值|MU max per window
    """
    frames = []
    for f in files:
        rows = []  # (chrom_or_None, line)
        current_chrom = None
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('//'):
                    tokens = line[2:].strip().split()
                    if tokens:
                        current_chrom = tokens[0]
                    continue
                rows.append((current_chrom, line))
        if not rows:
            continue
        n_cols = len(rows[0][1].split())
        if n_cols == 2:
            names = ['POS', 'MU']
        elif n_cols == 7:
            names = ['POS', 'START', 'END', 'VAR', 'SFS', 'LD', 'MU']
        else:
            # 静默跳过会让该染色体的MU数据无声消失,必须显式报错(合并层try/except降级为WARNING)
            # |silently skipping would drop this chromosome's MU data unnoticed;
            # raise instead (merger's try/except downgrades to a WARNING)
            raise ValueError(
                f"RAiSD 报告列数异常|Unexpected column count in {f}: {n_cols}")
        text = '\n'.join(line for _, line in rows)
        d = pd.read_csv(StringIO(text), sep='\t', header=None, names=names)
        # 染色体名:优先//分隔行,兜底文件名最后下划线后缀
        # |chrom: prefer // separator, fallback to last underscore suffix of filename
        fallback_chrom = Path(f).name.rsplit('_', 1)[-1]
        d['CHR'] = [rc if rc else fallback_chrom for rc, _ in rows]
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=['CHR', 'WIN', f'MU_{pop}'])
    df = pd.concat(frames, ignore_index=True)
    df['WIN'] = bin_pos(df['POS'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[['MU']].max()
    df = df.rename(columns={'MU': f'MU_{pop}'})
    return df


def load_sweed(files, pop, win) -> pd.DataFrame:
    """加载SweeD报告文件|Load SweeD report files

    - 每行格式|Line format: Position Likelihood(CLR) Alpha StartPos EndPos
    - 表头行(首列非数字)自动跳过|Non-numeric header line auto-skipped
    - `//` 分隔行(仅reports=0模式出现)取首token作染色体名,优先于文件名后缀
    - 染色体名:优先//分隔行,兜底文件名最后下划线后缀(与load_raisd一致)
    - 异常列数/无法解析的行抛错(由调用方捕获降级为WARNING)

    Args:
        files: 报告文件路径列表|Report file paths
        pop: 群体名|Population name
        win: 窗口大小|Window size

    Returns:
        DataFrame[CHR, WIN, CLR_{pop}],每窗口取CLR最大值|CLR max per window
    """
    frames = []
    for f in files:
        rows = []  # (chrom_or_None, pos, clr)
        current_chrom = None
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('//'):
                    tokens = line[2:].strip().split()
                    if tokens:
                        current_chrom = tokens[0]
                    continue
                parts = line.split()
                try:
                    # Position为%.4f浮点(如50200.0000),先转float再取整
                    # |Position is %.4f float (e.g. 50200.0000); float then int
                    pos = int(float(parts[0]))
                    clr = float(parts[1])
                except (ValueError, IndexError):
                    # 表头行|header line (Position Likelihood Alpha StartPos EndPos)
                    continue
                if len(parts) < 5:
                    # 静默跳过会让该染色体的CLR数据无声消失,必须显式报错
                    # |silently skipping would drop this chromosome's CLR data
                    # unnoticed; raise instead (caller downgrades to WARNING)
                    raise ValueError(
                        f"SweeD 报告列数异常|Unexpected column count in {f}: {len(parts)}")
                rows.append((current_chrom, pos, clr))
        if not rows:
            continue
        fallback_chrom = Path(f).name.rsplit('_', 1)[-1]
        d = pd.DataFrame(rows, columns=['CHR_SRC', 'POS', 'CLR'])
        d['CHR'] = [rc if rc else fallback_chrom for rc, _, _ in rows]
        d = d.drop(columns=['CHR_SRC'])
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=['CHR', 'WIN', f'CLR_{pop}'])
    df = pd.concat(frames, ignore_index=True)
    df['WIN'] = bin_pos(df['POS'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[['CLR']].max()
    df = df.rename(columns={'CLR': f'CLR_{pop}'})
    return df


def load_xpclr(files, label, win) -> pd.DataFrame:
    """加载XP-CLR报告文件(跨群体,像Fst)|Load XP-CLR report files (cross-pop, like Fst)

    - 每文件一条染色体(TSV,列含 chrom/pos_start/pos_stop/xpclr)
    - 表头行(首列非数字)跳过;空/NaN的xpclr行丢弃
    - 异常格式抛错(由调用方捕获降级)|malformed files raise (caller downgrades)

    Args:
        files: 每染色体TSV路径列表(同一群体对)|per-chrom TSV paths (same pop pair)
        label: 群体对标签(如 a_b)|population-pair label (e.g. a_b)
        win: 窗口大小|Window size

    Returns:
        DataFrame[CHR, WIN, XPCLR_{label}],每窗口取xpclr最大值|xpclr max per window
    """
    frames = []
    for f in files:
        d = pd.read_csv(f, sep='\t')
        if 'chrom' not in d.columns or 'xpclr' not in d.columns:
            raise ValueError(
                f"XP-CLR 报告缺关键列(chrom/xpclr)|XP-CLR report missing required "
                f"columns (chrom/xpclr): {f}")
        d = d.rename(columns={'chrom': 'CHR', 'pos_start': 'POS', 'xpclr': 'XPCLR'})
        # xpclr为2*(modelL-nullL);NaN(无有效SNP窗口)丢弃
        # |xpclr = 2*(modelL-nullL); drop NaN rows (windows with no valid SNPs)
        d = d[['CHR', 'POS', 'XPCLR']].dropna(subset=['XPCLR'])
        if not d.empty:
            frames.append(d)
    if not frames:
        return pd.DataFrame(columns=['CHR', 'WIN', f'XPCLR_{label}'])
    df = pd.concat(frames, ignore_index=True)
    df['WIN'] = bin_pos(df['POS'], win)
    df = df.groupby(['CHR', 'WIN'], as_index=False)[['XPCLR']].max()
    df = df.rename(columns={'XPCLR': f'XPCLR_{label}'})
    return df


class SweepMerger:
    """选择性扫荡统计合并器|Selective Sweep Statistics Merger"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.win = config.win

    # ---- 输入自动发现|input auto-discovery ----
    def discover_inputs(self) -> Dict[str, List[Tuple]]:
        """从stats_dir自动发现输入文件|Auto-discover input files from stats dir

        Returns:
            {'pi': [(pop, path)], 'tajd': [(pop, path)],
             'raisd': [(pop, [files])], 'sweed': [(pop, [files])],
             'fst': [(label, path)]}
        """
        stats_dir = Path(self.config.stats_dir)
        result = {'pi': [], 'tajd': [], 'raisd': [], 'sweed': [], 'xpclr': [], 'fst': []}

        for f in sorted(stats_dir.glob('*_windowed.pi')):
            result['pi'].append((f.name[: -len('_windowed.pi')], f))
        for f in sorted(stats_dir.glob('*_Tajima.D')):
            result['tajd'].append((f.name[: -len('_Tajima.D')], f))
        for f in sorted(stats_dir.glob('*_windowed.weir.fst')):
            result['fst'].append((f.name[: -len('_windowed.weir.fst')], f))
        # RAiSD_Report_<pop>_<chr>:前缀后、最后下划线后缀之前全部为pop
        # |pop = everything between "RAiSD_Report_" and last underscore suffix
        raisd_grouped = {}
        for f in sorted(stats_dir.glob('RAiSD_Report_*')):
            pop = f.name[len('RAiSD_Report_'):].rsplit('_', 1)[0]
            raisd_grouped.setdefault(pop, []).append(f)
        result['raisd'] = [(pop, files) for pop, files in sorted(raisd_grouped.items())]
        # SweeD_Report_<pop>_<chr>:与RAiSD相同的前缀后pop解析|same pop parsing
        sweed_grouped = {}
        for f in sorted(stats_dir.glob('SweeD_Report_*')):
            pop = f.name[len('SweeD_Report_'):].rsplit('_', 1)[0]
            sweed_grouped.setdefault(pop, []).append(f)
        result['sweed'] = [(pop, files) for pop, files in sorted(sweed_grouped.items())]
        # XPCLR_<a>_<b>_<chr>.tsv:跨群体对,前缀后、最后下划线后缀之前为群体对标签
        # |XPCLR_<a>_<b>_<chr>.tsv: cross-pop pair label between prefix and last underscore
        xpclr_grouped = {}
        for f in sorted(stats_dir.glob('XPCLR_*.tsv')):
            stem = f.name[: -len('.tsv')]  # XPCLR_<a>_<b>_<chr>
            label = stem[len('XPCLR_'):].rsplit('_', 1)[0]  # <a>_<b>
            xpclr_grouped.setdefault(label, []).append(f)
        result['xpclr'] = [(label, files) for label, files in sorted(xpclr_grouped.items())]

        for kind, label in [('pi', 'π'), ('tajd', "Tajima's D"),
                            ('raisd', 'RAiSD'), ('sweed', 'SweeD'),
                            ('xpclr', 'XP-CLR'), ('fst', 'Fst')]:
            if not result[kind]:
                self.logger.warning(f"未发现{label}输入文件|No {label} input files found in {stats_dir}")
        return result

    # ---- composite_score|composite score ----
    def build_composite_score(self, merged: pd.DataFrame,
                              pop_summary: Dict[str, int]):
        """构建composite_score|Build composite score

        Args:
            merged: 合并后的窗口表|Merged window table
            pop_summary: pop -> 样本数|pop -> sample count

        Returns:
            (score_df, merged): score_df为各分量rank列;merged追加composite_score/n_stats_supporting
        """
        score_components = []

        # π分量:双(多)群体取前两个PI列,比值极端性;单群体低π
        # |pi: first two PI cols ratio extremity; single pop low-pi
        pi_cols = [c for c in merged.columns if c.startswith('PI_')]
        if len(pi_cols) >= 2:
            a, b = pi_cols[0], pi_cols[1]
            merged['PI_ratio'] = merged[a] / merged[b].replace(0, np.nan)
            merged['PI_ratio_extremity'] = np.log2(merged['PI_ratio']).abs()
            score_components.append(rank_pct(merged['PI_ratio_extremity']))
            self.logger.info(f"π分量|Pi component: {a}/{b} 比值极端性|ratio extremity")
        elif len(pi_cols) == 1:
            score_components.append(rank_pct(merged[pi_cols[0]], ascending=False))
            self.logger.info(f"π分量|Pi component: 单群体低π|single-pop low pi")

        # Tajima's D:越负排名越高|more negative ranks higher
        tajd_cols = [c for c in merged.columns if c.startswith('TajD_')]
        for c in tajd_cols:
            score_components.append(rank_pct(merged[c], ascending=False))
        if tajd_cols:
            self.logger.info(f"Tajima's D分量|TajimaD components: {len(tajd_cols)} 个|count")

        # RAiSD μ:越高排名越高;低样本群体默认排除
        # |mu: higher ranks higher; low-n pops excluded by default
        mu_cols = [c for c in merged.columns if c.startswith('MU_')]
        for c in mu_cols:
            pop = c[len('MU_'):]
            if (pop in pop_summary and pop_summary[pop] < self.config.raisd_min_samples
                    and not self.config.include_mu_low_n):
                self.logger.warning(
                    f"群体{pop}样本量{pop_summary[pop]}<{self.config.raisd_min_samples},"
                    f"MU分量已从composite_score排除(--include-mu-low-n可强制加入)"
                    f"|Pop {pop} has {pop_summary[pop]} samples < "
                    f"{self.config.raisd_min_samples}; MU component excluded from "
                    f"composite score (--include-mu-low-n to force)")
                continue
            score_components.append(rank_pct(merged[c], ascending=True))
        if mu_cols:
            self.logger.info(f"RAiSD μ分量|RAiSD mu components: {len(mu_cols)} 个|count")

        # SweeD CLR:越高排名越高;低样本群体默认排除
        # |CLR: higher ranks higher; low-n pops excluded by default
        clr_cols = [c for c in merged.columns if c.startswith('CLR_')]
        for c in clr_cols:
            pop = c[len('CLR_'):]
            if (pop in pop_summary and pop_summary[pop] < self.config.sweed_min_samples
                    and not self.config.include_sweed_low_n):
                self.logger.warning(
                    f"群体{pop}样本量{pop_summary[pop]}<{self.config.sweed_min_samples},"
                    f"CLR分量已从composite_score排除(--include-sweed-low-n可强制加入)"
                    f"|Pop {pop} has {pop_summary[pop]} samples < "
                    f"{self.config.sweed_min_samples}; CLR component excluded from "
                    f"composite score (--include-sweed-low-n to force)")
                continue
            score_components.append(rank_pct(merged[c], ascending=True))
        if clr_cols:
            self.logger.info(f"SweeD CLR分量|SweeD CLR components: {len(clr_cols)} 个|count")

        # XP-CLR:跨群体CLR,越高排名越高;任一群体低样本默认排除
        # |XP-CLR: cross-pop CLR, higher ranks higher; excluded if either pop low-n
        xpclr_cols = [c for c in merged.columns if c.startswith('XPCLR_')]
        for c in xpclr_cols:
            label = c[len('XPCLR_'):]
            # 跨群体标签 a_b:两群体样本量任一低于阈值则排除
            # |cross-pop label a_b: exclude if either pop's n is below threshold
            pops = label.split('_')
            n_vals = [pop_summary.get(p) for p in pops]
            low = any(n is not None and n < self.config.xpclr_min_samples for n in n_vals)
            if low and not self.config.include_xpclr_low_n:
                self.logger.warning(
                    f"群体对{label}含低样本群体({pops}={[n for n in n_vals if n is not None]}),"
                    f"XP-CLR分量已从composite_score排除(--include-xpclr-low-n可强制加入)"
                    f"|Pair {label} includes a low-n pop ({pops}); "
                    f"XP-CLR component excluded from composite score "
                    f"(--include-xpclr-low-n to force)")
                continue
            score_components.append(rank_pct(merged[c], ascending=True))
        if xpclr_cols:
            self.logger.info(f"XP-CLR分量|XP-CLR components: {len(xpclr_cols)} 个|count")

        # Fst:越高排名越高|higher Fst ranks higher
        fst_cols = [c for c in merged.columns if c.startswith('Fst_')]
        for c in fst_cols:
            score_components.append(rank_pct(merged[c], ascending=True))
        if fst_cols:
            self.logger.info(f"Fst分量|Fst components: {len(fst_cols)} 个|count")

        if not score_components:
            raise RuntimeError(
                "没有可用的统计分量构建composite_score|No usable statistics to build composite score")

        score_df = pd.concat(score_components, axis=1)
        merged['composite_score'] = score_df.mean(axis=1, skipna=True)
        merged['n_stats_supporting'] = score_df.notna().sum(axis=1)
        return score_df, merged

    # ---- 候选区域|candidate regions ----
    def find_candidate_regions(self, merged: pd.DataFrame) -> pd.DataFrame:
        """候选区域:top分位数窗口合并|Candidate regions: top-quantile windows merged

        同染色体且相邻窗口间隔<=merge_gap时合并|Merge same-chr windows with gap <= merge_gap
        """
        thr = merged['composite_score'].quantile(1 - self.config.top_quantile)
        cand = merged[merged['composite_score'] >= thr].copy()
        cand = cand.sort_values(['CHR', 'WIN']).reset_index(drop=True)
        self.logger.info(
            f"composite_score阈值|Threshold (top {self.config.top_quantile*100:.1f}%): {thr:.4f}")
        if cand.empty:
            self.logger.warning("无候选窗口|No candidate windows above threshold")
            return pd.DataFrame(columns=['CHR', 'START', 'END', 'MAX_SCORE', 'N_WIN'])

        regions = []
        cur = None
        for _, row in cand.iterrows():
            if cur is None:
                cur = {'CHR': row['CHR'], 'START': row['WIN'],
                       'END': row['WIN'] + self.win,
                       'MAX_SCORE': row['composite_score'], 'N_WIN': 1}
            elif (row['CHR'] == cur['CHR']
                  and row['WIN'] - cur['END'] <= self.config.merge_gap):
                cur['END'] = row['WIN'] + self.win
                cur['MAX_SCORE'] = max(cur['MAX_SCORE'], row['composite_score'])
                cur['N_WIN'] += 1
            else:
                regions.append(cur)
                cur = {'CHR': row['CHR'], 'START': row['WIN'],
                       'END': row['WIN'] + self.win,
                       'MAX_SCORE': row['composite_score'], 'N_WIN': 1}
        if cur is not None:
            regions.append(cur)
        return pd.DataFrame(regions)

    # ---- 主入口|main entry ----
    def run(self) -> pd.DataFrame:
        """运行合并打分|Run merge and scoring"""
        self.logger.info("开始合并统计量|Starting statistics merge")

        inputs = self.discover_inputs()
        tables = []
        for pop, path in inputs['pi']:
            try:
                tables.append(load_pi(path, pop, self.win))
            except Exception as e:
                self.logger.warning(f"π文件解析失败,跳过|Failed to parse pi file, skipped: {path}: {e}")
        for pop, path in inputs['tajd']:
            try:
                tables.append(load_tajd(path, pop, self.win))
            except Exception as e:
                self.logger.warning(
                    f"TajimaD文件解析失败,跳过|Failed to parse TajimaD file, skipped: {path}: {e}")
        for pop, files in inputs['raisd']:
            try:
                tables.append(load_raisd(files, pop, self.win))
            except Exception as e:
                self.logger.warning(
                    f"RAiSD报告解析失败,跳过|Failed to parse RAiSD report, skipped: {pop}: {e}")
        for pop, files in inputs['sweed']:
            try:
                tables.append(load_sweed(files, pop, self.win))
            except Exception as e:
                self.logger.warning(
                    f"SweeD报告解析失败,跳过|Failed to parse SweeD report, skipped: {pop}: {e}")
        for label, files in inputs['xpclr']:
            try:
                tables.append(load_xpclr(files, label, self.win))
            except Exception as e:
                self.logger.warning(
                    f"XP-CLR报告解析失败,跳过|Failed to parse XP-CLR report, skipped: {label}: {e}")
        for label, path in inputs['fst']:
            try:
                tables.append(load_fst(path, label, self.win))
            except Exception as e:
                self.logger.warning(
                    f"Fst文件解析失败,跳过|Failed to parse Fst file, skipped: {path}: {e}")

        if not tables:
            raise RuntimeError(
                "02_stats/中无任何可用统计输入文件|No usable statistics input files in 02_stats/")

        merged = tables[0]
        for t in tables[1:]:
            merged = pd.merge(merged, t, on=['CHR', 'WIN'], how='outer')
        merged = merged.sort_values(['CHR', 'WIN']).reset_index(drop=True)
        self.logger.info(f"合并完成|Merged: {len(merged)} 个窗口|windows")

        # 读取群体样本量摘要(用于低样本MU排除)|read pop sample summary
        pop_summary = {}
        summary_file = Path(self.config.info_dir) / 'pop_summary.tsv'
        if summary_file.exists():
            try:
                sdf = pd.read_csv(summary_file, sep='\t')
                pop_summary = dict(zip(sdf['pop'], sdf['n_samples']))
            except Exception as e:
                self.logger.warning(f"读取pop_summary.tsv失败|Failed to read pop_summary.tsv: {e}")

        _, merged = self.build_composite_score(merged, pop_summary)

        # 全基因组表|genome-wide table
        genome_wide = Path(self.config.sweep_dir) / 'genome_wide_stats.tsv'
        merged.to_csv(genome_wide, sep='\t', index=False)
        self.logger.info(
            f"全基因组统计表|Genome-wide table -> {genome_wide} ({len(merged)} 窗口|windows)")

        # 候选区域|candidate regions
        regions = self.find_candidate_regions(merged)
        regions_file = Path(self.config.sweep_dir) / 'candidate_regions.tsv'
        regions.to_csv(regions_file, sep='\t', index=False)
        self.logger.info(
            f"候选区域|Candidate regions -> {regions_file} ({len(regions)} 个|regions)")

        return merged

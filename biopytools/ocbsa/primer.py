"""
OcBSA - BSA候选区域引物设计|OcBSA - BSA Candidate Region Primer Design

基于BSA分析结果设计分子标记引物，依赖primer3-py和BLAST。
"""

import os
import subprocess


class BsaPrimerDesigner:
    """BSA引物设计器|BSA Primer Designer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def _read_genome_chr(self, chr_name):
        """从基因组FASTA中读取指定染色体序列|Read chromosome sequence from genome FASTA"""
        chr_seq = ''
        flag = False
        with open(self.config.genome, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if line.strip()[1:] == chr_name:
                        flag = True
                    else:
                        flag = False
                elif flag:
                    chr_seq += line.strip()
        return chr_seq

    def pick_indel_sequences(self):
        """从基因组中提取INDEL侧翼序列|Extract INDEL flanking sequences from genome"""
        cfg = self.config
        chr_name = cfg.region.split(',')[0]
        start = int(cfg.region.split(',')[1])
        end = int(cfg.region.split(',')[2])
        flank = cfg.flank_length

        self.logger.info(f"提取INDEL侧翼序列|Extracting INDEL flanking sequences: {chr_name}:{start}-{end}")
        chr_seq = self._read_genome_chr(chr_name)
        if not chr_seq:
            self.logger.error(f"未找到染色体序列|Chromosome sequence not found: {chr_name}")
            return None, None

        out_dhhp_dict = {}
        indel_file = os.path.join(cfg.output_dir, f"{chr_name}_{start}_{end}_indel.fasta")
        with open(cfg.output.ocvalue_file, 'r') as dhhp_file, open(indel_file, 'w') as out_fasta:
            for line in dhhp_file:
                line_list = line.strip().split('\t')
                if line_list[0] == chr_name and start < int(line_list[1]) < end:
                    out_dhhp_dict[line_list[0] + '_' + line_list[1]] = line_list[:8] + [line_list[-2]]

                    if len(line_list[2]) > 4:
                        out_seq = (chr_seq[int(line_list[1]) - flank: int(line_list[1])]
                                  + line_list[2]
                                  + chr_seq[int(line_list[1]) + 1: int(line_list[1]) + flank + 1])
                    elif len(line_list[3]) > 4:
                        out_seq = (chr_seq[int(line_list[1]) - flank: int(line_list[1])]
                                  + line_list[3]
                                  + chr_seq[int(line_list[1]) + 1: int(line_list[1]) + flank + 1])
                    else:
                        continue
                    out_fasta.write(f">{chr_name}_{line_list[1]}\n{out_seq}\n")

        self.logger.info(f"提取到|Extracted {len(out_dhhp_dict)} 个INDEL位点")
        return indel_file, out_dhhp_dict

    def filter_by_blast(self, indel_fasta):
        """使用BLAST过滤非唯一序列|Filter non-unique sequences using BLAST"""
        cfg = self.config
        db_dir = os.path.join(cfg.output_dir, 'db')
        os.makedirs(db_dir, exist_ok=True)

        genome_db = os.path.join(db_dir, os.path.basename(cfg.genome))
        blast_out = os.path.join(cfg.output_dir, 'indel_genome.blast')

        # makeblastdb
        cmd_makedb = ['makeblastdb', '-in', cfg.genome, '-dbtype', 'nucl', '-parse_seqids', '-out', genome_db]
        self.logger.info(f"执行|Executing: 构建BLAST数据库|Build BLAST database")
        self.logger.info(f"命令|Command: {' '.join(cmd_makedb)}")
        subprocess.run(cmd_makedb, capture_output=True)

        # blastn
        cmd_blast = ['blastn', '-num_threads', '30', '-outfmt', '6', '-evalue', '1e-5',
                      '-db', genome_db, '-query', indel_fasta, '-out', blast_out]
        self.logger.info(f"执行|Executing: BLAST比对|BLAST alignment")
        self.logger.info(f"命令|Command: {' '.join(cmd_blast)}")
        subprocess.run(cmd_blast, capture_output=True)

        # 统计唯一比对|Count unique hits
        count_dict = {}
        with open(blast_out, 'r') as f:
            for line in f:
                line_list = line.strip().split('\t')
                if int(line_list[3]) < 60:
                    continue
                if line_list[0] not in count_dict:
                    count_dict[line_list[0]] = 0
                count_dict[line_list[0]] += 1

        uniq_list = [k for k, v in count_dict.items() if v == 1]
        self.logger.info(f"BLAST过滤|BLAST filter: {len(uniq_list)}/{len(count_dict)} 条唯一序列")

        # 输出过滤后的FASTA|Write filtered FASTA
        filtered_fasta = os.path.join(cfg.output_dir,
                                       f"{os.path.basename(indel_fasta).replace('.fasta', '')}_filter.fasta")
        out_dict = {}
        with open(indel_fasta, 'r') as f_in, open(filtered_fasta, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    seq_id = line.strip()[1:]
                    if seq_id in uniq_list:
                        f_out.write(line)
                        out_dict[seq_id] = ''
                elif seq_id in uniq_list:
                    out_dict[seq_id] += line.strip()
                    f_out.write(line)

        return out_dict, filtered_fasta

    def design_primers(self, seq_dict, dhhp_dict):
        """使用primer3设计引物|Design primers using primer3"""
        import primer3

        cfg = self.config
        global_args = {
            'PRIMER_NUM_RETURN': cfg.primer_num,
            'PRIMER_OPT_SIZE': cfg.primer_opt_size,
            'PRIMER_MIN_SIZE': cfg.primer_min_size,
            'PRIMER_MAX_SIZE': cfg.primer_max_size,
            'PRIMER_OPT_TM': 59.0,
            'PRIMER_MIN_TM': cfg.min_tm,
            'PRIMER_MAX_TM': cfg.max_tm,
            'PRIMER_MIN_GC': cfg.min_gc,
            'PRIMER_MAX_GC': cfg.max_gc,
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [cfg.product_min, cfg.product_max],
            'PRIMER_GC_CLAMP': 1,
            'PRIMER_PAIR_MAX_DIFF_TM': cfg.tm_diff
        }

        output_file = os.path.join(cfg.output_dir, 'primer_design_results.txt')
        flank = cfg.flank_length

        self.logger.info("开始引物设计|Starting primer design")
        with open(output_file, 'w') as out_f:
            header = ('Index\tChr\tPos\tRef\tAlt\tParent1\tParent2\t'
                      'Pool1\tPool2\tDhhp\tLeft_seq\tLeft_TM\tLeft_GC\t'
                      'Right_seq\tRight_TM\tRight_GC\tProduct_size\n')
            out_f.write(header)

            index = 0
            for key, value in seq_dict.items():
                seq_args = {
                    'SEQUENCE_ID': key,
                    'SEQUENCE_TEMPLATE': value,
                    'SEQUENCE_INCLUDED_REGION': [0, len(value) - 1],
                    'SEQUENCE_TARGET': [flank, len(value) - 2 * flank]
                }
                index += 1
                result = primer3.bindings.designPrimers(seq_args, global_args)

                for i in range(result['PRIMER_LEFT_NUM_RETURNED']):
                    out_f.write(f'primer_{index}_{i + 1}\t')
                    for x in dhhp_dict[key]:
                        out_f.write(f'{x}\t')
                    out_f.write(f"{result['PRIMER_LEFT_' + str(i) + '_SEQUENCE']}\t")
                    out_f.write(f"{result['PRIMER_LEFT_' + str(i) + '_TM']}\t")
                    out_f.write(f"{result['PRIMER_LEFT_' + str(i) + '_GC_PERCENT']}\t")
                    out_f.write(f"{result['PRIMER_RIGHT_' + str(i) + '_SEQUENCE']}\t")
                    out_f.write(f"{result['PRIMER_RIGHT_' + str(i) + '_TM']}\t")
                    out_f.write(f"{result['PRIMER_RIGHT_' + str(i) + '_GC_PERCENT']}\t")
                    out_f.write(f"{result['PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE']}\n")

        self.logger.info(f"引物设计结果已保存|Primer design results saved to: {output_file}")
        return True

    def run(self):
        """运行引物设计流程|Run primer design pipeline"""
        self.config.validate()

        self.logger.info("=" * 60)
        self.logger.info("BSA候选区域引物设计|BSA Candidate Region Primer Design")
        self.logger.info("=" * 60)
        self.logger.info(f"参考基因组|Reference genome: {self.config.genome}")
        self.logger.info(f"目标区间|Target region: {self.config.region}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        # 提取INDEL序列|Extract INDEL sequences
        indel_fasta, dhhp_dict = self.pick_indel_sequences()
        if not dhhp_dict:
            self.logger.error("未找到INDEL位点|No INDEL sites found")
            return False

        # BLAST过滤|BLAST filter
        seq_dict, filtered_fasta = self.filter_by_blast(indel_fasta)

        # 引物设计|Primer design
        self.design_primers(seq_dict, dhhp_dict)

        self.logger.info("引物设计完成|Primer design completed")
        return True

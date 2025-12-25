"""
🧬 祖先状态重建结果解析模块 | Ancestral State Reconstruction Result Parser Module
"""

import os
import pandas as pd
from pathlib import Path
from collections import defaultdict
try:
    from Bio import Phylo
    from io import StringIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

class AncestralStateParser:
    """🔬 祖先状态解析器 | Ancestral State Parser"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.state_file = None
        self.data = None
        
    def parse_ancestral_results(self):
        """解析祖先状态重建结果 | Parse ancestral state reconstruction results"""
        asr_prefix = self.config.output_path / f"{self.config.prefix}_ancestral"
        self.state_file = f"{asr_prefix}.state"
        
        if not os.path.exists(self.state_file):
            self.logger.warning(f"⚠️ 未找到祖先状态文件 | Ancestral state file not found: {self.state_file}")
            return False
        
        self.logger.info("=" * 60)
        self.logger.info("🔬 开始解析祖先状态重建结果 | Parsing ancestral state reconstruction results")
        self.logger.info("=" * 60)
        
        try:
            # 读取数据 | Read data
            self.logger.info(f"📖 读取文件 | Reading file: {self.state_file}")
            self.data = pd.read_csv(self.state_file, sep='\t', comment='#')
            
            # 生成各种分析 | Generate various analyses
            self._extract_sequences()
            self._find_uncertain_sites()
            self._analyze_evolution_hotspots()
            self._generate_summary_report()
            self._parse_tree_structure()
            self._trace_evolutionary_history()
            
            self.logger.info("✅ 祖先状态解析完成 | Ancestral state parsing completed")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 解析祖先状态文件时出错 | Error parsing ancestral state file: {e}")
            return False
    
    def _extract_sequences(self):
        """提取每个节点的完整序列 | Extract complete sequence for each node"""
        self.logger.info("📝 提取祖先序列 | Extracting ancestral sequences")
        
        # 按节点分组 | Group by node
        nodes = self.data['Node'].unique()
        self.logger.info(f"📊 共有 {len(nodes)} 个节点 | Total nodes: {len(nodes)}")
        
        # 输出文件 | Output file
        seq_file = self.config.output_path / f"{self.config.prefix}_ancestral_sequences.fasta"
        
        with open(seq_file, 'w') as f:
            for node in nodes:
                node_data = self.data[self.data['Node'] == node].sort_values('Site')
                sequence = ''.join(node_data['State'].values)
                
                # 写入FASTA格式 | Write in FASTA format
                f.write(f">{node}\n")
                # 每60个字符换行 | Line break every 60 characters
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')
        
        self.logger.info(f"✅ 祖先序列已保存 | Ancestral sequences saved: {seq_file}")
        
        # 显示前几个节点的序列长度 | Show sequence lengths for first few nodes
        self.logger.info("📏 序列长度统计 | Sequence length statistics:")
        for node in list(nodes)[:5]:
            node_data = self.data[self.data['Node'] == node]
            seq_len = len(node_data)
            self.logger.info(f"  • {node}: {seq_len} 个位点 | sites")
        
        if len(nodes) > 5:
            self.logger.info(f"  ... 还有 {len(nodes)-5} 个节点 | and {len(nodes)-5} more nodes")
    
    def _find_uncertain_sites(self, threshold=0.95):
        """找出不确定的位点 | Find uncertain sites"""
        self.logger.info(f"🔍 查找不确定位点 (概率 < {threshold}) | Finding uncertain sites (probability < {threshold})")
        
        # 获取所有概率列 | Get all probability columns
        prob_cols = [col for col in self.data.columns if col.startswith('p_')]
        
        # 找出最大概率 | Find maximum probability
        self.data['max_prob'] = self.data[prob_cols].max(axis=1)
        
        # 筛选不确定的位点 | Filter uncertain sites
        uncertain = self.data[self.data['max_prob'] < threshold]
        
        if len(uncertain) > 0:
            self.logger.info(f"⚠️ 发现 {len(uncertain)} 个不确定位点 | Found {len(uncertain)} uncertain sites")
            
            # 保存不确定位点 | Save uncertain sites
            uncertain_file = self.config.output_path / f"{self.config.prefix}_uncertain_sites.txt"
            uncertain[['Node', 'Site', 'State', 'max_prob']].to_csv(
                uncertain_file, sep='\t', index=False
            )
            self.logger.info(f"📄 不确定位点已保存 | Uncertain sites saved: {uncertain_file}")
            
            # 显示前10个不确定位点 | Show top 10 uncertain sites
            self.logger.info("🔝 最不确定的10个位点 | Top 10 most uncertain sites:")
            top_uncertain = uncertain.nsmallest(10, 'max_prob')[['Node', 'Site', 'State', 'max_prob']]
            for _, row in top_uncertain.iterrows():
                self.logger.info(f"  • {row['Node']} 位点{row['Site']}: {row['State']} (概率: {row['max_prob']:.4f})")
        else:
            self.logger.info(f"✅ 所有位点概率都 ≥ {threshold}，重建结果非常可靠 | All sites have probability ≥ {threshold}, very reliable")
    
    def _analyze_evolution_hotspots(self):
        """分析演化热点位点 | Analyze evolution hotspots"""
        self.logger.info("🌡️ 分析演化热点 | Analyzing evolution hotspots")
        
        # 统计每个位点有多少种不同的状态 | Count different states per site
        site_diversity = self.data.groupby('Site')['State'].apply(lambda x: len(set(x)))
        
        # 找出变化最多的位点 | Find most variable sites
        hotspots = site_diversity[site_diversity > 1].sort_values(ascending=False)
        
        if len(hotspots) > 0:
            self.logger.info(f"🔥 发现 {len(hotspots)} 个演化热点位点 | Found {len(hotspots)} evolution hotspot sites")
            
            # 保存热点位点 | Save hotspots
            hotspot_file = self.config.output_path / f"{self.config.prefix}_evolution_hotspots.txt"
            hotspots.to_csv(hotspot_file, sep='\t', header=['Num_States'])
            self.logger.info(f"📄 演化热点已保存 | Evolution hotspots saved: {hotspot_file}")
            
            # 显示前10个热点 | Show top 10 hotspots
            self.logger.info("🔝 变化最多的10个位点 | Top 10 most variable sites:")
            for site, count in list(hotspots.head(10).items()):
                states = self.data[self.data['Site'] == site]['State'].unique()
                self.logger.info(f"  • 位点 {site}: {count} 种状态 | states ({', '.join(states)})")
        else:
            self.logger.info("ℹ️ 所有位点在所有节点上都是保守的 | All sites are conserved across all nodes")
    
    def _generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        self.logger.info("📊 生成总结报告 | Generating summary report")
        
        report_file = self.config.output_path / f"{self.config.prefix}_ancestral_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("🧬 祖先状态重建分析总结 | Ancestral State Reconstruction Summary\n")
            f.write("=" * 70 + "\n\n")
            
            # 基本统计 | Basic statistics
            f.write("📊 基本统计信息 | Basic Statistics\n")
            f.write("-" * 70 + "\n")
            
            nodes = self.data['Node'].unique()
            sites = self.data['Site'].unique()
            
            f.write(f"节点数量 | Number of nodes: {len(nodes)}\n")
            f.write(f"位点数量 | Number of sites: {len(sites)}\n")
            f.write(f"序列类型 | Sequence type: ")
            
            # 判断序列类型 | Determine sequence type
            states = self.data['State'].unique()
            if len(states) == 4 and all(s in ['A', 'T', 'C', 'G'] for s in states):
                f.write("DNA\n")
            elif len(states) <= 20:
                f.write("蛋白质/Protein (氨基酸)\n")
            else:
                f.write(f"其他 | Other ({len(states)} 种状态)\n")
            
            f.write(f"\n观察到的状态类型 | Observed states: {', '.join(sorted(states))}\n")
            
            # 概率统计 | Probability statistics
            f.write("\n📈 概率分布统计 | Probability Distribution Statistics\n")
            f.write("-" * 70 + "\n")
            
            avg_prob = self.data['max_prob'].mean()
            min_prob = self.data['max_prob'].min()
            
            f.write(f"平均最大概率 | Average maximum probability: {avg_prob:.4f}\n")
            f.write(f"最小最大概率 | Minimum maximum probability: {min_prob:.4f}\n")
            
            # 可靠性评估 | Reliability assessment
            high_conf = len(self.data[self.data['max_prob'] >= 0.95])
            medium_conf = len(self.data[(self.data['max_prob'] >= 0.80) & (self.data['max_prob'] < 0.95)])
            low_conf = len(self.data[self.data['max_prob'] < 0.80])
            
            f.write(f"\n高置信度位点 (≥0.95) | High confidence (≥0.95): {high_conf} ({high_conf/len(self.data)*100:.2f}%)\n")
            f.write(f"中等置信度位点 (0.80-0.95) | Medium confidence (0.80-0.95): {medium_conf} ({medium_conf/len(self.data)*100:.2f}%)\n")
            f.write(f"低置信度位点 (<0.80) | Low confidence (<0.80): {low_conf} ({low_conf/len(self.data)*100:.2f}%)\n")
            
            # 状态频率统计 | State frequency statistics
            f.write("\n🔤 氨基酸/碱基频率 | Amino acid/Base frequencies\n")
            f.write("-" * 70 + "\n")
            
            state_counts = self.data['State'].value_counts()
            for state, count in state_counts.head(10).items():
                f.write(f"{state}: {count} ({count/len(self.data)*100:.2f}%)\n")
            
            # 输出文件列表 | Output files list
            f.write("\n📁 生成的文件 | Generated files\n")
            f.write("-" * 70 + "\n")
            f.write(f"1. 祖先序列 | Ancestral sequences: {self.config.prefix}_ancestral_sequences.fasta\n")
            
            uncertain_file = self.config.output_path / f"{self.config.prefix}_uncertain_sites.txt"
            if uncertain_file.exists():
                f.write(f"2. 不确定位点 | Uncertain sites: {self.config.prefix}_uncertain_sites.txt\n")
            
            hotspot_file = self.config.output_path / f"{self.config.prefix}_evolution_hotspots.txt"
            if hotspot_file.exists():
                f.write(f"3. 演化热点 | Evolution hotspots: {self.config.prefix}_evolution_hotspots.txt\n")
            
            f.write(f"4. 总结报告 | Summary report: {self.config.prefix}_ancestral_summary.txt\n")
            
            # 使用建议 | Usage recommendations
            f.write("\n💡 使用建议 | Usage Recommendations\n")
            f.write("-" * 70 + "\n")
            f.write("1. 查看 *_ancestral_sequences.fasta 获取完整的祖先序列\n")
            f.write("2. 关注 *_uncertain_sites.txt 中的低置信度位点\n")
            f.write("3. 分析 *_evolution_hotspots.txt 了解哪些位点经历了较多演化变化\n")
            f.write("4. 可以使用这些祖先序列进行进一步的功能分析或实验验证\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        self.logger.info(f"✅ 总结报告已保存 | Summary report saved: {report_file}")
        
        # 在日志中显示关键统计 | Show key statistics in log
        nodes = self.data['Node'].unique()
        sites = self.data['Site'].unique()
        avg_prob = self.data['max_prob'].mean()
        
        self.logger.info(f"📊 关键统计 | Key Statistics:")
        self.logger.info(f"  • 节点数 | Nodes: {len(nodes)}")
        self.logger.info(f"  • 位点数 | Sites: {len(sites)}")
        self.logger.info(f"  • 平均置信度 | Average confidence: {avg_prob:.4f}")

    def _parse_tree_structure(self):
        """解析树结构，生成CSV格式的节点-物种映射"""
        self.logger.info("🌳 解析树拓扑结构 | Parsing tree topology")
        
        asr_prefix = self.config.output_path / f"{self.config.prefix}_ancestral"
        tree_file = f"{asr_prefix}.treefile"
        
        if not os.path.exists(tree_file):
            self.logger.warning(f"⚠️ 未找到树文件 | Tree file not found: {tree_file}")
            return
        
        if not HAS_BIOPYTHON:
            self.logger.warning("⚠️ 未安装BioPython，无法解析树结构 | BioPython not installed")
            self.logger.info("💡 安装方法：pip install biopython")
            return
        
        try:
            tree = Phylo.read(tree_file, "newick")
            
            # 获取所有叶子节点（物种）
            all_leaves = tree.get_terminals()
            
            # 先找出最大深度
            max_depth = 0
            species_paths = []
            
            for leaf in all_leaves:
                path = tree.get_path(leaf)
                node_names = [n.name for n in path if n.name and not n.is_terminal()]
                max_depth = max(max_depth, len(node_names))
                species_paths.append({
                    'species': leaf.name,
                    'path': node_names  # 从直接父节点到根节点
                })
            
            # 构建数据框
            rows = []
            for item in species_paths:
                row = {}
                # 填充节点列（从Level1到LevelN）
                for i in range(max_depth):
                    if i < len(item['path']):
                        row[f'Level{i+1}'] = item['path'][i]
                    else:
                        row[f'Level{i+1}'] = ''  # 如果深度不够，填空
                row['Species'] = item['species']
                rows.append(row)
            
            # 保存为CSV
            df = pd.DataFrame(rows)
            # 确保列的顺序：Level1, Level2, ..., Species
            columns = [f'Level{i+1}' for i in range(max_depth)] + ['Species']
            df = df[columns]
            
            lineage_file = self.config.output_path / f"{self.config.prefix}_species_lineage.csv"
            df.to_csv(lineage_file, index=False)
            
            self.logger.info(f"✅ 物种路径表已保存 | Species lineage saved: {lineage_file}")
            self.logger.info(f"📊 树的最大深度 | Maximum tree depth: {max_depth}")
            self.logger.info(f"📊 总物种数 | Total species: {len(all_leaves)}")
            
            # 显示前几行示例
            self.logger.info("📋 前3个物种的路径示例:")
            for i in range(min(3, len(rows))):
                path_str = ' > '.join([rows[i][f'Level{j+1}'] for j in range(max_depth) if rows[i][f'Level{j+1}']])
                self.logger.info(f"  • {rows[i]['Species']}: {path_str}")
                    
        except Exception as e:
            self.logger.error(f"❌ 解析树结构时出错 | Error parsing tree structure: {e}")
    
    # def _trace_evolutionary_history(self):
    #     """追溯每个物种的演化历史，输出CSV格式"""
    #     self.logger.info("🔬 追溯演化历史 | Tracing evolutionary history")
        
    #     lineage_file = self.config.output_path / f"{self.config.prefix}_species_lineage.csv"
    #     seq_file = self.config.output_path / f"{self.config.prefix}_ancestral_sequences.fasta"
        
    #     if not lineage_file.exists() or not seq_file.exists():
    #         self.logger.warning("⚠️ 缺少必需文件，跳过演化历史追溯")
    #         return
        
    #     try:
    #         from Bio import SeqIO
            
    #         # 读取数据
    #         df_lineage = pd.read_csv(lineage_file)
    #         ancestral_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(seq_file, "fasta")}
    #         modern_seqs = {rec.id: str(rec.seq).replace('-', '') 
    #                     for rec in SeqIO.parse(self.config.input_file, "fasta")}
            
    #         self.logger.info(f"📊 开始追溯 {len(df_lineage)} 个物种的演化历史")
            
    #         all_records = []
            
    #         for idx, row in df_lineage.iterrows():
    #             species = row['Species']
                
    #             # 获取路径
    #             path = [row[col] for col in df_lineage.columns 
    #                 if col.startswith('Level') and pd.notna(row[col]) and row[col] != '']
                
    #             if len(path) == 0 or species not in modern_seqs:
    #                 all_records.append({
    #                     'Species': species,
    #                     'Level': 0,
    #                     'Node': 'Root',
    #                     'Sequence_Length': 0,
    #                     'Num_Mutations': 0,
    #                     'Position': '',
    #                     'From_AA': '',
    #                     'To_AA': '',
    #                     'Note': 'Basal branch'
    #                 })
    #                 continue
                
    #             prev_seq = None
    #             for i, node in enumerate(path):
    #                 if node not in ancestral_seqs:
    #                     continue
                    
    #                 curr_seq = ancestral_seqs[node]
                    
    #                 if prev_seq:
    #                     mutations = self._compare_sequences(prev_seq, curr_seq, node, node)
                        
    #                     if mutations:
    #                         for mut in mutations:
    #                             all_records.append({
    #                                 'Species': species,
    #                                 'Level': i + 1,
    #                                 'Node': node,
    #                                 'Sequence_Length': len(curr_seq),
    #                                 'Num_Mutations': len(mutations),
    #                                 'Position': mut['Position'],
    #                                 'From_AA': mut['From'],
    #                                 'To_AA': mut['To'],
    #                                 'Note': ''
    #                             })
    #                     else:
    #                         # 无突变也记录
    #                         all_records.append({
    #                             'Species': species,
    #                             'Level': i + 1,
    #                             'Node': node,
    #                             'Sequence_Length': len(curr_seq),
    #                             'Num_Mutations': 0,
    #                             'Position': '',
    #                             'From_AA': '',
    #                             'To_AA': '',
    #                             'Note': 'Conserved'
    #                         })
    #                 else:
    #                     # 最原始祖先
    #                     all_records.append({
    #                         'Species': species,
    #                         'Level': i + 1,
    #                         'Node': node,
    #                         'Sequence_Length': len(curr_seq),
    #                         'Num_Mutations': 0,
    #                         'Position': '',
    #                         'From_AA': '',
    #                         'To_AA': '',
    #                         'Note': 'Ancestral'
    #                     })
                    
    #                 prev_seq = curr_seq
                
    #             # 现代序列
    #             modern_seq = modern_seqs[species]
    #             if prev_seq and len(prev_seq) == len(modern_seq):
    #                 mutations = self._compare_sequences(prev_seq, modern_seq, path[-1], species)
                    
    #                 if mutations:
    #                     for mut in mutations:
    #                         all_records.append({
    #                             'Species': species,
    #                             'Level': len(path) + 1,
    #                             'Node': species,
    #                             'Sequence_Length': len(modern_seq),
    #                             'Num_Mutations': len(mutations),
    #                             'Position': mut['Position'],
    #                             'From_AA': mut['From'],
    #                             'To_AA': mut['To'],
    #                             'Note': 'Modern'
    #                         })
            
    #         # 保存CSV
    #         df_trace = pd.DataFrame(all_records)
    #         trace_file = self.config.output_path / f"{self.config.prefix}_evolutionary_trace.csv"
    #         df_trace.to_csv(trace_file, index=False)
            
    #         self.logger.info(f"✅ 演化历史表已保存: {trace_file}")
    #         self.logger.info(f"📊 共记录 {len(all_records)} 条演化事件")
            
    #         # 统计信息
    #         total_mutations = df_trace[df_trace['Position'] != ''].shape[0]
    #         self.logger.info(f"📊 共记录 {total_mutations} 个突变事件")
            
    #     except Exception as e:
    #         self.logger.error(f"❌ 追溯演化历史时出错: {e}")

    def _trace_evolutionary_history(self):
        """追溯每个物种的演化历史，输出CSV和树状结构"""
        self.logger.info("🔬 追溯演化历史 | Tracing evolutionary history")
        
        lineage_file = self.config.output_path / f"{self.config.prefix}_species_lineage.csv"
        seq_file = self.config.output_path / f"{self.config.prefix}_ancestral_sequences.fasta"
        
        if not lineage_file.exists() or not seq_file.exists():
            self.logger.warning("⚠️ 缺少必需文件，跳过演化历史追溯")
            return
        
        try:
            from Bio import SeqIO
            
            # 读取数据
            df_lineage = pd.read_csv(lineage_file)
            ancestral_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(seq_file, "fasta")}
            modern_seqs = {rec.id: str(rec.seq).replace('-', '') 
                        for rec in SeqIO.parse(self.config.input_file, "fasta")}
            
            self.logger.info(f"📊 开始追溯 {len(df_lineage)} 个物种的演化历史")
            
            all_records = []
            
            for idx, row in df_lineage.iterrows():
                species = row['Species']
                
                # 获取路径
                path = [row[col] for col in df_lineage.columns 
                    if col.startswith('Level') and pd.notna(row[col]) and row[col] != '']
                
                if len(path) == 0 or species not in modern_seqs:
                    all_records.append({
                        'Species': species,
                        'Level': 0,
                        'Node': 'Root',
                        'Sequence_Length': 0,
                        'Num_Mutations': 0,
                        'Position': '',
                        'From_AA': '',
                        'To_AA': '',
                        'Note': 'Basal branch'
                    })
                    continue
                
                prev_seq = None
                for i, node in enumerate(path):
                    if node not in ancestral_seqs:
                        continue
                    
                    curr_seq = ancestral_seqs[node]
                    
                    if prev_seq:
                        mutations = self._compare_sequences(prev_seq, curr_seq, node, node)
                        
                        if mutations:
                            for mut in mutations:
                                all_records.append({
                                    'Species': species,
                                    'Level': i + 1,
                                    'Node': node,
                                    'Sequence_Length': len(curr_seq),
                                    'Num_Mutations': len(mutations),
                                    'Position': mut['Position'],
                                    'From_AA': mut['From'],
                                    'To_AA': mut['To'],
                                    'Note': ''
                                })
                        else:
                            all_records.append({
                                'Species': species,
                                'Level': i + 1,
                                'Node': node,
                                'Sequence_Length': len(curr_seq),
                                'Num_Mutations': 0,
                                'Position': '',
                                'From_AA': '',
                                'To_AA': '',
                                'Note': 'Conserved'
                            })
                    else:
                        all_records.append({
                            'Species': species,
                            'Level': i + 1,
                            'Node': node,
                            'Sequence_Length': len(curr_seq),
                            'Num_Mutations': 0,
                            'Position': '',
                            'From_AA': '',
                            'To_AA': '',
                            'Note': 'Ancestral'
                        })
                    
                    prev_seq = curr_seq
                
                # 现代序列
                modern_seq = modern_seqs[species]
                if prev_seq and len(prev_seq) == len(modern_seq):
                    mutations = self._compare_sequences(prev_seq, modern_seq, path[-1], species)
                    
                    if mutations:
                        for mut in mutations:
                            all_records.append({
                                'Species': species,
                                'Level': len(path) + 1,
                                'Node': species,
                                'Sequence_Length': len(modern_seq),
                                'Num_Mutations': len(mutations),
                                'Position': mut['Position'],
                                'From_AA': mut['From'],
                                'To_AA': mut['To'],
                                'Note': 'Modern'
                            })
            
            # 保存CSV
            df_trace = pd.DataFrame(all_records)
            trace_file = self.config.output_path / f"{self.config.prefix}_evolutionary_trace.csv"
            df_trace.to_csv(trace_file, index=False)
            
            self.logger.info(f"✅ 演化历史CSV已保存: {trace_file}")
            
            # 生成树状结构文件
            tree_viz_file = self.config.output_path / f"{self.config.prefix}_mutation_tree.txt"
            
            with open(tree_viz_file, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("🌳 物种演化突变树 | Species Evolutionary Mutation Tree\n")
                f.write("=" * 80 + "\n\n")
                
                for species in df_lineage['Species'].unique():
                    species_data = df_trace[df_trace['Species'] == species].sort_values('Level')
                    
                    if species_data.empty:
                        continue
                    
                    f.write("\n" + "=" * 80 + "\n")
                    f.write(f"物种: {species}\n")
                    f.write("=" * 80 + "\n\n")
                    
                    prev_level = None
                    for idx, row in species_data.iterrows():
                        level = row['Level']
                        node = row['Node']
                        
                        # 只在level变化时输出节点信息
                        if level != prev_level:
                            if prev_level is not None:
                                f.write("    │\n")
                            
                            # 节点名称
                            if row['Note'] == 'Ancestral':
                                f.write(f"{node} (Level {level}) - 最原始祖先\n")
                            elif row['Note'] == 'Modern':
                                f.write(f"{node} - 现代序列\n")
                            elif row['Note'] == 'Basal branch':
                                f.write(f"{node} - 基部分支，无中间祖先\n")
                            else:
                                f.write(f"{node} (Level {level})\n")
                            
                            # 序列（前60个字符）
                            if node in ancestral_seqs:
                                seq = ancestral_seqs[node][:60]
                                if len(ancestral_seqs[node]) > 60:
                                    seq += '...'
                                f.write(f"    序列: {seq}\n")
                            elif node in modern_seqs:
                                seq = modern_seqs[node][:60]
                                if len(modern_seqs[node]) > 60:
                                    seq += '...'
                                f.write(f"    序列: {seq}\n")
                            
                            prev_level = level
                        
                        # 突变信息
                        if row['Position'] != '' and pd.notna(row['Position']):
                            f.write(f"    ├─ 位点 {int(row['Position'])}: {row['From_AA']} → {row['To_AA']}\n")
                        elif row['Num_Mutations'] == 0 and row['Note'] == 'Conserved':
                            f.write(f"    ├─ 无突变（保守）\n")
                    
                    f.write("\n")
            
            self.logger.info(f"✅ 突变树已保存: {tree_viz_file}")
            
            # 统计信息
            total_mutations = df_trace[df_trace['Position'] != ''].shape[0]
            self.logger.info(f"📊 共记录 {len(all_records)} 条演化事件")
            self.logger.info(f"📊 共记录 {total_mutations} 个突变事件")
            
        except Exception as e:
            self.logger.error(f"❌ 追溯演化历史时出错: {e}")
    
    def _compare_sequences(self, seq1, seq2, node1, node2):
        """比较两个序列，找出突变位点"""
        mutations = []
        
        min_len = min(len(seq1), len(seq2))
        
        for i in range(min_len):
            if seq1[i] != seq2[i]:
                mutations.append({
                    'Position': i + 1,
                    'From': seq1[i],
                    'To': seq2[i]
                })
        
        return mutations

"""
基因组组装质量控制模块 | Genome Assembly Quality Control Module 🏆
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from .utils import CommandRunner, get_file_stats, format_size

class AssemblyQualityAssessor:
    """组装质量评估器 | Assembly Quality Assessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.qc_dir = Path(config.output_dir) / "quality_assessment"
        self.qc_dir.mkdir(parents=True, exist_ok=True)
    
    def assess_all_results(self, assembly_files: Dict[str, List[str]], scaffold_results: Dict[str, Dict[str, str]] = None) -> Dict[str, Dict[str, any]]:
        """评估所有组装和挂载结果 | Assess all assembly and scaffolding results"""
        self.logger.info("🏆 开始组装质量评估 | Starting assembly quality assessment")
        
        quality_results = {}
        
        # 评估原始组装质量
        quality_results['assemblies'] = self._assess_assemblies(assembly_files)
        
        # 评估Hi-C挂载质量 (如果有)
        if scaffold_results:
            quality_results['scaffolds'] = self._assess_scaffolds(scaffold_results)
        
        # 运行BUSCO评估 (如果可用)
        quality_results['busco'] = self._run_busco_assessment(assembly_files, scaffold_results)
        
        # 生成综合质量报告
        self._generate_comprehensive_report(quality_results)
        
        # 质量等级评定
        quality_results['grades'] = self._grade_assembly_quality(quality_results)
        
        self.logger.info("✅ 组装质量评估完成 | Assembly quality assessment completed")
        return quality_results
    
    def _assess_assemblies(self, assembly_files: Dict[str, List[str]]) -> Dict[str, Dict[str, any]]:
        """评估原始组装质量 | Assess original assembly quality"""
        self.logger.info("📊 评估原始组装质量 | Assessing original assembly quality")
        
        assembly_results = {}
        
        for assembly_type, fasta_list in assembly_files.items():
            assembly_results[assembly_type] = {}
            
            for fasta_file in fasta_list:
                assembly_name = Path(fasta_file).stem
                self.logger.info(f"🔍 分析{assembly_type}组装: {assembly_name}")
                
                # 基本统计信息
                basic_stats = get_file_stats(fasta_file, self.logger)
                
                # 详细连续性分析
                continuity_stats = self._calculate_continuity_metrics(fasta_file)
                
                # 组装完整性检查
                completeness_stats = self._check_assembly_completeness(fasta_file, assembly_name)
                
                assembly_results[assembly_type][assembly_name] = {
                    'file_path': fasta_file,
                    'basic_stats': basic_stats,
                    'continuity': continuity_stats,
                    'completeness': completeness_stats
                }
        
        return assembly_results
    
    def _assess_scaffolds(self, scaffold_results: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, any]]:
        """评估Hi-C挂载质量 | Assess Hi-C scaffolding quality"""
        self.logger.info("🔗 评估Hi-C挂载质量 | Assessing Hi-C scaffolding quality")
        
        scaffold_assessment = {}
        
        for assembly_type, results in scaffold_results.items():
            scaffold_assessment[assembly_type] = {}
            
            for assembly_name, files in results.items():
                self.logger.info(f"🎯 评估{assembly_type}挂载: {assembly_name}")
                
                # 查找最终scaffold文件
                scaffold_file = self._find_final_scaffold_file(files)
                if not scaffold_file:
                    self.logger.warning(f"⚠️ 未找到{assembly_name}的最终scaffold文件")
                    continue
                
                # 基本统计
                basic_stats = get_file_stats(scaffold_file, self.logger)
                
                # 挂载质量分析
                scaffolding_metrics = self._calculate_scaffolding_metrics(scaffold_file)
                
                # Hi-C一致性评估 (如果有.hic文件)
                hic_consistency = self._assess_hic_consistency(files, assembly_name)
                
                scaffold_assessment[assembly_type][assembly_name] = {
                    'scaffold_file': scaffold_file,
                    'basic_stats': basic_stats,
                    'scaffolding_metrics': scaffolding_metrics,
                    'hic_consistency': hic_consistency,
                    'output_files': files
                }
        
        return scaffold_assessment
    
    def _calculate_continuity_metrics(self, fasta_file: str) -> Dict[str, any]:
        """计算连续性指标 | Calculate continuity metrics"""
        try:
            # 获取所有序列长度
            cmd = f"seqkit fx2tab -l {fasta_file} | cut -f2 | sort -nr"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            lengths = [int(x) for x in result.stdout.strip().split('\n') if x.strip()]
            
            if not lengths:
                return {}
            
            total_length = sum(lengths)
            num_sequences = len(lengths)
            
            # 计算Nx指标
            metrics = {
                'total_length': total_length,
                'num_sequences': num_sequences,
                'longest_sequence': lengths[0] if lengths else 0,
                'shortest_sequence': lengths[-1] if lengths else 0,
                'mean_length': total_length // num_sequences if num_sequences > 0 else 0
            }
            
            # 计算N50, N90等
            for x in [10, 25, 50, 75, 90]:
                nx_value = self._calculate_nx(lengths, x)
                metrics[f'N{x}'] = nx_value
                
                # 计算Lx (达到Nx需要的序列数)
                cumulative = 0
                target = total_length * (x / 100)
                for i, length in enumerate(lengths):
                    cumulative += length
                    if cumulative >= target:
                        metrics[f'L{x}'] = i + 1
                        break
            
            return metrics
            
        except Exception as e:
            self.logger.warning(f"⚠️ 连续性计算失败: {e}")
            return {}
    
    def _calculate_nx(self, lengths: List[int], x: int) -> int:
        """计算Nx值 | Calculate Nx value"""
        total = sum(lengths)
        target = total * (x / 100)
        cumulative = 0
        
        for length in sorted(lengths, reverse=True):
            cumulative += length
            if cumulative >= target:
                return length
        return 0
    
    def _check_assembly_completeness(self, fasta_file: str, assembly_name: str) -> Dict[str, any]:
        """检查组装完整性 | Check assembly completeness"""
        completeness = {}
        
        try:
            # 基本完整性指标
            stats = get_file_stats(fasta_file, self.logger)
            if stats and 'sum_len' in stats:
                assembly_size = int(stats['sum_len'].replace(',', ''))
                expected_size = self._parse_genome_size_bp(self.config.genome_size)
                
                if expected_size > 0:
                    size_ratio = assembly_size / expected_size
                    completeness['size_ratio'] = size_ratio
                    completeness['expected_size'] = expected_size
                    completeness['actual_size'] = assembly_size
                    
                    # 大小合理性评估
                    if 0.9 <= size_ratio <= 1.1:
                        completeness['size_assessment'] = 'good'
                    elif 0.8 <= size_ratio < 0.9 or 1.1 < size_ratio <= 1.3:
                        completeness['size_assessment'] = 'acceptable'
                    else:
                        completeness['size_assessment'] = 'concerning'
            
            # 检查N含量
            n_content = self._calculate_n_content(fasta_file)
            completeness['n_content'] = n_content
            
        except Exception as e:
            self.logger.warning(f"⚠️ 完整性检查失败: {e}")
        
        return completeness
    
    def _calculate_scaffolding_metrics(self, scaffold_file: str) -> Dict[str, any]:
        """计算挂载质量指标 | Calculate scaffolding quality metrics"""
        metrics = {}
        
        try:
            # 基本连续性指标
            basic_continuity = self._calculate_continuity_metrics(scaffold_file)
            metrics.update(basic_continuity)
            
            # 染色体级别分析
            scaffold_analysis = self._analyze_chromosome_scaffolds(scaffold_file)
            metrics['chromosome_analysis'] = scaffold_analysis
            
            # Gap分析
            gap_analysis = self._analyze_gaps(scaffold_file)
            metrics['gap_analysis'] = gap_analysis
            
        except Exception as e:
            self.logger.warning(f"⚠️ 挂载指标计算失败: {e}")
        
        return metrics
    
    def _analyze_chromosome_scaffolds(self, scaffold_file: str) -> Dict[str, any]:
        """分析染色体级scaffold | Analyze chromosome-level scaffolds"""
        analysis = {}
        
        try:
            # 获取序列长度信息
            cmd = f"seqkit fx2tab -l {scaffold_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            sequences = []
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        name, length = parts[0], int(parts[1])
                        sequences.append((name, length))
            
            # 按长度排序
            sequences.sort(key=lambda x: x[1], reverse=True)
            
            # 分析大scaffold (可能的染色体)
            large_scaffolds = [seq for seq in sequences if seq[1] >= 10000000]  # >10Mb
            
            analysis['total_scaffolds'] = len(sequences)
            analysis['large_scaffolds'] = len(large_scaffolds)
            analysis['large_scaffold_lengths'] = [seq[1] for seq in large_scaffolds]
            
            if large_scaffolds:
                analysis['largest_scaffold'] = large_scaffolds[0][1]
                analysis['large_scaffolds_total_length'] = sum(seq[1] for seq in large_scaffolds)
                
                # 估算染色体数量 (基于大scaffold数量)
                analysis['estimated_chromosomes'] = len(large_scaffolds)
            
        except Exception as e:
            self.logger.warning(f"⚠️ 染色体分析失败: {e}")
        
        return analysis
    
    def _analyze_gaps(self, scaffold_file: str) -> Dict[str, any]:
        """分析Gap情况 | Analyze gaps"""
        gap_info = {}
        
        try:
            # 计算N的数量和分布
            cmd = f"seqkit fx2tab {scaffold_file} | cut -f2 | tr -cd 'N' | wc -c"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            total_ns = int(result.stdout.strip()) if result.stdout.strip() else 0
            
            gap_info['total_ns'] = total_ns
            
            # 计算N区域的数量
            cmd = f"seqkit fx2tab {scaffold_file} | cut -f2 | grep -o 'N\\+' | wc -l"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            gap_regions = int(result.stdout.strip()) if result.stdout.strip() else 0
            
            gap_info['gap_regions'] = gap_regions
            gap_info['average_gap_size'] = total_ns // gap_regions if gap_regions > 0 else 0
            
        except Exception as e:
            self.logger.warning(f"⚠️ Gap分析失败: {e}")
        
        return gap_info
    
    def _assess_hic_consistency(self, files: Dict[str, str], assembly_name: str) -> Dict[str, any]:
        """评估Hi-C一致性 | Assess Hi-C consistency"""
        consistency = {}
        
        try:
            # 如果有.hic文件，分析Hi-C map质量
            hic_file = files.get('hic_file')
            if hic_file and os.path.exists(hic_file):
                # 这里可以使用juicer_tools来分析Hi-C质量
                # 暂时记录文件存在状态
                consistency['hic_file_available'] = True
                consistency['hic_file_path'] = hic_file
                
                # 可以添加更多Hi-C质量分析
                # 例如: 距离衰减曲线、分辨率分析等
                
            else:
                consistency['hic_file_available'] = False
                
        except Exception as e:
            self.logger.warning(f"⚠️ Hi-C一致性评估失败: {e}")
        
        return consistency
    
    def _run_busco_assessment(self, assembly_files: Dict[str, List[str]], scaffold_results: Dict[str, Dict[str, str]] = None) -> Dict[str, Dict[str, any]]:
        """运行BUSCO评估 | Run BUSCO assessment"""
        self.logger.info("🧬 运行BUSCO完整性评估 | Running BUSCO completeness assessment")
        
        busco_results = {}
        
        # 检查BUSCO是否可用
        try:
            result = subprocess.run(['busco', '--version'], capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                self.logger.warning("⚠️ BUSCO不可用，跳过BUSCO评估")
                return {}
        except:
            self.logger.warning("⚠️ BUSCO不可用，跳过BUSCO评估")
            return {}
        
        # 为主要组装文件运行BUSCO
        files_to_assess = []
        
        # 收集组装文件
        for assembly_type, fasta_list in assembly_files.items():
            if assembly_type == 'primary':  # 主要关注primary组装
                files_to_assess.extend([(f, f"{Path(f).stem}_assembly") for f in fasta_list])
        
        # 收集挂载文件
        if scaffold_results:
            for assembly_type, results in scaffold_results.items():
                if assembly_type == 'primary':  # 主要关注primary组装的挂载结果
                    for assembly_name, files in results.items():
                        scaffold_file = self._find_final_scaffold_file(files)
                        if scaffold_file:
                            files_to_assess.append((scaffold_file, f"{assembly_name}_scaffolded"))
        
        # 运行BUSCO评估
        for fasta_file, name in files_to_assess[:3]:  # 限制评估数量以节省时间
            busco_result = self._run_single_busco(fasta_file, name)
            if busco_result:
                busco_results[name] = busco_result
        
        return busco_results
    
    def _run_single_busco(self, fasta_file: str, name: str) -> Optional[Dict[str, any]]:
        """运行单个BUSCO评估 | Run single BUSCO assessment"""
        self.logger.info(f"🔬 BUSCO评估: {name}")
        
        busco_dir = self.qc_dir / "busco" / name
        busco_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建BUSCO命令
        cmd_parts = [
            'busco',
            f'-i {fasta_file}',
            f'-o {name}',
            '-m genome',
            f'--out_path {busco_dir.parent}',
            f'-c {min(self.config.threads, 16)}'  # 限制BUSCO线程数
        ]
        
        # 添加lineage参数
        if self.config.busco_lineage and self.config.busco_lineage != "auto":
            cmd_parts.append(f'-l {self.config.busco_lineage}')
        
        cmd = " ".join(cmd_parts)
        
        # 运行BUSCO (设置较长超时时间)
        success = self.cmd_runner.run(cmd, f"BUSCO评估 - {name}", timeout=7200)
        
        if not success:
            self.logger.warning(f"⚠️ BUSCO评估失败: {name}")
            return None
        
        # 解析BUSCO结果
        return self._parse_busco_results(busco_dir / name)
    
    def _parse_busco_results(self, busco_output_dir: Path) -> Optional[Dict[str, any]]:
        """解析BUSCO结果 | Parse BUSCO results"""
        try:
            # 查找BUSCO短摘要文件
            summary_files = list(busco_output_dir.glob("short_summary*.txt"))
            if not summary_files:
                return None
            
            summary_file = summary_files[0]
            
            with open(summary_file, 'r') as f:
                content = f.read()
            
            # 解析BUSCO统计
            busco_stats = {}
            for line in content.split('\n'):
                line = line.strip()
                if 'Complete BUSCOs' in line:
                    busco_stats['complete'] = self._extract_busco_number(line)
                elif 'Complete and single-copy BUSCOs' in line:
                    busco_stats['single_copy'] = self._extract_busco_number(line)
                elif 'Complete and duplicated BUSCOs' in line:
                    busco_stats['duplicated'] = self._extract_busco_number(line)
                elif 'Fragmented BUSCOs' in line:
                    busco_stats['fragmented'] = self._extract_busco_number(line)
                elif 'Missing BUSCOs' in line:
                    busco_stats['missing'] = self._extract_busco_number(line)
                elif 'Total BUSCO groups searched' in line:
                    busco_stats['total'] = self._extract_busco_number(line)
            
            # 计算百分比
            if 'total' in busco_stats and busco_stats['total'] > 0:
                total = busco_stats['total']
                for key in ['complete', 'single_copy', 'duplicated', 'fragmented', 'missing']:
                    if key in busco_stats:
                        busco_stats[f'{key}_percent'] = (busco_stats[key] / total) * 100
            
            return busco_stats
            
        except Exception as e:
            self.logger.warning(f"⚠️ BUSCO结果解析失败: {e}")
            return None
    
    def _extract_busco_number(self, line: str) -> int:
        """提取BUSCO数字 | Extract BUSCO number"""
        # 例如: "C:98.9%[S:98.1%,D:0.8%],F:0.6%,M:0.5%,n:2586"
        # 或者: "1234	Complete BUSCOs (C)"
        import re
        numbers = re.findall(r'\d+', line)
        return int(numbers[0]) if numbers else 0
    
    def _find_final_scaffold_file(self, files: Dict[str, str]) -> Optional[str]:
        """查找最终scaffold文件 | Find final scaffold file"""
        # 优先顺序: scaffolds_final -> final_fasta -> scaffolds_FINAL
        priority_keys = ['scaffolds_final', 'final_fasta', 'scaffolds_FINAL']
        
        for key in priority_keys:
            if key in files and os.path.exists(files[key]):
                return files[key]
        
        # 如果没找到，返回任意存在的文件
        for file_path in files.values():
            if os.path.exists(file_path) and file_path.endswith('.fa'):
                return file_path
        
        return None
    
    def _calculate_n_content(self, fasta_file: str) -> float:
        """计算N含量百分比 | Calculate N content percentage"""
        try:
            # 获取总长度和N的数量
            cmd = f"seqkit stats {fasta_file} -T | tail -n1 | cut -f5"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            total_length = int(result.stdout.strip().replace(',', ''))
            
            cmd = f"seqkit fx2tab {fasta_file} | cut -f2 | tr -cd 'N' | wc -c"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            n_count = int(result.stdout.strip())
            
            return (n_count / total_length) * 100 if total_length > 0 else 0
            
        except Exception as e:
            self.logger.warning(f"⚠️ N含量计算失败: {e}")
            return 0.0
    
    def _parse_genome_size_bp(self, genome_size: str) -> int:
        """解析基因组大小为bp | Parse genome size to bp"""
        try:
            size = genome_size.lower()
            if size.endswith('k'):
                return int(float(size[:-1]) * 1000)
            elif size.endswith('m'):
                return int(float(size[:-1]) * 1000000)
            elif size.endswith('g'):
                return int(float(size[:-1]) * 1000000000)
            else:
                return int(size)
        except:
            return 0
    
    def _generate_comprehensive_report(self, quality_results: Dict[str, Dict[str, any]]):
        """生成综合质量报告 | Generate comprehensive quality report"""
        report_file = self.qc_dir / "comprehensive_quality_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🏆 基因组组装综合质量报告 | Comprehensive Genome Assembly Quality Report\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"项目名称 | Project: {self.config.project_name}\n")
            f.write(f"组装策略 | Assembly Strategy: {self.config.assembly_strategy}\n")
            f.write(f"数据类型 | Data Types: {', '.join(self.config.detected_data_types)}\n")
            f.write(f"Hi-C策略 | Hi-C Strategy: {self.config.hic_strategy}\n\n")
            
            # 组装质量报告
            if 'assemblies' in quality_results:
                f.write("📊 组装质量评估 | Assembly Quality Assessment\n")
                f.write("-" * 60 + "\n")
                
                for assembly_type, results in quality_results['assemblies'].items():
                    f.write(f"\n{assembly_type.upper()} 组装:\n")
                    
                    for name, data in results.items():
                        f.write(f"  组装名称 | Assembly: {name}\n")
                        
                        if 'basic_stats' in data:
                            stats = data['basic_stats']
                            f.write(f"    序列数 | Sequences: {stats.get('num', 'N/A')}\n")
                            f.write(f"    总长度 | Total Length: {stats.get('sum_len', 'N/A')}\n")
                            f.write(f"    最大长度 | Max Length: {stats.get('max_len', 'N/A')}\n")
                        
                        if 'continuity' in data and data['continuity']:
                            cont = data['continuity']
                            f.write(f"    N50: {cont.get('N50', 'N/A')}\n")
                            f.write(f"    N90: {cont.get('N90', 'N/A')}\n")
                            f.write(f"    L50: {cont.get('L50', 'N/A')}\n")
                        
                        if 'completeness' in data and data['completeness']:
                            comp = data['completeness']
                            if 'size_ratio' in comp:
                                f.write(f"    大小比例 | Size Ratio: {comp['size_ratio']:.2f}\n")
                            if 'n_content' in comp:
                                f.write(f"    N含量 | N Content: {comp['n_content']:.2f}%\n")
                        
                        f.write("\n")
            
            # 挂载质量报告
            if 'scaffolds' in quality_results:
                f.write("🔗 挂载质量评估 | Scaffolding Quality Assessment\n")
                f.write("-" * 60 + "\n")
                
                for assembly_type, results in quality_results['scaffolds'].items():
                    f.write(f"\n{assembly_type.upper()} 挂载:\n")
                    
                    for name, data in results.items():
                        f.write(f"  挂载名称 | Scaffold: {name}\n")
                        
                        if 'basic_stats' in data:
                            stats = data['basic_stats']
                            f.write(f"    序列数 | Sequences: {stats.get('num', 'N/A')}\n")
                            f.write(f"    总长度 | Total Length: {stats.get('sum_len', 'N/A')}\n")
                        
                        if 'scaffolding_metrics' in data:
                            metrics = data['scaffolding_metrics']
                            if 'N50' in metrics:
                                f.write(f"    Scaffold N50: {metrics['N50']}\n")
                            
                            if 'chromosome_analysis' in metrics:
                                chrom = metrics['chromosome_analysis']
                                f.write(f"    大scaffold数 | Large Scaffolds: {chrom.get('large_scaffolds', 'N/A')}\n")
                                f.write(f"    预估染色体数 | Est. Chromosomes: {chrom.get('estimated_chromosomes', 'N/A')}\n")
                        
                        f.write("\n")
            
            # BUSCO评估报告
            if 'busco' in quality_results and quality_results['busco']:
                f.write("🧬 BUSCO完整性评估 | BUSCO Completeness Assessment\n")
                f.write("-" * 60 + "\n")
                
                for name, busco_data in quality_results['busco'].items():
                    f.write(f"\n  {name}:\n")
                    if 'complete_percent' in busco_data:
                        f.write(f"    完整BUSCO | Complete: {busco_data['complete_percent']:.1f}%\n")
                    if 'single_copy_percent' in busco_data:
                        f.write(f"    单拷贝 | Single-copy: {busco_data['single_copy_percent']:.1f}%\n")
                    if 'duplicated_percent' in busco_data:
                        f.write(f"    重复 | Duplicated: {busco_data['duplicated_percent']:.1f}%\n")
                    if 'missing_percent' in busco_data:
                        f.write(f"    缺失 | Missing: {busco_data['missing_percent']:.1f}%\n")
        
        self.logger.info(f"📋 综合质量报告已生成: {report_file}")
    
    def _grade_assembly_quality(self, quality_results: Dict[str, Dict[str, any]]) -> Dict[str, str]:
        """评定组装质量等级 | Grade assembly quality"""
        grades = {}
        
        try:
            # 基于关键指标进行评级
            for assembly_type, results in quality_results.get('assemblies', {}).items():
                for name, data in results.items():
                    grade = self._calculate_single_grade(data, quality_results.get('busco', {}), name)
                    grades[f"{assembly_type}_{name}"] = grade
            
            # 挂载结果评级
            for assembly_type, results in quality_results.get('scaffolds', {}).items():
                for name, data in results.items():
                    grade = self._calculate_scaffold_grade(data, quality_results.get('busco', {}), name)
                    grades[f"scaffold_{assembly_type}_{name}"] = grade
        
        except Exception as e:
            self.logger.warning(f"⚠️ 质量评级失败: {e}")
        
        return grades
    
    def _calculate_single_grade(self, assembly_data: Dict[str, any], busco_data: Dict[str, any], name: str) -> str:
        """计算单个组装的质量等级 | Calculate quality grade for single assembly"""
        score = 0
        max_score = 0
        
        # N50评分 (30分)
        max_score += 30
        continuity = assembly_data.get('continuity', {})
        n50 = continuity.get('N50', 0)
        if n50 >= 10000000:  # >10Mb
            score += 30
        elif n50 >= 1000000:  # >1Mb
            score += 20
        elif n50 >= 100000:  # >100kb
            score += 10
        
        # BUSCO评分 (40分)
        max_score += 40
        busco_key = f"{name}_assembly"
        if busco_key in busco_data:
            complete_percent = busco_data[busco_key].get('complete_percent', 0)
            if complete_percent >= 95:
                score += 40
            elif complete_percent >= 90:
                score += 30
            elif complete_percent >= 80:
                score += 20
            elif complete_percent >= 70:
                score += 10
        
        # 大小合理性评分 (20分)
        max_score += 20
        completeness = assembly_data.get('completeness', {})
        if 'size_assessment' in completeness:
            if completeness['size_assessment'] == 'good':
                score += 20
            elif completeness['size_assessment'] == 'acceptable':
                score += 15
            elif completeness['size_assessment'] == 'concerning':
                score += 5
        
        # N含量评分 (10分)
        max_score += 10
        n_content = completeness.get('n_content', 100)  # 默认高N含量
        if n_content <= 1:
            score += 10
        elif n_content <= 5:
            score += 7
        elif n_content <= 10:
            score += 5
        
        # 计算最终等级
        if max_score > 0:
            percentage = (score / max_score) * 100
            if percentage >= 90:
                return "A (优秀)"
            elif percentage >= 80:
                return "B (良好)"
            elif percentage >= 70:
                return "C (中等)"
            elif percentage >= 60:
                return "D (及格)"
            else:
                return "F (不及格)"
        
        return "未评级"
    
    def _calculate_scaffold_grade(self, scaffold_data: Dict[str, any], busco_data: Dict[str, any], name: str) -> str:
        """计算挂载结果的质量等级 | Calculate quality grade for scaffolds"""
        # 挂载结果通常质量更高，评分标准更严格
        score = 0
        max_score = 0
        
        # Scaffold N50评分 (40分)
        max_score += 40
        metrics = scaffold_data.get('scaffolding_metrics', {})
        n50 = metrics.get('N50', 0)
        if n50 >= 50000000:  # >50Mb (染色体级别)
            score += 40
        elif n50 >= 10000000:  # >10Mb
            score += 30
        elif n50 >= 1000000:  # >1Mb
            score += 15
        
        # 染色体数量评分 (30分)
        max_score += 30
        chrom_analysis = metrics.get('chromosome_analysis', {})
        large_scaffolds = chrom_analysis.get('large_scaffolds', 0)
        if large_scaffolds > 0:
            if large_scaffolds <= 50:  # 合理的染色体数量
                score += 30
            elif large_scaffolds <= 100:
                score += 20
            else:
                score += 10
        
        # BUSCO评分 (20分)
        max_score += 20
        scaffold_busco = f"{name}_scaffolded"
        if scaffold_busco in busco_data:
            complete_percent = busco_data[scaffold_busco].get('complete_percent', 0)
            if complete_percent >= 95:
                score += 20
            elif complete_percent >= 90:
                score += 15
            elif complete_percent >= 85:
                score += 10
        
        # Gap情况评分 (10分)
        max_score += 10
        gap_analysis = metrics.get('gap_analysis', {})
        gap_regions = gap_analysis.get('gap_regions', 0)
        if gap_regions == 0:
            score += 10
        elif gap_regions <= 100:
            score += 7
        elif gap_regions <= 1000:
            score += 5
        
        # 计算最终等级
        if max_score > 0:
            percentage = (score / max_score) * 100
            if percentage >= 90:
                return "A+ (参考级别)"
            elif percentage >= 80:
                return "A (优秀)"
            elif percentage >= 70:
                return "B (良好)"
            elif percentage >= 60:
                return "C (中等)"
            else:
                return "D (需改进)"
        
        return "未评级"

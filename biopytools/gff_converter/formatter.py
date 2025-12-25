import re
from typing import List, Dict
from .gff_parser import GFFFeature
from .id_generator import IDGenerator, HierarchicalIDGenerator

class GFFFormatter:
    """GFF格式化器 | GFF Formatter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self._feature_counters = {}  # 添加这行，初始化计数器
    
    def write_output(self, header_lines: List[str], features: List[GFFFeature]):
        """写入输出文件 | Write output file"""
        self.logger.info(f"💾 开始写入输出文件 | Starting to write output file: {self.config.output_file}")
        
        with open(self.config.output_file, 'w') as f:
            # 写入头部 | Write header
            for header in header_lines:
                f.write(header + '\n')
            
            # 写入特征 | Write features
            for feature in features:
                f.write(feature.to_gff_line() + '\n')
        
        self.logger.info(f"✅ 输出文件写入完成 | Output file writing completed")
        self.logger.info(f"📁 输出文件路径 | Output file path: {self.config.output_file}")

    def format_features(self, features: List[GFFFeature], id_mapping: Dict[str, str]) -> List[GFFFeature]:
        """格式化特征列表 | Format feature list"""
        self.logger.info("🎨 开始格式化GFF特征 | Starting to format GFF features")
        
        # 1. 直接按染色体和位置排序所有特征 | Sort all features by chromosome and position
        self.logger.info("📏 按位置排序所有特征 | Sorting all features by position")
        
        def natural_sort_key(feature):
            chromosome = feature.seqid
            position = feature.start
            
            # 提取染色体编号进行自然排序 | Extract chromosome number for natural sorting
            import re
            chr_match = re.search(r'(\d+)', chromosome)
            chr_num = int(chr_match.group(1)) if chr_match else float('inf')
            
            return (chr_num, position)
        
        # 排序所有特征 | Sort all features
        sorted_features = sorted(features, key=natural_sort_key)
        
        # 2. 检测并标记需要合并的gene和非编码RNA | Detect and mark genes and ncRNAs to be merged
        self.logger.info("🔍 检测重复的基因和非编码RNA | Detecting duplicate genes and ncRNAs")
        gene_notes, features_to_skip = self._detect_gene_ncrna_duplicates(sorted_features)
        
        # 3. 建立ID映射表 | Build ID mapping table
        self.logger.info("🗂️ 建立ID映射表 | Building ID mapping table")
        complete_id_mapping = self._build_complete_mapping(sorted_features, id_mapping, features_to_skip)
        
        # 4. 格式化所有特征（跳过标记的特征）| Format all features (skip marked features)
        self.logger.info("✨ 格式化特征 | Formatting features")
        formatted_features = []
        
        for feature in sorted_features:
            feature_key = (feature.seqid, feature.start, feature.end, feature.type)
            
            # 跳过需要合并的非编码RNA特征 | Skip ncRNA features to be merged
            if feature_key in features_to_skip:
                continue
                
            formatted_feature = self._format_single_feature_simple(feature, complete_id_mapping, gene_notes)
            formatted_features.append(formatted_feature)
        
        self.logger.info(f"✅ 特征格式化完成 | Feature formatting completed: {len(formatted_features)} features")
        return formatted_features

    def _detect_gene_ncrna_duplicates(self, sorted_features: List[GFFFeature]) -> tuple:
        """检测需要合并的基因和非编码RNA | Detect genes and ncRNAs to be merged"""
        gene_notes = {}  # 基因ID -> Note内容
        features_to_skip = set()  # 需要跳过的特征
        
        # 按坐标分组特征 | Group features by coordinates
        coord_groups = {}
        for feature in sorted_features:
            coord_key = (feature.seqid, feature.start, feature.end)
            if coord_key not in coord_groups:
                coord_groups[coord_key] = []
            coord_groups[coord_key].append(feature)
        
        # 检查每个坐标组 | Check each coordinate group
        for coord_key, group in coord_groups.items():
            gene_features = [f for f in group if f.type == 'gene']
            ncrna_features = [f for f in group if f.type in ['lncRNA', 'tRNA', 'rRNA', 'miRNA', 'snRNA', 'snoRNA']]
            
            # 如果有基因和非编码RNA在相同位置 | If gene and ncRNA at same position
            if len(gene_features) == 1 and len(ncrna_features) == 1:
                gene = gene_features[0]
                ncrna = ncrna_features[0]
                
                # 检查是否真的是父子关系 | Check if they are really parent-child
                gene_id = gene.get_attribute('ID', '').replace('gene-', '')
                ncrna_parent = ncrna.get_attribute('Parent', '').replace('gene-', '')
                
                if gene_id == ncrna_parent:
                    # 标记基因需要添加Note | Mark gene to add Note
                    gene_notes[gene_id] = ncrna.type
                    
                    # 标记非编码RNA需要跳过 | Mark ncRNA to skip
                    feature_key = (ncrna.seqid, ncrna.start, ncrna.end, ncrna.type)
                    features_to_skip.add(feature_key)
                    
                    self.logger.info(f"🔗 合并 {ncrna.type} 到基因: {gene_id}")
        
        return gene_notes, features_to_skip

    def _build_complete_mapping(self, sorted_features: List[GFFFeature], id_mapping: Dict[str, str], features_to_skip: set) -> Dict[str, str]:
        """建立完整的ID映射表 | Build complete ID mapping table"""
        complete_mapping = {}
        complete_mapping.update(id_mapping)  # 基因ID映射
        
        # 收集转录本并去重编号 | Collect transcripts and deduplicate numbering
        transcript_counters = {}  # 每个基因的转录本计数器
        
        for feature in sorted_features:
            # 跳过需要删除的特征 | Skip features to be deleted
            feature_key = (feature.seqid, feature.start, feature.end, feature.type)
            if feature_key in features_to_skip:
                continue
                
            if feature.type in ['mRNA', 'lncRNA', 'tRNA', 'rRNA', 'transcript']:
                parent_id = feature.get_attribute('Parent', '').replace('gene-', '')
                new_gene_id = id_mapping.get(parent_id, parent_id)
                
                # 转录本计数 | Transcript counting
                if new_gene_id not in transcript_counters:
                    transcript_counters[new_gene_id] = {}
                
                transcript_type = feature.type
                if transcript_type not in transcript_counters[new_gene_id]:
                    transcript_counters[new_gene_id][transcript_type] = 0
                
                transcript_counters[new_gene_id][transcript_type] += 1
                transcript_num = transcript_counters[new_gene_id][transcript_type]
                
                original_id = feature.get_attribute('ID', '').replace('rna-', '')
                new_transcript_id = f"{new_gene_id}.{transcript_type}{transcript_num}"
                complete_mapping[original_id] = new_transcript_id
        
        return complete_mapping

    def _format_single_feature_simple(self, feature: GFFFeature, complete_mapping: Dict[str, str], gene_notes: Dict[str, str]) -> GFFFeature:
        """简化的单特征格式化 | Simplified single feature formatting"""
        new_feature = GFFFeature(
            seqid=feature.seqid, source=feature.source, type=feature.type,
            start=feature.start, end=feature.end, score=feature.score,
            strand=feature.strand, phase=feature.phase, attributes={}
        )
        
        original_id = feature.get_attribute('ID', '').replace('gene-', '').replace('rna-', '').replace('exon-', '').replace('cds-', '')
        
        if feature.type == 'gene':
            # 基因特征 | Gene feature
            new_id = complete_mapping.get(original_id, original_id)
            new_feature.set_attribute('ID', new_id)
            new_feature.set_attribute('Name', new_id)
            
            # 如果需要添加Note | If need to add Note
            if original_id in gene_notes:
                new_feature.set_attribute('Note', gene_notes[original_id])
        
        elif feature.type in ['mRNA', 'lncRNA', 'tRNA', 'rRNA', 'transcript']:
            # 转录本特征 | Transcript features
            new_id = complete_mapping.get(original_id, original_id)
            parent_id = feature.get_attribute('Parent', '').replace('gene-', '')
            new_parent_id = complete_mapping.get(parent_id, parent_id)
            
            new_feature.set_attribute('ID', new_id)
            new_feature.set_attribute('Name', new_id)
            new_feature.set_attribute('Parent', new_parent_id)
        
        else:
            # 其他特征 | Other features
            parent_id = feature.get_attribute('Parent', '').replace('gene-', '').replace('rna-', '')
            new_parent_id = complete_mapping.get(parent_id, parent_id)
            
            # 生成特征编号（简化版：按出现顺序）| Generate feature number (simplified: by appearance order)
            counter_key = f"{new_parent_id}_{feature.type}"
            if counter_key not in self._feature_counters:
                self._feature_counters[counter_key] = 0
            self._feature_counters[counter_key] += 1
            
            feature_num = self._feature_counters[counter_key]
            
            # 生成新ID | Generate new ID
            if feature.type == 'exon':
                new_id = f"{new_parent_id}.exon{feature_num}"
            elif feature.type == 'CDS':
                new_id = f"{new_parent_id}.cds{feature_num}"
            elif feature.type == 'five_prime_UTR':
                new_id = f"{new_parent_id}.5utrp{feature_num}"
            elif feature.type == 'three_prime_UTR':
                new_id = f"{new_parent_id}.3utrp{feature_num}"
            else:
                new_id = f"{new_parent_id}.{feature.type.lower()}{feature_num}"
            
            new_feature.set_attribute('ID', new_id)
            new_feature.set_attribute('Parent', new_parent_id)
        
        return new_feature
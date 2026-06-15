"""
InterProScan结果解析模块|InterProScan Result Parser Module
解析XML/TSV/JSON格式的InterProScan输出，提取GO和Pathway信息
Parse XML/TSV/JSON format InterProScan output, extract GO and Pathway information
"""

import xml.etree.ElementTree as ET
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict
from .go_database import GODatabase


@dataclass
class ProteinMatch:
    """蛋白质匹配结果|Protein match result"""
    protein_id: str
    md5: str
    length: int
    database: str
    signature_id: str
    signature_description: str
    start: int
    end: int
    score: float
    interpro_id: str = ""
    interpro_description: str = ""
    go_terms: List[str] = field(default_factory=list)
    pathways: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """转换为字典|Convert to dictionary"""
        return {
            'Protein ID': self.protein_id,
            'MD5': self.md5,
            'Length': self.length,
            'Database': self.database,
            'Signature ID': self.signature_id,
            'Description': self.signature_description,
            'Start': self.start,
            'End': self.end,
            'Score': self.score,
            'InterPro ID': self.interpro_id,
            'InterPro Description': self.interpro_description,
            'GO Terms': ';'.join(self.go_terms) if self.go_terms else '',
            'Pathways': ';'.join(self.pathways) if self.pathways else ''
        }


@dataclass
class ProteinSummary:
    """蛋白质注释汇总|Protein annotation summary"""
    protein_id: str
    length: int
    total_matches: int
    databases: Set[str] = field(default_factory=set)
    interpro_ids: Set[str] = field(default_factory=set)
    go_terms: Set[str] = field(default_factory=set)
    pathways: Set[str] = field(default_factory=set)

    def to_dict(self) -> dict:
        """转换为字典|Convert to dictionary"""
        return {
            'Protein ID': self.protein_id,
            'Length': self.length,
            'Total Matches': self.total_matches,
            'Database Count': len(self.databases),
            'Databases': ';'.join(sorted(self.databases)),
            'InterPro Count': len(self.interpro_ids),
            'GO Terms Count': len(self.go_terms),
            'GO Terms': ';'.join(sorted(self.go_terms)) if self.go_terms else '',
            'Pathways Count': len(self.pathways),
            'Pathways': ';'.join(sorted(self.pathways)) if self.pathways else ''
        }


class InterProScanParser:
    """InterProScan结果解析器|InterProScan Result Parser"""

    def __init__(self, logger: logging.Logger, parse_goterms: bool = True, parse_pathways: bool = True,
                 go_database_path: Optional[str] = None):
        self.logger = logger
        self.parse_goterms = parse_goterms
        self.parse_pathways = parse_pathways
        self.matches: List[ProteinMatch] = []
        self.summaries: Dict[str, ProteinSummary] = {}

        # 初始化GO数据库|Initialize GO database
        self.go_database = None
        try:
            # 默认使用内置数据库，如果指定了外部JSON则使用外部数据库
            # Use built-in database by default, or use external database if specified
            self.go_database = GODatabase(
                go_json_path=go_database_path,
                logger=logger,
                use_builtin=True  # 默认使用内置数据库|Default use built-in database
            )
            self.logger.info(f"GO数据库已加载 ({self.go_database.get_source()}|GO database loaded ({self.go_database.get_source()}): {self.go_database.get_stats()}")
        except Exception as e:
            self.logger.warning(f"GO数据库加载失败，继续但不填充GO名称|GO database load failed, continuing without GO names: {str(e)}")

    def parse_xml(self, xml_file: str) -> List[ProteinMatch]:
        """解析XML格式输出|Parse XML format output"""
        self.logger.info(f"解析XML文件|Parsing XML file: {xml_file}")

        if not Path(xml_file).exists():
            self.logger.error(f"XML文件不存在|XML file does not exist: {xml_file}")
            return []

        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            matches = []
            protein_count = 0

            # 定义匹配类型标签|Define match type tags
            match_types = ['hmmer3-match', 'coils-match', 'mobidblite-match', 'panther-match',
                          'superfamily-match', 'blastprodom-match', 'printsmatch',
                          'smart-match', 'tigrfam-match', 'gene3d-match', 'pfam-match']

            for protein in root.findall('protein'):
                # 获取蛋白质ID|Get protein ID
                xref = protein.find('xref')
                protein_id = xref.get('id') if xref is not None else ""
                name = xref.get('name') if xref is not None else ""

                # 获取序列信息|Get sequence info
                sequence = protein.find('sequence')
                md5 = sequence.get('md5') if sequence is not None else ""
                seq_text = sequence.text if sequence is not None else ""
                length = len(seq_text) if seq_text else 0

                protein_count += 1

                # 遍历所有类型的匹配|Iterate all match types
                for match_type in match_types:
                    for match in protein.findall(f'matches/{match_type}'):
                        # 获取signature信息|Get signature info
                        signature = match.find('signature')
                        if signature is None:
                            continue

                        signature_id = signature.get('ac', signature.get('name', ''))
                        signature_desc = signature.get('desc', '')
                        signature_name = signature.get('name', '')

                        # 获取数据库信息|Get database info
                        sig_lib = signature.find('signature-library-release')
                        database = sig_lib.get('library', '') if sig_lib is not None else match_type.replace('-match', '').upper()

                        # 获取InterPro信息|Get InterPro info
                        entry = signature.find('entry')
                        interpro_id = entry.get('ac', '') if entry is not None else ''
                        interpro_description = entry.get('desc', '') if entry is not None else ''

                        # 获取位置信息|Get location info
                        locations = match.find('.//locations')
                        if locations is None:
                            continue

                        for location in locations:
                            # 处理不同类型的location|Handle different location types
                            loc_tag = location.tag
                            start = int(location.get('start', 0))
                            end = int(location.get('end', 0))
                            score = float(location.get('score', 0)) if location.get('score') else 0.0
                            evalue = location.get('evalue', '')

                            # 解析GO Terms|Parse GO Terms
                            go_terms = []
                            if self.parse_goterms:
                                for go in signature.findall('.//go-xref'):
                                    go_id = go.get('id')
                                    go_term = go.get('name')
                                    go_category = go.get('category')
                                    if go_id:
                                        go_terms.append(f"{go_id}|{go_term}|{go_category}" if go_term else go_id)
                                # 也在entry下查找GO|Also find GO under entry
                                if entry is not None:
                                    for go in entry.findall('.//go-xref'):
                                        go_id = go.get('id')
                                        go_term = go.get('name')
                                        go_category = go.get('category')
                                        if go_id:
                                            go_terms.append(f"{go_id}|{go_term}|{go_category}" if go_term else go_id)

                            # 解析Pathways|Parse Pathways
                            pathways = []
                            if self.parse_pathways:
                                for pathway in signature.findall('.//pathway-xref'):
                                    pathway_db = pathway.get('db')
                                    pathway_id = pathway.get('id')
                                    pathway_name = pathway.get('name')
                                    if pathway_id:
                                        pathways.append(f"{pathway_db}:{pathway_id}|{pathway_name}" if pathway_name else f"{pathway_db}:{pathway_id}")
                                # 也在entry下查找pathway|Also find pathway under entry
                                if entry is not None:
                                    for pathway in entry.findall('.//pathway-xref'):
                                        pathway_db = pathway.get('db')
                                        pathway_id = pathway.get('id')
                                        pathway_name = pathway.get('name')
                                        if pathway_id:
                                            pathways.append(f"{pathway_db}:{pathway_id}|{pathway_name}" if pathway_name else f"{pathway_db}:{pathway_id}")

                            # 构建描述|Build description
                            description = signature_desc if signature_desc else signature_name

                            protein_match = ProteinMatch(
                                protein_id=protein_id,
                                md5=md5,
                                length=length,
                                database=database,
                                signature_id=signature_id,
                                signature_description=description,
                                start=start,
                                end=end,
                                score=score,
                                interpro_id=interpro_id,
                                interpro_description=interpro_description,
                                go_terms=go_terms,
                                pathways=pathways
                            )

                            matches.append(protein_match)

                # 更新汇总信息|Update summary
                protein_matches = [m for m in matches if m.protein_id == protein_id]
                self._update_summary(protein_id, length, protein_matches)

            self.logger.info(f"XML解析完成|XML parsing completed: {protein_count} proteins, {len(matches)} matches")
            self.matches = matches
            return matches

        except Exception as e:
            self.logger.error(f"XML解析失败|XML parsing failed: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return []

    def parse_tsv(self, tsv_file: str) -> List[ProteinMatch]:
        """解析TSV格式输出|Parse TSV format output"""
        self.logger.info(f"解析TSV文件|Parsing TSV file: {tsv_file}")

        if not Path(tsv_file).exists():
            self.logger.error(f"TSV文件不存在|TSV file does not exist: {tsv_file}")
            return []

        matches = []

        try:
            with open(tsv_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) < 14:
                        continue

                    # TSV格式：0:ProteinID 1:MD5 2:Length 3:DB 4:SignatureID 5:Desc 6:Start 7:End 8:Score 9:Status 10:Date 11:InterProID 12:InterProDesc 13:GO 14:Pathway
                    protein_id = parts[0]
                    md5 = parts[1]
                    length = int(parts[2])
                    database = parts[3]
                    signature_id = parts[4]
                    signature_description = parts[5]
                    start = int(parts[6])
                    end = int(parts[7])
                    score = self._parse_score(parts[8])
                    interpro_id = parts[11] if len(parts) > 11 else ""
                    interpro_description = parts[12] if len(parts) > 12 else ""

                    # 解析GO Terms (列13)|Parse GO Terms (column 13)
                    go_terms = []
                    if self.parse_goterms and len(parts) > 13 and parts[13] != '-':
                        go_list = parts[13].split('|')
                        for g in go_list:
                            g = g.strip()
                            if g and g != '-':
                                # 如果有GO数据库，添加名称和ontology信息|If GO database exists, add name and ontology info
                                if self.go_database:
                                    go_id = g.split('(')[0].strip()
                                    go_info = self.go_database.get_go_info(go_id)
                                    go_name = go_info.get('name', '')
                                    go_ontology = go_info.get('ontology', '')
                                    # 格式：GO:ID|Name|Ontology(Source) 或 GO:ID(Source)
                                    if go_name and go_ontology:
                                        source = g.split('(')[-1].rstrip(')') if '(' in g else ''
                                        go_terms.append(f"{go_id}|{go_name}|{go_ontology}({source})")
                                    else:
                                        go_terms.append(g)
                                else:
                                    go_terms.append(g)

                    # 解析Pathways (列14)|Parse Pathways (column 14)
                    pathways = []
                    if self.parse_pathways and len(parts) > 14 and parts[14] != '-':
                        pathway_list = parts[14].split('|')
                        pathways = [p.strip() for p in pathway_list if p.strip() and p.strip() != '-']

                    match = ProteinMatch(
                        protein_id=protein_id,
                        md5=md5,
                        length=length,
                        database=database,
                        signature_id=signature_id,
                        signature_description=signature_description,
                        start=start,
                        end=end,
                        score=score,
                        interpro_id=interpro_id,
                        interpro_description=interpro_description,
                        go_terms=go_terms,
                        pathways=pathways
                    )

                    matches.append(match)

            self.logger.info(f"TSV解析完成|TSV parsing completed: {len(matches)} matches")

            # 更新汇总|Update summaries
            self._update_summaries_from_matches(matches)

            self.matches = matches
            return matches

        except Exception as e:
            self.logger.error(f"TSV解析失败|TSV parsing failed: {str(e)}")
            return []

    def parse_json(self, json_file: str) -> List[ProteinMatch]:
        """解析JSON格式输出|Parse JSON format output"""
        self.logger.info(f"解析JSON文件|Parsing JSON file: {json_file}")

        if not Path(json_file).exists():
            self.logger.error(f"JSON文件不存在|JSON file does not exist: {json_file}")
            return []

        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            matches = []

            for protein in data.get('results', []):
                protein_id = protein.get('xrefs', [{}])[0].get('name', '')
                md5 = protein.get('md5', '')
                length = protein.get('sequence', '').get('length', 0)

                for match in protein.get('matches', []):
                    db = match.get('dbname', '')
                    signature_id = match.get('id', '')
                    signature_description = match.get('signature-library-release', {}).get('library', '')

                    location = match.get('locations', [{}])[0]
                    start = location.get('start', 0)
                    end = location.get('end', 0)
                    score = location.get('score', 0.0)

                    ipr = match.get('ipr', {})
                    interpro_id = ipr.get('ipr_id', '')
                    interpro_description = ipr.get('description', '')

                    # GO Terms
                    go_terms = []
                    if self.parse_goterms:
                        for go in match.get('go-terms', []):
                            go_id = go.get('identifier')
                            go_name = go.get('name')
                            go_category = go.get('category')
                            if go_id:
                                go_terms.append(f"{go_id}|{go_name}|{go_category}" if go_name else go_id)

                    # Pathways
                    pathways = []
                    if self.parse_pathways:
                        for pathway in match.get('pathways', []):
                            pathway_db = pathway.get('db')
                            pathway_id = pathway.get('id')
                            pathway_name = pathway.get('name')
                            if pathway_id:
                                pathways.append(f"{pathway_db}:{pathway_id}|{pathway_name}" if pathway_name else f"{pathway_db}:{pathway_id}")

                    protein_match = ProteinMatch(
                        protein_id=protein_id,
                        md5=md5,
                        length=length,
                        database=db,
                        signature_id=signature_id,
                        signature_description=signature_description,
                        start=start,
                        end=end,
                        score=score,
                        interpro_id=interpro_id,
                        interpro_description=interpro_description,
                        go_terms=go_terms,
                        pathways=pathways
                    )

                    matches.append(protein_match)

            self.logger.info(f"JSON解析完成|JSON parsing completed: {len(matches)} matches")
            self.matches = matches

            # 更新汇总|Update summaries
            self._update_summaries_from_matches(matches)

            return matches

        except Exception as e:
            self.logger.error(f"JSON解析失败|JSON parsing failed: {str(e)}")
            return []

    def _parse_score(self, score_str: str) -> float:
        """解析分数|Parse score"""
        try:
            return float(score_str)
        except (ValueError, TypeError):
            return 0.0

    def _update_summary(self, protein_id: str, length: int, matches: List[ProteinMatch]):
        """更新蛋白质汇总信息|Update protein summary"""
        if protein_id not in self.summaries:
            self.summaries[protein_id] = ProteinSummary(
                protein_id=protein_id,
                length=length,
                total_matches=0
            )

        summary = self.summaries[protein_id]
        summary.total_matches += len(matches)

        for match in matches:
            if match.database:
                summary.databases.add(match.database)
            if match.interpro_id and match.interpro_id != '-':
                summary.interpro_ids.add(match.interpro_id)
            summary.go_terms.update(match.go_terms)
            summary.pathways.update(match.pathways)

    def _update_summaries_from_matches(self, matches: List[ProteinMatch]):
        """从matches更新汇总信息|Update summaries from matches"""
        for match in matches:
            if match.protein_id not in self.summaries:
                self.summaries[match.protein_id] = ProteinSummary(
                    protein_id=match.protein_id,
                    length=match.length,
                    total_matches=0
                )

            summary = self.summaries[match.protein_id]
            summary.total_matches += 1

            if match.database:
                summary.databases.add(match.database)
            if match.interpro_id and match.interpro_id != '-':
                summary.interpro_ids.add(match.interpro_id)
            summary.go_terms.update(match.go_terms)
            summary.pathways.update(match.pathways)

    def get_matches_dataframe(self):
        """获取matches的pandas DataFrame|Get pandas DataFrame of matches"""
        import pandas as pd
        if not self.matches:
            return None
        return pd.DataFrame([m.to_dict() for m in self.matches])

    def get_summaries_dataframe(self):
        """获取summaries的pandas DataFrame|Get pandas DataFrame of summaries"""
        import pandas as pd
        if not self.summaries:
            return None
        return pd.DataFrame([s.to_dict() for s in self.summaries.values()])

    def get_database_statistics(self) -> Dict[str, int]:
        """获取数据库统计信息|Get database statistics"""
        stats = defaultdict(int)
        for match in self.matches:
            stats[match.database] += 1
        return dict(sorted(stats.items(), key=lambda x: x[1], reverse=True))

    def get_go_statistics(self) -> Dict[str, int]:
        """获取GO术语统计|Get GO term statistics"""
        if not self.parse_goterms:
            return {}

        stats = defaultdict(int)
        for match in self.matches:
            for go in match.go_terms:
                # 提取GO类别|Extract GO category
                if '|BP|' in go:
                    stats['Biological Process'] += 1
                elif '|CC|' in go:
                    stats['Cellular Component'] += 1
                elif '|MF|' in go:
                    stats['Molecular Function'] += 1
        return dict(stats)

    def get_expanded_go_dataframe(self):
        """获取展开的GO注释DataFrame（每个GO一行）|Get expanded GO annotation DataFrame (one GO per row)"""
        import pandas as pd

        if not self.parse_goterms:
            return None

        go_records = []
        for match in self.matches:
            if not match.go_terms:
                continue

            for go in match.go_terms:
                # 解析GO字符串格式:
                # 格式1: GO:ID|Name|Ontology(Source)
                # 格式2: GO:ID|Name|Category (旧格式)
                # 格式3: GO:ID(Source)
                parts = go.split('|')
                go_id_full = parts[0].strip() if parts else ''

                # 提取纯GO ID（移除来源后缀）|Extract pure GO ID (remove source suffix)
                go_id = go_id_full.split('(')[0].strip()

                go_name = ''
                go_category = ''

                # 检查是否已经有名称和ontology（来自GO数据库）|Check if name and ontology already present (from GO database)
                if len(parts) >= 3:
                    go_name = parts[1].strip()
                    go_category_raw = parts[2].strip()
                    # 移除来源后缀|Remove source suffix
                    go_category = go_category_raw.split('(')[0].strip()
                elif len(parts) == 2:
                    # 可能是旧格式 GO:ID|Category|Source 或其他
                    go_name = parts[1].strip()

                # 如果有GO数据库且名称为空，尝试从数据库获取|If GO database exists and name is empty, try to get from database
                if not go_name and self.go_database:
                    go_info = self.go_database.get_go_info(go_id)
                    go_name = go_info.get('name', '')
                    go_category = go_info.get('ontology', '')

                go_records.append({
                    'Protein ID': match.protein_id,
                    'Protein Length': match.length,
                    'Database': match.database,
                    'Domain ID': match.signature_id,
                    'Domain Description': match.signature_description,
                    'Domain Start': match.start,
                    'Domain End': match.end,
                    'InterPro ID': match.interpro_id,
                    'GO ID': go_id,
                    'GO Name': go_name,
                    'GO Category': go_category
                })

        if not go_records:
            return None

        return pd.DataFrame(go_records)

    def get_expanded_pathway_dataframe(self):
        """获取展开的Pathway注释DataFrame（每个Pathway一行）|Get expanded Pathway annotation DataFrame (one Pathway per row)"""
        import pandas as pd

        if not self.parse_pathways:
            return None

        pathway_records = []
        for match in self.matches:
            if not match.pathways:
                continue

            for pathway in match.pathways:
                # 解析Pathway字符串格式: DB:ID|Name 或 DB:ID
                parts = pathway.split('|')
                pathway_full_id = parts[0] if parts else ''
                pathway_name = parts[1] if len(parts) > 1 else ''

                # 分离数据库和ID|Split database and ID
                if ':' in pathway_full_id:
                    pathway_db, pathway_id = pathway_full_id.split(':', 1)
                else:
                    pathway_db = ''
                    pathway_id = pathway_full_id

                pathway_records.append({
                    'Protein ID': match.protein_id,
                    'Protein Length': match.length,
                    'Database': match.database,
                    'Domain ID': match.signature_id,
                    'Domain Description': match.signature_description,
                    'Domain Start': match.start,
                    'Domain End': match.end,
                    'InterPro ID': match.interpro_id,
                    'Pathway Database': pathway_db,
                    'Pathway ID': pathway_id,
                    'Pathway Name': pathway_name
                })

        if not pathway_records:
            return None

        return pd.DataFrame(pathway_records)

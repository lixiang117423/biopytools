"""
HMMsearch domtblout文件解析模块|HMMsearch domtblout File Parser Module
"""

import pandas as pd
from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class DomainHit:
    """Domain命中数据类|Domain Hit Data Class"""
    target_name: str
    target_accession: str
    target_length: int
    query_name: str
    query_accession: str
    query_length: int
    full_evalue: float
    full_score: float
    full_bias: float
    domain_number: int
    domain_total: int
    domain_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    accuracy: float
    description: str


class DomtbloutParser:
    """domtblout文件解析器|domtblout File Parser"""

    def __init__(self, logger):
        self.logger = logger

    def parse(self, domtblout_file: str) -> List[DomainHit]:
        """
        解析domtblout文件|Parse domtblout file

        Returns:
            List[DomainHit]: Domain命中列表|List of domain hits
        """
        self.logger.info(f"开始解析domtblout文件|Parsing domtblout file: {domtblout_file}")

        hits = []
        with open(domtblout_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # 跳过注释行和空行|Skip comment lines and empty lines
                if line.startswith('#') or not line:
                    continue

                try:
                    hit = self._parse_line(line)
                    if hit:
                        hits.append(hit)
                except Exception as e:
                    self.logger.warning(f"行{line_num}解析失败，跳过|Failed to parse line {line_num}, skipping: {e}")
                    continue

        self.logger.info(f"解析完成，共{len(hits)}条domain命中|Parsing completed, {len(hits)} domain hits found")
        return hits

    def _parse_line(self, line: str) -> Optional[DomainHit]:
        """
        解析单行数据|Parse single line

        domtblout格式(23列)|domtblout format (23 columns):
        1. target name
        2. target accession
        3. target length
        4. query name
        5. query accession
        6. query length
        7. full sequence E-value
        8. full sequence score
        9. full sequence bias
        10. domain number
        11. domain total (of this query)
        12. domain E-value (c-Evalue)
        13. domain i-Evalue
        14. domain score
        15. domain bias
        16. hmm from
        17. hmm to
        18. alignment from
        19. alignment to
        20. env from
        21. env to
        22. accuracy
        23. description
        """
        fields = line.split()

        # 至少需要23列|Need at least 23 columns
        if len(fields) < 23:
            return None

        return DomainHit(
            target_name=fields[0],
            target_accession=fields[1] if fields[1] != '-' else '',
            target_length=int(fields[2]),
            query_name=fields[3],
            query_accession=fields[4] if fields[4] != '-' else '',
            query_length=int(fields[5]),
            full_evalue=float(fields[6]),
            full_score=float(fields[7]),
            full_bias=float(fields[8]),
            domain_number=int(fields[9]),
            domain_total=int(fields[10]),
            domain_evalue=float(fields[11]),  # c-Evalue
            domain_score=float(fields[13]),  # 注意：跳过i-Evalue (fields[12]) | Skip i-Evalue
            domain_bias=float(fields[14]),   # 注意：偏移1列 | Offset by 1
            hmm_from=int(fields[15]),        # 偏移1列 | Offset by 1
            hmm_to=int(fields[16]),
            ali_from=int(fields[17]),
            ali_to=int(fields[18]),
            env_from=int(fields[19]),
            env_to=int(fields[20]),
            accuracy=float(fields[21]),      # 偏移1列 | Offset by 1
            description=' '.join(fields[22:]) if len(fields) > 22 else ''
        )

    def hits_to_dataframe(self, hits: List[DomainHit]) -> pd.DataFrame:
        """
        将Domain命中列表转换为DataFrame|Convert domain hits list to DataFrame

        Args:
            hits: Domain命中列表|List of domain hits

        Returns:
            pd.DataFrame: 结果数据框|Result dataframe
        """
        self.logger.info("转换为DataFrame格式|Converting to DataFrame format")

        data = []
        for hit in hits:
            data.append({
                'target_name': hit.target_name,
                'target_accession': hit.target_accession,
                'target_length': hit.target_length,
                'query_name': hit.query_name,
                'query_accession': hit.query_accession,
                'query_length': hit.query_length,
                'full_evalue': hit.full_evalue,
                'full_score': hit.full_score,
                'full_bias': hit.full_bias,
                'domain_number': hit.domain_number,
                'domain_total': hit.domain_total,
                'domain_evalue': hit.domain_evalue,
                'domain_score': hit.domain_score,
                'domain_bias': hit.domain_bias,
                'hmm_from': hit.hmm_from,
                'hmm_to': hit.hmm_to,
                'align_from': hit.ali_from,
                'align_to': hit.ali_to,
                'env_from': hit.env_from,
                'env_to': hit.env_to,
                'accuracy': hit.accuracy,
                'description': hit.description
            })

        df = pd.DataFrame(data)
        self.logger.info(f"DataFrame创建完成，共{len(df)}行{len(df.columns)}列|DataFrame created with {len(df)} rows and {len(df.columns)} columns")
        return df

    def filter_hits(self, hits: List[DomainHit],
                    evalue_threshold: float = None,
                    score_threshold: float = None) -> List[DomainHit]:
        """
        过滤Domain命中|Filter domain hits

        Args:
            hits: Domain命中列表|List of domain hits
            evalue_threshold: E-value阈值(保留小于该值的)|E-value threshold (keep values less than this)
            score_threshold: 分数阈值(保留大于该值的)|Score threshold (keep values greater than this)

        Returns:
            List[DomainHit]: 过滤后的命中列表|Filtered list of domain hits
        """
        filtered = hits

        if evalue_threshold is not None:
            before = len(filtered)
            filtered = [h for h in filtered if h.domain_evalue <= evalue_threshold]
            after = len(filtered)
            self.logger.info(f"E-value过滤|E-value filter (<= {evalue_threshold}): {before} -> {after}")

        if score_threshold is not None:
            before = len(filtered)
            filtered = [h for h in filtered if h.domain_score >= score_threshold]
            after = len(filtered)
            self.logger.info(f"分数过滤|Score filter (>= {score_threshold}): {before} -> {after}")

        return filtered

    def generate_summary(self, hits: List[DomainHit]) -> Dict:
        """
        生成统计摘要|Generate statistical summary

        Args:
            hits: Domain命中列表|List of domain hits

        Returns:
            Dict: 统计摘要|Statistical summary
        """
        self.logger.info("生成统计摘要|Generating statistical summary")

        # 统计每个目标基因的domain数量|Count domains per target gene
        target_domains = {}
        for hit in hits:
            target = hit.target_name
            if target not in target_domains:
                target_domains[target] = 0
            target_domains[target] += 1

        # 统计domain数量分布|Count domain number distribution
        domain_count_dist = {}
        for count in target_domains.values():
            if count not in domain_count_dist:
                domain_count_dist[count] = 0
            domain_count_dist[count] += 1

        summary = {
            'total_hits': len(hits),
            'unique_targets': len(target_domains),
            'domain_count_distribution': domain_count_dist,
            'targets_by_domain_count': target_domains
        }

        self.logger.info(f"统计摘要生成完成|Summary generated: {summary['total_hits']} hits, {summary['unique_targets']} unique targets")
        return summary

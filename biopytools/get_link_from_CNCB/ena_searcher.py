"""
ENA回退链接搜索器|ENA Fallback Link Searcher

当CNCB FTP上找不到Run ID时,通过ENA Portal API批量查询下载链接
|Batch-query ENA download links via the ENA Portal API when CNCB FTP misses a Run ID
"""

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict, List, Optional

from .utils import get_requests, create_retry_session

# 默认ENA Portal API地址|Default ENA Portal API endpoint
DEFAULT_ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# 单次请求的最大accession数量,超过则分片|Max accessions per request, split beyond this
DEFAULT_BATCH_SIZE = 100

# 逐个查询阶段的并发线程数|Concurrency for the individual query phase
DEFAULT_MAX_WORKERS = 8

# filereport API查询参数|filereport API query params
QUERY_FIELDS = 'run_accession,fastq_ftp,sra_ftp'


class ENALinkSearcher:
    """ENA下载链接搜索器|ENA Download Link Searcher

    两阶段查询|Two-phase query:
    1. 批量: 逗号分隔一次查多个ID(快速路径)|Batch: one request for many IDs (fast path)
    2. 逐个: ENA批量查询会被任一无效ID整体污染(实测只返回表头),批量未命中的ID并发逐个查询
       |Individual: ENA batch queries are poisoned wholesale by any invalid ID (returns header
       only in practice), so IDs missed by the batch are re-queried individually in parallel

    requests未安装时不抛异常,查询优雅降级为空结果
    |Without requests installed, queries degrade gracefully to empty results
    """

    def __init__(self, logger=None, api_url: str = DEFAULT_ENA_API,
                 timeout: int = 60, retry_attempts: int = 3,
                 session=None,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 max_workers: int = DEFAULT_MAX_WORKERS,
                 session_factory: Optional[Callable] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.api_url = api_url
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.batch_size = batch_size
        self.max_workers = max_workers
        self.requests = get_requests()
        if self.requests is not None:
            self.session = session or self._create_session(retry_attempts)
            # 逐个查询阶段每次新建会话,避免多线程共享Session
            # |Fresh session per individual query to avoid cross-thread sharing
            self.session_factory = session_factory or (
                lambda: self._create_session(self.retry_attempts)
            )
        else:
            self.session = None
            self.session_factory = None

    def _create_session(self, retry_attempts: int):
        """创建带重试的HTTP会话|Create HTTP session with retries"""
        if self.requests is None:
            raise ImportError("requests未安装|requests not installed")
        return create_retry_session(self.requests, retry_attempts)

    def search_links(self, run_ids: List[str]) -> Dict[str, List[str]]:
        """
        查询ENA下载链接|Query ENA download links

        Args:
            run_ids: Run ID列表|List of Run IDs

        Returns:
            Run ID到下载链接列表的映射|Mapping of Run ID to download link list
        """
        # 去重保序|Dedupe preserving order
        run_ids = list(dict.fromkeys(run_ids))
        if not run_ids:
            return {}

        # requests未安装时优雅降级|Degrade gracefully without requests
        if self.requests is None:
            self.logger.error(
                "未安装requests,ENA回退不可用|requests not installed, ENA fallback unavailable"
            )
            return {}

        results: Dict[str, List[str]] = {}
        missing: List[str] = []

        # 阶段1: 批量查询|Phase 1: batch query
        for i in range(0, len(run_ids), self.batch_size):
            batch = run_ids[i:i + self.batch_size]
            if len(batch) == 1:
                missing.append(batch[0])
                continue

            batch_links = self._query_batch(batch)
            for run_id in batch:
                if run_id in batch_links:
                    results[run_id] = batch_links[run_id]
                else:
                    missing.append(run_id)

        # 阶段2: 批量未命中的ID并发逐个查询|Phase 2: individual queries for batch misses
        if missing:
            self.logger.info(
                f"ENA逐个查询|ENA individual queries: {len(missing)} 个IDs|IDs"
            )
            with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                for run_id, links in zip(missing, pool.map(self._query_single, missing)):
                    if links:
                        results[run_id] = links

        return results

    def _build_params(self, run_ids: List[str]) -> dict:
        """构造API查询参数|Build API query params"""
        return {
            'accession': ','.join(run_ids),
            'result': 'read_run',
            'format': 'tsv',
            'fields': QUERY_FIELDS,
        }

    def _query_batch(self, run_ids: List[str]) -> Dict[str, List[str]]:
        """查询单批accession|Query a single batch of accessions"""
        requests = self.requests
        if requests is None:
            return {}

        params = self._build_params(run_ids)

        try:
            self.logger.info(
                f"ENA批量查询|ENA batch query: {len(run_ids)} 个IDs|IDs -> {self.api_url}"
            )
            response = self.session.get(self.api_url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return self._parse_filereport(response.text, run_ids)
        except requests.exceptions.RequestException as e:
            self.logger.warning(f"ENA批量查询失败|ENA batch query failed: {e}")
            return {}

    def _query_single(self, run_id: str) -> List[str]:
        """查询单个accession|Query a single accession"""
        requests = self.requests
        if requests is None:
            return []

        params = self._build_params([run_id])

        try:
            self.logger.debug(f"ENA单个查询|ENA single query: {run_id}")
            session = self.session_factory()
            response = session.get(self.api_url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return self._parse_filereport(response.text, [run_id]).get(run_id, [])
        except requests.exceptions.RequestException as e:
            self.logger.warning(f"ENA单个查询失败|ENA single query failed: {run_id}: {e}")
            return []

    def _parse_filereport(self, content: str, requested_ids: List[str]) -> Dict[str, List[str]]:
        """
        解析filereport TSV内容|Parse filereport TSV content

        每个Run优先取fastq_ftp链接,没有时回退到sra_ftp
        |Prefer fastq_ftp links per run, fall back to sra_ftp

        Args:
            content: API返回的TSV文本|TSV text returned by the API
            requested_ids: 请求的Run ID列表|Requested Run IDs

        Returns:
            Run ID到下载链接列表的映射|Mapping of Run ID to download link list
        """
        lines = [line for line in content.splitlines() if line.strip()]
        if not lines:
            return {}

        header = lines[0].split('\t')
        if 'run_accession' not in header:
            return {}

        idx_run = header.index('run_accession')
        idx_fastq = header.index('fastq_ftp') if 'fastq_ftp' in header else None
        idx_sra = header.index('sra_ftp') if 'sra_ftp' in header else None

        # 只保留请求范围内的Run,防止API返回study级结果时混入无关行
        # |Keep only requested runs; guards against study-level rows leaking in
        requested_set = set(requested_ids)

        links_by_run: Dict[str, List[str]] = {}
        for line in lines[1:]:
            cols = line.split('\t')
            if idx_run >= len(cols) or not cols[idx_run].strip():
                continue

            run_id = cols[idx_run].strip()
            if run_id not in requested_set:
                continue
            urls = self._extract_urls(cols, idx_fastq)
            if not urls:
                urls = self._extract_urls(cols, idx_sra)

            if urls:
                links_by_run[run_id] = urls

        return links_by_run

    @staticmethod
    def _extract_urls(cols: List[str], idx: Optional[int]) -> List[str]:
        """提取某列中的URL列表(分号分隔)|Extract semicolon-separated URLs from a column"""
        if idx is None or idx >= len(cols):
            return []

        urls = []
        for url in cols[idx].split(';'):
            url = url.strip()
            if not url:
                continue
            # ENA API返回的路径不带scheme,wget需要完整URL,补ftp://前缀
            # |ENA API omits the scheme; wget needs a full URL, prefix ftp://
            if not url.startswith(('ftp://', 'http://', 'https://')):
                url = f"ftp://{url}"
            urls.append(url)
        return urls

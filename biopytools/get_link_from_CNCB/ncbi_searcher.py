"""
NCBI SDL数据定位回退搜索器|NCBI SRA Data Locator Fallback Searcher

CNCB FTP镜像与ENA都没有链接的INSDC Run(如仅存NCBI S3的PacBio source BAM,
ENA filereport有元数据但fastq_ftp/sra_ftp为空),通过NCBI SDL API定位真实文件位置
|For INSDC runs missed by both the CNCB FTP mirror and ENA (e.g. NCBI-only PacBio
source BAMs, where the ENA filereport lists metadata with empty link columns),
locate the real files via the NCBI SDL API
"""

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Dict, List, Optional

from .utils import get_requests

# 默认SDL API地址|Default SDL API endpoint
DEFAULT_SDL_API = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

# 并发查询线程数|Concurrency for per-run queries
DEFAULT_MAX_WORKERS = 8


class NCBISDLSearcher:
    """NCBI SDL链接搜索器|NCBI SDL Link Searcher

    每个Run一次请求: GET {api}?acc={run}&accept-alternate-locations=yes,
    解析JSON中 files[].locations[].link 得到可下载URL
    |One request per run: GET {api}?acc={run}&accept-alternate-locations=yes,
    parsing files[].locations[].link from the JSON response

    requests未安装时不抛异常,查询优雅降级为空结果
    |Without requests installed, queries degrade gracefully to empty results
    """

    def __init__(self, logger=None, api_url: str = DEFAULT_SDL_API,
                 timeout: int = 60, retry_attempts: int = 3,
                 session=None,
                 max_workers: int = DEFAULT_MAX_WORKERS,
                 session_factory: Optional[Callable] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.api_url = api_url
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.max_workers = max_workers
        self.requests = get_requests()
        if self.requests is not None:
            self.session = session or self._create_session(retry_attempts)
            # 并发查询阶段每次新建会话,避免多线程共享Session
            # |Fresh session per query to avoid cross-thread sharing
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
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry

        session = self.requests.Session()
        retries = Retry(
            total=retry_attempts,
            backoff_factor=0.5,
            status_forcelist=[500, 502, 503, 504]
        )
        session.mount('https://', HTTPAdapter(max_retries=retries))
        return session

    def search_links(self, run_ids: List[str]) -> Dict[str, List[str]]:
        """
        查询NCBI SDL下载链接|Query NCBI SDL download links

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
                "未安装requests,NCBI SDL回退不可用|requests not installed, NCBI SDL fallback unavailable"
            )
            return {}

        self.logger.info(f"NCBI SDL查询|NCBI SDL query: {len(run_ids)} 个IDs|IDs")
        results: Dict[str, List[str]] = {}
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            for run_id, links in zip(run_ids, pool.map(self._query_single, run_ids)):
                if links:
                    results[run_id] = links
        return results

    def _query_single(self, run_id: str) -> List[str]:
        """查询单个Run的SDL记录|Query the SDL record of a single run"""
        requests = self.requests
        if requests is None:
            return []

        params = {"acc": run_id, "accept-alternate-locations": "yes"}
        try:
            self.logger.debug(f"NCBI SDL单个查询|NCBI SDL single query: {run_id}")
            session = self.session_factory()
            response = session.get(self.api_url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return self._parse_sdl(response.json(), run_id)
        except requests.exceptions.RequestException as e:
            self.logger.warning(f"NCBI SDL查询失败|NCBI SDL query failed: {run_id}: {e}")
            return []
        except ValueError as e:
            # SDL返回非JSON时降级为空|Non-JSON responses degrade to empty
            self.logger.warning(f"NCBI SDL响应解析失败|NCBI SDL response parse failed: {run_id}: {e}")
            return []

    def _parse_sdl(self, payload: dict, run_id: str) -> List[str]:
        """
        解析SDL JSON提取文件位置链接|Parse SDL JSON for file location links

        Args:
            payload: SDL API返回的JSON|JSON returned by the SDL API
            run_id: 请求的Run ID|Requested Run ID

        Returns:
            下载链接列表|List of download links, or empty if no files located
        """
        links = []
        for result in payload.get("result") or []:
            for file_entry in result.get("files") or []:
                for location in file_entry.get("locations") or []:
                    link = location.get("link")
                    if link:
                        links.append(link)
        return links

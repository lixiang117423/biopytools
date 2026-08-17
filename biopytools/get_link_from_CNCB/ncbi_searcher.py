"""
NCBI SDL数据定位回退搜索器|NCBI SRA Data Locator Fallback Searcher

CNCB FTP镜像与ENA都没有链接的INSDC Run(如仅存NCBI S3的PacBio source BAM,
ENA filereport有元数据但fastq_ftp/sra_ftp为空),通过NCBI SDL API定位真实文件位置
|For INSDC runs missed by both the CNCB FTP mirror and ENA (e.g. NCBI-only PacBio
source BAMs, where the ENA filereport lists metadata with empty link columns),
locate the real files via the NCBI SDL API
"""

import logging
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional

from .utils import get_requests, create_retry_session

# 默认SDL API地址|Default SDL API endpoint
DEFAULT_SDL_API = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

# 并发查询线程数|Concurrency for per-run queries
DEFAULT_MAX_WORKERS = 8


class NCBISDLSearcher:
    """NCBI SDL链接搜索器|NCBI SDL Link Searcher

    每个Run一次请求: GET {api}?acc={run}&accept-alternate-locations=yes,
    从JSON中为每个文件条目选一个可用位置链接
    |One request per run: GET {api}?acc={run}&accept-alternate-locations=yes,
    picking one usable location link per file entry from the JSON

    会话管理:注入session=时所有查询直接复用该会话;未注入时每个工作线程
    惰性自建一个带重试的会话(线程本地,复用TLS连接),search_links结束后统一关闭
    |Sessions: an injected session= is reused for all queries; otherwise each
    worker thread lazily builds one retry session (thread-local, TLS reuse),
    all closed together when search_links finishes

    requests未安装时不抛异常,查询优雅降级为空结果
    |Without requests installed, queries degrade gracefully to empty results
    """

    def __init__(self, logger=None, api_url: str = DEFAULT_SDL_API,
                 timeout: int = 60, retry_attempts: int = 3,
                 session=None,
                 max_workers: int = DEFAULT_MAX_WORKERS):
        self.logger = logger or logging.getLogger(__name__)
        self.api_url = api_url
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.max_workers = max_workers
        self.requests = get_requests()
        self._injected_session = session
        self._thread_local = threading.local()
        # 自建会话登记表,查询结束后关闭|Owned sessions, closed after queries
        self._owned_sessions: List = []

    def _get_session(self):
        """取当前线程的会话:注入优先,否则线程本地自建|Thread session: injected first, else thread-local"""
        if self._injected_session is not None:
            return self._injected_session
        session = getattr(self._thread_local, "session", None)
        if session is None:
            session = self._create_session(self.retry_attempts)
            self._thread_local.session = session
            self._owned_sessions.append(session)
        return session

    def _close_owned_sessions(self):
        """关闭自建会话;注入的会话归调用方,不动|Close owned sessions; injected ones belong to the caller"""
        for session in self._owned_sessions:
            try:
                session.close()
            except Exception as e:
                self.logger.debug(f"关闭会话失败|Failed to close session: {e}")
        self._owned_sessions.clear()

    def _create_session(self, retry_attempts: int):
        """创建带重试的HTTP会话|Create HTTP session with retries"""
        if self.requests is None:
            raise ImportError("requests未安装|requests not installed")
        return create_retry_session(self.requests, retry_attempts)

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
        try:
            with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                for run_id, links in zip(run_ids, pool.map(self._query_single, run_ids)):
                    if links:
                        results[run_id] = links
        finally:
            self._close_owned_sessions()
        return results

    def _query_single(self, run_id: str) -> List[str]:
        """查询单个Run的SDL记录|Query the SDL record of a single run"""
        requests = self.requests
        if requests is None:
            return []

        params = {"acc": run_id, "accept-alternate-locations": "yes"}
        try:
            self.logger.debug(f"NCBI SDL单个查询|NCBI SDL single query: {run_id}")
            session = self._get_session()
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
        解析SDL JSON,每个文件条目只选一个可下载链接|Parse SDL JSON, one usable link per file entry

        同一文件的locations[]是同一字节的云镜像(S3/GCS/NCBI SOS),全部收集会让下载脚本
        重复下载同一文件;payRequired位置需请求方付费签名,纯wget不可用;noqual(.lite)
        是无质量精简副本,与完整文件数据重复——均跳过
        |A file's locations[] are cloud mirrors of the same bytes; collecting all of them
        makes the download script fetch the same file repeatedly. payRequired locations
        need requester-pays signing (unusable with plain wget), and noqual (.lite) entries
        are quality-less duplicates of the full file — all skipped
        """
        # SDL异常时可能返回null或数组,防御性降级|SDL may return null/list on errors; degrade
        if not isinstance(payload, dict):
            self.logger.warning(
                f"NCBI SDL响应非预期格式|NCBI SDL unexpected payload format: {run_id}"
            )
            return []

        links = []
        for result in payload.get("result") or []:
            if not isinstance(result, dict):
                continue
            for file_entry in result.get("files") or []:
                if not isinstance(file_entry, dict):
                    continue
                # noqual(.lite)与完整文件重复|noqual (.lite) duplicates the full file
                if file_entry.get("noqual"):
                    self.logger.debug(
                        f"跳过noqual精简副本|Skipping noqual lite copy: "
                        f"{run_id}/{file_entry.get('name')}"
                    )
                    continue
                link = self._pick_file_link(file_entry, run_id)
                if link:
                    links.append(link)
        return links

    def _pick_file_link(self, file_entry: dict, run_id: str) -> Optional[str]:
        """
        为单个文件选择第一个免费可用位置|Pick the first free usable location for one file

        Args:
            file_entry: SDL的files[]条目|A files[] entry from SDL
            run_id: 请求的Run ID|Requested Run ID

        Returns:
            下载链接;全部位置付费或无link时返回None(并告警)
            |Download link; None (with a warning) when every location is pay-required or linkless
        """
        file_name = file_entry.get("name") or run_id
        for location in file_entry.get("locations") or []:
            if not isinstance(location, dict):
                continue
            link = location.get("link")
            # 无link=待取回(rehydrationRequired);payRequired=请求方付费,wget不可用
            # |No link = pending rehydration; payRequired = requester-pays, unusable via wget
            if link and not location.get("payRequired"):
                return link
        self.logger.warning(
            f"文件无可用的免费下载位置(仅付费镜像或待取回)|File has no free downloadable "
            f"location (pay-required mirrors or pending rehydration only): {run_id}/{file_name}"
        )
        return None

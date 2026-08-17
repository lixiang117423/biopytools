"""
GSA HTTP下载链接搜索器|GSA HTTP Download Link Searcher

CRR是CNCB GSA的原生Run ID,不属于INSDC镜像(FTP /INSDC*)和ENA。
GSA数据通过HTTPS目录列表公开: https://download.cncb.ac.cn/gsa2/{CRA}/{CRR}/
旧项目在 /gsa/{CRA}/{CRR}/(gsa2查不到时回退)
|CRR is GSA's native run ID, absent from the INSDC FTP mirror and ENA.
GSA data is exposed via HTTPS autoindex: https://download.cncb.ac.cn/gsa2/{CRA}/{CRR}/
with fallback to /gsa/{CRA}/{CRR}/ for older projects
"""

import logging
import re
from typing import List, Optional

from .utils import get_requests

# 默认GSA下载主机|Default GSA download host
DEFAULT_GSA_BASE = "https://download.cncb.ac.cn"

# NGDC GSA搜索页(服务端渲染结果,可反查CRR所属CRA)
# |NGDC GSA search page (server-rendered results, resolves a CRR's parent CRA)
GSA_SEARCH_URL = "https://ngdc.cncb.ac.cn/gsa/search"

# NGDC拒绝python-requests默认UA(直接断连,实测),须伪装浏览器UA
# |NGDC drops the default python-requests UA (verified live); send a browser UA
BROWSER_UA = "Mozilla/5.0 (X11; Linux x86_64) biopytools/get-link-from-CNCB"

# 依次尝试的目录布局|Directory layouts to try in order
GSA_LAYOUTS = ("/gsa2/{cra}/{run}/", "/gsa/{cra}/{run}/")

# GSA数据文件名模式: {CRR}.fq.gz / {CRR}_f1.fq.gz / {CRR}_r2.fq.gz 等
# |GSA file naming patterns: {CRR}.fq.gz / {CRR}_f1.fq.gz / {CRR}_r2.fq.gz etc.
FILE_PATTERN = r"^{run}(_f1|_f2|_r1|_r2|_1|_2)?\.(fq|fastq)\.gz$"


class GSAHTTPSearcher:
    """GSA下载链接搜索器|GSA Download Link Searcher"""

    def __init__(self, logger=None, base_url: str = DEFAULT_GSA_BASE,
                 timeout: int = 60, retry_attempts: int = 3, session=None):
        self.logger = logger or logging.getLogger(__name__)
        self.base_url = base_url
        self.timeout = timeout
        self.requests = get_requests()
        # 默认会话带重试适配器(镜像ENA搜索器);未安装requests时session为None
        # |Default session mounts a retry adapter (mirrors the ENA searcher);
        # session stays None without requests installed
        if session is not None:
            self.session = session
        elif self.requests is not None:
            self.session = self._create_session(retry_attempts)
        else:
            self.session = None
        # 统一补浏览器UA(NGDC拒绝默认requests UA);调用方显式设置的UA优先
        # |Apply the browser UA to whichever session won (NGDC drops the default
        # requests UA); a caller-provided UA takes precedence
        if self.session is not None:
            # 注入的测试假会话可能没有headers属性,防御性跳过
            # |Injected test fakes may lack a headers attribute; skip defensively
            headers = getattr(self.session, "headers", None)
            if headers is not None:
                headers.setdefault("User-Agent", BROWSER_UA)
        # CRR→CRA解析结果缓存|Cache of CRR→CRA resolutions
        self._cra_cache = {}

    def resolve_cra(self, run_id: str) -> Optional[str]:
        """
        反查CRR所属的CRA项目编号|Resolve the parent CRA accession of a CRR run

        通过NGDC GSA搜索页(服务端渲染)解析,供项目列填了PRJCA等BioProject编号时自动转换
        |Uses the server-rendered NGDC GSA search page; auto-converts when the project
        column holds a BioProject ID like PRJCA instead of the CRA accession

        Args:
            run_id: GSA Run ID(CRR开头)|GSA Run ID (CRR prefix)

        Returns:
            CRA编号,查不到或网络失败返回None|CRA accession, or None if not found / network error
        """
        if run_id in self._cra_cache:
            return self._cra_cache[run_id]

        if self.requests is None:
            self.logger.error(
                "未安装requests,无法反查CRA|requests not installed, cannot resolve CRA"
            )
            return None

        try:
            self.logger.info(f"NGDC反查CRA|Resolving CRA via NGDC: {run_id}")
            response = self.session.get(
                GSA_SEARCH_URL, params={"searchTerm": run_id}, timeout=self.timeout
            )
            response.raise_for_status()
            # 结果行为 href="browse/{CRA}/{run}",run必须精确匹配防前缀撞名
            # |Result rows look like href="browse/{CRA}/{run}"; exact match avoids prefix collisions
            pattern = re.compile(r'href="browse/(CRA\d+)/%s"' % re.escape(run_id))
            match = pattern.search(response.text)
            cra = match.group(1) if match else None
        except Exception as e:
            # 反查是逐Run的可选步骤,任何异常(含response.text解码错误)都降级为失败,
            # 不得中断整批提取|Resolution is a per-run optional step; any exception
            # (incl. response.text decode errors) degrades to a miss and must not
            # abort the whole batch
            self.logger.warning(f"NGDC反查失败|CRA resolution failed: {run_id}: {e}")
            cra = None

        # 仅缓存成功结果;失败(网络抖动/未收录)不缓存,同一Run再现时可重试
        # |Cache only successes; failures (network blips/unindexed) are not cached
        # so the same run is retried when it reappears
        if cra is not None:
            self._cra_cache[run_id] = cra
        return cra

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

    def search_files(self, cra_accession: str, run_id: str) -> List[str]:
        """
        查找CRR Run的下载文件URL|Find download file URLs for a CRR run

        Args:
            cra_accession: GSA项目编号(CRA开头)|GSA project accession (CRA prefix)
            run_id: GSA Run ID(CRR开头)|GSA Run ID (CRR prefix)

        Returns:
            文件URL列表|List of file URLs, or empty if not found
        """
        if self.requests is None:
            self.logger.error(
                "未安装requests,GSA查询不可用|requests not installed, GSA lookup unavailable"
            )
            return []

        for layout in GSA_LAYOUTS:
            dir_path = layout.format(cra=cra_accession, run=run_id)
            url = f"{self.base_url}{dir_path}"
            self.logger.info(f"GSA查询|GSA query: {url}")
            content = self._fetch(url)
            if content is None:
                continue

            names = self._parse_autoindex(content, run_id)
            if names:
                self.logger.info(
                    f"GSA找到 {len(names)} 个文件|GSA found {len(names)} files: {run_id}"
                )
                return [f"{self.base_url}{dir_path}{name}" for name in sorted(names)]

        self.logger.warning(f"GSA未找到文件|No GSA files found for: {run_id}")
        return []

    def _fetch(self, url: str) -> Optional[str]:
        """GET目录列表,404或异常返回None|GET a directory listing; None on 404 or error"""
        requests = self.requests
        try:
            response = self.session.get(url, timeout=self.timeout)
            if response.status_code == 404:
                return None
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            self.logger.warning(f"GSA查询失败|GSA query failed: {url}: {e}")
            return None

    def _parse_autoindex(self, content: str, run_id: str) -> List[str]:
        """
        解析Apache autoindex HTML,提取该Run的数据文件|Parse autoindex HTML for the run's data files

        Args:
            content: 目录列表HTML|Directory listing HTML
            run_id: Run ID|Run ID

        Returns:
            匹配的文件名列表|List of matching filenames
        """
        pattern = re.compile(FILE_PATTERN.format(run=re.escape(run_id)))
        names = []
        for match in re.finditer(r'href="([^"]+)"', content):
            name = match.group(1)
            if pattern.fullmatch(name):
                names.append(name)
        return names

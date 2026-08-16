#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conda 环境合并分析|Conda env merge analysis
比对|Compare: 超算实际环境(conda_envs_backup/*.yaml) x 代码引用清单(docs/conda_env_audit.csv)
纯标准库|Pure stdlib, 不依赖 PyYAML
用法|Usage: python3 scripts/analyze_conda_envs_merge.py
"""

import csv
import os
import re
import sys
from collections import Counter, defaultdict

BACKUP = "conda_envs_backup"
AUDIT_CSV = "docs/conda_env_audit.csv"

# 常见底层依赖前缀,不计入"工具"|Common dependency prefixes, not counted as tools
NOISE_PREFIXES = (
    "lib", "_libgcc", "_openmp", "alsa-", "brotli", "bzip2", "c-ares", "ca-",
    "cairo", "contourpy", "cycler", "cyrus-", "dbus", "double-", "font-",
    "fontconfig", "fonts-", "freetype", "fribidi", "gdk-", "gettext", "giflib",
    "glib", "graphite2", "gstreamer", "gtk", "harfbuzz", "icu", "jpeg", "lcms2",
    "lerc", "libdeflate", "libffi", "libgcc", "libglib", "libiconv", "libidn2",
    "libllvm", "libnghttp2", "libpng", "libssh2", "libstdcxx", "libtiff",
    "libtool", "libunistring", "libuuid", "libwebp", "libxcb", "libxml2",
    "libxslt", "libzlib", "lz4-", "lzo", "ncurses", "nspr", "nss", "openjpeg",
    "openssl", "p11-kit", "pango", "pcre2", "pixman", "poppler", "pyparsing",
    "readline", "snappy", "sqlite", "tk", "tzdata", "xz", "zlib", "zstd",
    "krb5", "keyutils", "lame", "libaec", "libcap", "libcups", "libevent",
    "libexpat", "libgfortran", "libgomp", "libnsl", "libpq", "libsqlite",
    "libssh", "libsystemd0", "libtirpc", "libudev", "libuv", "libvpx",
    "libxkbcommon", "libzopfli", "libblas", "libcblas", "liblapack", "libopenblas",
    "mpfr", "gmp", "gsl", "fftw", "hdf5", "netcdf4", "curl", "certifi",
    "wheel", "setuptools", "pip", "packaging", "platformdirs", "pluggy",
    "pytest", "six", "tzlocal", "tornado", "yaml", "pyyaml", "expat",
    "colorama", "decorator", "exceptiongroup", "execnet", "filelock", "fsspec",
    "importlib", "iniconfig", "jmespath", "markupsafe", "mdurl", "more-itertools",
    "msgpack", "munkres", "networkx", "numpy-base", "olefile", "pathlib2",
    "psutil", "pycparser", "pygments", "pyproject", "pysocks", "python-dateutil",
    "python-editor", "python-json-logger", "python-tzdata", "python_abi", "pyyaml",
    "regex", "s3transfer", "sacremoses", "sed", "shellingham", "sqlalchemy",
    "structlog", "tabulate", "text-unidecode", "threadpoolctl", "toolz",
    "tomli", "toml", "tqdm", "typing", "unicodedata2", "urllib3", "wcwidth",
    "win-unicode-console", "wrapt", "zipp", "zstandard", "zope", "appdirs",
    "attrs", "automat", "backcall", "backports", "bcrypt", "bitarray", "bleach",
    "boltons", "bottleneck", "bs4", "bsdiff4", "bwidget", "cached_property",
    "cachetools", "charset-normalizer", "click", "cloudpickle", "colorlog",
    "commonmark", "constantly", "cpp-htslib", "cramjam", "crcmod", "cryptography",
    "cssselect", "cyclonedx", "dill", "docutils", "entrypoints", "et_xmlfile",
    "et_xmlfile", "faker", "fastobo", "flit", "flit-core", "fonttools",
    "formulaic", "frozendict", "fsspec", "future", "gast", "glob2",
    "google-auth", "graphviz", "grpcio", "hatch", "h11", "h5py", "heapdict",
    "httpretty", "humanfriendly", "hyperlink", "idna", "imagesize", "inflection",
    "iniconfig", "isodate", "itemadapter", "jaraco", "jedi", "jeepney",
    "jinja2", "joblib", "jupyter", "kiwisolver", "ld_impl_linux-64", "llvmlite",
    "locket", "lxml", "mako", "marshmallow", "matplotlib", "matplotlib-base",
    "mccabe", "mdit-py", "mistune", "mkl", "more-itertools", "mpmath", "multidict",
    "natsort", "nbconvert", "notebook", "numba", "numpy", "openjdk", "openpyxl",
    "orjson", "pandas", "paramiko", "parso", "partd", "pathspec", "patsy",
    "pendulum", "pexpect", "phonenumbers", "pickleshare", "pillow", "pint",
    "pip-audit", "pkginfo", "plotly", "ply", "prompt_toolkit", "propcache",
    "ptyprocess", "pure_eval", "py-cpuinfo", "pyasn1", "pycosat", "pycurl",
    "pydantic", "pygments", "pyjwt", "pynacl", "pyopenssl", "pyparsing",
    "pyqt", "pysam", "pytest-cov", "python-picard", "pytz", "pyu2f", "pyzmq",
    "qtconsole", "qtpy", "rapidfuzz", "requests", "rich", "rpds-py", "ruamel",
    "s2n", "scikit-learn", "scipy", "seaborn", "secretstorage", "selenium",
    "send2trash", "setproctitle", "setuptools_scm", "sip", "soupsieve",
    "stack_data", "statsmodels", "sympy", "syrupy", "testpath", "tinycss2",
    "traitlets", "typer", "uc-micro-py", "ujson", "unidecode", "uri-template",
    "validate-pyproject", "virtualenv", "wcwidth", "webencodings", "widgetsnbextension",
)


def parse_env_yaml(path):
    """解析 conda env export yaml|Parse conda env export yaml (no PyYAML)"""
    name = None
    conda_deps = {}
    pip_deps = {}
    in_deps = False
    in_pip = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("name:") and name is None:
                name = line.split(":", 1)[1].strip()
                continue
            if line == "dependencies:":
                in_deps = True
                continue
            if in_deps:
                if line == "prefix:":
                    break
                if line.startswith("- pip:"):
                    in_pip = True
                    continue
                m = re.match(r"\s*-\s+(\S+?)(?:==?|=)([^= ]+)", line)
                if in_pip:
                    if m:
                        pip_deps[m.group(1).lower()] = m.group(2)
                else:
                    if m:
                        conda_deps[m.group(1).lower()] = m.group(2)
    return name, conda_deps, pip_deps


def dep_version(deps, name):
    v = deps.get(name)
    if v is not None:
        return v
    for k, val in deps.items():
        if k.startswith(name + "="):
            return val
    return None


def main():
    # 1. 解析全部 yaml|Parse all yamls
    envs = {}
    for fn in sorted(os.listdir(BACKUP)):
        if not fn.endswith(".yaml"):
            continue
        name, conda_deps, pip_deps = parse_env_yaml(os.path.join(BACKUP, fn))
        if name is None:
            name = fn[:-5]
        envs[name] = {
            "file": fn, "conda": conda_deps, "pip": pip_deps,
            "n_conda": len(conda_deps), "n_pip": len(pip_deps),
            "python": dep_version(conda_deps, "python"),
            "r": dep_version(conda_deps, "r-base"),
            "openjdk": dep_version(conda_deps, "openjdk"),
            "perl": dep_version(conda_deps, "perl"),
        }

    # 2. 代码引用清单|Code-cited canonical envs
    cited = set()
    cited_code = set()
    cited_docs = set()
    with open(AUDIT_CSV, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            c = row["canonical_env"]
            cited.add(c)
            (cited_code if row["scope"] == "code" else cited_docs).add(c)

    # 3. 比对|Compare (名字大小写不敏感)
    yaml_lower = {k.lower(): k for k in envs}
    matched = {}
    dead = []
    for c in sorted(cited):
        hit = yaml_lower.get(c.lower())
        if hit:
            matched[c] = hit
        else:
            dead.append(c)
    orphans = []
    matched_names = {v for v in matched.values()}
    for k in sorted(envs):
        if k not in matched_names:
            orphans.append(k)

    print(f"=== 超算实际环境|HPC actual envs: {len(envs)} ===")
    print(f"=== 代码引用环境|Code-cited canonical envs: {len(cited)} (code:{len(cited_code)}, docs-only:{len(cited_docs)}) ===")
    print(f"=== 匹配|Matched: {len(matched)} ===")
    print(f"=== 死依赖(代码引用但超算不存在)|Dead refs: {len(dead)} ===")
    for c in dead:
        print(f"  {c}")
    print(f"=== 孤儿环境(超算有但代码未引用)|Orphan envs: {len(orphans)} ===")
    print("  " + ", ".join(orphans))

    # 4. python/R 版本分布|Interpreter version distribution
    print()
    print("=== Python 版本分布|Python version distribution ===")
    pyc = Counter(e["python"] for e in envs.values())
    for v, n in sorted(pyc.items(), key=lambda kv: str(kv[0])):
        print(f"  python {v}: {n} 个环境")
    print("=== R 版本分布|R version distribution ===")
    rc = Counter(e["r"] for e in envs.values())
    for v, n in sorted(rc.items(), key=lambda kv: str(kv[0])):
        print(f"  r-base {v}: {n} 个环境")
    print("=== 含 openjdk 的环境|Envs with openjdk ===")
    for k, e in envs.items():
        if e["openjdk"]:
            print(f"  {k}: openjdk {e['openjdk']}")

    # 5. 每个环境的"特征工具"|Signature tools per env (top-level tools, noise filtered)
    print()
    print("=== 特征工具(前20个)|Signature tools (top 20) ===")
    for k, e in envs.items():
        tools = []
        for pkg in sorted(e["conda"]):
            if pkg.startswith(NOISE_PREFIXES):
                continue
            if re.search(r"^(py|r-|bioconductor-)", pkg) and not re.match(
                    r"^(python|r-base|perl|openjdk)$", pkg):
                continue
            tools.append(pkg)
        sig = tools[:25]
        print(f"{k} [{e['n_conda']}+{e['n_pip']} pkgs]: {', '.join(sig)}")


if __name__ == "__main__":
    main()

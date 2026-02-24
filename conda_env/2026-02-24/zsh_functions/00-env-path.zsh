# =============================================================================
#  00-env-path.zsh - 环境变量与路径设置模块
#  Environment Variables & Path Configuration Module
# =============================================================================

# 加载统一路径配置（如果还未加载）
if [[ -z "$ZSH_PATH_CONFIG_LOADED" ]]; then
  source "${ZSH_FUNCTIONS_DIR:-$HOME/zsh/functions}/00-path-config.zsh"
fi

# -----------------------------------------------------------------------------
#  环境变量与路径设置 (Environment & Path)
# -----------------------------------------------------------------------------

# 本地用户目录优先于系统路径
export PATH="${LOCAL_BIN_DIR:-$HOME/.local/bin}:$PATH"

# 特定软件路径 - 使用统一配置
if [[ -n "${SOFTWARE_DIR}" ]]; then
  export PATH="${SOFTWARE_DIR}/metaWRAP/bin:$PATH"
  export PATH="${SOFTWARE_DIR}/HiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin:$PATH"
  export PATH="${SOFTWARE_DIR}/HiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin/utils:$PATH"
fi

# 数据库路径 - 使用统一配置
if [[ -n "${DATABASE_DIR}" ]]; then
  export BUSCO_DATASETS_PATH="${DATABASE_DIR}/busco"
fi

# 模块加载成功标记
export ZSH_MODULE_ENV_PATH_LOADED=1
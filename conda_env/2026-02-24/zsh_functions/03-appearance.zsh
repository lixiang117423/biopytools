# =============================================================================
#  03-appearance.zsh - 颜色与外观配置模块
#  Colors & Appearance Configuration Module
# =============================================================================

# -----------------------------------------------------------------------------
#  详细的文件类型颜色配置 (Detailed File Type Colors)
# -----------------------------------------------------------------------------
export LS_COLORS='rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=00:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32'

# 文本文件
export LS_COLORS="$LS_COLORS:*.txt=00;33:*.log=00;33:*.md=01;33"

# 图片文件
export LS_COLORS="$LS_COLORS:*.jpg=01;35:*.jpeg=01;35:*.png=01;35:*.gif=01;35"

# 音频文件
export LS_COLORS="$LS_COLORS:*.mp3=01;36:*.wav=01;36:*.flac=01;36"

# 视频文件
export LS_COLORS="$LS_COLORS:*.mp4=01;31:*.avi=01;31:*.mkv=01;31"

# 压缩文件
export LS_COLORS="$LS_COLORS:*.zip=01;31:*.tar=01;31:*.gz=01;31:*.bz2=01;31"

# 代码文件
export LS_COLORS="$LS_COLORS:*.py=01;32:*.sh=01;32:*.js=01;32:*.html=01;32"

# 文档文件
export LS_COLORS="$LS_COLORS:*.pdf=00;32:*.doc=00;32:*.docx=00;32"

# 模块加载成功标记
export ZSH_MODULE_APPEARANCE_LOADED=1
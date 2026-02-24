# =============================================================================
#  04-aliases.zsh - 别名定义模块
#  Aliases Definition Module
# =============================================================================

# -----------------------------------------------------------------------------
#  系统与导航别名 (System & Navigation Aliases)
# -----------------------------------------------------------------------------
alias ls='ls --color=auto'
alias l='ls -CF --color=auto'
alias ll='ls -lFh --color=auto'
alias la='ls -A --color=auto'
alias lh='ls -lh --color=auto'
alias lk='ls -lhtr'          # 按修改时间逆序排序
alias lw='ll | wc -l'         # 计算文件和目录总数
alias lz="eza -l --icons --group-directories-first --time-style=long-iso --sort=size"
alias lt="eza -l --icons --group-directories-first --time-style=long-iso --sort=time"
alias less="less -S"

alias ..='cd ..'
alias ...='cd .. && cd ..'
alias ....='cd ../../..'
alias .....='cd ../../../..'
alias grep='grep --color=auto'
alias fgrep='fgrep --color=auto'
alias egrep='egrep --color=auto'
alias hh="head"
alias tt='tail'
alias np='nohup'
alias pp='python'
alias cp='cp -r'
alias rr='rm'
alias xx='cat'
alias xz="cat ~/submitted_jobs.txt"
alias md="mkdir -p"

# -----------------------------------------------------------------------------
#  系统工具别名 (System Tools Aliases)
# -----------------------------------------------------------------------------
alias top='${MINIFORGE3_DIR:-$HOME/miniforge3}/envs/BioinfTools/bin/btop'

# -----------------------------------------------------------------------------
#  安全删除配置 (Safe Delete Configuration)
# -----------------------------------------------------------------------------
alias trash-put="${MINIFORGE3_DIR:-$HOME/miniforge3}/envs/Augustus_v.3.5.0/bin/trash-put"
alias truerm='command rm'     # 原始rm命令
alias rf='rm -rf'             # 危险：强制删除

# -----------------------------------------------------------------------------
#  Git 版本控制别名 (Git Version Control Aliases)
# -----------------------------------------------------------------------------
alias gs='git status'
alias ga='git add'
alias gc='git commit'
alias gp='git push'
alias gl='git log --oneline --graph'

# -----------------------------------------------------------------------------
#  集群管理别名 (Cluster Management Aliases)
# -----------------------------------------------------------------------------
alias cs="csub"
alias cv="cview"
# alias cj="cjobs"
# alias cj='cj_formatted'
alias cj='cjj'
alias cc='scancel -n'

# -----------------------------------------------------------------------------
#  Conda/Mamba 包管理别名 (Conda/Mamba Package Management Aliases)
# -----------------------------------------------------------------------------
alias mm='mamba'
alias mmi="mamba install --freeze-installed"
alias mmb="mamba install --freeze-installed -c bioconda"
alias mml='mamba env list'
alias mmd='conda deactivate'
alias mmt='mamba activate bioinftools'

# -----------------------------------------------------------------------------
#  下载与网络工具别名 (Download & Network Tools Aliases)
# -----------------------------------------------------------------------------
alias get='axel -n 30'
alias d='axel -n 30'
alias getaws="aws s3 cp --no-sign-request"

# -----------------------------------------------------------------------------
#  生物信息学工具别名 (Bioinformatics Tools Aliases)
# -----------------------------------------------------------------------------
alias gg="gunzip *.gz"
alias ss="sed -i -e 's/ //g' -e 's/dna:.\{1,1000\}//g' *.fa"
alias edge="${SOFTWARE_DIR:-$HOME/software}/edgeturbo/edgeturbo-client/edgeturbo"
alias bb="biopytools"
alias singularity='${MINIFORGE3_DIR:-$HOME/miniforge3}/envs/singularity_v.3.8.7/bin/singularity'
alias vcf2pcacluster='${SOFTWARE_DIR:-$HOME/software}/VCF2PCACluster-1.42/bin/VCF2PCACluster'
alias cpbackup='rsync -avc --no-t --delete ${BIOPYTOOLS_DIR:-$HOME/software/biopytools}/ ${BIOPYTOOLS_DIR:-$HOME/software/biopytools_backup}/'
# -----------------------------------------------------------------------------
#  复杂别名（通过函数实现） (Complex Aliases via Functions)
# -----------------------------------------------------------------------------
alias pg='_pg(){ ps -aux | grep "$1" };_pg'
alias git2hub='_f() { git add . && git commit -m "$1" && git push origin; }; _f'
alias gh='_gh(){ grep ">" "$1" | head; };_gh'
alias ww='_ww(){ wget "$1" -O "$2"; };_ww'
alias dl='wget -c'

# -----------------------------------------------------------------------------
#  备份系统别名 (Backup System Aliases)
# -----------------------------------------------------------------------------
alias backup-envs='backup_envs'
alias backup-conda='backup_envs'
alias smart-backup='smart_backup_envs'
alias check-changes='check_env_changes'
alias backup-status='backup_status'
alias setup-backup='setup_weekly_backup'
alias backup-history='backup_history'
alias analyze-history='analyze_history'
alias search-history='search_history'
alias debug-envs='debug_conda_envs'
alias test-backup='test_backup_envs'
alias health-check='health_check'

# -----------------------------------------------------------------------------
#  其他别名
# -----------------------------------------------------------------------------
alias uz='x'
alias extract='x'
alias unpack='x'
alias xk='x --keep'
alias xr='x --remove'
alias ccc='cd ~/project'


# 模块加载成功标记
export ZSH_MODULE_ALIASES_LOADED=1

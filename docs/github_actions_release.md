# GitHub Actions Release 管理配置说明

本配置实现了智能的 Release 创建和自动清理功能，避免 Release 数量和存储空间无限增长。

## 功能特性

### 1. 智能 Release 创建
- ✅ **仅在打标签时触发**：避免每次提交都生成 Release
- ✅ **自动识别版本类型**：区分正式版本和测试版本（alpha/beta/rc）
- ✅ **自动生成 Release Notes**：使用 GitHub 的自动生成功能
- ✅ **支持附件上传**：可配置上传构建产物

### 2. 自动清理旧 Release
- ✅ **保留策略**：
  - 正式版本：保留最新 10 个
  - 测试版本（alpha/beta/rc）：保留最新 5 个
- ✅ **时间限制**：删除超过 30 天的旧 Release
- ✅ **智能分类**：自动区分正式版本和测试版本
- ✅ **详细日志**：输出清理记录和释放的存储空间
- ✅ **Dry Run 模式**：支持测试运行，不实际删除

## 文件结构

```
.github/
├── workflows/
│   └── release.yml           # 主 workflow 文件
└── scripts/
    └── cleanup_releases.py    # 清理脚本
```

## 使用方法

### 创建 Release

**方式一：推送标签（推荐）**

```bash
# 创建正式版本
git tag v1.0.0
git push origin v1.0.0

# 创建测试版本
git tag v1.0.0-beta.1
git push origin v1.0.0-beta.1
```

标签命名规范：
- 正式版本：`v1.0.0`, `v2.1.3` 等
- 测试版本：包含 `alpha`, `beta`, `rc`, `pre` 等关键词
  - `v1.0.0-alpha.1`
  - `v1.0.0-beta.1`
  - `v1.0.0-rc.1`

**方式二：手动触发清理**

在 GitHub Actions 页面手动运行 workflow，选择 "cleanup-only" 模式。

### 配置参数

在 `.github/workflows/release.yml` 中可以调整以下参数：

```yaml
env:
  GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
  GITHUB_REPOSITORY: ${{ github.repository }}
  # 清理参数
  KEEP_RELEASES: 10           # 正式版本保留数量
  KEEP_PRERELEASES: 5         # 测试版本保留数量
  DELETE_AFTER_DAYS: 30       # 删除超过天数的 Release
  DRY_RUN: false              # 设为 true 仅测试不实际删除
```

### 自定义 Release 内容

编辑 `release.yml` 中的 "Create Release" 步骤，添加你的构建命令和上传文件：

```yaml
- name: Build project
  run: |
    # 示例：Python 项目
    pip install build
    python -m build

    # 示例：Node.js 项目
    npm run build

- name: Create Release
  uses: softprops/action-gh-release@v2
  with:
    files: |
      dist/*.whl
      dist/*.tar.gz
      # 添加更多文件...
```

## 清理规则详解

### 保留策略

1. **正式版本**（如 v1.0.0, v2.0.0）
   - 保留最新 10 个
   - 即使超过 10 个，如果未超过 30 天也会保留

2. **测试版本**（如 v1.0.0-beta.1, v2.0.0-alpha.1）
   - 保留最新 5 个
   - 更严格的保留策略

### 删除条件

Release 需要同时满足以下条件才会被删除：
1. 超过保留数量（正式版本 > 10，测试版本 > 5）
2. 创建时间超过 30 天

示例：

```
假设有 15 个正式版本：

v1.0.0 (创建于 100 天前)  ❌ 删除（超过保留数量 + 超过30天）
v1.1.0 (创建于 80 天前)   ❌ 删除（超过保留数量 + 超过30天）
v1.2.0 (创建于 60 天前)   ❌ 删除（超过保留数量 + 超过30天）
v1.3.0 (创建于 40 天前)   ❌ 删除（超过保留数量 + 超过30天）
v1.4.0 (创建于 35 天前)   ❌ 删除（超过保留数量 + 超过30天）
v1.5.0 (创建于 25 天前)   ✅ 保留（超过保留数量但未超过30天）
v1.6.0 (创建于 20 天前)   ✅ 保留（超过保留数量但未超过30天）
v1.7.0 ~ v2.0.0          ✅ 保留（最新10个版本）
```

## 查看清理日志

清理任务执行后，可以在以下位置查看日志：

1. **GitHub Actions 运行日志**
   - 进入 GitHub 仓库 → Actions → Release Management → 选择运行记录

2. **清理日志 Artifact**
   - Actions 运行页面 → Artifacts 部分
   - 下载 `cleanup-log-*.zip` 查看详细日志

日志内容包括：
- 当前 Release 统计
- 每个被删除 Release 的详细信息
- 附件列表和下载次数
- 释放的存储空间统计

## Dry Run 模式（测试运行）

在正式运行前，可以先测试清理效果：

```yaml
env:
  DRY_RUN: true  # 仅测试，不实际删除
```

Dry Run 模式会：
- ✅ 显示将要删除的 Release
- ✅ 计算释放的存储空间
- ✅ 输出完整的清理日志
- ❌ 不实际删除任何内容

## 手动触发清理

如果需要立即执行清理而不创建新 Release：

1. 进入 GitHub 仓库 → Actions
2. 选择 "Release Management" workflow
3. 点击 "Run workflow"
4. 勾选 "cleanup_only" 选项
5. 点击 "Run workflow"

## 高级配置

### 修改构建步骤

根据项目类型，修改 "Build project" 步骤：

**Python 项目：**
```yaml
- name: Build Python package
  run: |
    pip install build wheel
    python -m build
```

**Node.js 项目：**
```yaml
- name: Build Node.js package
  run: |
    npm install
    npm run build
```

**Go 项目：**
```yaml
- name: Build Go binaries
  run: |
    GOOS=linux GOARCH=amd64 go build -o dist/app-linux-amd64
    GOOS=darwin GOARCH=amd64 go build -o dist/app-darwin-amd64
    GOOS=windows GOARCH=amd64 go build -o dist/app-windows-amd64.exe
```

### 添加多个附件

```yaml
- name: Create Release
  uses: softprops/action-gh-release@v2
  with:
    files: |
      dist/*.whl
      dist/*.tar.gz
      README.md
      LICENSE
      docs/*.pdf
```

### 自定义 Release 说明

可以在 workflow 中添加生成 Release Notes 的步骤：

```yaml
- name: Generate Release Notes
  id: release_notes
  run: |
    echo "## 更新内容" > release_notes.md
    echo "" >> release_notes.md
    git log --pretty=format:"- %s" $(git describe --tags --abbrev=0 HEAD^)..HEAD >> release_notes.md

- name: Create Release
  uses: softprops/action-gh-release@v2
  with:
    body_path: release_notes.md
```

## 注意事项

1. **权限要求**
   - workflow 需要写入权限
   - 已在配置中添加 `permissions: contents: write`

2. **标签推送**
   - 只有推送标签才会创建 Release
   - 推送到 main 分支不会创建 Release

3. **首次运行建议**
   - 先设置 `DRY_RUN: true` 测试
   - 检查清理日志确认无误
   - 再改为 `DRY_RUN: false` 正式运行

4. **重要版本保护**
   - 重要的 LTS 版本可以手动标记为 "Keep"
   - 或调整 `KEEP_RELEASES` 参数增加保留数量

## 故障排查

### Release 未创建
- 检查标签格式是否符合 `v*.*.*` 模式
- 查看 Actions 日志确认是否有错误

### 清理未执行
- 检查 `GITHUB_TOKEN` 权限是否正确
- 确认 Release 数量是否超过保留阈值
- 查看清理日志了解详细信息

### 删除了重要版本
- 从 Git 标签恢复：`git tag <tag_name> <commit_sha>`
- 重新推送标签触发 Release 创建
- 考虑调整保留参数

## 最佳实践

1. **定期检查清理日志**：确保没有误删重要版本
2. **使用语义化版本**：方便识别正式版本和测试版本
3. **备份重要附件**：重要文件建议单独存档
4. **测试后再正式运行**：先用 Dry Run 模式测试
5. **监控存储空间**：在 GitHub Settings 查看存储使用情况

## 相关链接

- [GitHub Actions 文档](https://docs.github.com/en/actions)
- [softprops/action-gh-release](https://github.com/softprops/action-gh-release)
- [PyGithub 文档](https://github.com/PyGithub/PyGithub)
- [语义化版本规范](https://semver.org/lang/zh-CN/)

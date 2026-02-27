# GitHub Actions Release 管理

智能的 Release 创建和自动清理解决方案，避免 GitHub Release 数量和存储空间无限增长。

## 快速开始

### 1. 创建 Release

```bash
# 创建正式版本
git tag v1.0.0
git push origin v1.0.0

# 创建测试版本
git tag v1.0.0-beta.1
git push origin v1.0.0-beta.1
```

### 2. 自动清理规则

- **正式版本**：保留最新 10 个
- **测试版本**（alpha/beta/rc）：保留最新 5 个
- **时间限制**：删除超过 30 天的旧 Release

### 3. 手动清理

进入 GitHub Actions → "Release Management" → "Run workflow" → 勾选 "cleanup_only"

## 功能特性

✅ 仅在打标签时创建 Release，避免每次提交都生成
✅ 自动区分正式版本和测试版本
✅ 智能清理策略，保留重要版本
✅ 详细的清理日志和存储空间统计
✅ 支持 Dry Run 模式测试

## 文件说明

- `.github/workflows/release.yml` - 主 workflow 配置
- `.github/scripts/cleanup_releases.py` - 清理脚本
- `docs/github_actions_release.md` - 详细文档

## 配置参数

在 `release.yml` 中修改：

```yaml
KEEP_RELEASES: 10           # 正式版本保留数量
KEEP_PRERELEASES: 5         # 测试版本保留数量
DELETE_AFTER_DAYS: 30       # 删除超过天数的 Release
DRY_RUN: false              # 测试模式（true 不实际删除）
```

## 查看完整文档

详见 [docs/github_actions_release.md](../docs/github_actions_release.md)

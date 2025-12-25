# Git 开发工作流程文档

## 📋 分支说明

- **main 分支**：公开的稳定版本，只包含已完成并测试通过的功能
- **development 分支**：开发分支，包含所有代码（包括未完成的功能）

## 🚀 日常开发流程

### 1. 开始开发前

```bash
# 切换到 development 分支
git checkout development

# 拉取最新代码（如果有多设备开发）
git pull origin development

# 查看当前状态
git status
```

### 2. 开发过程中

```bash
# 编写代码...

# 查看修改了哪些文件
git status

# 查看具体改动内容
git diff

# 添加文件到暂存区
git add .                    # 添加所有文件
# 或
git add 文件名               # 添加指定文件

# 提交到本地仓库
git commit -m "描述你的改动内容"

# 推送到 GitHub（重要！防止丢失）
git push origin development
```

### 3. 每天结束时

```bash
# 确保所有改动都已提交并推送
git add .
git commit -m "YYYYMMDD 今日开发内容描述"
git push origin development
```

## ✅ 功能完成后合并到 main

### 步骤 1：查看提交历史

```bash
# 确保在 development 分支
git checkout development

# 查看提交历史（找到要合并的功能）
git log --oneline

# 输出示例：
# a1b2c3d 完成用户登录功能    ← 这个要合并
# e4f5g6h 修复登录页面bug     ← 这个也要合并
# i7j8k9l 开发中的支付功能     ← 这个还没完成，不合并
# m1n2o3p 20251225_commit_for_backup
```

**记下要合并的 commit ID**（例如：`a1b2c3d` 和 `e4f5g6h`）

### 步骤 2：切换到 main 分支

```bash
# 切换到 main 分支
git checkout main

# 拉取最新的 main 分支（如果有多人协作）
git pull origin main
```

### 步骤 3：Cherry-pick 选择性合并

```bash
# 合并指定的提交（按顺序，从旧到新）
git cherry-pick e4f5g6h
git cherry-pick a1b2c3d

# 如果有多个连续的提交，可以用范围
# git cherry-pick e4f5g6h^..a1b2c3d
```

### 步骤 4：推送到 GitHub

```bash
# 推送 main 分支
git push origin main

# 切回 development 继续开发
git checkout development
```

## 🔧 常用命令速查

### 查看状态

```bash
git status              # 查看当前状态
git branch              # 查看本地分支（当前分支有*标记）
git branch -a           # 查看所有分支（包括远程）
git log --oneline       # 查看提交历史
git log --oneline -10   # 查看最近10条提交
```

### 撤销操作

```bash
# 撤销工作区的修改（还未 add）
git checkout -- 文件名

# 撤销已添加到暂存区的文件（已 add 但未 commit）
git reset HEAD 文件名

# 撤销最后一次提交但保留改动（已 commit 但未 push）
git reset --soft HEAD~1

# 查看某个提交的详细信息
git show commit-id
```

### 分支切换

```bash
git checkout development    # 切换到 development 分支
git checkout main          # 切换到 main 分支
```

## ⚠️ 注意事项

### 1. 每天必做：推送到 GitHub
```bash
git push origin development
```
**防止本地代码丢失！**

### 2. 提交信息要清晰
```bash
# ✅ 好的提交信息
git commit -m "添加用户登录功能"
git commit -m "修复登录页面样式bug"
git commit -m "20251225 完成用户模块开发"

# ❌ 不好的提交信息
git commit -m "update"
git commit -m "fix"
git commit -m "aaa"
```

### 3. Cherry-pick 的顺序很重要
- 按照**时间顺序**（从旧到新）进行 cherry-pick
- 如果提交之间有依赖关系，必须按顺序合并

### 4. 合并前先测试
- 在 development 分支充分测试功能
- 确认没有 bug 后再合并到 main

### 5. 不要直接在 main 分支开发
- main 只用于接收 development 的完成功能
- 所有开发都在 development 分支进行

## 🆘 常见问题

### Q1: 忘记在哪个分支了？
```bash
git branch
# 当前分支前面有 * 号
```

### Q2: 提交了但忘记推送？
```bash
git push origin development
```

### Q3: Cherry-pick 时遇到冲突？
```bash
# 1. 手动解决冲突文件
# 2. 添加解决后的文件
git add .
# 3. 继续 cherry-pick
git cherry-pick --continue
```

### Q4: 想取消 cherry-pick？
```bash
git cherry-pick --abort
```

### Q5: 想查看 main 和 development 的差异？
```bash
git checkout main
git log development --oneline --not main
# 显示 development 有但 main 没有的提交
```

## 📝 工作流程总结

```
开始开发
   ↓
切换到 development 分支
   ↓
编写代码
   ↓
git add . → git commit → git push
   ↓
功能完成并测试 OK？
   ↓ 是
切换到 main 分支
   ↓
git cherry-pick 选择提交
   ↓
git push origin main
   ↓
切回 development 继续开发
```

## 🎯 快速操作模板

**每天开发：**
```bash
git checkout development
# ...写代码...
git add .
git commit -m "YYYYMMDD 开发内容"
git push origin development
```

**功能完成合并到 main：**
```bash
git checkout development
git log --oneline           # 找到 commit-id
git checkout main
git cherry-pick <commit-id>
git push origin main
git checkout development
```

---

**记住：development 随便折腾，main 保持干净！**

**每天必推送 development 分支到 GitHub！**

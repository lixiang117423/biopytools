#!/usr/bin/env python3
"""
GitHub Releases 清理脚本

功能：
- 自动删除超过指定天数的旧 Release
- 保留指定数量的最新版本
- 区分正式版本和测试版本（alpha/beta/rc）
- 输出详细的清理日志和存储空间估算
"""

import os
import sys
import logging
from datetime import datetime, timedelta
from typing import List, Dict, Optional
from dataclasses import dataclass

# 尝试导入 PyGithub，如果没有则安装
try:
    from github import Github, GithubException
except ImportError:
    print("正在安装 PyGithub...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "PyGithub"])
    from github import Github, GithubException


@dataclass
class ReleaseInfo:
    """Release 信息"""
    id: int
    tag_name: str
    name: str
    created_at: datetime
    prerelease: bool
    draft: bool
    assets: List[Dict]
    total_size: int

    def is_prerelease(self) -> bool:
        """判断是否为测试版本"""
        prerelease_keywords = ['alpha', 'beta', 'rc', 'pre', 'dev']
        tag_lower = self.tag_name.lower()
        return self.prerelease or any(keyword in tag_lower for keyword in prerelease_keywords)


class ReleaseCleaner:
    """Release 清理器"""

    def __init__(self):
        self.gh_token = os.getenv('GITHUB_TOKEN')
        self.repo_name = os.getenv('GITHUB_REPOSITORY')
        self.keep_releases = int(os.getenv('KEEP_RELEASES', '10'))
        self.keep_prereleases = int(os.getenv('KEEP_PRERELEASES', '5'))
        self.delete_after_days = int(os.getenv('DELETE_AFTER_DAYS', '30'))
        self.dry_run = os.getenv('DRY_RUN', 'false').lower() == 'true'

        # 设置日志
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('cleanup_log.txt', mode='w', encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

        # 统计信息
        self.deleted_releases = []
        self.total_freed_space = 0

    def run(self):
        """执行清理任务"""
        self.logger.info("=" * 60)
        self.logger.info("GitHub Releases 清理任务开始")
        self.logger.info("=" * 60)

        if self.dry_run:
            self.logger.warning("⚠️  DRY RUN 模式：不会实际删除任何内容")

        try:
            # 初始化 GitHub 客户端
            g = Github(self.gh_token)
            repo = g.get_repo(self.repo_name)

            # 获取所有 releases
            self.logger.info(f"正在获取 {self.repo_name} 的所有 releases...")
            releases = []
            for release in repo.get_releases():
                # 计算附件总大小
                total_size = sum(asset.size for asset in release.get_assets())

                release_info = ReleaseInfo(
                    id=release.id,
                    tag_name=release.tag_name,
                    name=release.name or release.tag_name,
                    created_at=release.created_at,
                    prerelease=release.prerelease,
                    draft=release.draft,
                    assets=[{
                        'name': asset.name,
                        'size': asset.size,
                        'download_count': asset.download_count
                    } for asset in release.get_assets()],
                    total_size=total_size
                )
                releases.append(release_info)

            if not releases:
                self.logger.info("未找到任何 releases")
                return

            # 分类 releases
            stable_releases = [r for r in releases if not r.is_prerelease()]
            prereleases = [r for r in releases if r.is_prerelease()]

            # 按创建时间排序（新的在前）
            stable_releases.sort(key=lambda x: x.created_at, reverse=True)
            prereleases.sort(key=lambda x: x.created_at, reverse=True)

            self.logger.info(f"\n📊 当前 Release 统计：")
            self.logger.info(f"  - 正式版本：{len(stable_releases)} 个")
            self.logger.info(f"  - 测试版本：{len(prereleases)} 个")
            self.logger.info(f"  - 总计：{len(releases)} 个")

            # 显示当前配置
            self.logger.info(f"\n⚙️  清理配置：")
            self.logger.info(f"  - 保留正式版本：{self.keep_releases} 个")
            self.logger.info(f"  - 保留测试版本：{self.keep_prereleases} 个")
            self.logger.info(f"  - 删除超过：{self.delete_after_days} 天")
            self.logger.info(f"  - Dry Run：{self.dry_run}")

            # 执行清理
            self._cleanup_releases(repo, stable_releases, is_prerelease=False)
            self._cleanup_releases(repo, prereleases, is_prerelease=True)

            # 输出清理报告
            self._print_summary()

            self.logger.info("\n✅ 清理任务完成！")

        except GithubException as e:
            self.logger.error(f"GitHub API 错误: {e}")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"未知错误: {e}")
            sys.exit(1)

    def _cleanup_releases(self, repo, releases: List[ReleaseInfo], is_prerelease: bool):
        """清理指定类型的 releases"""
        release_type = "测试版本" if is_prerelease else "正式版本"
        keep_count = self.keep_prereleases if is_prerelease else self.keep_releases

        self.logger.info(f"\n{'=' * 60}")
        self.logger.info(f"正在清理 {release_type}（保留最新 {keep_count} 个）")
        self.logger.info(f"{'=' * 60}")

        if len(releases) <= keep_count:
            self.logger.info(f"当前 {release_type} 数量（{len(releases)}）未超过保留数量（{keep_count}），无需清理")
            return

        # 确定需要删除的 releases
        cutoff_date = datetime.now() - timedelta(days=self.delete_after_days)
        to_delete = []

        # 跳过需要保留的最新版本
        for i, release in enumerate(releases[keep_count:], start=keep_count):
            if release.created_at < cutoff_date:
                to_delete.append(release)
            else:
                # 即使超过保留数量，如果未超过天数限制，也暂时保留
                age_days = (datetime.now() - release.created_at).days
                self.logger.info(
                    f"⏭️  保留 {release.tag_name} "
                    f"(创建于 {age_days} 天前，未超过 {self.delete_after_days} 天限制)"
                )

        if not to_delete:
            self.logger.info(f"没有需要删除的 {release_type}")
            return

        self.logger.warning(f"\n🗑️  计划删除 {len(to_delete)} 个 {release_type}：")

        # 执行删除
        for release in to_delete:
            age_days = (datetime.now() - release.created_at).days
            size_mb = release.total_size / (1024 * 1024)

            self.logger.warning(
                f"\n  - {release.tag_name}"
                f"\n    名称: {release.name}"
                f"\n    创建时间: {release.created_at.strftime('%Y-%m-%d %H:%M:%S')} ({age_days} 天前)"
                f"\n    附件数量: {len(release.assets)}"
                f"\n    总大小: {size_mb:.2f} MB"
                f"\n    附件列表:"
            )

            for asset in release.assets:
                asset_size_kb = asset['size'] / 1024
                self.logger.warning(
                    f"      • {asset['name']} "
                    f"({asset_size_kb:.2f} KB, "
                    f"下载 {asset['download_count']} 次)"
                )

            # 实际删除
            if not self.dry_run:
                try:
                    gh_release = repo.get_release(release.id)
                    gh_release.delete()
                    self.logger.info(f"    ✅ 已删除")
                except GithubException as e:
                    self.logger.error(f"    ❌ 删除失败: {e}")
                    continue
            else:
                self.logger.info(f"    🔍 [DRY RUN] 将删除")

            # 记录统计信息
            self.deleted_releases.append({
                'tag': release.tag_name,
                'type': release_type,
                'created_at': release.created_at,
                'size': release.total_size
            })
            self.total_freed_space += release.total_size

    def _print_summary(self):
        """打印清理摘要"""
        self.logger.info(f"\n{'=' * 60}")
        self.logger.info("📋 清理摘要")
        self.logger.info(f"{'=' * 60}")

        if not self.deleted_releases:
            self.logger.info("没有删除任何 Release")
            return

        # 按类型分组
        stable_deleted = [r for r in self.deleted_releases if r['type'] == '正式版本']
        prerelease_deleted = [r for r in self.deleted_releases if r['type'] == '测试版本']

        self.logger.info(f"\n删除的正式版本：{len(stable_deleted)} 个")
        for r in stable_deleted:
            size_mb = r['size'] / (1024 * 1024)
            self.logger.info(f"  - {r['tag']} ({size_mb:.2f} MB)")

        self.logger.info(f"\n删除的测试版本：{len(prerelease_deleted)} 个")
        for r in prerelease_deleted:
            size_mb = r['size'] / (1024 * 1024)
            self.logger.info(f"  - {r['tag']} ({size_mb:.2f} MB)")

        # 存储空间统计
        total_mb = self.total_freed_space / (1024 * 1024)
        total_gb = total_mb / 1024

        self.logger.info(f"\n💾 释放的存储空间：")
        self.logger.info(f"  - {self.total_freed_space:,} bytes")
        self.logger.info(f"  - {total_mb:.2f} MB")
        self.logger.info(f"  - {total_gb:.2f} GB")


if __name__ == '__main__':
    cleaner = ReleaseCleaner()
    cleaner.run()

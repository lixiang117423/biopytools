#!/usr/bin/env python3
"""
版本更新脚本 | Version Bump Script
用于自动更新项目版本号 | Automatically update project version numbers

用法 | Usage:
    python scripts/bump_version.py [major|minor|patch] [options]

示例 | Examples:
    python scripts/bump_version.py patch              # 1.0.0 -> 1.0.1
    python scripts/bump_version.py minor              # 1.0.0 -> 1.1.0
    python scripts/bump_version.py major              # 1.0.0 -> 2.0.0
    python scripts/bump_version.py patch --dry-run    # 模拟运行
    python scripts/bump_version.py minor --status beta # 设置为beta版本
"""

import re
import sys
import argparse
import subprocess
from pathlib import Path
from typing import Tuple, Optional, List

class VersionBumper:
    """版本更新器"""
    
    def __init__(self, project_root: Path = None):
        self.project_root = project_root or Path(__file__).parent.parent
        self.version_file = self.project_root / "biopytools" / "_version.py"
        
    def get_current_version(self) -> str:
        """获取当前版本号"""
        if not self.version_file.exists():
            raise FileNotFoundError(f"版本文件不存在: {self.version_file}")
        
        content = self.version_file.read_text(encoding='utf-8')
        version_match = re.search(r'__version__ = ["\']([^"\']*)["\']', content)
        
        if not version_match:
            raise ValueError("无法解析当前版本号")
        
        return version_match.group(1)
    
    def parse_version(self, version: str) -> Tuple[int, int, int, Optional[str]]:
        """解析版本号"""
        # 匹配 1.2.3-alpha.1 格式
        pattern = r'^(\d+)\.(\d+)\.(\d+)(?:-([a-zA-Z]+(?:\.\d+)?))?$'
        match = re.match(pattern, version)
        
        if not match:
            raise ValueError(f"无效的版本格式: {version}")
        
        major, minor, patch, prerelease = match.groups()
        return int(major), int(minor), int(patch), prerelease
    
    def bump_version(self, current_version: str, bump_type: str, 
                    status: Optional[str] = None) -> str:
        """更新版本号"""
        major, minor, patch, _ = self.parse_version(current_version)
        
        if bump_type == 'major':
            new_version = f"{major + 1}.0.0"
        elif bump_type == 'minor':
            new_version = f"{major}.{minor + 1}.0"
        elif bump_type == 'patch':
            new_version = f"{major}.{minor}.{patch + 1}"
        else:
            raise ValueError(f"无效的更新类型: {bump_type}")
        
        # 添加状态标识
        if status and status != 'stable':
            new_version += f"-{status}"
        
        return new_version
    
    def update_version_file(self, new_version: str, status: str = "stable") -> None:
        """更新版本文件"""
        content = self.version_file.read_text(encoding='utf-8')
        
        # 解析纯版本号（去除预发布标识）
        clean_version = new_version.split('-')[0]
        major, minor, patch, _ = self.parse_version(clean_version)
        
        # 更新 __version__
        content = re.sub(
            r'__version__ = ["\'][^"\']*["\']',
            f'__version__ = "{new_version}"',
            content
        )
        
        # 更新 __version_info__
        content = re.sub(
            r'__version_info__ = \([^)]*\)',
            f'__version_info__ = ({major}, {minor}, {patch})',
            content
        )
        
        # 更新 VERSION_STATUS
        content = re.sub(
            r'VERSION_STATUS = ["\'][^"\']*["\']',
            f'VERSION_STATUS = "{status}"',
            content
        )
        
        self.version_file.write_text(content, encoding='utf-8')
    
    def update_pyproject_toml(self, new_version: str) -> None:
        """更新 pyproject.toml 中的版本号（如果使用静态版本）"""
        pyproject_file = self.project_root / "pyproject.toml"
        
        if not pyproject_file.exists():
            return
        
        content = pyproject_file.read_text(encoding='utf-8')
        
        # 检查是否使用动态版本
        if 'dynamic = ["version"]' in content:
            print("检测到动态版本配置，跳过 pyproject.toml 更新")
            return
        
        # 更新静态版本
        clean_version = new_version.split('-')[0]  # 移除预发布标识
        content = re.sub(
            r'version = ["\'][^"\']*["\']',
            f'version = "{clean_version}"',
            content
        )
        
        pyproject_file.write_text(content, encoding='utf-8')
        print(f"已更新 pyproject.toml 版本为: {clean_version}")
    
    def update_init_file(self) -> None:
        """更新 __init__.py 文件中的导入"""
        init_file = self.project_root / "biopytools" / "__init__.py"
        
        if not init_file.exists():
            return
        
        content = init_file.read_text(encoding='utf-8')
        
        # 确保有正确的导入语句
        import_line = "from ._version import __version__, __version_info__, get_version_string, get_build_info"
        
        if import_line not in content:
            # 在适当位置添加导入
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if line.startswith('"""') and i > 0:
                    # 在文档字符串后添加导入
                    end_docstring = -1
                    for j in range(i + 1, len(lines)):
                        if lines[j].strip().endswith('"""'):
                            end_docstring = j
                            break
                    
                    if end_docstring > 0:
                        lines.insert(end_docstring + 2, import_line)
                        break
            
            init_file.write_text('\n'.join(lines), encoding='utf-8')
            print("已更新 __init__.py 导入语句")
    
    def validate_version(self, version: str) -> bool:
        """验证版本号格式"""
        try:
            self.parse_version(version.split('-')[0])  # 验证基础版本号
            return True
        except ValueError:
            return False
    
    def get_git_status(self) -> Tuple[bool, List[str]]:
        """获取Git状态"""
        try:
            # 检查是否有未提交的更改
            result = subprocess.run(
                ['git', 'status', '--porcelain'],
                capture_output=True,
                text=True,
                cwd=self.project_root
            )
            
            if result.returncode != 0:
                return False, ["Git 命令执行失败"]
            
            changes = result.stdout.strip().split('\n') if result.stdout.strip() else []
            clean = len(changes) == 0 or (len(changes) == 1 and not changes[0])
            
            return clean, changes
        
        except FileNotFoundError:
            return False, ["Git 未安装或不在 PATH 中"]

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='版本更新脚本 | Version Bump Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s patch              # 补丁版本更新 1.0.0 -> 1.0.1
  %(prog)s minor              # 次要版本更新 1.0.0 -> 1.1.0
  %(prog)s major              # 主要版本更新 1.0.0 -> 2.0.0
  %(prog)s patch --dry-run    # 模拟运行，不实际更改
  %(prog)s minor --status beta # 创建beta版本
        """
    )
    
    parser.add_argument(
        'bump_type',
        choices=['major', 'minor', 'patch'],
        help='版本更新类型 | Version bump type'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='模拟运行，不实际更改文件 | Dry run without making changes'
    )
    
    parser.add_argument(
        '--status',
        choices=['alpha', 'beta', 'rc', 'stable'],
        default='stable',
        help='版本状态 | Version status (default: stable)'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='强制更新，忽略Git状态检查 | Force update, ignore Git status'
    )
    
    parser.add_argument(
        '--check-git',
        action='store_true',
        default=True,
        help='检查Git工作目录状态 | Check Git working directory status'
    )
    
    args = parser.parse_args()
    
    try:
        bumper = VersionBumper()
        
        # 检查Git状态
        if args.check_git and not args.force:
            is_clean, changes = bumper.get_git_status()
            if not is_clean:
                print("警告：工作目录有未提交的更改:")
                for change in changes[:10]:  # 只显示前10个
                    if change.strip():
                        print(f"  {change}")
                if len(changes) > 10:
                    print(f"  ... 还有 {len(changes) - 10} 个更改")
                
                if not args.force:
                    response = input("是否继续？(y/N): ")
                    if response.lower() not in ['y', 'yes']:
                        print("操作已取消")
                        sys.exit(1)
        
        # 获取当前版本
        try:
            current_version = bumper.get_current_version()
        except (FileNotFoundError, ValueError) as e:
            print(f"错误：{e}")
            sys.exit(1)
        
        # 计算新版本
        try:
            new_version = bumper.bump_version(current_version, args.bump_type, 
                                            args.status if args.status != 'stable' else None)
        except ValueError as e:
            print(f"错误：{e}")
            sys.exit(1)
        
        # 显示版本变更
        print(f"当前版本: {current_version}")
        print(f"新版本:   {new_version}")
        print(f"更新类型: {args.bump_type}")
        print(f"版本状态: {args.status}")
        
        if args.dry_run:
            print("\n模拟运行完成（未实际更改文件）")
            return
        
        # 确认更新
        if not args.force:
            response = input("\n确认更新版本？(y/N): ")
            if response.lower() not in ['y', 'yes']:
                print("操作已取消")
                sys.exit(0)
        
        # 执行更新
        print("\n开始更新版本文件...")
        
        bumper.update_version_file(new_version, args.status)
        print(f"✓ 已更新 {bumper.version_file.relative_to(bumper.project_root)}")
        
        bumper.update_pyproject_toml(new_version)
        bumper.update_init_file()
        
        print(f"\n版本已成功更新到 {new_version}")
        print("\n下一步操作:")
        print("1. 更新 CHANGELOG.md")
        print("2. 运行测试确保一切正常")
        print("3. 提交更改并推送到远程: git add . && git commit -m 'bump version to {}'".format(new_version), " && git push origin main")
        
        # if args.status == 'stable':
        #     print("4. 创建标签: git tag v{}".format(new_version))
        #     print("5. 推送到远程: git push origin main --tags")
        
    except KeyboardInterrupt:
        print("\n操作被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"未预期的错误：{e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
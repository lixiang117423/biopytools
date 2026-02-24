import os
import json
import subprocess
import sys
from datetime import datetime
import shutil  # <-- 新增导入 shutil 模块

# 主备份目录
BASE_OUTPUT_DIR = "conda_env_backups"

def get_current_date_str():
    """获取 'YYYY-MM-DD' 格式的当前日期字符串。"""
    return datetime.now().strftime('%Y-%m-%d')

def get_conda_envs():
    """获取所有 conda 环境的名称和路径。"""
    try:
        result = subprocess.run(
            ['conda', 'env', 'list', '--json'],
            capture_output=True, text=True, check=True
        )
        data = json.loads(result.stdout)
        return data.get('envs', [])
    except (subprocess.CalledProcessError, FileNotFoundError, json.JSONDecodeError) as e:
        print(f"❌ 错误：无法获取 Conda 环境列表。请确保 Conda 已正确安装并配置在 PATH 中。")
        print(f"具体错误: {e}")
        sys.exit(1)

def export_env(env_path, output_dir):
    """导出一个指定路径的 conda 环境。"""
    if not env_path:
        return
    env_name = os.path.basename(env_path)
    if 'envs' not in os.path.normpath(env_path).split(os.sep):
        print(f"ℹ️  跳过 'base' 环境 ({env_name})...")
        return

    print(f"▶️  正在导出环境: {env_name} ...")
    output_file = os.path.join(output_dir, f"{env_name}.yml")
    try:
        command = ['conda', 'env', 'export', '-n', env_name, '--no-builds']
        with open(output_file, 'w') as f:
            subprocess.run(
                command, check=True, stdout=f, stderr=subprocess.PIPE, text=True
            )
        print(f"✅ 成功导出到 {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ 导出环境 '{env_name}' 失败！")
        print(f"   错误信息: {e.stderr.strip()}")

# def backup_zshrc(destination_dir):
#     """备份 ~/.zshrc 文件到指定目录。"""
#     # os.path.expanduser('~') 会自动将 '~' 转换成你的主目录路径
#     zshrc_path = os.path.join(os.path.expanduser('~'), '.zshrc')
    
#     print("-------------------------------------")
#     print(f"▶️  正在尝试备份 .zshrc 文件...")

#     if os.path.exists(zshrc_path):
#         try:
#             # shutil.copy2 会保留文件的元数据（如修改时间）
#             shutil.copy2(zshrc_path, destination_dir)
#             print(f"✅ 成功备份 .zshrc 到 {destination_dir}")
#         except Exception as e:
#             print(f"❌ 备份 .zshrc 失败！")
#             print(f"   错误信息: {e}")
#     else:
#         print("ℹ️  未找到 ~/.zshrc 文件，跳过备份。")

def backup_zshrc(destination_dir):
    """备份 ~/.zshrc 文件到指定目录。"""
    # os.path.expanduser('~') 会自动将 '~' 转换成你的主目录路径
    zshrc_path = os.path.join(os.path.expanduser('~'), '.zshrc')
    
    print("-------------------------------------")
    print(f"▶️  正在尝试备份 .zshrc 文件...")

    if os.path.exists(zshrc_path):
        try:
            # 构建目标文件路径，重命名为 zshrc（去掉前面的点）
            destination_file = os.path.join(destination_dir, 'zshrc')
            
            # shutil.copy2 会保留文件的元数据（如修改时间）
            shutil.copy2(zshrc_path, destination_file)
            print(f"✅ 成功备份 .zshrc 到 {destination_file}")
        except Exception as e:
            print(f"❌ 备份 .zshrc 失败！")
            print(f"   错误信息: {e}")
    else:
        print("ℹ️  未找到 ~/.zshrc 文件，跳过备份。")

def main():
    """主函数"""
    today_str = get_current_date_str()
    dated_output_dir = os.path.join(BASE_OUTPUT_DIR, today_str)

    if not os.path.exists(dated_output_dir):
        os.makedirs(dated_output_dir)
        print(f"📂 创建归档目录: {dated_output_dir}")

    env_paths = get_conda_envs()
    if env_paths:
        print(f"\n发现 {len(env_paths)} 个环境。开始导出...")
        print("-------------------------------------")
        for path in env_paths:
            export_env(path, dated_output_dir)
            print("-------------------------------------")
    else:
        print("未找到任何 Conda 环境。")
    
    # 在所有环境导出完成后，执行 .zshrc 的备份
    backup_zshrc(dated_output_dir)

    print("\n🎉 所有任务完成！")
    print(f"备份文件已保存在 '{dated_output_dir}' 目录中。")

if __name__ == "__main__":
    main()

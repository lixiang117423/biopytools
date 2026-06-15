import sys
import os


def resource_path(relative_path, parent_path=False):
    '''返回资源绝对路径。'''
    if hasattr(sys, '_MEIPASS'):
        # PyInstaller会创建临时文件夹temp
        # 并把路径存储在_MEIPASS中
        base_path = sys._MEIPASS
    else:
        base_path = os.path.abspath('.')
    if not parent_path:
        return os.path.join(base_path, relative_path)
    else:
        return os.path.abspath(os.path.join(base_path, relative_path, ".."))


def make_dirs():
    dirs = ["Excel_Files", "None_Pretreated_Files", "Pretreated_Files"]
    root_path = resource_path("", True)
    for dir in dirs:
        if not os.path.exists(os.path.join(root_path, dir)):
            os.mkdir(os.path.join(root_path, dir))

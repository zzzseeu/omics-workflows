"""Access packaged workflow resources without depending on a Git checkout."""

from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path


@contextmanager
def package_directory():
    """以上下文方式提供框架包根目录。

    返回安装资源的Path；如资源需要临时解包，保持上下文直到分析进程结束再清理。
    """
    with as_file(files("omics_workflows")) as directory:
        yield directory


def source_boundary(directory):
    """根据包目录返回代码保护边界Path：开发安装保护整个源码仓库，wheel安装保护包目录。"""
    directory = Path(directory).resolve()
    if directory.parent.name == "src" and (directory.parent.parent / "pyproject.toml").is_file():
        return directory.parent.parent
    return directory


def workflow_defaults(name):
    """根据流程名称读取包内defaults.yaml，返回UTF-8字符串，不依赖当前工作目录。"""
    return files(f"omics_workflows.workflows.{name}").joinpath("defaults.yaml").read_text(encoding="utf-8")

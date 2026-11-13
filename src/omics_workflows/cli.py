"""Unified command for installed omics analysis workflows."""

import importlib
import sys

from omics_workflows import __version__
from omics_workflows.common.arguments import argument_parser

# New built-in workflows register their entry module here; imports stay lazy.
WORKFLOWS = {
    "assembly": ("omics_workflows.workflows.assembly.runner", "基因组组装与质量评估"),
}


def main(argv=None, prog=None):
    """解析统一命令并延迟导入选中的流程，返回模块退出码；版本、列表和帮助不导入分析依赖。"""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argument_parser(
        description="统一生信分析流程命令；选择流程后可通过该流程的 --help 查看参数。",
        prog=prog or "omics-workflows",
        epilog=(
            "用法示例：\n"
            "  omics-workflows list\n"
            "  omics-workflows assembly --help\n"
            "  omics-workflows --version"
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="类型：开关，无需赋值。显示omics-workflows版本后退出。",
    )
    parser.add_argument(
        "workflow",
        choices=[*WORKFLOWS, "list"],
        help="类型：字符串枚举；必填。assembly运行基因组组装；list列出已注册流程。流程参数通过对应子命令的 --help 查看。",
    )
    if not argv:
        parser.print_help()
        return 0
    if argv[0].startswith("-"):
        parser.parse_args(argv)
        return 0
    selected = parser.parse_args(argv[:1])
    if selected.workflow == "list":
        parser.parse_args(argv)
        for name, (_, description) in WORKFLOWS.items():
            print(f"{name}\t{description}")
        return 0
    module = importlib.import_module(WORKFLOWS[selected.workflow][0])
    return module.main(argv[1:], prog=f"{parser.prog} {selected.workflow}")

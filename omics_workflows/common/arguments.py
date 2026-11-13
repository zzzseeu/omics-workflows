"""Shared help presentation for workflow commands and standalone analysis steps."""

import argparse


def argument_parser(description, prog=None, epilog=None):
    """建立统一的中文命令行解析器。

    输入为命令说明、可选程序名和补充说明；返回带有显式--help开关的ArgumentParser。
    调用者逐项定义参数类型、用途和默认行为；本函数不改变各参数的解析规则。
    """
    parser = argparse.ArgumentParser(
        description=description,
        prog=prog,
        epilog=epilog,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="类型：开关，无需赋值。显示参数含义、类型和默认行为后退出。",
    )
    return parser

#!/usr/bin/env python3
"""`python -m locuskit` —— 启动网页 GUI（打包成可执行文件时也用这个入口）。

    python -m locuskit                     # 起服务并打开浏览器
    python -m locuskit --port 9000         # 换端口
    python -m locuskit --root /path/work   # 指定工作目录（含 code/ 与 data/）
    python -m locuskit --no-browser        # 不自动开浏览器
"""
from __future__ import annotations

import argparse
import os
import sys
import threading
import webbrowser


def default_root():
    """打包后取可执行文件所在目录；源码运行时取包的上层目录。"""
    if getattr(sys, "frozen", False):
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="locuskit", description="locuskit GUI")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8088)
    ap.add_argument("--root", default=None,
                    help="工作目录（含 code/ 与 data/）；默认取可执行文件/包所在目录")
    ap.add_argument("--no-browser", action="store_true", help="不自动打开浏览器")
    a = ap.parse_args(argv)

    root = os.path.abspath(a.root or default_root())
    os.chdir(root)

    # 绝对导入：打包成可执行文件时本文件被当顶层脚本运行，相对导入会失败
    try:
        from locuskit import gui
    except ImportError:                      # 作为包内模块运行时
        from . import gui
    gui.set_root(root)

    url = f"http://{a.host}:{a.port}/"
    print(f"locuskit GUI\n  工作目录: {root}\n  {url}\n（Ctrl+C 退出）", flush=True)
    if not a.no_browser:
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    gui.APP.run(host=a.host, port=a.port, threaded=True, debug=False)


if __name__ == "__main__":
    main()

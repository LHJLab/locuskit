#!/usr/bin/env python3
"""打包成单文件可执行程序（Linux/macOS → dist/locuskit-gui，Windows → .exe）。

用法（先激活一个**干净**的 venv，装好依赖）：

    python -m pip install numpy matplotlib flask pyyaml pyinstaller
    python build_exe.py                # 产物在 dist/
    python build_exe.py --debug         # 保留控制台输出与控制台窗口
    python build_exe.py --dist out/     # 换输出目录

为什么单独写这个脚本：本仓库的目录名就是包名，`import locuskit` 需要把**包的父目录**
放进 `--paths`，入口是包内的 `__main__.py`。手动敲 PyInstaller 很容易把这两个路径写错。
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent          # 仓库根（也是包目录）
PKG_PARENT = HERE.parent                        # 让 `import locuskit` 成立

EXCLUDES = [
    # 环境里常见的重依赖与 GUI 后端，本工具用不到
    "tkinter", "pandas", "scipy", "IPython", "pytest", "notebook",
    "PyQt5", "PyQt6", "PySide2", "PySide6", "shiboken2", "shiboken6",
    "matplotlib.backends.backend_qt", "matplotlib.backends.backend_qt5agg",
    "matplotlib.backends.backend_webagg",
]


def check_env():
    problems = []
    try:
        import PyInstaller  # noqa: F401
    except ImportError:
        problems.append("缺少 PyInstaller：python -m pip install pyinstaller")
    for mod in ("numpy", "matplotlib", "flask", "yaml"):
        try:
            __import__(mod)
        except ImportError:
            problems.append(f"缺少 {mod}：python -m pip install {mod}")
    if problems:
        print("环境不完整：", file=sys.stderr)
        for p in problems:
            print("  -", p, file=sys.stderr)
        raise SystemExit(1)

    # conda 的 numpy 外链 MKL，PyInstaller 收集不到 → 运行时报 libmkl_core.so 缺失
    try:
        import numpy
        np_dir = Path(numpy.__file__).resolve().parent
        if not (np_dir.parent / "numpy.libs").is_dir():
            print("警告：numpy 没有自带 numpy.libs（可能来自 conda）。"
                  "打包后若报 libmkl_core.so 缺失，请改用干净 venv 重装 numpy。",
                  file=sys.stderr)
    except Exception:
        pass


def main(argv=None):
    ap = argparse.ArgumentParser(description="打包 locuskit 为单文件可执行程序")
    ap.add_argument("--dist", default=str(HERE / "dist"), help="输出目录（默认 dist/）")
    ap.add_argument("--work", default=str(HERE / "build"), help="中间目录（默认 build/）")
    ap.add_argument("--name", default="locuskit-gui", help="可执行文件名")
    ap.add_argument("--debug", action="store_true", help="保留控制台窗口与调试输出")
    ap.add_argument("--hooks", default=str(HERE / "_pyinstaller_hooks"),
                    help="额外 hook 目录（内置 hook 加载失败时用）")
    ap.add_argument("--no-hooks", action="store_true", help="不使用额外 hook")
    a = ap.parse_args(argv)

    check_env()

    args = [
        "--onefile",
        "--name", a.name,
        "--clean",
        "--noconfirm",
        "--paths", str(PKG_PARENT),
        "--collect-submodules", "locuskit",
        "--distpath", a.dist,
        "--workpath", a.work,
        "--specpath", a.work,
    ]
    if not a.debug:
        args += ["--windowed"] if os.name == "nt" else []
    if not a.no_hooks and Path(a.hooks).is_dir():
        args += ["--additional-hooks-dir", a.hooks]
    for mod in EXCLUDES:
        args += ["--exclude-module", mod]
    args.append(str(HERE / "__main__.py"))

    print("[1/2] PyInstaller 参数：")
    print("      " + " ".join(args))
    from PyInstaller.__main__ import run as pyi_run
    pyi_run(args)

    exe = Path(a.dist) / (a.name + (".exe" if os.name == "nt" else ""))
    print("\n[2/2] 完成：" + str(exe))
    if exe.exists():
        print(f"      大小 {exe.stat().st_size / 1048576:.0f} MB")
        print("\n自检：")
        print(f"  {exe} --help")
        print(f"  {exe} --no-browser --port 8088 &")
        print("  curl -s -o /dev/null -w '%{http_code}\\n' http://127.0.0.1:8088/")
        print("\n分发：把可执行文件放到含 code/ 与 data/ 的目录后双击即可。")
    else:
        print("      未找到产物，请检查上面的 PyInstaller 输出。", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

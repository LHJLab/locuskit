PyInstaller 额外 hook（仅在打包时报内置 hook 加载失败时需要）

用法：
    --additional-hooks-dir locuskit/_pyinstaller_hooks

- hook-numpy.py       显式收集 numpy.libs 下的 OpenBLAS/gfortran 动态库
- hook-matplotlib.py  收集 mpl-data（字体与样式表），后端限定 Agg/PDF/SVG

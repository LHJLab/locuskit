"""最小 matplotlib hook：字体与样式表必须带上，后端用 Agg。"""
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = collect_data_files('matplotlib')
hiddenimports = ['matplotlib.backends.backend_agg',
                 'matplotlib.backends.backend_pdf',
                 'matplotlib.backends.backend_svg']
try:
    hiddenimports += collect_submodules('matplotlib.backends')
except Exception:
    pass

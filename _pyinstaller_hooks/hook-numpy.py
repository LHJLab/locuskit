"""numpy hook：显式收集 numpy.libs 下的 OpenBLAS/gfortran 动态库。

anaconda 的 numpy 外链 MKL（PyInstaller 收集不到，运行时报 libmkl_core.so 缺失），
PyPI 的 wheel 把 OpenBLAS 放在 site-packages/numpy.libs/，这里确保一起打包。
"""
import glob
import os

from PyInstaller.utils.hooks import collect_dynamic_libs, collect_submodules, get_package_paths

hiddenimports = ['numpy.core._multiarray_umath', 'numpy.linalg._umath_linalg',
                 'numpy.fft._pocketfft_internal']
try:
    hiddenimports += collect_submodules('numpy.core')
except Exception:
    pass

binaries = collect_dynamic_libs('numpy') or []
try:
    _pkg, _ = get_package_paths('numpy')
    _libs = os.path.join(os.path.dirname(_pkg), 'numpy.libs')
    if os.path.isdir(_libs):
        binaries += [(p, '.' if 'openblas' in os.path.basename(p).lower() else 'numpy.libs')
                     for p in glob.glob(os.path.join(_libs, '*.so*'))]
except Exception:
    pass

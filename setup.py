from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
cymodule = 'breakfit_cython_ext'

setup(
  name='breakfit_cython',
  ext_modules=[Extension(cymodule, [cymodule + '.pyx'], libraries=['m'],compiler_directives={'linetrace': True}, extra_compile_args = ['-ffast-math'], define_macros=[('CYTHON_TRACE', '1')])],
  cmdclass={'build_ext': build_ext},
)

# cython: profile=True
# cython: linetrace=True

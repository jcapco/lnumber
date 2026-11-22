#for linux, macOS and windows builds
#the library is a plain shared library loaded with ctypes (see lnumber/__init__.py),
#not a CPython extension module; compat/ replaces GMP and the tr1 headers
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from sys import version_info, platform
import os
try:
  from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:
  from wheel.bdist_wheel import bdist_wheel

def rem_opt(l,keys):
  return [y for y in l if y not in keys]
  
class BuildExt(build_ext):
  def build_extensions(self):
    if self.compiler.compiler_type == "msvc":
      #link the C/C++ runtime statically (/MT instead of /MD)
      if not getattr(self.compiler, "initialized", True):
        self.compiler.initialize(self.plat_name)
      for opts in ("compile_options", "compile_options_debug"):
        if hasattr(self.compiler, opts):
          setattr(self.compiler, opts, [{"/MD": "/MT", "/MDd": "/MTd"}.get(y, y) for y in getattr(self.compiler, opts)])
    else:
      keys = ["-Wstrict-prototypes","-fwrapv","-g","-O1","-O2","-fstack-protector-strong","-fno-strict-aliasing"]
      self.compiler.compiler_so = rem_opt(self.compiler.compiler_so, keys)
    build_ext.build_extensions(self)

  #no PyInit_... symbol: the library is loaded with ctypes
  def get_export_symbols(self, ext):
    return []

  #file names expected by lnumber/__init__.py (setuptools adds the package folder)
  def get_ext_filename(self, ext_name):
    return "lib_lnumber.dll" if os.name == "nt" else "lnumber.so"

#the library does not use the Python C API: one wheel per platform for all Python 3 versions
class BdistWheel(bdist_wheel):
  def get_tag(self):
    python, abi, plat = bdist_wheel.get_tag(self)
    return "py3", "none", plat

#python 2.X support
if version_info.major==2: 
  name='lnumber'
  version='0.20'
  packages=['lnumber']
  author = 'Jose Capco'
  #url='https://github.com/jcapco/lnumber' #deprecated in the new pipy
  data_files = [("", ["license.txt"])]
  os.environ["CC"] = "gcc"
  os.environ["CXX"] = "gcc"

if platform == 'win32':
  compile_args = ['/O2','/DNDEBUG','/DLIBLNUMBER_EXPORTS','/openmp','/FIcmath']
  link_args = []
else:
  compile_args = ['-std=c++11','-O3','-DNDEBUG','-DLIBLNUMBER_EXPORTS','-Wall','-Wextra','-Wno-unknown-pragmas','-Wno-sign-compare','-fpic']
  if platform == 'darwin':
    #libomp is provided by the build environment (see the CI workflow)
    compile_args += ['-Xpreprocessor','-fopenmp']
    link_args = ['-lomp']
  else:
    compile_args += ['-flto','-fopenmp']
    link_args = ['-shared','-lstdc++','-lm','-lgomp','-fopenmp']

module = Extension('lnumber.lnumber',
  sources = ['lnumber/src/laman_number.cpp','lnumber/src/lib.cpp'],
  include_dirs = ['compat','lnumber/inc','lnumber'],
  language = 'c++',
  extra_compile_args=compile_args,
  extra_link_args=link_args)

setup(
  name='lnumber',
  version='0.20',
  packages=['lnumber'],
  #python_tag=python_tag,    
  cmdclass={'build_ext': BuildExt, 'bdist_wheel': BdistWheel},
  include_package_data=True,
  ext_modules=[module],
  author = 'Jose Capco',
  #url='https://github.com/jcapco/lnumber',
  data_files = [("", ["license.txt"])]
)

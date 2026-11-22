#for linux builds
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from sys import version_info
import os

def rem_opt(l,keys):
  return [y for y in l if y not in keys]
  
class BuildExt(build_ext):
  def build_extensions(self):
    keys = ["-Wstrict-prototypes","-fwrapv","-g","-O1","-O2","-fstack-protector-strong","-fno-strict-aliasing"]
    self.compiler.compiler_so = rem_opt(self.compiler.compiler_so, keys)
    build_ext.build_extensions(self)

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

module = Extension('lnumber.lnumber',
  sources = ['lnumber/src/laman_number.cpp','lnumber/src/lib.cpp'],
  include_dirs = ['lnumber/inc','lnumber'],
  extra_compile_args=['-std=c++11','-O3','-s','-DNDEBUG','-DLIBLNUMBER_EXPORTS','-flto','-fopenmp','-m64','-Wall','-Wextra','-Wno-unknown-pragmas','-Wno-sign-compare','-fpic'],
  extra_link_args=['-shared','-lstdc++','-lm','-lgmp','-lgmpxx','-lgomp','-fopenmp'])

setup(
  name='lnumber',
  version='0.20',
  packages=['lnumber'],
  #python_tag=python_tag,    
  cmdclass={'build_ext': BuildExt},
  include_package_data=True,
  ext_modules=[module],
  author = 'Jose Capco',
  #url='https://github.com/jcapco/lnumber',
  data_files = [("", ["license.txt"])]
)

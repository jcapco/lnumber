from ctypes import cdll, c_char, c_size_t
from os import path, name
from sys import version_info

if name == "nt":    
  #window specific issue, cannot directly load with absolute path
  if version_info.major==3: 
    lib = cdll.LoadLibrary(path.join(path.dirname(__file__),"./lib_lnumber.dll"))
  else: #nt and py2 specific issue, cannot directly load dll with absolute path
    from os import environ
    dir = path.abspath(path.dirname(__file__))
    if dir not in environ['PATH']:
      environ['PATH'] = dir + ';' + environ['PATH']
    lib = cdll.LoadLibrary("lib_lnumber.dll")
else:
  dir = path.abspath(path.dirname(__file__))
  
  if version_info.major==2: 
    lib = cdll.LoadLibrary(dir+"/lnumber.so")
  else: #gcc gives platform suffix in Python 3.x
    from os import listdir
    for file in listdir(dir):
      if file.endswith(".so"): break  
    lib = cdll.LoadLibrary(dir+"/"+file)

#print lnumber(252590061719913632,12)
def lnumber(graph):
  global lib
  return lib.laman_number(str(graph).encode("utf-8"))

#print lnumbers(252590061719913632)
def lnumbers(graph):
  global lib
  return lib.laman_number_spherical(str(graph).encode("utf-8"))

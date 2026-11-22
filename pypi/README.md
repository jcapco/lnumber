This is a python source and wheel for the lnumber module. Although this can be installed using

`pip install lnumber`

In Linux or MacOS:

`pip` will install by compiling from source.
For best performance, the .so file is compiled:

`gcc -c ../src/laman_number.cpp -I../inc -I../ -std=c++11 -O3 -s -DNDEBUG -flto -fopenmp -fpic -m64 -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -fpic -o ../laman_number.o`

`gcc -c ../src/lib.cpp -I../inc -I../ -std=c++11 -O3 -s -DNDEBUG -flto -fopenmp -fpic -m64 -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -fpic -o ../lib.o`

`gcc -shared -lstdc++ -lm -lgmp -lgmpxx -lgomp -o ./lnumber.so ../laman_number.o ../lib.o`

and renaming the `__init\__.py` as `lnumber.py` and placing this and the generated .so files into a searchable python path (e.g. `site-packages` or a path in `PYTHONPATH`) 

In Windows:
`pip` can only install using the binary distribution.
This installation requires Microsoft Visual C++ 2008 Redistributable. If you cannot run in windows due to missing dll, please install the MSVC 2008 redistributable.

More information and documentation about usage, limitations and the C++ program and compilation can be found here:
https://github.com/jcapco/lnumber

For further question please email me (Jose Capco) at jcapco-at-risc-dot-jku-dot-at


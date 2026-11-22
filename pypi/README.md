The `lnumber` package computes the Laman number of a Laman graph, i.e. its number of realizations, both on the plane and on the sphere.

```
pip install lnumber
```

Prebuilt wheels are available for Python 3 on the following 64-bit platforms, so `pip` installs `lnumber` without compiling anything and without any additional dependency:

* Linux (glibc and musl): x86_64, aarch64
* macOS: arm64 (Apple Silicon), x86_64 (Intel)
* Windows: AMD64, ARM64

On other platforms (e.g. 32-bit systems, or Linux on ppc64le, s390x, riscv64), `pip` builds the package from source. This requires a C++11 compiler with OpenMP support; GMP is not needed. On macOS, OpenMP can be installed with Homebrew, which keeps it outside the default search paths:
```
brew install libomp
CPPFLAGS="-I$(brew --prefix libomp)/include" LDFLAGS="-L$(brew --prefix libomp)/lib" pip install lnumber
```

Python 2 and 32-bit Windows are no longer supported; the last release supporting them is 0.20.

**Note for macOS users:** the macOS wheels include their own copy of the OpenMP runtime (`libomp`). If another package that bundles its own `libomp` (e.g. PyTorch) is imported in the same Python process, the program may stop with `OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized`. In that case, run `lnumber` in a separate process.

Example: a graph is encoded as an integer whose bits are the upper triangle (without diagonal) of its adjacency matrix, read row by row.
```
from lnumber import *

# 808 number of planar realizations of the graph 1813573113164
lnumber(1813573113164)

# 1280 number of spherical realizations of the graph 1813573113164
lnumbers(1813573113164)
```

More information and documentation about usage, limitations and the C++ program and compilation can be found here:
https://github.com/jcapco/lnumber

For further question please email me (Jose Capco) at jcapco-at-risc-dot-jku-dot-at

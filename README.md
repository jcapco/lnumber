# Toolkit for Computing the Laman Number
The Laman number library `lnumber` is meant to be used as an external dynamic library by programs such as the Nauty Laman plugin. The library implements a fast algorithm, original executable implemented by Christoph Koutschan, that computes the Laman number of a Laman graph based on our paper (with M. Gallet, G. Grasegger, C. Koutschan, N. Lubbes and J. Schicho) see [Capco et. al. The number of Realization of a Laman Graph](http://www.koutschan.de/data/laman/). 

The `nauty` plugin is a fork of the [Nauty-Laman plugin](https://github.com/martinkjlarsson/nauty-laman-plugin) originally written by Martin Larrson for the utility 
program `geng` provided by [Nauty](http://pallini.di.uniroma1.it/). The plugin counts the number of non-isomorphic Laman graphs of a given number of 
vertices. The quick generation of the Laman graphs was implemented by Martin Larrson. I integrated parallel computing of Laman numbers into the plugin via dependency on the `lnumber` library. Parallel computing is done via `openmp`.

## Installation with compiling

**Disclaimer:** All MSVC project files are for VS2008 and setup for Release. You have to manually setup other build configurations if you want them. The project files should be forward-compatible to newer MSVC. `nauty` and `lnumber` dependencies for `geng` are setup to link dynamically.

**For the `lnumber` library:** 
1. Clone this repo.
1. Build the `lnumber` library. 
    * You can build with gcc like this (in the `lnumber` folder)
    ```makefile
    gcc -c ./src/laman_number.cpp -I./inc -std=c++11 -O3 -s \
      -DNDEBUG -flto -fopenmp -fpic -m64 -Wall -Wextra -Wno-unknown-pragmas \
      -Wno-sign-compare -fpic -o ./laman_number.o

    gcc -c ./src/lib.cpp -I./inc -std=c++11 -O3 -s -DNDEBUG \
      -flto -fopenmp -fpic -m64 -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare \
      -fpic -o ./lib.o
    
    gcc -shared -lstdc++ -lm -lgmp -lgmpxx -lgomp -o ./liblnumber.a ./lib.o -fopenmp 
    ```
    * You can use the VS2008 project file `./vs2008/lnumber.vcproj` to build with MSVC. I recommend to link with the fork of `gmp` called `mpir`. For VS2008, I provided the `mpir` library in the `./vs2008/lib` folder. For other MSVC versions, you have to download and build [`mpir`](https://github.com/wbhart/mpir).

**For the `nauty` plugin:**
1. Build the library `lnumber` as described above
1. Download and save [`nauty`](http://pallini.di.uniroma1.it/) source files in the `./nauty` folder. 
1. Build all objects from the `nauty` makefile (if you are building with gcc). You can use the VS2008 project file, `./vs2008/nauty.vcproj` to build nauty with MSVC. There are no dependencies.
1. `geng` has to be rebuilt 
in C++ with the flag `-D'PLUGIN="<path to this repo>/prunelaman.h"'` and linked with `lnumber`. It is also important to activate `openmp` when building, otherwise the program may still use multiple threads but there will be no speedup for parallel computing and multiple core computers will have roughly the same runtime as a single core one.
    * If you build with gcc, you can build `geng` as follows (in the `geng` folder)
    ```makefile
    g++  -Wno-write-strings -std=c++98 -I./ -o geng -fopenmp -O3 -s \ 
      -mpopcnt -march=native -D'PLUGIN="<path to this repo>/prunelaman.h"' \
      -DMAXN=WORDSIZE -DWORDSIZE=32 geng.c gtoolsW.o nautyW1.o nautilW1.o \
      naugraphW1.o schreier.o naurng.o  -L. -lstdc++ -lm -llnumber \
      -lgmpxx -lgomp -lgmp
    ```
    * If you are using MSVC you can compile and link using the VS2008 project file `./vs2008/geng.vcproj`. You will need to link with `nauty` and `mpir`. For VS2008, these dependencies (`lib` and `dll`) in  the `./vs2008/lib/` and `./bin` folders. 
    
5. Once `geng` is compiled you can run it with the `-M` parameter to compute the maximum Laman numbers (see usage)

## Installation with prebuilt binaries (Windows)

1. The executable files for Windows can be copied from the `./bin` folder
1. Run the file `geng` (see usage).

## Usage
**For the `nauty` plugin:**

The plugin was originally developed by Martin Larrson who added the following parameters to `geng`:
* `-K#`: generate (k,l)-tight graphs where l = k(k+1)/2. Minimum degree and number of edges will default to k and kn-l, respectively. Sparse graphs can be generated by manually providing the minimum and maximum number of edges (e.g. `0:999`). In that case, the minimum degree will default to zero.
* `-L#`: provides the l when generating (k,l)-sparse or (k,l)-tight graphs.
* `-H`: generate (k,l)-tight graphs constructible by [Henneberg type I moves](https://en.wikipedia.org/wiki/Laman_graph#Henneberg_construction). k defaults to 2 but can be set using `-K#`. l is always k(k+1)/2.
* `-N#`: all (complete graps) graphs with this number of nodes or fewer are considered (tight) sparse. The default value is max(⌊k⌋,2) or the highest n such that a complete graph on n vertices is (k,l)-sparse.

Both `-K` and `-L` accept rational numbers making it possible to generate, e.g., (3/2,2)-tight graphs (see results below). Note, however, that denominators equal to their numerator are ignored, e.g., `-K2/2` is equivalent to `-K2`. If rational arguments are not needed, define the macro `INT_KL` before compiling for a small increase in performance.

After forking the plugin, I added the following parameter:
* `-Mm:n`: computes the Laman number while each (nonisomorphic) Laman graph is generated (`-K2` must be given) and keeps the graph with maximum Laman number. The output will parse the maximum Laman number for the given number of vertices. Computation is done in parallel using openmp. If n is 0 then it will compute the maximum of all the Laman graphs it generates. If n>0, then it computes the maximum of every n+m graph it generates. This is useful when computing the maximum Laman number with multiple processes. One runs `-Mm:n` in n different processes with m=0,1,..., n-1. The maximum Laman number is the maximum of all the numbers from the output of these n processes.

**Note:** The `geng` plugin will parse the *cpu clock* and the *wall clock*. Since we are computing in parallel, we use the *wall clock* (`omp_get_wtime`) in the benchmarks below.

## Algorithm
The pebble game algorithm (see [Lee and Streinu (2008) Pebble game algorithms and sparse graphs](https://www.sciencedirect.com/science/article/pii/S0012365X07005602)) is used to generated the Laman graphs. While generating the graphs our algorithm (see [Capco et. al. The number of Realization of a Laman Graph](http://www.koutschan.de/data/laman/)), implemented in the `liblnumber` library, is used in parallel to compute the Laman numbers. 

<hr>

**What follows are for realizations of Laman graph on the plane!**

## Results and execution times
The tables below show the execution time when generating graphs while computing the maximum laman numbers on the plane. The setup used are
* laptop: Windows 10 Pro, 4-core Intel® Core™ i5-8350U CPU @1.70GHz, 16GiB 
* ippo: Debian GNU/Linux 11, 4-core Intel® Core™ i7-2600 CPU @3.40GHz, 16GiB 
* qft1: Debian GNU/Linux 11, total 16-core, 8-core Intel Xeon® CPU E5-2670 @2.6GHz, 386GiB 
* qft10: Debian GNU/Linux 11, total 16-core, AMD EPYC 7262 8-core, 2.9GHz, 2036GiB 
* leo5-64: [LEO5](https://www.uibk.ac.at/zid/systeme/hpc-systeme/leo5/) UIBK supercomputer, allocating one node with 64-cores.
* leo5-64xn: [LEO5](https://www.uibk.ac.at/zid/systeme/hpc-systeme/leo5/) UIBK supercomputer, n processes, each process in 1 node with 64-cores.

### Laman graphs with maximum Laman numbers on the plane
OEIS entry for number of Laman Graphs: [A227117](https://oeis.org/A227117 "Number of minimally rigid graphs in 2D on n vertices.")<br>
OEIS entry for maximum Laman numbers: [A306420](https://oeis.org/A306420)

Command: `geng $v -K2 -u -M$m:$n`

[Laman graphs](https://en.wikipedia.org/wiki/Laman_graph) are exactly the (2,3)-tight graphs. When increasing n by one, for large n, the number of graphs increases by a factor of approximately 30 while the execution time (for both computing Laman numbers on the plane and generating the graphs) increases by a factor of approximately 60. `geng` parallelizes well over physical cores but poorly over logical cores. When computing Laman numbers and employing [embarassingly parallel](https://en.wikipedia.org/wiki/Embarrassingly_parallel) methods with 4x multiprocessing with the same total number of cores, the execution time increases by a factor of approximately 1.2.


v             |     9    |   10    |    11      |    12      |      13       |     14         |
--------------|:--------:|:-------:|:----------:|:----------:|:-------------:|:--------------:|
Laman graphs  | 7 222    | 110 132 |  2 039 273 | 44 176 717 | 1 092 493 042 | 30 322 994 747 |
Max. Laman No.(plane)| 344      | 880     | 2 288      | 6 180      | 15 536        | 42 780 |
laptop    |   0.9 s  | 28 s    | 25.4 min   | \*| \*  | \*  |
ippo      |   0.2 s  | 6 s     | 3.7 min    | 3.48 hrs      |   \*  | \*  |
qft1      |   0.1 s  | 2.8 s   | 1.6 min    | 1.2 hrs       |  2.72 days  | \* |
qft10     |   0.07s  | 1.7 s   | 1 min      | 0.8 hrs       |  1.84 days  | \*  |
leo5-64   | \*       | \*      | \*         | 19 min        |   16.78 hrs  | \*  |
leo5-64x4 | \*  | \*  | \*  | 5.63 min  | 4.7hrs  | \* |
leo5-64x11 | \*  | \*  | \*  | \*  | \*  | 5.55days |

\* Not Measured

## Acknowledgement

Many thanks to my co-authors in the first paper proposing the algorithm and other people who motivated me to do this:
Matteo Gallet, Georg Grasegger, Christoph Koutschan, Jan Legersky, Niels Lubbes and Josef Schicho.

Some of the benchmarks have been achieved using the LEO HPC infrastructure of the University of Innsbruck. 
I would also like to thank [RISC](https://risc.jku.at/) for allowing me to use the RISC-DESY server cluster for some of my benchmarks.

## Citing
 
[![DOI](https://zenodo.org/badge/683425893.svg)](https://zenodo.org/badge/latestdoi/683425893)


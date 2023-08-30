# Maximum Laman Number
The Laman number library `liblnumber` is meant to be used as an external dynamic library by programs such as the Nauty Laman plugin. The library implements a fast algorithm that computes Laman number of a Laman graph based on our paper (with M. Gallet, G. Grasseger, C. Koutschan, N. Lubbes and J. Schicho) see [Capco et. al. The number of Realization of a Laman Graph](http://www.koutschan.de/data/laman/).

The `nauty` plugin is a fork of the [Nauty-Laman plugin](https://github.com/martinkjlarsson/nauty-laman-plugin) originally written by Martin Larrson for the utility 
program `geng` provided by [Nauty](http://pallini.di.uniroma1.it/). The plugin counts the number of non-isomorphic Laman graphs of a given number of 
vertices. The quick generation of the Laman graphs was implemented by Martin Larrson. I integrated parallel computing of Laman numbers into the plugin via dependency on the `liblnumber` library. Parallel computing is done via openmp.

## Installation with compiling

**Disclaimer:** All MSVC project files are setup for Release. You have to manually setup other build configurations that you want.

1. Clone this repo.
1. Download and save [Nauty](http://pallini.di.uniroma1.it/) source files in the `./nauty` folder. 
1. Build all objects from the `nauty` makefile (if you are building with gcc). You can use the VS2008 project file, `./vs2008/nauty.vcproj` to build nauty with MSVC. There are no dependencies.
1. Build the `liblnumber` library. 
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
1. `geng` has to be rebuilt 
in C++ with the flag `-D'PLUGIN="<path to this repo>/prunelaman.h"'` and linked with `lnumber`. It is also important to activate `openmp` when building, otherwise the program may still use multiple threads but there will be no speedup for parallel computing and multiple core computers will have roughly the same runtime as a single core one.
    * If you bild with gcc, you can build `geng` as follows 
    ```makefile
    g++  -Wno-write-strings -std=c++98 -I./ -o geng -fopenmp -O3 -s \ 
      -mpopcnt -march=native -D'PLUGIN="<path to this repo>/prunelaman.h"' \
      -DMAXN=WORDSIZE -DWORDSIZE=32 geng.c gtoolsW.o nautyW1.o nautilW1.o \
      naugraphW1.o schreier.o naurng.o  -L. -lstdc++ -lm -llnumber \
      -lgmpxx -lgomp -lgmp
    ```
    * If you
or if you are using Windows with MSVC you can compile and link using the VS2008 project file that I provided in VS2008 folder. The project file should be forward-compatible to newer MSVC. You will need to change the directory of the dependencies in the VS2008 project file

6. Once `geng` is compiled you can run it with the parameters to compute the maximum Laman numbers (see usage)

## Installation with prebuilt binaries

**Disclaimer:** I cannot gaurantee the prebuilt binaries for any Linux distribution. I know it is tricky to prepare prebuilt binaries for Linux. Nevertheless, I share my Debian builds. 

1. The executable files for Windows and Debian can be copied from the `./bin` folder
1. Run the file `geng` (see usage).

## Usage
The plugin was originally developed by Martin Larrson with the following additional parameters to `geng`:
* `-K#`: generate (k,l)-tight graphs where l = k(k+1)/2. Minimum degree and number of edges will default to k and kn-l, respectively. Sparse graphs can be generated by manually providing the minimum and maximum number of edges (e.g. `0:999`). In that case, the minimum degree will default to zero.
* `-L#`: provides the l when generating (k,l)-sparse or (k,l)-tight graphs.
* `-H`: generate (k,l)-tight graphs constructible by [Henneberg type I moves](https://en.wikipedia.org/wiki/Laman_graph#Henneberg_construction). k defaults to 2 but can be set using `-K#`. l is always k(k+1)/2.
* `-N#`: all (complete graps) graphs with this number of nodes or fewer are considered (tight) sparse. The default value is max(⌊k⌋,2) or the highest n such that a complete graph on n vertices is (k,l)-sparse.

Both `-K` and `-L` accept rational numbers making it possible to generate, e.g., (3/2,2)-tight graphs (see results below). Note, however, that denominators equal to their numerator are ignored, e.g., `-K2/2` is equivalent to `-K2`. If rational arguments are not needed, define the macro `INT_KL` before compiling for a small increase in performance.

After forking the plugin, I added the following parameter:
* `-M`: computes the Laman number while each (nonisomorphic) Laman graph is generated (`-K2` must be given) and keeps the graph with maximum Laman number. The output will parse the maximum Laman number for the given number of vertices. Computation is done in parallel using openmp.

**Note:** `geng` is going to parse the cpu clock. Since we are computing in parallel, we use the wall clock (parsed by `omp`) in the benchmark below.

## Algorithm
The pebble game algorithm (see [Lee and Streinu (2008) Pebble game algorithms and sparse graphs](https://www.sciencedirect.com/science/article/pii/S0012365X07005602)) is used to generated the Laman graphs. While generating the graphs our algorithm (see [Capco et. al. The number of Realization of a Laman Graph](http://www.koutschan.de/data/laman/)), implemented in the `liblnumber` library, is used in parallel to compute the Laman numbers. 

## Results and execution times
The tables below show the execution time when generating graphs while computing the maximum laman number. The setup used are
* laptop: Windows 10 Pro, 4-core/8-thread Intel® Core™ i5-8350U CPU @1.70 GHz, 16GiB 
* ippo: Debian GNU/Linux 11, 4-core/8-thread Intel® Core™ i7-2600 CPU @3.40 GHz, 16GiB 
* qft1: Debian GNU/Linux 11, 32-core/64-thread Intel Xeon® CPU E5-2670 @ 2.6GHz, 386GiB 
* qft10: Debian GNU/Linux 11, 32-core/64-thread AMD EPYC 7262 8-Core Processor, 2055.545 MHz, 2036GiB 

### Laman graphs with maximum Laman numbers
OEIS entry for number of Laman Graphs: [A227117](https://oeis.org/A227117 "Number of minimally rigid graphs in 2D on n vertices.")<br>
OEIS entry for Laman numbers: [A306420](https://oeis.org/A306420)

Command: `geng $n -K2 -u -M`

[Laman graphs](https://en.wikipedia.org/wiki/Laman_graph), minimally rigid graphs in 2D, are exactly the (2,3)-tight graphs. When increasing n by one, for large n, the number of graphs increases by a factor of approximately 30 while the execution time increases by a factor of approximately 60. `geng` parallelizes well over physical cores but poorly over logical cores. 


n             |     9    |   10    |    11      |    12      |      13       |     14         |
--------------|:--------:|:-------:|:----------:|:----------:|:-------------:|:--------------:|
Laman graphs  | 7 222    | 110 132 |  2 039 273 | 44 176 717 | 1 092 493 042 | 30 322 994 747 |
Max. Laman No.| 344      | 880     | 2 288      | 6 180      | 15 536        | *Not measured* |
laptop        |   0.9 s  | 28 s    | 25.4 min   | *Not measured*| *Not measured*  | *Not measured*  |
ippo         |   0.2 s  | 6 s     | 3.7 min    | 3.48 hrs   |   *Not measured*  | *Not measured*  |
qft1      |   0.1 s  | 2.8 s   | 1.6 min    | 1.2 hrs    |   *Not measured*  | *Not measured*  |
qft10     |   0.07s  | 1.7 s   | 1 min      | 0.8 hrs    |   *Not measured*  | *Not measured*  |



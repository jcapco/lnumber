#pragma once

#if defined(_MSC_VER) && defined(_WIN32)
  #ifdef LIBLNUMBER_EXPORTS
    #define LIBLNUMBER_LIBRARY_INTERFACE __declspec(dllexport)
  #else
    #define LIBLNUMBER_LIBRARY_INTERFACE __declspec(dllimport)
  #endif
#else
  #define LIBLNUMBER_LIBRARY_INTERFACE extern
  #include <cstddef>
#endif

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable: 4146 4244 4800)
#include <mpirxx.h>
#pragma warning(pop)
#else
  #ifdef __cplusplus
    #include <gmpxx.h>
  #else
    #include <gmp.h>
  #endif
#endif

#ifdef __cplusplus
extern "C"
{ 
#endif
  //convert symmetric matrix a[i,j] to flat symmetric no diagonal upper triangular 
  //nxn-matrix by getting the index of [i,j] for i>j
  inline int idx_flat(int i, int j)
  {
    return int(i*(i-1)/2 + j);
  }

  //stops MSVC from mangling
  size_t LIBLNUMBER_LIBRARY_INTERFACE laman_number(char* graph, size_t verts);
  //compute laman number of a nauty dense graph
  size_t LIBLNUMBER_LIBRARY_INTERFACE laman_number_nauty(unsigned int* ngraph, int verts);
  size_t LIBLNUMBER_LIBRARY_INTERFACE laman_numbern(mpz_ptr nptr, size_t verts);

#ifdef __cplusplus
}
#endif


// Compatibility header used only by the Python package build (pypi/setup.py).
//
// lnumber/inc/lib.h and lnumber/inc/laman_number.h include <gmpxx.h> (or
// <mpirxx.h> on Windows). The Python build does not link GMP: this header
// replaces the small subset of GMP used by lnumber/src/lib.cpp, namely
// building a number from a decimal string and reading its bits.
//
// Graphs have at most 26 vertices, i.e. at most 325 bits, so a plain
// vector of bits is enough. Characters other than decimal digits are
// ignored; as with GMP, the input is not validated.
//
// If lib.cpp starts using further GMP functions, the Python build fails to
// compile and this header has to be extended.
#pragma once

#include <cstddef>
#include <vector>

struct lnumber_mpz
{
  std::vector<unsigned char> bits; // least significant bit first, no trailing zeros
};

typedef lnumber_mpz* mpz_ptr;

// Number of bits (base 2 only); 1 for zero, as in GMP.
inline size_t mpz_sizeinbase(const lnumber_mpz* n, int /*base*/)
{
  return n->bits.empty() ? 1 : n->bits.size();
}

inline int mpz_tstbit(const lnumber_mpz* n, size_t i)
{
  return i < n->bits.size() ? n->bits[i] : 0;
}

class mpz_class
{
  lnumber_mpz value;

public:
  // Decimal string only (base 10), converted by repeated division by 2.
  mpz_class(const char* s, int /*base*/)
  {
    std::vector<unsigned char> digits; // most significant digit first
    for (; *s; ++s)
      if (*s >= '0' && *s <= '9') digits.push_back(*s - '0');

    while (!digits.empty())
    {
      std::vector<unsigned char> quotient;
      int remainder = 0;
      for (size_t i = 0; i < digits.size(); ++i)
      {
        int current = remainder * 10 + digits[i];
        remainder = current % 2;
        if (!quotient.empty() || current / 2) quotient.push_back(current / 2);
      }
      value.bits.push_back(remainder);
      digits.swap(quotient);
    }
    while (!value.bits.empty() && !value.bits.back()) value.bits.pop_back();
  }

  mpz_ptr get_mpz_t() { return &value; }
};

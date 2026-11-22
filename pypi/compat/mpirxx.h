// Compatibility header used only by the Python package build (pypi/setup.py).
//
// On Windows, lnumber/inc/lib.h and lnumber/inc/laman_number.h include
// <mpirxx.h> instead of <gmpxx.h>; both are served by the same replacement.
#pragma once

#include "gmpxx.h"

///////////////////////////////////////////////////////////////////////
//
// minkbinder.cc
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Mon Jul  9 13:59:45 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MINKBINDER_CC_DEFINED
#define MINKBINDER_CC_DEFINED

#include "cmodel/minkbinder.h"

//#define USE_ASM

namespace
{
  double square(double x) { return x * x; }
}

double MinkowskiBinder::minkDist2(mtx_const_iter x1) const
{
#ifndef USE_ASM

  double wt_sum = 0.0;
  mtx_const_iter wt = itsAttWeights;
  mtx_const_iter x2 = itsX2;

  for (; wt.has_more(); ++wt, ++x1, ++x2)
    {
      wt_sum += (*wt) * square((*x1) - (*x2));
    }
  return sqrt(wt_sum);

#else

// layout of *this
//      12 bytes -- wts
// +0        4 bytes -- data pointer
// +4        4 bytes -- stride
// +8        4 bytes -- stop pointer
//      12 bytes -- x2
// +12       4 bytes -- data pointer
// +16       4 bytes -- stride
// +20       4 bytes -- stop pointer
// +24  8 bytes  -- itsR
// +32  8 bytes  -- itsRinv

asm(
    // "pushl %ebp" assumed
    // "movl %esp,%ebp" assumed

    "subl $4,%esp;"             // {1c} 4b on stack
    "movl 8(%ebp),%edx;"        // {1c} copy the "this" pointer to edx

    "fldz;"                     // {2c} initialize the result to 0.0

    "pushl %edi;"               // {1c}
    "pushl %esi;"               // {1c}
    "pushl %ebx;"               // {1c}

#define WTS_REG "%eax"
    "movl (%edx),"WTS_REG";"    // {1c} copy wts.data to WTS_REG

    "movl 8(%edx),%ebx;"        // {1c} copy wts.stop to ebx

#define X2_REG "%ecx"
    "movl 12(%edx),"X2_REG";"   // {1c} copy x2.data to X2_REG

    "cmpl %ebx,"WTS_REG";"      // {1c}
    "jae cleanup;"              // {3/1c} jump if (wts.data >= wts.stop)

    "movl 4(%edx),%edi;"        // {1c} copy wts.stride to edi
    "sall $3,%edi;"             // {2c} wts.stride *= 8 (which is sizeof(double))
    "movl %edi,-4(%ebp);"       // {1c} copy wts.stride to stack-4

    "movl 16(%ebp),%edi;"       // {1c} copy x1.stride to edi
    "sall $3,%edi;"             // {2c} x1.stride *= 8

    "movl 16(%edx),%esi;"       // {1c} copy x2.stride to esi

    "sall $3,%esi;"             // {2c} x2.stride *= 8
    ".p2align 4,,7;"

#define X1_REG "%edx"

    "movl 12(%ebp),"X1_REG";"   // {1c} copy x1.data to edx

"sum_loop:;"
    "fldl ("X1_REG");"          // {2c} top = *x1.data
    "fsubl ("X2_REG");"         // {3/1c} top -= *x2.data
    "fldl ("WTS_REG");"         // {2c} top = *wts.data
    "fmul %st(1),%st;"          // {3/1c} top *= (*x1.data - *x2.data)

    "addl -4(%ebp),"WTS_REG";"  // {2c} wts.data += wts.stride
    "addl %edi,"X1_REG";"       // {1c} x1.data += x1.stride
    "addl %esi,"X2_REG";"       // {1c} x2.data += x2.stride

    "fmulp %st,%st(1);"         // {3/1c} top *= (*x1.data - *x2.data) and pop st(0)
    "faddp %st,%st(1);"         // {3/1c} add term to wt_sum and pop st(0)

    "cmpl %ebx,"WTS_REG";"      // {1c}
    "jb sum_loop;"              // {3/1c} loop if (wts.data < wts.stop)

    "fsqrt;"                    // {70c} st(0) = sqrt(wt_sum)

"cleanup:;"
    "leal -16(%ebp),%esp;"      // {1c}
    "popl %ebx;"                // {4c}
    "popl %esi;"                // {4c}
    "popl %edi;"                // {4c}
    // "movl %ebp,%esp;" assumed
    // "popl %ebp;" assumed
    // "ret;" assumed
    ); // end asm
#endif
}

static const char vcid_minkbinder_cc[] = "$Id$ $URL$";
#endif // !MINKBINDER_CC_DEFINED

///////////////////////////////////////////////////////////////////////
//
// cmodelutil.h
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Wed Jul 31 14:52:08 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELUTIL_H_DEFINED
#define CMODELUTIL_H_DEFINED

class Mtx;
class slice;

namespace CModelUtil
{
  void clampRows(Mtx& src, int firstrow, int nrows, const Mtx& hilo);

  Mtx getStoredExemplars(const slice& otherParams, int nstored,
                         const Mtx& hilo0, const Mtx& hilo1);

  Mtx getHiLo(const Mtx& src);
}

static const char vcid_cmodelutil_h[] = "$Header$";
#endif // !CMODELUTIL_H_DEFINED

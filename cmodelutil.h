///////////////////////////////////////////////////////////////////////
//
// cmodelutil.h
//
// Copyright (c) 2002-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Wed Jul 31 14:52:08 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELUTIL_H_DEFINED
#define CMODELUTIL_H_DEFINED

class mtx;
class slice;

namespace CModelUtil
{
  void clampRows(mtx& src, int firstrow, int nrows, const mtx& hilo);

  mtx getStoredExemplars(const slice& otherParams, int nstored,
                         const mtx& hilo0, const mtx& hilo1);

  mtx getHiLo(const mtx& src);
}

static const char vcid_cmodelutil_h[] = "$Id$ $URL$";
#endif // !CMODELUTIL_H_DEFINED

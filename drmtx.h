///////////////////////////////////////////////////////////////////////
//
// drmtx.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Tue Apr 17 11:21:42 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef DRMTX_H_DEFINED
#define DRMTX_H_DEFINED

#include <matrix.h>

///////////////////////////////////////////////////////////////////////
//
// DualRepMtx is a helper class for transitioning code from using
// mxArray*'s to using mtx's. It implements a dual-representation
// object which can be accessed and assigned as either a mtx or an
// mxArray*. The alternate representations are updated on demand
// (lazily).
//
///////////////////////////////////////////////////////////////////////


class DualRepMtx {
private:
  mutable mtx itsMtx;
  mutable bool itsMtxIsValid;

  mutable mxArray* itsArray;
  mutable bool itsArrayIsValid;

  void updateMtx() const
  {
    itsMtx = mtx(itsArray);
    itsMtxIsValid = true;
  }

  void updateArray() const
  {
    mlfAssign(&itsArray, make_mxarray(itsMtx));
    itsArrayIsValid = true;
  }

public:
  DualRepMtx() :
    itsMtx(0,0),
    itsMtxIsValid(true),
    itsArray(mclGetUninitializedArray()),
    itsArrayIsValid(false)
  {}

  DualRepMtx(const mtx& other) :
    itsMtx(other),
    itsMtxIsValid(true),
    itsArray(mclGetUninitializedArray()),
    itsArrayIsValid(false)
  {}

  ~DualRepMtx()
  {
    mxDestroyArray(itsArray);
  }

  void assignArray(mxArray* rhs)
  {
    mlfAssign(&itsArray, rhs);
    itsArrayIsValid = true;
    itsMtxIsValid = false;
  }

  void assignMtx(const mtx& rhs)
  {
    itsMtx = rhs;
    itsMtxIsValid = true;
    itsArrayIsValid = false;
  }

  mxArray* asArray() const
  {
    if (!itsArrayIsValid) updateArray();
    return itsArray;
  }

  const mtx& asMtx() const
  {
    if (!itsMtxIsValid) updateMtx();
    return itsMtx;
  }

  mtx& ncMtx()
  {
    if (!itsMtxIsValid) updateMtx();
    itsArrayIsValid = false;
    return itsMtx;
  }
};

static const char vcid_drmtx_h[] = "$Header$";
#endif // !DRMTX_H_DEFINED

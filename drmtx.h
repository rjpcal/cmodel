///////////////////////////////////////////////////////////////////////
//
// drmtx.h
//
// Copyright (c) 1998-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 17 11:21:42 2001
// written: Mon Feb  4 18:12:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef DRMTX_H_DEFINED
#define DRMTX_H_DEFINED

#include "libmatlbm.h"

///////////////////////////////////////////////////////////////////////
//
// DualRepMtx is a helper class for transitioning code from using
// mxArray*'s to using Mtx's. It implements a dual-representation
// object which can be accessed and assigned as either a Mtx or an
// mxArray*. The alternate representations are updated on demand
// (lazily).
//
///////////////////////////////////////////////////////////////////////


class DualRepMtx {
private:
  mutable Mtx itsMtx;
  mutable bool itsMtxIsValid;

  mutable mxArray* itsArray;
  mutable bool itsArrayIsValid;

  void updateMtx() const
  {
         itsMtx = Mtx(itsArray);
         itsMtxIsValid = true;
  }

  void updateArray() const
  {
         mlfAssign(&itsArray, itsMtx.makeMxArray());
         itsArrayIsValid = true;
  }

public:
  DualRepMtx() :
         itsMtx(0,0),
         itsMtxIsValid(true),
         itsArray(mclGetUninitializedArray()),
         itsArrayIsValid(false)
  {}

  DualRepMtx(const Mtx& other) :
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

  void assignMtx(const Mtx& rhs)
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

  const Mtx& asMtx() const
  {
         if (!itsMtxIsValid) updateMtx();
         return itsMtx;
  }

  Mtx& ncMtx()
  {
         if (!itsMtxIsValid) updateMtx();
         itsArrayIsValid = false;
         return itsMtx;
  }
};

static const char vcid_drmtx_h[] = "$Header$";
#endif // !DRMTX_H_DEFINED

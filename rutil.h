///////////////////////////////////////////////////////////////////////
//
// rutil.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 19 16:20:06 2001
// written: Wed Feb 21 14:14:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef RUTIL_H_DEFINED
#define RUTIL_H_DEFINED

#define validateInput(x) (x) = mlfAssign(&(x), mclVa((x), #x))

#define validateVariable(x) (x) = mlfAssign(&(x), mclVv((x), #x))

#include "libmatlb.h"

class Rat {
public:
  Rat(mxArray* a) :
	 mrows(mxGetM(a)),
	 ncols(mxGetN(a)),
	 data(mxGetPr(a)),
	 array(a)
  {}

  void print() const
    {
		mexPrintf("array = %p, mrows = %d, ncols = %d\n", array, mrows, ncols);
		for(int i = 0; i < mrows; ++i)
		  {
			 for(int j = 0; j < ncols; ++j)
				mexPrintf("%7.4f   ", at(i,j));
			 mexPrintf("\n");
		  }
		mexPrintf("\n");
	 }

  int index(int row, int col) const { return row + (col*mrows); }

  double* address(int row, int col) { return data + index(row, col); }

  double& at(int row, int col) { return data[index(row, col)]; }

  double at(int row, int col) const { return data[index(row, col)]; }

  double& at(int elem) { return data[elem]; }

  double at(int elem) const { return data[elem]; }

  int length() const { return (mrows > ncols) ? mrows : ncols; }

  int nelems() const { return mrows*ncols; }

  int mrows;
  int ncols;
  double* data;
  void* array;
};

static const char vcid_rutil_h[] = "$Header$";
#endif // !RUTIL_H_DEFINED

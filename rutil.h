///////////////////////////////////////////////////////////////////////
//
// rutil.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 19 16:20:06 2001
// written: Tue Feb 20 14:21:59 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef RUTIL_H_DEFINED
#define RUTIL_H_DEFINED

#define validateInput(x) (x) = mclVa((x), #x)

#define validateVariable(x) (x) = mclVv((x), #x)

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
		int i;
		mexPrintf("array = %p, mrows = %d, ncols = %d\n", array, mrows, ncols);
		for(i=0; i < mrows * ncols; ++i) {
		  mexPrintf("%7.4f\n", data[i]);
		}
		mexPrintf("\n");
	 }

  int length() const { return (mrows > ncols) ? mrows : ncols; }

  int nelems() const { return mrows*ncols; }

  int mrows;
  int ncols;
  double* data;
  void* array;
};

static const char vcid_rutil_h[] = "$Header$";
#endif // !RUTIL_H_DEFINED

///////////////////////////////////////////////////////////////////////
//
// rutil.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 19 16:20:06 2001
// written: Tue Mar  6 12:12:22 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef RUTIL_H_DEFINED
#define RUTIL_H_DEFINED

#define validateInput(x) (x) = mlfAssign(&(x), mclVa((x), #x))

#define validateVariable(x) (x) = mlfAssign(&(x), mclVv((x), #x))

#include "libmatlb.h"

class Slice {
  const double* data;
  int const stride;
public:
  Slice(const double* d, int s) : data(d), stride(s) {}

  double operator[](int i) const { return *(data + stride*i); }

  double operator*() const { return *data; }

  void bump() { data += stride; }
};

class Rat {
private:
  enum StorageType { BORROWED, OWNED };

  Rat& operator=(const Rat&); // not allowed

  Rat(const Rat& other, int column /* zero-based */) :
	 mrows_(other.mrows_),
	 ncols_(1),
	 data_(other.data_ + (other.mrows_*column)),
	 storage_(BORROWED)
  {}

public:
  Rat(mxArray* a) :
	 mrows_(mxGetM(a)),
	 ncols_(mxGetN(a)),
	 data_(mxGetPr(a)),
	 storage_(BORROWED)
  {}

  Rat(const Rat& other) :
	 mrows_(other.mrows_),
	 ncols_(other.ncols_),
	 data_(0),
	 storage_(other.storage_)
  {
	 switch(storage_) {
	 case BORROWED:
		data_ = other.data_;
		break;
	 case OWNED:
		data_ = new double[nelems()];
		for (int i = 0; i < nelems(); ++i)
		  data_[i] = other.data_[i];
		break;
	 }
  }

  ~Rat()
  {
	 switch(storage_) {
	 case OWNED:
		delete [] data_;
		break;
	 case BORROWED:
		break;
	 }
  }

  void print() const
    {
		mexPrintf("mrows = %d, ncols = %d\n", mrows_, ncols_);
		for(int i = 0; i < mrows_; ++i)
		  {
			 for(int j = 0; j < ncols_; ++j)
				mexPrintf("%7.4f   ", at(i,j));
			 mexPrintf("\n");
		  }
		mexPrintf("\n");
	 }

  int index(int row, int col) const { return row + (col*mrows_); }

  double* address(int row, int col) { return data_ + index(row, col); }

  const double* address(int row, int col) const { return data_ + index(row, col); }

  double& at(int row, int col) { return data_[index(row, col)]; }

  double at(int row, int col) const { return data_[index(row, col)]; }

  double& at(int elem) { return data_[elem]; }

  double at(int elem) const { return data_[elem]; }

  int length() const { return (mrows_ > ncols_) ? mrows_ : ncols_; }

  int nelems() const { return mrows_*ncols_; }

  int mrows() const { return mrows_; }

  int ncols() const { return ncols_; }

  double* data() { return data_; }

  const double* data() const { return data_; }

  Rat columnSlice(int column) const { return Rat(*this, column); }

private:
  int mrows_;
  int ncols_;
  double* data_;
  StorageType storage_;
};

static const char vcid_rutil_h[] = "$Header$";
#endif // !RUTIL_H_DEFINED

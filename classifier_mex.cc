///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Thu Mar  8 10:25:22 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_MEX_CC_DEFINED
#define CLASSIFIER_MEX_CC_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "classifier.h"

#include "libmatlb.h"

///////////////////////////////////////////////////////////////////////
//
// The function "mexLibrary" is a Compiler-generated mex wrapper,
// suitable for building a MEX-function. It initializes any persistent
// variables as well as a function table for use by the feval
// function. It then calls the function "mlxClassifier". Finally, it
// clears the feval table and exits.
//
///////////////////////////////////////////////////////////////////////

extern "C"
mex_information mexLibrary() {

  static mexFunctionTableEntry function_table[1] = {
	 { "classifier", mlxClassifier, 4, 1, &_local_function_table_classifier }
  };

  static const char * path_list_[1] = { "/matlab_r.12/toolbox/matlab/specfun" };

  static _mexInitTermTableEntry init_term_table[1] = {
	 { InitializeModule_classifier, TerminateModule_classifier },
  };

  static _mex_information _mex_info
	 = { 1, 1, function_table, 0, NULL, 1, path_list_, 1, init_term_table };

  return &_mex_info;
}

static const char vcid_classifier_mex_cc[] = "$Header$";
#endif // !CLASSIFIER_MEX_CC_DEFINED

///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Thu Mar  8 10:50:05 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_MEX_CC_DEFINED
#define CLASSIFIER_MEX_CC_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "classifier.h"
#include "classifier_mex.h"

#include "rutil.h"
#include "trace.h"

#include "libmatlb.h"

void InitializeModule_classifier(void) {
}

void TerminateModule_classifier(void) {
}

_mexLocalFunctionTable _local_function_table_classifier
  = { 0, (mexFunctionTableEntry *)NULL };

static mxArray * Mclassifier(int /* nargout_ */,
									  mxArray * modelParams_mx,
									  mxArray * objParams_mx,
									  mxArray * observedIncidence_mx,
									  mxArray * numStoredExemplars_mx)
{
DOTRACE("Mclassifier");

#ifdef LOCAL_PROF
  if (mxGetScalar(numStoredExemplars_mx) == -1) {
	 ofstream ofs("profdata.out");
	 Util::Prof::printAllProfData(ofs);
	 return mxCreateScalarDouble(-1.0);
  }

  if (mxGetScalar(numStoredExemplars_mx) == -2) {
	 Util::Prof::resetAllProfData();
	 return mxCreateScalarDouble(-2.0);
  }
#endif

  //---------------------------------------------------------------------
  //
  // Set up
  //
  //---------------------------------------------------------------------

  mexLocalFunctionTable save_local_function_table_ =
    mclSetCurrentLocalFunctionTable(&_local_function_table_classifier);

  mclCopyArray(&modelParams_mx); // THIS IS REQUIRED
  validateInput(modelParams_mx);

  validateInput(objParams_mx);
  const Rat objParams(objParams_mx);

  validateInput(observedIncidence_mx);
  const Rat observedIncidence(observedIncidence_mx);

  int numStoredExemplars = int(mxGetScalar(numStoredExemplars_mx));

  //---------------------------------------------------------------------
  //
  // Call the real computational function for each set of model params
  //
  //---------------------------------------------------------------------

  Rat allModelParams(modelParams_mx);

  mxArray* ll_mx = mxCreateDoubleMatrix(allModelParams.ncols(), 1, mxREAL);
  Rat ll(ll_mx);

  ModelCssm model(objParams,
						observedIncidence,
						numStoredExemplars);

  for (int i = 0; i < allModelParams.ncols(); ++i)
	 {
		Rat modelParams(allModelParams.columnSlice(i));
		ll.at(i) = model.loglikelihoodFor(modelParams);
	 }

  //---------------------------------------------------------------------
  //
  // Clean up
  //
  //---------------------------------------------------------------------

  mxDestroyArray(modelParams_mx);

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  return ll_mx;
}

/*
 * The function "mlfClassifier" contains the normal interface for the
 * "classifier" function. This function processes any input arguments
 * and passes them to the implementation version of the function,
 * appearing above.
*/
mxArray* mlfClassifier(mxArray* modelParams_mx,
								mxArray* objParams_mx,
								mxArray* observedIncidence_mx,
								mxArray* numStoredExemplars_mx)
{
DOTRACE("mlfClassifier");

  int nargout = 1;

  mlfEnterNewContext(0, 4,
							modelParams_mx, objParams_mx,
							observedIncidence_mx, numStoredExemplars_mx);

  mxArray * ll = Mclassifier(nargout, modelParams_mx, objParams_mx,
									  observedIncidence_mx, numStoredExemplars_mx);

  mlfRestorePreviousContext(0, 4,
									 modelParams_mx, objParams_mx,
									 observedIncidence_mx, numStoredExemplars_mx);

  return mlfReturnValue(ll);
}

/*
 * The function "mlxClassifier" contains the feval interface for the
 * "classifier" function. The feval function calls the implementation
 * version of classifier through this function. This function
 * processes any input arguments and passes them to the implementation
 * version of the function, appearing above.  */
void mlxClassifier(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
DOTRACE("mlxClassifier");

  mxArray * mprhs[4];
  mxArray * mplhs[1];
  int i;
  if (nlhs > 1) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of outputs (1).");
  }
  if (nrhs > 4) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of inputs (4).");
  }
  for (i = 0; i < 1; ++i) {
	 mplhs[i] = mclGetUninitializedArray();
  }
  for (i = 0; i < 4 && i < nrhs; ++i) {
	 mprhs[i] = prhs[i];
  }
  for (; i < 4; ++i) {
	 mprhs[i] = NULL;
  }

  mlfEnterNewContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  mplhs[0] = Mclassifier(nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  mlfRestorePreviousContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  plhs[0] = mplhs[0];
}

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

///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Thu Mar  8 15:38:30 2001
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

#include "cmodelcssm.h"
#include "cmodelgcm.h"
#include "cmodelpbi.h"
#include "cmodelwpsm.h"

#include "error.h"
#include "rutil.h"
#include "strings.h"
#include "trace.h"

#include "util/pointers.h"

#include <fstream.h>
#include "libmatlb.h"

void InitializeModule_classifier() {
#ifdef LOCAL_DEBUG
  mexPrintf("loading 'classifier' mex file\n");
#endif
}

void TerminateModule_classifier() {
#ifdef LOCAL_DEBUG
  mexPrintf("unloading 'classifier' mex file\n");
#endif
}

_mexLocalFunctionTable _local_function_table_classifier
  = { 0, (mexFunctionTableEntry *)NULL };

namespace Local {
  fixed_string getString(const mxArray* arr)
  {
	 fixed_string str(mxGetM(arr) * mxGetN(arr) * sizeof(mxChar));

	 mxGetString(arr, str.data(), str.length() + 1);

	 return str;
  }
}

shared_ptr<Classifier> makeClassifier(const fixed_string& whichType,
												  const Rat& objParams,
												  const Rat& observedIncidence,
												  mxArray* optArgs_mx)
{
DOTRACE("makeClassifier");
  if (whichType == "cssm")
	 {
		int numStoredExemplars = optArgs_mx ? int(mxGetScalar(optArgs_mx)) : 0;
		return shared_ptr<Classifier>(
		  new CModelCssm(objParams, observedIncidence, numStoredExemplars));
	 }
  else if (whichType == "gcm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelGcm(objParams, observedIncidence));
	 }
  else if (whichType == "pbi")
	 {
	   return shared_ptr<Classifier>(
		  new CModelPbi(objParams, observedIncidence));
	 }
  else if (whichType == "wpsm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelWpsm(objParams, observedIncidence));
	 }
  else
	 {
		ErrorWithMsg err("unknown classifier type: ");
		err.appendMsg(whichType.c_str());
		throw err;
	 }
}

static mxArray* Mclassifier(int /* nargout_ */,
									 mxArray* modelParams_mx,
									 mxArray* modelName_mx,
									 mxArray* actionRequest_mx,
									 mxArray* objParams_mx,
									 mxArray* observedIncidence_mx,
									 mxArray* optArgs_mx)
{
DOTRACE("Mclassifier");

  try {

	 fixed_string modelName = Local::getString(modelName_mx);

	 fixed_string actionRequest = Local::getString(actionRequest_mx);


#ifdef LOCAL_PROF
	 if (optArgs_mx && mxGetScalar(optArgs_mx) == -1) {
		ofstream ofs("profdata.out");
		Util::Prof::printAllProfData(ofs);
		return mxCreateScalarDouble(-1.0);
	 }

	 if (optArgs_mx && mxGetScalar(optArgs_mx) == -2) {
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

	 Rat allModelParams(modelParams_mx);

	 mxArray* result_mx = mxCreateDoubleMatrix(allModelParams.ncols(), 1, mxREAL);
	 Rat result(result_mx);

	 shared_ptr<Classifier> model = makeClassifier(modelName,
																  objParams,
																  observedIncidence,
																  optArgs_mx);

	 int multiplier = 1;
	 // check for minus sign
	 if (actionRequest.c_str()[0] == '-')
		{
		  multiplier = -1;
		  fixed_string trimmed = actionRequest.c_str() + 1;
		  actionRequest = trimmed;
		}

	 //---------------------------------------------------------------------
	 //
	 // Call the requested computational function for each set of model
	 // params
	 //
	 //---------------------------------------------------------------------

	 if ( actionRequest == "ll" || actionRequest == "llc" )
		{
		  for (int i = 0; i < allModelParams.ncols(); ++i)
			 {
				Rat modelParams(allModelParams.columnSlice(i));
				result.at(i) = multiplier * model->currentLogL(modelParams);
			 }
		}
	 else if ( actionRequest == "llf" )
		{
		  for (int i = 0; i < allModelParams.ncols(); ++i)
			 {
				Rat modelParams(allModelParams.columnSlice(i));
				result.at(i) = multiplier * model->fullLogL();
			 }
		}
	 else if ( actionRequest == "dev" )
		{
		  for (int i = 0; i < allModelParams.ncols(); ++i)
			 {
				Rat modelParams(allModelParams.columnSlice(i));
				result.at(i) = multiplier * model->deviance(modelParams);
			 }
		}
	 else
		{
		  ErrorWithMsg err("unknown model action: ");
		  err.appendMsg(actionRequest.c_str());
		  throw err;
 		}

	 //---------------------------------------------------------------------
	 //
	 // Clean up
	 //
	 //---------------------------------------------------------------------

	 mxDestroyArray(modelParams_mx);

	 mclSetCurrentLocalFunctionTable(save_local_function_table_);

	 return result_mx;
  }
  catch (ErrorWithMsg& err) {
	 mexErrMsgTxt(err.msg_cstr());
  }
  catch (...) {
	 mexErrMsgTxt("an unknown C++ exception occurred.");
  }
}




///////////////////////////////////////////////////////////////////////
//
// The function "mlfClassifier" contains the normal interface for the
// "classifier" function. This function processes any input arguments
// and passes them to the implementation version of the function,
// appearing above.
//
///////////////////////////////////////////////////////////////////////

extern "C"
mxArray* mlfClassifier(mxArray* modelParams_mx,
							  mxArray* modelName_mx,
							  mxArray* actionRequest_mx,
							  mxArray* objParams_mx,
							  mxArray* observedIncidence_mx,
							  mxArray* optArgs_mx)
{
DOTRACE("mlfClassifier");

  int nargout = 1;

  mlfEnterNewContext(0, CLASSIFIER_NARGIN,
							modelParams_mx, modelName_mx, actionRequest_mx,
							objParams_mx, observedIncidence_mx, optArgs_mx);

  mxArray* result = Mclassifier(nargout,
										  modelParams_mx, modelName_mx, actionRequest_mx,
										  objParams_mx,
										  observedIncidence_mx, optArgs_mx);

  mlfRestorePreviousContext(0, CLASSIFIER_NARGIN,
									 modelParams_mx, modelName_mx, actionRequest_mx, 
									 objParams_mx, observedIncidence_mx,
									 optArgs_mx);

  return mlfReturnValue(result);
}




///////////////////////////////////////////////////////////////////////
//
// The function "mlxClassifier" contains the feval interface for the
// "classifier" function. The feval function calls the implementation
// version of classifier through this function. This function
// processes any input arguments and passes them to the implementation
// version of the function, appearing above.
//
///////////////////////////////////////////////////////////////////////

extern "C"
void mlxClassifier(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
DOTRACE("mlxClassifier");

  const int NUM_OUTPUTS = 1;

  mxArray * mprhs[CLASSIFIER_NARGIN];
  mxArray * mplhs[NUM_OUTPUTS];
  int i;
  if (nlhs > NUM_OUTPUTS) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of outputs (1).");
  }
  if (nrhs > CLASSIFIER_NARGIN) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of inputs.");
  }
  for (i = 0; i < NUM_OUTPUTS; ++i) {
	 mplhs[i] = mclGetUninitializedArray();
  }
  for (i = 0; i < CLASSIFIER_NARGIN && i < nrhs; ++i) {
	 mprhs[i] = prhs[i];
  }
  for (; i < CLASSIFIER_NARGIN; ++i) {
	 mprhs[i] = NULL;
  }

  mplhs[0] = mlfClassifier(mprhs[0], mprhs[1], mprhs[2],
									mprhs[3], mprhs[4], mprhs[5]);

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
#ifndef LOCAL_DEBUG
	 { "classifier", mlxClassifier, 5, 1, &_local_function_table_classifier }
#else
	 { "dclassifier", mlxClassifier, 5, 1, &_local_function_table_classifier }	 
#endif
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

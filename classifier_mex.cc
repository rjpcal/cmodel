///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Fri May 11 16:30:10 2001
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
#include "cmodelcuevalidity.h"
#include "cmodelgcm.h"
#include "cmodelpbi.h"
#include "cmodelwpsm.h"

#include "util/error.h"
#include "mexbuf.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "rutil.h"
#include "util/strings.h"
#include "util/pointers.h"

#include <iostream.h>
#include "libmatlb.h"

#include "util/trace.h"

namespace {
  shared_ptr<CModelCssm>* recentModel = 0;
  Mtx* recentObjParams = 0;
  int recentNumStored = -1;

  MexBuf* mexBuf = 0;
}

void InitializeModule_classifier()
{
DOTRACE("InitializeModule_classifier");
  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  cout = mexBuf;
  cerr = mexBuf;
#else
  cout.rdbuf(mexBuf);
  cerr.rdbuf(mexBuf);
#endif

#ifdef LOCAL_DEBUG
  mexPrintf("loading 'classifier' mex file\n");
#endif
  recentModel = new shared_ptr<CModelCssm>(0);
  recentObjParams = new Mtx(0,0);
}

void TerminateModule_classifier()
{
DOTRACE("TerminateModule_classifier");
#ifdef LOCAL_DEBUG
  mexPrintf("unloading 'classifier' mex file\n");
#endif

  delete recentObjParams;
  delete recentModel;

  delete mexBuf;
}

static mexFunctionTableEntry classifierFunctionTable[1] = {
#ifndef LOCAL_DEBUG
  { "classifier", mlxClassifier, 6, 1, &_local_function_table_classifier }
#else
  { "dclassifier", mlxClassifier, 6, 1, &_local_function_table_classifier }	 
#endif
};

_mexLocalFunctionTable _local_function_table_classifier
  = { 1, classifierFunctionTable };

shared_ptr<Classifier> makeClassifier(const fixed_string& whichType,
												  const Mtx& objParams,
												  mxArray* extraArgs_mx)
{
DOTRACE("makeClassifier");
  if (whichType == "cssm")
	 {

		int numStoredExemplars = 0;
		if (extraArgs_mx && mxIsStruct(extraArgs_mx))
		  {
			 mxArray* ns_mx = mxGetField(extraArgs_mx, 0, "numStoredExemplars");
			 if (ns_mx)
				numStoredExemplars = int(mxGetScalar(ns_mx));
		  }

		if ( (numStoredExemplars == recentNumStored) &&
			  (objParams == *recentObjParams) )
		  {
			 DOTRACE("use old");

			 return *recentModel;
		  }
		else
		  {
			 DOTRACE("make new");

			 *recentObjParams = objParams;
			 recentObjParams->makeUnique();

			 recentNumStored = numStoredExemplars;

          recentModel->reset(
              new CModelCssm(*recentObjParams,
                             CModelExemplar::EXP_DECAY, recentNumStored));

			 return *recentModel;
        }
	 }
  else if (whichType == "gcm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelGcm(objParams, CModelExemplar::EXP_DECAY));
	 }
  else if (whichType == "adm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelGcm(objParams, CModelExemplar::LINEAR_DECAY));
	 }
  else if (whichType == "pbi")
	 {
	   return shared_ptr<Classifier>(
		  new CModelPbi(objParams));
	 }
  else if (whichType == "wpsm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelWpsm(objParams, CModelExemplar::EXP_DECAY));
	 }
  else if (whichType == "wpm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelWpsm(objParams, CModelExemplar::LINEAR_DECAY));
	 }
  else if (whichType == "wcvm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelCueValidity(objParams, CModelCueValidity::NO_FREQ_WEIGHT));
	 }
  else if (whichType == "wfcvm")
	 {
	   return shared_ptr<Classifier>(
		  new CModelCueValidity(objParams, CModelCueValidity::FREQ_WEIGHT));
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
									 mxArray* extraArgs_mx)
{
DOTRACE("Mclassifier");

  try {

	 fixed_string modelName = MxWrapper::extractString(modelName_mx);

	 fixed_string actionRequest = MxWrapper::extractString(actionRequest_mx);

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
	 if (extraArgs_mx && mxIsStruct(extraArgs_mx))
		{
		  mxArray* debugFlag = mxGetField(extraArgs_mx, 0, "debugFlag");

		  if (debugFlag)
			 {
				if (mxGetScalar(debugFlag) == -1) {
				  Util::Prof::printAllProfData(cerr);
				  return mxCreateScalarDouble(-1.0);
				}

				if (mxGetScalar(debugFlag) == -2) {
				  Util::Prof::resetAllProfData();
				  return mxCreateScalarDouble(-2.0);
				}
			 }
		}
#endif

	 //---------------------------------------------------------------------
	 //
	 // Set up
	 //
	 //---------------------------------------------------------------------

	 mexLocalFunctionTable save_local_function_table_ =
		mclSetCurrentLocalFunctionTable(&_local_function_table_classifier);

	 validateInput(modelParams_mx);
	 // This Mtx will copy the data leaving the original mxArray untouched
	 Mtx allModelParams(modelParams_mx);

	 const Mtx objParams(MxWrapper::extractStructField(extraArgs_mx,
																		"objParams"));

	 shared_ptr<Classifier> model =
		makeClassifier(modelName, objParams, extraArgs_mx);

	 Classifier::RequestResult res = model->handleRequest(actionRequest,
																			allModelParams,
																			extraArgs_mx);

	 if ( !res.requestHandled )
		{
		  dynamic_string msg("unknown model action: ");
		  msg.append(actionRequest.c_str());
		  mexWarnMsgTxt(msg.c_str());
 		}

	 //---------------------------------------------------------------------
	 //
	 // Clean up
	 //
	 //---------------------------------------------------------------------

	 mclSetCurrentLocalFunctionTable(save_local_function_table_);

	 return res.result.release();
  }
  catch (ErrorWithMsg& err) {
	 mexErrMsgTxt(err.msg_cstr());
  }
  catch (...) {
	 mexErrMsgTxt("an unknown C++ exception occurred.");
  }

  return (mxArray*) 0; // can't happen, but placate compiler
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
							  mxArray* extraArgs_mx)
{
DOTRACE("mlfClassifier");

  int nargout = 1;

  mlfEnterNewContext(0, CLASSIFIER_NARGIN,
                     modelParams_mx, modelName_mx, actionRequest_mx,
							extraArgs_mx);

  mxArray* result = Mclassifier(nargout,
                                modelParams_mx, modelName_mx, actionRequest_mx,
										  extraArgs_mx);

  mlfRestorePreviousContext(0, CLASSIFIER_NARGIN,
                            modelParams_mx, modelName_mx, actionRequest_mx, 
									 extraArgs_mx);

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
void mlxClassifier(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[])
{
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

  mplhs[0] = mlfClassifier(mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  plhs[0] = mplhs[0];
}




///////////////////////////////////////////////////////////////////////
//
// The function "mexLibrary" is a Compiler-generated mex wrapper,
// which (apparently, since I can't find any docs for this API)
// returns information about the functions provided in the containing
// MEX file.
//
///////////////////////////////////////////////////////////////////////

extern "C"
mex_information mexLibrary()
{
DOTRACE("mexLibrary");

  static _mexInitTermTableEntry init_term_table[1] = {
	 { InitializeModule_classifier, TerminateModule_classifier },
  };

  static _mex_information _mex_info
	 = { 1, 1, classifierFunctionTable, 0, NULL, 0, NULL, 1, init_term_table };

  return &_mex_info;
}

static const char vcid_classifier_mex_cc[] = "$Header$";
#endif // !CLASSIFIER_MEX_CC_DEFINED

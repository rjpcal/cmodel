///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Mon Feb 25 14:20:36 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_MEX_CC_DEFINED
#define CLASSIFIER_MEX_CC_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "classifier.h"

#include "cmodelcssm.h"
#include "cmodelcuevalidity.h"
#include "cmodelgcm.h"
#include "cmodelpbi.h"
#include "cmodelspc.h"
#include "cmodelwpsm.h"

#include "mtx/mtx.h"

#include "mx/mexpkg.h"
#include "mx/mx.h"

#include "util/error.h"
#include "util/pointers.h"
#include "util/strings.h"

#include <exception>
#include <iostream>
#include <libmatlb.h>

#include "util/trace.h"

#ifdef LOCAL_DEBUG
#define MEXFUNCNAME "dclassifier"
#else
#define MEXFUNCNAME "classifier"
#endif

namespace
{
  shared_ptr<CModelCssm>* recentModel = 0;
  Mtx* recentObjParams = 0;
  int recentNumStored = -1;

  MexPkg* mexPkg = 0;
}

void InitializeModule_classifier()
{
  DOTRACE("InitializeModule_classifier");

  mexPrintf("loading '" MEXFUNCNAME "' mex file\n");

  recentModel = new shared_ptr<CModelCssm>(0);
  recentObjParams = new Mtx(0,0);
}

void TerminateModule_classifier()
{
  DOTRACE("TerminateModule_classifier");

  Util::Prof::printAtExit(false);

  mexPrintf("unloading '" MEXFUNCNAME "' mex file...\n");

  mexPrintf("\tdeleting recentObjParams...\n");
  delete recentObjParams;

  mexPrintf("\tdeleting recentModel...\n");
  delete recentModel;

  mexPrintf("\tdeleting mexPkg...\n");
  delete mexPkg;
  mexPrintf("\tdone.\n");
}

shared_ptr<Classifier> makeClassifier(const fstring& whichType,
                                      const Mtx& objParams,
                                      mxArray* extraArgs_mx)
{
  DOTRACE("makeClassifier");
  if (whichType == "cssm")
    {

      int numStoredExemplars =
        Mx::getIntField(extraArgs_mx, "numStoredExemplars");

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

          recentModel->reset
            (new CModelCssm(*recentObjParams,
                            CModelExemplar::EXP_DECAY, recentNumStored));

          return *recentModel;
        }
    }
  else if (whichType == "gcm")
    {
      return shared_ptr<Classifier>
        (new CModelGcm(objParams, CModelExemplar::EXP_DECAY));
    }
  else if (whichType == "adm")
    {
      return shared_ptr<Classifier>
        (new CModelGcm(objParams, CModelExemplar::LINEAR_DECAY));
    }
  else if (whichType == "pbi")
    {
      return shared_ptr<Classifier>(new CModelPbi(objParams));
    }
  else if (whichType == "wpsm")
    {
      return shared_ptr<Classifier>
        (new CModelWpsm(objParams, CModelExemplar::EXP_DECAY));
    }
  else if (whichType == "wpm")
    {
      return shared_ptr<Classifier>
        (new CModelWpsm(objParams, CModelExemplar::LINEAR_DECAY));
    }
  else if (whichType == "wcvm")
    {
      return shared_ptr<Classifier>
        (new CModelCueValidity(objParams, CModelCueValidity::NO_FREQ_WEIGHT));
    }
  else if (whichType == "wfcvm")
    {
      return shared_ptr<Classifier>
        (new CModelCueValidity(objParams, CModelCueValidity::FREQ_WEIGHT));
    }
  else if (whichType == "spc")
    {
      int numStoredExemplars =
        Mx::getIntField(extraArgs_mx, "numStoredExemplars");

      return shared_ptr<Classifier>
        (new CModelSPC(objParams, numStoredExemplars));
    }
  else
    {
      throw Util::Error(fstring("unknown classifier type: ", whichType));
    }

  // Can't get here, but placate compiler with a return value
  return shared_ptr<Classifier>(0);
}

static mxArray* Mclassifier(const Mtx& modelParams,
                            const fstring& modelName,
                            const fstring& actionRequest,
                            mxArray* extraArgs_mx)
{
  DOTRACE("Mclassifier");

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
  if (Mx::hasField(extraArgs_mx, "debugFlag"))
    {
      int debugFlag = Mx::getIntField(extraArgs_mx, "debugFlag");

      if (debugFlag == -1)       Util::Prof::printAllProfData(std::cerr);
      else if (debugFlag == -2)  Util::Prof::resetAllProfData();

      return mxCreateScalarDouble(debugFlag);
    }
#endif

  const Mtx objParams(Mx::getField(extraArgs_mx, "objParams"));

  shared_ptr<Classifier> model =
    makeClassifier(modelName, objParams, extraArgs_mx);

  Classifier::RequestResult res = model->handleRequest(actionRequest,
                                                       modelParams,
                                                       extraArgs_mx);

  if ( !res.requestHandled )
    {
      fstring msg("unknown model action: ", actionRequest);
      mexWarnMsgTxt(msg.c_str());
    }

  return res.result.release();
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

namespace
{
  const int NARGIN = 4;
  const int NARGOUT = 1;
}

extern "C"
void mlxClassifier(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[])
{
  DOTRACE("mlxClassifier");

  mxArray* mprhs[NARGIN];

  if (nlhs != NARGOUT)
    {
      mexErrMsgTxt("Error: classifier was called with the wrong "
                   "number of outputs (should be 1).");
    }

  if (nrhs != NARGIN)
    {
      mexErrMsgTxt("Error: classifier was called with the wrong "
                   "number of inputs (should be 4).");
    }

  for (int i = 0; i < NARGIN; ++i)
    {
      mprhs[i] = prhs[i];
    }


  mlfEnterNewContext(0, NARGIN, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  try
    {
      plhs[0] = Mclassifier(Mtx(mprhs[0], Mtx::COPY), // model params
                            Mx::getString(mprhs[1]), // model name
                            Mx::getString(mprhs[2]), // action request
                            mprhs[3]); // extra args struct
    }
  catch (Util::Error& err)
    {
      mexErrMsgTxt(err.msg_cstr());
    }
  catch (std::exception& err)
    {
      mexErrMsgTxt(err.what());
    }
  catch (...)
    {
      mexErrMsgTxt("an unknown C++ exception occurred.");
    }

  mlfRestorePreviousContext(0, NARGIN, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);
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

  if (mexPkg == 0)
    mexPkg = new MexPkg;

  static _mexInitTermTableEntry init_term_table[1] =
    { { InitializeModule_classifier, TerminateModule_classifier }, };

  static mexFunctionTableEntry function_table[1] =
    {
      { MEXFUNCNAME, mlxClassifier, NARGIN, NARGOUT,
        (_mexLocalFunctionTable*)0 }
    };

  static _mex_information _mex_info
    = { 1, 1, function_table, 0, NULL, 0, NULL, 1, init_term_table };

  return &_mex_info;
}

static const char vcid_classifier_mex_cc[] = "$Header$";
#endif // !CLASSIFIER_MEX_CC_DEFINED

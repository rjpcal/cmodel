///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Thu Mar  8 09:49:21 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_MEX_CC_DEFINED
#define CLASSIFIER_MEX_CC_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "cmodel/classifier.h"
#include "cmodel/cmodelcssm.h"
#include "cmodel/cmodelcuevalidity.h"
#include "cmodel/cmodelgcm.h"
#include "cmodel/cmodelpbi.h"
#include "cmodel/cmodelrxm.h"
#include "cmodel/cmodelspc.h"
#include "cmodel/cmodelwpsm.h"

#include "mtx/mtx.h"

#include "mx/mexpkg.h"
#include "mx/mx.h"

#include "util/error.h"
#include "util/pointers.h"
#include "util/strings.h"

#include <mex.h>

#include "util/trace.h"

#ifdef LOCAL_DEBUG
#define MEXFUNCNAME "dclassifier"
#else
#define MEXFUNCNAME "classifier"
#endif

namespace
{
  class MyMexPkg : public MexPkg
  {
  public:
    MyMexPkg(ExitFcn f) :
      MexPkg(MEXFUNCNAME, f),
      recentCssm(0),
      recentRxm(0)
    {}

    virtual ~MyMexPkg() {}

    shared_ptr<CModelCssm> recentCssm;
    shared_ptr<CModelRxm> recentRxm;
  };

  MyMexPkg* mexPkg = 0;
}

shared_ptr<Classifier> makeClassifier(const fstring& whichType,
                                      const mtx& objParams,
                                      const mxArray* extraArgs_mx)
{
  DOTRACE("<classifier_mex.cc>::makeClassifier");
  if (whichType == "cssm" || whichType == "rxmlin")
    {
      const CModelExemplar::TransferFunction tfunc =
        whichType == "cssm"
        ? CModelExemplar::EXP_DECAY
        : CModelExemplar::LINEAR_DECAY;

      const int numStoredExemplars =
        Mx::get_int_field(extraArgs_mx, "numStoredExemplars");

      if ( mexPkg->recentCssm.get() != 0 &&
           numStoredExemplars == mexPkg->recentCssm->numStoredExemplars() &&
           objParams == mexPkg->recentCssm->objParams() &&
           tfunc == mexPkg->recentCssm->transferFunction() )
        {
          DOTRACE("<classifier_mex.cc>::makeClassifier-use old cssm");

          return mexPkg->recentCssm;
        }
      else
        {
          DOTRACE("<classifier_mex.cc>::makeClassifier-make new cssm");

          // To avoid relying on transient matlab storage:
          mtx uniqObjParams = objParams;
          uniqObjParams.make_unique();

          mexPkg->recentCssm.reset
            (new CModelCssm(uniqObjParams,
                            tfunc,
                            numStoredExemplars));

          return mexPkg->recentCssm;
        }
    }
  else if (whichType == "rxm" || whichType == "rxme")
    {
      const CModelExemplar::TransferFunction tfunc =
        whichType == "rxm"
        ? CModelExemplar::LINEAR_DECAY
        : CModelExemplar::EXP_DECAY;

      const int numStoredExemplars =
        Mx::get_int_field(extraArgs_mx, "numStoredExemplars");

      if ( mexPkg->recentRxm.get() != 0 &&
           numStoredExemplars == mexPkg->recentRxm->numStoredExemplars() &&
           objParams == mexPkg->recentRxm->objParams() &&
           tfunc == mexPkg->recentRxm->transferFunction() )
        {
          DOTRACE("<classifier_mex.cc>::makeClassifier-use old rxm");

          return mexPkg->recentRxm;
        }
      else
        {
          DOTRACE("<classifier_mex.cc>::makeClassifier-make new rxm");

          // To avoid relying on transient matlab storage:
          mtx uniqObjParams = objParams;
          uniqObjParams.make_unique();

          mexPkg->recentRxm.reset
            (new CModelRxm(uniqObjParams,
                           tfunc,
                           numStoredExemplars));

          return mexPkg->recentRxm;
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
        Mx::get_int_field(extraArgs_mx, "numStoredExemplars");

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


///////////////////////////////////////////////////////////////////////
//
// classifier()
//
///////////////////////////////////////////////////////////////////////

namespace
{
  const int NARGIN = 4;
  const int NARGOUT = 1;
}

void classifier(int nlhs, mxArray* plhs[],
                int nrhs, const mxArray* prhs[])
{
DOTRACE("<classifier_mex.cc>::classifier");

  const mtx      modelParams   (prhs[0], mtx::COPY);
  const fstring  modelName     (Mx::as_string(prhs[1]));
  const fstring  actionRequest (Mx::as_string(prhs[2]));
  const mxArray* extraArgs_mx  (prhs[3]);

  const mtx      objParams     (Mx::get_field(extraArgs_mx, "objParams"));

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

  plhs[0] = res.result.release();
}


///////////////////////////////////////////////////////////////////////
//
// mexLibrary()
//
///////////////////////////////////////////////////////////////////////

namespace
{
  void terminateModule()
  {
    delete mexPkg;
    mexPkg = 0;
  }
}

extern "C"
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
DOTRACE("mexFunction");

  if (mexPkg == 0)
    {
      mexPkg = new MyMexPkg(&terminateModule);
      mexPkg->addFcn(MEXFUNCNAME, &classifier, NARGIN, NARGOUT);
    }

  mexPkg->invokeFcn(nlhs, plhs, nrhs, prhs);
}

static const char vcid_classifier_mex_cc[] = "$Header$";
#endif // !CLASSIFIER_MEX_CC_DEFINED

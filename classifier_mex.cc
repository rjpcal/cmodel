///////////////////////////////////////////////////////////////////////
//
// classifier_mex.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:49:21 2001
// written: Wed Jul 31 15:21:59 2002
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

#include "mx/mexfcn.h"
#include "mx/mexpkg.h"
#include "mx/mx.h"

#include "util/error.h"
#include "util/pointers.h"
#include "util/strings.h"

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
                                      const Mtx& objParams,
                                      const mxArray* extraArgs_mx)
{
  DOTRACE("makeClassifier");
  if (whichType == "cssm" || whichType == "rxmlin")
    {
      const CModelExemplar::TransferFunction tfunc =
        whichType == "cssm"
        ? CModelExemplar::EXP_DECAY
        : CModelExemplar::LINEAR_DECAY;

      const int numStoredExemplars =
        Mx::getIntField(extraArgs_mx, "numStoredExemplars");

      if ( mexPkg->recentCssm.get() != 0 &&
           numStoredExemplars == mexPkg->recentCssm->numStoredExemplars() &&
           objParams == mexPkg->recentCssm->objParams() &&
           tfunc == mexPkg->recentCssm->transferFunction() )
        {
          DOTRACE("use old cssm");

          return mexPkg->recentCssm;
        }
      else
        {
          DOTRACE("make new cssm");

          // To avoid relying on transient matlab storage:
          Mtx uniqObjParams = objParams;
          uniqObjParams.makeUnique();

          mexPkg->recentCssm.reset
            (new CModelCssm(uniqObjParams,
                            tfunc,
                            numStoredExemplars));

          return mexPkg->recentCssm;
        }
    }
  else if (whichType == "rxm")
    {
      const CModelExemplar::TransferFunction tfunc =
        CModelExemplar::LINEAR_DECAY;

      const int numStoredExemplars =
        Mx::getIntField(extraArgs_mx, "numStoredExemplars");

      if ( mexPkg->recentRxm.get() != 0 &&
           numStoredExemplars == mexPkg->recentRxm->numStoredExemplars() &&
           objParams == mexPkg->recentRxm->objParams() &&
           tfunc == mexPkg->recentRxm->transferFunction() )
        {
          DOTRACE("use old rxm");

          return mexPkg->recentRxm;
        }
      else
        {
          DOTRACE("make new rxm");

          // To avoid relying on transient matlab storage:
          Mtx uniqObjParams = objParams;
          uniqObjParams.makeUnique();

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

void classifier(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[])
{
  const Mtx      modelParams   (prhs[0], Mtx::COPY);
  const fstring  modelName     (Mx::getString(prhs[1]));
  const fstring  actionRequest (Mx::getString(prhs[2]));
  const mxArray* extraArgs_mx  (prhs[3]);

  const Mtx      objParams     (Mx::getField(extraArgs_mx, "objParams"));

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
  void terminateModule() { delete mexPkg; }
}

extern "C"
mex_information mexLibrary()
{
  DOTRACE("mexLibrary");

  if (mexPkg == 0)
    {
      mexPkg = new MyMexPkg(terminateModule);
      mexPkg->addFcn(MexFcn(MEXFUNCNAME, classifier, NARGIN, NARGOUT));
    }

  return mexPkg->mexInfo();
}

static const char vcid_classifier_mex_cc[] = "$Header$";
#endif // !CLASSIFIER_MEX_CC_DEFINED

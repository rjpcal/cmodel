///////////////////////////////////////////////////////////////////////
//
// annealVisitParameters.cc
//
// Copyright (c) 1998-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb  4 18:14:46 2002
// written: Thu Feb  7 13:54:49 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALVISITPARAMETERS_CC_DEFINED
#define ANNEALVISITPARAMETERS_CC_DEFINED

/*
 * MATLAB Compiler: 2.1
 * Date: Fri Mar 23 17:17:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "-h" "annealVisitParameters"
 */
#include "annealVisitParameters.h"

#include "util/error.h"
#include "mexbuf.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "rutil.h"
#include "util/strings.h"

#include <iostream.h>
#include "libmatlbm.h"

#include "util/trace.h"
#include "util/debug.h"

static mxArray* _mxarray20_;
static mxArray* _mxarray21_;
static mxArray* _mxarray22_;
static mxArray* _mxarray23_;

static mxArray* _mxarray26_;

namespace
{
  MexBuf* mexBuf = 0;
}

void InitializeModule_annealVisitParameters()
{
  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  cout = mexBuf;
  cerr = mexBuf;
#else
  cout.rdbuf(mexBuf);
  cerr.rdbuf(mexBuf);
#endif

  _mxarray20_ = mclInitializeDouble(0.0);
  _mxarray21_ = mclInitializeDouble(1.0);
  _mxarray22_ = mclInitializeDouble(2.0);
  _mxarray23_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
  _mxarray26_ = mclInitializeDouble(708.3964185322641);
}

void TerminateModule_annealVisitParameters()
{
  delete mexBuf;

  mxDestroyArray(_mxarray26_);
  mxDestroyArray(_mxarray23_);
  mxDestroyArray(_mxarray22_);
  mxDestroyArray(_mxarray21_);
  mxDestroyArray(_mxarray20_);
}

mxArray* MannealVisitParameters(int nargout_,
                                mxArray* bestModel,
                                mxArray* valueScalingRange,
                                mxArray* deltas,
                                mxArray* bounds,
                                mxArray* canUseMatrix,
                                mxArray* FUN,
                                mxArray* temp,
                                mxArray* varargin);

Mtx makeTestModels(int x_zerobased,
                   const Mtx& bestModel,
                   const Mtx& valueScalingRange,
                   const Mtx& deltas,
                   const Mtx& bounds);

mxArray* doFuncEvals(bool canUseMatrix,
                     const Mtx& models,
                     mxArray* func,
                     mxArray* varargin);

int sampleFromPdf_zerobased(mxArray* temp, mxArray* costs);


mxArray* makePDF(mxArray* temp,
                 mxArray* costs);

mxArray* eprob(mxArray* temp,
               mxArray* costs);

static mexFunctionTableEntry local_function_table_[5]
  = { { "makeTestModels",
        NULL /*mlxmakeTestModels*/, 5, 1, NULL },
      { "doFuncEvals", NULL/*mlxdoFuncEvals*/, 2, 1, NULL },
      { "sampleFromPdf", NULL/*mlxsampleFromPdf*/, 2, 1, NULL },
      { "makePDF", NULL/*mlxmakePDF*/, 2, 1, NULL },
      { "eprob", NULL/*mlxeprob*/, 2, 1, NULL } };

_mexLocalFunctionTable _local_function_table_annealVisitParameters
  = { 5, local_function_table_ };

/*
 * The function "mlfAnnealVisitParameters" contains the normal interface for
 * the "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray* mlfAnnealVisitParameters(mxArray* bestModel,
                                   mxArray* valueScalingRange,
                                   mxArray* deltas,
                                   mxArray* bounds,
                                   mxArray* canUseMatrix,
                                   mxArray* FUN,
                                   mxArray* temp,
                                   ...)
{
DOTRACE("mlfAnnealVisitParameters");
    mxArray* varargin = NULL;
    int nargout = 1;
    mxArray* S = mclGetUninitializedArray();
    mlfVarargin(&varargin, temp, 0);
    mlfEnterNewContext(
      0,
      -8,
      bestModel,
      valueScalingRange,
      deltas,
      bounds,
      canUseMatrix,
      FUN,
      temp,
      varargin);
    S
      = MannealVisitParameters(
          nargout,
          bestModel,
          valueScalingRange,
          deltas,
          bounds,
          canUseMatrix,
          FUN,
          temp,
          varargin);
    mlfRestorePreviousContext(
      0,
      7,
      bestModel,
      valueScalingRange,
      deltas,
      bounds,
      canUseMatrix,
      FUN,
      temp);
    mxDestroyArray(varargin);
    return mlfReturnValue(S);
}

/*
 * The function "mlxAnnealVisitParameters" contains the feval interface for the
 * "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). The feval function calls the implementation version of
 * annealVisitParameters through this function. This function processes any
 * input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
void mlxAnnealVisitParameters(int nlhs,
                              mxArray* plhs[],
                              int nrhs,
                              mxArray* prhs[])
{
DOTRACE("mlxAnnealVisitParameters");
    mxArray* mprhs[8];
    mxArray* mplhs[1];
    int i;
    if (nlhs > 1) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters Line: 1 "
                     "Column: 1 The function \"annealVisitParameters\" was "
                     "called with more than the declared number of outputs "
                     "(1).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 7 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 7; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0,
      7,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6]);
    mprhs[7] = NULL;
    mlfAssign(&mprhs[7], mclCreateVararginCell(nrhs - 7, prhs + 7));
    mplhs[0]
      = MannealVisitParameters(
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4],
          mprhs[5],
          mprhs[6],
          mprhs[7]);

    mlfRestorePreviousContext(
      0,
      7,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6]);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[7]);
}

/*
 * The function "MannealVisitParameters" is the implementation version of the
 * "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function S = visitAllParameters(...
     * bestModel, valueScalingRange, deltas, bounds, ...
     * canUseMatrix, FUN, ...
     * temp, ...
     * varargin)
 */

mxArray* MannealVisitParameters(int nargout_,
                                mxArray* bestModel_mx,
                                mxArray* valueScalingRange_mx,
                                mxArray* deltas_mx,
                                mxArray* bounds_mx,
                                mxArray* canUseMatrix_mx,
                                mxArray* FUN_mx,
                                mxArray* temp_mx,
                                mxArray* varargin_mx)
{
DOTRACE("MannealVisitParameters");

  try
    {

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
      if (varargin_mx && mxGetScalar(mxGetCell(varargin_mx,0)) == -1)
        {
          Util::Prof::printAllProfData(cerr);
          return mxCreateScalarDouble(-1.0);
        }

      if (varargin_mx && mxGetScalar(mxGetCell(varargin_mx,0)) == -2)
        {
          Util::Prof::resetAllProfData();
          return mxCreateScalarDouble(-2.0);
        }
#endif

      mexLocalFunctionTable save_local_function_table_ =
        mclSetCurrentLocalFunctionTable(&_local_function_table_annealVisitParameters);

      mxArray* costs_mx = mclGetUninitializedArray();
      mxArray* x_mx = mclGetUninitializedArray();

      validateInput(bestModel_mx);
      mclCopyArray(&bestModel_mx);

      Mtx bestModel(bestModel_mx, Mtx::REFER);

      validateInput(valueScalingRange_mx);
      mclCopyArray(&valueScalingRange_mx);

      validateInput(deltas_mx);
      mclCopyArray(&deltas_mx);

      validateInput(bounds_mx);
      mclCopyArray(&bounds_mx);

      validateInput(canUseMatrix_mx);
      const bool canUseMatrix = (mxGetPr(canUseMatrix_mx)[0] != 0.0);

      validateInput(FUN_mx);
      mclCopyArray(&FUN_mx);

      validateInput(temp_mx);
      mclCopyArray(&temp_mx);

      validateInput(varargin_mx);
      mclCopyArray(&varargin_mx);

      int nevals = 0;
      int s_zerobased = 0;

      // for x = find(deltas' ~= 0)
      {
        mclForLoopIterator viter__;
        for (mclForStart(
                         &viter__,
                         mlfFind(
                                 NULL,
                                 NULL,
                                 mclNe(mlfCtranspose(deltas_mx), _mxarray20_)),
                         NULL,
                         NULL);
             mclForNext(&viter__, &x_mx);
             )
          {

            int x_zerobased = int(mxGetScalar(x_mx)) - 1;

            // modelmatrix = makeTestModels(x, bestModel, valueScalingRange, ...
            // deltas, bounds);
            const Mtx valueScalingRange(valueScalingRange_mx, Mtx::BORROW);
            const Mtx deltas(deltas_mx, Mtx::BORROW);
            const Mtx bounds(bounds_mx, Mtx::BORROW);
            Mtx modelmatrix = makeTestModels(x_zerobased,
                                             bestModel,
                                             valueScalingRange,
                                             deltas,
                                             bounds);

            // costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
            mlfAssign(&costs_mx,
                      doFuncEvals(canUseMatrix,
                                  modelmatrix,
                                  FUN_mx,
                                  mlfIndexRef(varargin_mx,
                                              "{?}",
                                              mlfCreateColonIndex())
                                  ));

            // S.nevals = S.nevals + length(costs);
            nevals += ( mxGetM(costs_mx) > mxGetN(costs_mx) ?
                        mxGetM(costs_mx) : mxGetN(costs_mx) );

            // Sample from probability distribution
            // s = sampleFromPdf(temp, costs);
            s_zerobased = sampleFromPdf_zerobased(temp_mx, costs_mx);

            // bestModel(x, 1) = modelmatrix(x, s);
            bestModel.at(x_zerobased, 0) = modelmatrix.at(x_zerobased, s_zerobased);
          }

        mclDestroyForLoopIterator(viter__);
      }

      const char* fieldNames[] = { "nevals", "newModel", "cost" };
      mxArray* output = mxCreateStructMatrix(1,1,3,fieldNames);

      // S.nevals = 0;
      mxSetField(output, 0, "nevals", mxCreateScalarDouble(nevals));

      // S.newModel = bestModel;
      mxSetField(output, 0, "newModel", bestModel_mx);

      // S.cost = costs(s);
      mxSetField(output, 0, "cost",
                 mxCreateScalarDouble(mxGetPr(costs_mx)[s_zerobased]));

      mxDestroyArray(x_mx);
      mxDestroyArray(costs_mx);
      mxDestroyArray(varargin_mx);
      mxDestroyArray(temp_mx);
      mxDestroyArray(FUN_mx);
      mxDestroyArray(bounds_mx);
      mxDestroyArray(deltas_mx);
      mxDestroyArray(valueScalingRange_mx);

      mclSetCurrentLocalFunctionTable(save_local_function_table_);

      return output;
    }
  catch (Util::Error& err)
    {
      mexErrMsgTxt(err.msg_cstr());
    }
  catch (...)
    {
      mexErrMsgTxt("an unknown C++ exception occurred.");
    }

  return (mxArray*) 0; // can't happen, but placate compiler
}

/*
 * The function "MmakeTestModels" is the implementation
 * version of the "annealVisitParameters/makeTestModels" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 28-35). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function models = makeTestModels(x, bestModel, valueScalingRange, deltas, bounds)
 */
Mtx makeTestModels(int x_zerobased,
                   const Mtx& bestModel,
                   const Mtx& valueScalingRange,
                   const Mtx& deltas,
                   const Mtx& bounds)
{
DOTRACE("makeTestModels");

  const double current_x_val = bestModel.at(x_zerobased);

  const double delta_x = deltas.at(x_zerobased);

  const double lowerbound_x = bounds.at(x_zerobased,0);
  const double upperbound_x = bounds.at(x_zerobased,1);

  int num_testmodels = 0;
  {for (MtxConstIter iter = valueScalingRange.rowIter(0); iter.hasMore(); ++iter)
    {
      double value = current_x_val + (*iter) * delta_x;
      if ( (value > lowerbound_x) && (value < upperbound_x) )
        ++num_testmodels;
    }
  }

  Mtx newModels(bestModel.mrows(), num_testmodels);

  const Slice bestModelCol = bestModel.column(0);

  for (int i = 0; i < num_testmodels; ++i)
    {
      newModels.column(i) = bestModelCol;
    }

  {
    int col = 0;
    for (MtxConstIter iter = valueScalingRange.rowIter(0); iter.hasMore(); ++iter)
    {
      double value = current_x_val + (*iter) * delta_x;
      if ( (value > lowerbound_x) && (value < upperbound_x) )
        newModels.at(x_zerobased, col++) = value;
    }
  }

  return newModels;
}

/*
 * The function "doFuncEvals" is the implementation
 * version of the "annealVisitParameters/doFuncEvals" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 35-48). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function costs = doFuncEvals(canUseMatrix, models, func, varargin)
 */

mxArray* doFuncEvals(bool canUseMatrix,
                     const Mtx& models,
                     mxArray* func,
                     mxArray* varargin_mx)
{
DOTRACE("doFuncEvals");

  mxArray* costs_mx = mclGetUninitializedArray();
  mxArray* e = mclGetUninitializedArray();
  mxArray* NM = mclGetUninitializedArray();
  mxArray* models_mx = mclGetUninitializedArray();
  mlfAssign(&models_mx, models.makeMxArray());
  mclCopyArray(&func);
  mclCopyArray(&varargin_mx);

  if (canUseMatrix)
    {
      DOTRACE("eval w/ matrix");
      // costs = feval(func, models, varargin{:});
#if 0
      mlfAssign(&costs_mx,
                mlfFeval(mclValueVarargout(),
                         func,
                         models_mx,
                         mlfIndexRef(varargin_mx, "{?}", mlfCreateColonIndex()),
                         NULL));

#else

      mxArray* plhs[1] = { 0 };

      const int MAX_NRHS = 32;
      mxArray* prhs[MAX_NRHS];

      prhs[0] = mxDuplicateArray(models_mx);

      int nrhs = 1;

      const int nvararg = mxGetNumberOfElements(varargin_mx);

      for (int i = 0; i < nvararg && nrhs < MAX_NRHS; ++i)
        {
          DebugEvalNL(i);
          prhs[nrhs++] = mxDuplicateArray(mxGetCell(varargin_mx, i));
        }

        fstring cmd_name = MxWrapper::extractString(func);

        int result = mexCallMATLAB(1, plhs, nrhs, prhs, cmd_name.c_str());

      if (result != 0) mexErrMsgTxt("mexCallMATLAB failed in doFuncEvals");

      mlfAssign(&costs_mx, plhs[0]);
#endif
    }
  else
    {
      DOTRACE("eval w/o matrix");
      // NM = size(models, 2);
      mlfAssign(&NM, mlfSize(mclValueVarargout(), models_mx, _mxarray22_));

      // costs = zeros(NM, 1);
      mlfAssign(&costs_mx, mlfZeros(NM, _mxarray21_, NULL));

      // for e = 1:NM
      {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(NM);
        if (v_ > e_)
          {
            mlfAssign(&e, _mxarray23_);
          }
        else
          {
            //costs(e) = feval(func, models(:,e), varargin{:});
            for (; ; )
              {
                mclIntArrayAssign1(
                                   &costs_mx,
                                   mlfFeval(
                                            mclValueVarargout(),
                                            func,
                                            mclArrayRef2(
                                                         models_mx,
                                                         mlfCreateColonIndex(),
                                                         mlfScalar(v_)),
                                            mlfIndexRef(
                                                        varargin_mx,
                                                        "{?}",
                                                        mlfCreateColonIndex()),
                                            NULL),
                                   v_);
                if (v_ == e_)
                  {
                    break;
                  }
                ++v_;
              }
            mlfAssign(&e, mlfScalar(v_));
          }
      }
    }

  mclValidateOutput(costs_mx, 1, 1, "costs", "annealVisitParameters/doFuncEvals");

  mxDestroyArray(NM);
  mxDestroyArray(e);
  mxDestroyArray(varargin_mx);
  mxDestroyArray(func);
  mxDestroyArray(models_mx);

  return costs_mx;
}

/*
 * The function "sampleFromPdf" is the implementation
 * version of the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function s = sampleFromPdf(temp, costs)
 */

int sampleFromPdf_zerobased(mxArray* temp_mx, mxArray* costs_mx)
{
DOTRACE("sampleFromPdf_zerobased");

  mexLocalFunctionTable save_local_function_table_ =
    mclSetCurrentLocalFunctionTable(&_local_function_table_annealVisitParameters);

  mxArray* s_mx = mclGetUninitializedArray();
  mxArray* cutoff = mclGetUninitializedArray();
  mxArray* dist = mclGetUninitializedArray();
  mclCopyArray(&temp_mx);
  mclCopyArray(&costs_mx);

  // dist = makePDF(temp, costs);
  mlfAssign(&dist, makePDF(temp_mx, costs_mx));

  // cutoff = rand;
  mlfAssign(&cutoff, mlfNRand(1, NULL));

  // s = find(cumsum(dist) >= cutoff);
  mlfAssign(&s_mx,
            mlfFind(NULL,
                    NULL,
                    mclGe(mlfCumsum(dist, NULL), cutoff)));

  // s = s(1);
  int s_zerobased = int(mxGetPr(s_mx)[0]) - 1;

  mxDestroyArray(dist);
  mxDestroyArray(cutoff);
  mxDestroyArray(costs_mx);
  mxDestroyArray(temp_mx);

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  mxDestroyArray(s_mx);

  return s_zerobased;
}

/*
 * The function "makePDF" is the implementation version
 * of the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = makePDF(temp, costs)
 */
mxArray* makePDF(mxArray* temp_mx, mxArray* costs_mx)
{
DOTRACE("makePDF");

  mexLocalFunctionTable save_local_function_table_ =
    mclSetCurrentLocalFunctionTable(&_local_function_table_annealVisitParameters);

  mxArray* pdf = mclGetUninitializedArray();
  mxArray* ans = mclGetUninitializedArray();
  mxArray* w = mclGetUninitializedArray();
  mxArray* good = mclGetUninitializedArray();
  mxArray* bad = mclGetUninitializedArray();
  mclCopyArray(&temp_mx);
  mclCopyArray(&costs_mx);

  // Forms exponential probability distribution given a temperature and vector of
  // costs.  Internal function for simulated annealing algorithm.

  // bad = find(isnan(costs));
  mlfAssign(&bad, mlfFind(NULL, NULL, mlfIsnan(costs_mx)));

  // if isempty(bad)
  if (mlfTobool(mlfIsempty(bad)))
    {
      mlfAssign(&pdf, eprob(temp_mx, costs_mx));
    }
  else
    {
      // good = find(~isnan(costs));
      mlfAssign(&good, mlfFind(NULL, NULL, mclNot(mlfIsnan(costs_mx))));

      // w = costs(good);
      mlfAssign(&w, mclArrayRef1(costs_mx, good));

      mlfAssign(&pdf, eprob(temp_mx, w));

      // pdf(good) = pdf;
      mclArrayAssign1(&pdf, pdf, good);

      // pdf(bad) = 0;
      mclArrayAssign1(&pdf, _mxarray20_, bad);

      cerr << "Warning: the cost function generated a NaN "
           << "for one or more models.";
    }

  mclValidateOutput(pdf, 1, 1, "pdf", "annealVisitParameters/makePDF");

  mxDestroyArray(bad);
  mxDestroyArray(good);
  mxDestroyArray(w);
  mxDestroyArray(ans);
  mxDestroyArray(costs_mx);
  mxDestroyArray(temp_mx);

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  return pdf;
}

/*
 * The function "eprob" is the implementation version of
 * the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = eprob(temp, costs)
 */
mxArray* eprob(mxArray* temp_mx, mxArray* costs_mx)
{
DOTRACE("eprob");

  mexLocalFunctionTable save_local_function_table_ =
    mclSetCurrentLocalFunctionTable(&_local_function_table_annealVisitParameters);

  mxArray* pdf = mclGetUninitializedArray();
  mxArray* scale = mclGetUninitializedArray();
  mxArray* mpdf = mclGetUninitializedArray();
  mxArray* toobig = mclGetUninitializedArray();
  mclCopyArray(&temp_mx);
  mclCopyArray(&costs_mx);

  // Scales cost vector and calculates exponential probability distribution.  The
  // scaling is necessary to permit wide ranges in temperature.  Internal
  // function for simulated annealing algorithm.

  // toobig = 708.3964185322641;
  mlfAssign(&toobig, _mxarray26_);

  // pdf = costs/temp;
  mlfAssign(&pdf, mclMrdivide(costs_mx, temp_mx));

  // mpdf = max(pdf);
  mlfAssign(&mpdf, mlfMax(NULL, pdf, NULL, NULL));

  // if mpdf>toobig
  if (mclGtBool(mpdf, toobig))
    {
      // scale = mpdf/toobig;
      mlfAssign(&scale, mclMrdivide(mpdf, toobig));

      // pdf = exp(-pdf/scale);
      mlfAssign(&pdf, mlfExp(mclMrdivide(mclUminus(pdf), scale)));

      // pdf = pdf/max(pdf);
      mlfAssign(&pdf, mclMrdivide(pdf, mlfMax(NULL, pdf, NULL, NULL)));

      // pdf = pdf.^scale;
      mlfAssign(&pdf, mlfPower(pdf, scale));
    }
  else
    {
      // pdf = exp(-pdf);
      mlfAssign(&pdf, mlfExp(mclUminus(pdf)));

      // pdf = pdf/max(pdf);
      mlfAssign(&pdf, mclMrdivide(pdf, mlfMax(NULL, pdf, NULL, NULL)));
    }

  // pdf = pdf/sum(pdf);
  mlfAssign(&pdf, mclMrdivide(pdf, mlfSum(pdf, NULL)));

  mclValidateOutput(pdf, 1, 1, "pdf", "annealVisitParameters/eprob");

  mxDestroyArray(toobig);
  mxDestroyArray(mpdf);
  mxDestroyArray(scale);
  mxDestroyArray(costs_mx);
  mxDestroyArray(temp_mx);

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  return pdf;
}

static const char vcid_annealVisitParameters_cc[] = "$Header$";
#endif // !ANNEALVISITPARAMETERS_CC_DEFINED

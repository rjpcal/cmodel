/*
 * MATLAB Compiler: 2.1
 * Date: Fri Mar 23 17:17:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "-h" "annealVisitParameters" 
 */
#include "annealVisitParameters.h"
#include "libmatlbm.h"

static mxArray * _mxarray20_;
static mxArray * _mxarray21_;
static mxArray * _mxarray22_;
static mxArray * _mxarray23_;

static mxArray * _mxarray26_;

void InitializeModule_annealVisitParameters(void) {
    _mxarray20_ = mclInitializeDouble(0.0);
    _mxarray21_ = mclInitializeDouble(1.0);
    _mxarray22_ = mclInitializeDouble(2.0);
    _mxarray23_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray26_ = mclInitializeDouble(708.3964185322641);
}

void TerminateModule_annealVisitParameters(void) {
    mxDestroyArray(_mxarray26_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray20_);
}

static mxArray * mlfAnnealVisitParameters_makeTestModels(mxArray * x,
                                                         mxArray * bestModel,
                                                         mxArray * valueScalingRange,
                                                         mxArray * deltas,
                                                         mxArray * bounds);
static void mlxAnnealVisitParameters_makeTestModels(int nlhs,
                                                    mxArray * plhs[],
                                                    int nrhs,
                                                    mxArray * prhs[]);
static mxArray * mlfAnnealVisitParameters_doFuncEvals(mxArray * canUseMatrix,
                                                      mxArray * models,
                                                      mxArray * func,
                                                      ...);
static void mlxAnnealVisitParameters_doFuncEvals(int nlhs,
                                                 mxArray * plhs[],
                                                 int nrhs,
                                                 mxArray * prhs[]);
static mxArray * mlfAnnealVisitParameters_sampleFromPdf(mxArray * temp,
                                                        mxArray * costs);
static void mlxAnnealVisitParameters_sampleFromPdf(int nlhs,
                                                   mxArray * plhs[],
                                                   int nrhs,
                                                   mxArray * prhs[]);
static mxArray * mlfAnnealVisitParameters_makePDF(mxArray * temp,
                                                  mxArray * costs);
static void mlxAnnealVisitParameters_makePDF(int nlhs,
                                             mxArray * plhs[],
                                             int nrhs,
                                             mxArray * prhs[]);
static mxArray * mlfAnnealVisitParameters_eprob(mxArray * temp,
                                                mxArray * costs);
static void mlxAnnealVisitParameters_eprob(int nlhs,
                                           mxArray * plhs[],
                                           int nrhs,
                                           mxArray * prhs[]);
static mxArray * MannealVisitParameters(int nargout_,
                                        mxArray * bestModel,
                                        mxArray * valueScalingRange,
                                        mxArray * deltas,
                                        mxArray * bounds,
                                        mxArray * canUseMatrix,
                                        mxArray * FUN,
                                        mxArray * temp,
                                        mxArray * varargin);
static mxArray * MannealVisitParameters_makeTestModels(int nargout_,
                                                       mxArray * x,
                                                       mxArray * bestModel,
                                                       mxArray * valueScalingRange,
                                                       mxArray * deltas,
                                                       mxArray * bounds);
static mxArray * MannealVisitParameters_doFuncEvals(int nargout_,
                                                    mxArray * canUseMatrix,
                                                    mxArray * models,
                                                    mxArray * func,
                                                    mxArray * varargin);
static mxArray * MannealVisitParameters_sampleFromPdf(int nargout_,
                                                      mxArray * temp,
                                                      mxArray * costs);
static mxArray * MannealVisitParameters_makePDF(int nargout_,
                                                mxArray * temp,
                                                mxArray * costs);
static mxArray * MannealVisitParameters_eprob(int nargout_,
                                              mxArray * temp,
                                              mxArray * costs);

static mexFunctionTableEntry local_function_table_[5]
  = { { "makeTestModels",
        mlxAnnealVisitParameters_makeTestModels, 5, 1, NULL },
      { "doFuncEvals", mlxAnnealVisitParameters_doFuncEvals, -4, 1, NULL },
      { "sampleFromPdf", mlxAnnealVisitParameters_sampleFromPdf, 2, 1, NULL },
      { "makePDF", mlxAnnealVisitParameters_makePDF, 2, 1, NULL },
      { "eprob", mlxAnnealVisitParameters_eprob, 2, 1, NULL } };

_mexLocalFunctionTable _local_function_table_annealVisitParameters
  = { 5, local_function_table_ };

/*
 * The function "mlfAnnealVisitParameters" contains the normal interface for
 * the "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). This function processes any input arguments and passes them to
 * the implementation version of the function, appearing above.
 */
mxArray * mlfAnnealVisitParameters(mxArray * bestModel,
                                   mxArray * valueScalingRange,
                                   mxArray * deltas,
                                   mxArray * bounds,
                                   mxArray * canUseMatrix,
                                   mxArray * FUN,
                                   mxArray * temp,
                                   ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * S = mclGetUninitializedArray();
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
                              mxArray * plhs[],
                              int nrhs,
                              mxArray * prhs[]) {
    mxArray * mprhs[8];
    mxArray * mplhs[1];
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
 * The function "mlfAnnealVisitParameters_makeTestModels" contains the normal
 * interface for the "annealVisitParameters/makeTestModels" M-function from
 * file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 28-35). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfAnnealVisitParameters_makeTestModels(mxArray * x,
                                                         mxArray * bestModel,
                                                         mxArray * valueScalingRange,
                                                         mxArray * deltas,
                                                         mxArray * bounds) {
    int nargout = 1;
    mxArray * models = mclGetUninitializedArray();
    mlfEnterNewContext(0, 5, x, bestModel, valueScalingRange, deltas, bounds);
    models
      = MannealVisitParameters_makeTestModels(
          nargout, x, bestModel, valueScalingRange, deltas, bounds);
    mlfRestorePreviousContext(
      0, 5, x, bestModel, valueScalingRange, deltas, bounds);
    return mlfReturnValue(models);
}

/*
 * The function "mlxAnnealVisitParameters_makeTestModels" contains the feval
 * interface for the "annealVisitParameters/makeTestModels" M-function from
 * file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 28-35). The feval function calls the implementation version of
 * annealVisitParameters/makeTestModels through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
static void mlxAnnealVisitParameters_makeTestModels(int nlhs,
                                                    mxArray * plhs[],
                                                    int nrhs,
                                                    mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/makeTestModels "
							"Line: 28 Column: 1 The function "
							"\"annealVisitParameters/makeTestModels\" was called with "
							"more than the declared number of outputs (1).");
    }
    if (nrhs > 5) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/makeTestModels "
							"Line: 28 Column: 1 The function "
							"\"annealVisitParameters/makeTestModels\" was called with "
							"more than the declared number of inputs (5).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 5 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 5; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mplhs[0]
      = MannealVisitParameters_makeTestModels(
          nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mlfRestorePreviousContext(
      0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfAnnealVisitParameters_doFuncEvals" contains the normal
 * interface for the "annealVisitParameters/doFuncEvals" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 35-48). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfAnnealVisitParameters_doFuncEvals(mxArray * canUseMatrix,
                                                      mxArray * models,
                                                      mxArray * func,
                                                      ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * costs = mclGetUninitializedArray();
    mlfVarargin(&varargin, func, 0);
    mlfEnterNewContext(0, -4, canUseMatrix, models, func, varargin);
    costs
      = MannealVisitParameters_doFuncEvals(
          nargout, canUseMatrix, models, func, varargin);
    mlfRestorePreviousContext(0, 3, canUseMatrix, models, func);
    mxDestroyArray(varargin);
    return mlfReturnValue(costs);
}

/*
 * The function "mlxAnnealVisitParameters_doFuncEvals" contains the feval
 * interface for the "annealVisitParameters/doFuncEvals" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 35-48). The feval function calls the implementation version of
 * annealVisitParameters/doFuncEvals through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
static void mlxAnnealVisitParameters_doFuncEvals(int nlhs,
                                                 mxArray * plhs[],
                                                 int nrhs,
                                                 mxArray * prhs[]) {
    mxArray * mprhs[4];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
		mexErrMsgTxt("Run-time Error: File: annealVisitParameters/doFuncEvals "
						 "Line: 35 Column: 1 The function "
						 "\"annealVisitParameters/doFuncEvals\" was called with "
						 "more than the declared number of outputs (1).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mprhs[3] = NULL;
    mlfAssign(&mprhs[3], mclCreateVararginCell(nrhs - 3, prhs + 3));
    mplhs[0]
      = MannealVisitParameters_doFuncEvals(
          nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
    mxDestroyArray(mprhs[3]);
}

/*
 * The function "mlfAnnealVisitParameters_sampleFromPdf" contains the normal
 * interface for the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfAnnealVisitParameters_sampleFromPdf(mxArray * temp,
                                                        mxArray * costs) {
    int nargout = 1;
    mxArray * s = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    s = MannealVisitParameters_sampleFromPdf(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(s);
}

/*
 * The function "mlxAnnealVisitParameters_sampleFromPdf" contains the feval
 * interface for the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). The feval function calls the implementation version of
 * annealVisitParameters/sampleFromPdf through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
static void mlxAnnealVisitParameters_sampleFromPdf(int nlhs,
                                                   mxArray * plhs[],
                                                   int nrhs,
                                                   mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/sampleFromPdf "
							"Line: 48 Column: 1 The function "
							"\"annealVisitParameters/sampleFromPdf\" was called with "
							"more than the declared number of outputs (1).");
    }
    if (nrhs > 2) {
        mexErrMsgTxt("Run-time  Error: File: annealVisitParameters/sampleFromPdf "
							"Line: 48 Column: 1 The function "
							"\"annealVisitParameters/sampleFromPdf\" was called with "
							"more than the declared number of inputs (2).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = MannealVisitParameters_sampleFromPdf(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfAnnealVisitParameters_makePDF" contains the normal
 * interface for the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfAnnealVisitParameters_makePDF(mxArray * temp,
                                                  mxArray * costs) {
    int nargout = 1;
    mxArray * pdf = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    pdf = MannealVisitParameters_makePDF(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(pdf);
}

/*
 * The function "mlxAnnealVisitParameters_makePDF" contains the feval interface
 * for the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). The feval function calls the implementation version of
 * annealVisitParameters/makePDF through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
static void mlxAnnealVisitParameters_makePDF(int nlhs,
                                             mxArray * plhs[],
                                             int nrhs,
                                             mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/makePDF "
							"Line: 55 Column: 1 The function "
							"\"annealVisitParameters/makePDF\" was called with "
							"more than the declared number of outputs (1).");
    }
    if (nrhs > 2) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/makePDF "
							"Line: 55 Column: 1 The function "
							"\"annealVisitParameters/makePDF\" was called with "
							"more than the declared number of inputs (2).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = MannealVisitParameters_makePDF(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfAnnealVisitParameters_eprob" contains the normal interface
 * for the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfAnnealVisitParameters_eprob(mxArray * temp,
                                                mxArray * costs) {
    int nargout = 1;
    mxArray * pdf = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    pdf = MannealVisitParameters_eprob(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(pdf);
}

/*
 * The function "mlxAnnealVisitParameters_eprob" contains the feval interface
 * for the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). The feval function calls the implementation version of
 * annealVisitParameters/eprob through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
static void mlxAnnealVisitParameters_eprob(int nlhs,
                                           mxArray * plhs[],
                                           int nrhs,
                                           mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mexErrMsgTxt("Run-time  Error: File: annealVisitParameters/eprob "
							"Line: 72 Column: 1 The function "
							"\"annealVisitParameters/eprob\" was called with "
							"more than the declared number of outputs (1).");
    }
    if (nrhs > 2) {
        mexErrMsgTxt("Run-time Error: File: annealVisitParameters/eprob "
							"Line: 72 Column: 1 The function "
							"\"annealVisitParameters/eprob\" was called with "
							"more than the declared number of inputs (2).");
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = MannealVisitParameters_eprob(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
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
 */
static mxArray * MannealVisitParameters(int nargout_,
                                        mxArray * bestModel,
                                        mxArray * valueScalingRange,
                                        mxArray * deltas,
                                        mxArray * bounds,
                                        mxArray * canUseMatrix,
                                        mxArray * FUN,
                                        mxArray * temp,
                                        mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * S = mclGetUninitializedArray();
    mxArray * s = mclGetUninitializedArray();
    mxArray * costs = mclGetUninitializedArray();
    mxArray * modelmatrix = mclGetUninitializedArray();
    mxArray * x = mclGetUninitializedArray();
    mclCopyArray(&bestModel);
    mclCopyArray(&valueScalingRange);
    mclCopyArray(&deltas);
    mclCopyArray(&bounds);
    mclCopyArray(&canUseMatrix);
    mclCopyArray(&FUN);
    mclCopyArray(&temp);
    mclCopyArray(&varargin);
    /*
     * bestModel, valueScalingRange, deltas, bounds, ...
     * canUseMatrix, FUN, ...
     * temp, ...
     * varargin)
     * 
     * S.nevals = 0;
     */
    mlfIndexAssign(&S, ".nevals", _mxarray20_);
    /*
     * 
     * for x = find(deltas' ~= 0)
     */
    {
        mclForLoopIterator viter__;
        for (mclForStart(
               &viter__,
               mclVe(
                 mlfFind(
                   NULL,
                   NULL,
                   mclNe(mlfCtranspose(mclVa(deltas, "deltas")), _mxarray20_))),
               NULL,
               NULL);
             mclForNext(&viter__, &x);
             ) {
            /*
             * 
             * modelmatrix = makeTestModels(x, bestModel, valueScalingRange, ...
             */
            mlfAssign(
              &modelmatrix,
              mlfAnnealVisitParameters_makeTestModels(
                mclVv(x, "x"),
                mclVa(bestModel, "bestModel"),
                mclVa(valueScalingRange, "valueScalingRange"),
                mclVa(deltas, "deltas"),
                mclVa(bounds, "bounds")));
            /*
             * deltas, bounds);
             * 
             * costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
             */
            mlfAssign(
              &costs,
              mlfAnnealVisitParameters_doFuncEvals(
                mclVa(canUseMatrix, "canUseMatrix"),
                mclVv(modelmatrix, "modelmatrix"),
                mclVa(FUN, "FUN"),
                mclVe(
                  mlfIndexRef(
                    mclVsa(varargin, "varargin"),
                    "{?}",
                    mlfCreateColonIndex())),
                NULL));
            /*
             * S.nevals = S.nevals + length(costs);
             */
            mlfIndexAssign(
              &S,
              ".nevals",
              mclFeval(
                mclValueVarargout(),
                mlxPlus,
                mclVe(mlfIndexRef(mclVsv(S, "S"), ".nevals")),
                mlfScalar(mclLengthInt(mclVv(costs, "costs"))),
                NULL));
            /*
             * 
             * % Sample from probability distribution
             * s = sampleFromPdf(temp, costs);
             */
            mlfAssign(
              &s,
              mlfAnnealVisitParameters_sampleFromPdf(
                mclVa(temp, "temp"), mclVv(costs, "costs")));
            /*
             * 
             * bestModel(x, 1) = modelmatrix(x, s);
             */
            mclArrayAssign2(
              &bestModel,
              mclArrayRef2(
                mclVsv(modelmatrix, "modelmatrix"),
                mclVsv(x, "x"),
                mclVsv(s, "s")),
              mclVsv(x, "x"),
              _mxarray21_);
        /*
         * end
         */
        }
        mclDestroyForLoopIterator(viter__);
    }
    /*
     * 
     * S.newModel = bestModel;
     */
    mlfIndexAssign(&S, ".newModel", mclVsa(bestModel, "bestModel"));
    /*
     * 
     * S.cost = costs(s);
     */
    mlfIndexAssign(
      &S, ".cost", mclArrayRef1(mclVsv(costs, "costs"), mclVsv(s, "s")));
    mclValidateOutput(S, 1, nargout_, "S", "annealVisitParameters");
    mxDestroyArray(x);
    mxDestroyArray(modelmatrix);
    mxDestroyArray(costs);
    mxDestroyArray(s);
    mxDestroyArray(varargin);
    mxDestroyArray(temp);
    mxDestroyArray(FUN);
    mxDestroyArray(canUseMatrix);
    mxDestroyArray(bounds);
    mxDestroyArray(deltas);
    mxDestroyArray(valueScalingRange);
    mxDestroyArray(bestModel);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return S;
    /*
     * 
     * 
     */
}

/*
 * The function "MannealVisitParameters_makeTestModels" is the implementation
 * version of the "annealVisitParameters/makeTestModels" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 28-35). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function models = makeTestModels(x, bestModel, valueScalingRange, deltas, bounds)
 */
static mxArray * MannealVisitParameters_makeTestModels(int nargout_,
                                                       mxArray * x,
                                                       mxArray * bestModel,
                                                       mxArray * valueScalingRange,
                                                       mxArray * deltas,
                                                       mxArray * bounds) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * models = mclGetUninitializedArray();
    mxArray * xv = mclGetUninitializedArray();
    mclCopyArray(&x);
    mclCopyArray(&bestModel);
    mclCopyArray(&valueScalingRange);
    mclCopyArray(&deltas);
    mclCopyArray(&bounds);
    /*
     * xv = bestModel(x) + valueScalingRange*deltas(x);
     */
    mlfAssign(
      &xv,
      mclPlus(
        mclVe(mclArrayRef1(mclVsa(bestModel, "bestModel"), mclVsa(x, "x"))),
        mclMtimes(
          mclVa(valueScalingRange, "valueScalingRange"),
          mclVe(mclArrayRef1(mclVsa(deltas, "deltas"), mclVsa(x, "x"))))));
    /*
     * xv = xv(find( (xv<=bounds(x,2)) & (xv>=bounds(x,1)) ));
     */
    mlfAssign(
      &xv,
      mclArrayRef1(
        mclVsv(xv, "xv"),
        mlfFind(
          NULL,
          NULL,
          mclAnd(
            mclLe(
              mclVv(xv, "xv"),
              mclVe(
                mclArrayRef2(
                  mclVsa(bounds, "bounds"), mclVsa(x, "x"), _mxarray22_))),
            mclGe(
              mclVv(xv, "xv"),
              mclVe(
                mclArrayRef2(
                  mclVsa(bounds, "bounds"), mclVsa(x, "x"), _mxarray21_)))))));
    /*
     * models = bestModel*ones(1, length(xv));
     */
    mlfAssign(
      &models,
      mclMtimes(
        mclVa(bestModel, "bestModel"),
        mclVe(
          mlfOnes(
            _mxarray21_, mlfScalar(mclLengthInt(mclVv(xv, "xv"))), NULL))));
    /*
     * models(x,:) = xv;
     */
    mclArrayAssign2(
      &models, mclVsv(xv, "xv"), mclVsa(x, "x"), mlfCreateColonIndex());
    mclValidateOutput(
      models, 1, nargout_, "models", "annealVisitParameters/makeTestModels");
    mxDestroyArray(xv);
    mxDestroyArray(bounds);
    mxDestroyArray(deltas);
    mxDestroyArray(valueScalingRange);
    mxDestroyArray(bestModel);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return models;
    /*
     * 
     * 
     */
}

/*
 * The function "MannealVisitParameters_doFuncEvals" is the implementation
 * version of the "annealVisitParameters/doFuncEvals" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 35-48). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function costs = doFuncEvals(canUseMatrix, models, func, varargin)
 */
static mxArray * MannealVisitParameters_doFuncEvals(int nargout_,
                                                    mxArray * canUseMatrix,
                                                    mxArray * models,
                                                    mxArray * func,
                                                    mxArray * varargin) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * costs = mclGetUninitializedArray();
    mxArray * e = mclGetUninitializedArray();
    mxArray * NM = mclGetUninitializedArray();
    mclCopyArray(&canUseMatrix);
    mclCopyArray(&models);
    mclCopyArray(&func);
    mclCopyArray(&varargin);
    /*
     * 
     * if canUseMatrix
     */
    if (mlfTobool(mclVa(canUseMatrix, "canUseMatrix"))) {
        /*
         * costs = feval(func, models, varargin{:});
         */
        mlfAssign(
          &costs,
          mlfFeval(
            mclValueVarargout(),
            mclVa(func, "func"),
            mclVa(models, "models"),
            mclVe(
              mlfIndexRef(
                mclVsa(varargin, "varargin"), "{?}", mlfCreateColonIndex())),
            NULL));
    /*
     * else
     */
    } else {
        /*
         * NM = size(models, 2);
         */
        mlfAssign(
          &NM,
          mlfSize(mclValueVarargout(), mclVa(models, "models"), _mxarray22_));
        /*
         * costs = zeros(NM, 1);
         */
        mlfAssign(&costs, mlfZeros(mclVv(NM, "NM"), _mxarray21_, NULL));
        /*
         * for e = 1:NM
         */
        {
            int v_ = mclForIntStart(1);
            int e_ = mclForIntEnd(mclVv(NM, "NM"));
            if (v_ > e_) {
                mlfAssign(&e, _mxarray23_);
            } else {
                /*
                 * costs(e) = feval(func, models(:,e), varargin{:});
                 * end
                 */
                for (; ; ) {
                    mclIntArrayAssign1(
                      &costs,
                      mlfFeval(
                        mclValueVarargout(),
                        mclVa(func, "func"),
                        mclVe(
                          mclArrayRef2(
                            mclVsa(models, "models"),
                            mlfCreateColonIndex(),
                            mlfScalar(v_))),
                        mclVe(
                          mlfIndexRef(
                            mclVsa(varargin, "varargin"),
                            "{?}",
                            mlfCreateColonIndex())),
                        NULL),
                      v_);
                    if (v_ == e_) {
                        break;
                    }
                    ++v_;
                }
                mlfAssign(&e, mlfScalar(v_));
            }
        }
    /*
     * end
     */
    }
    mclValidateOutput(
      costs, 1, nargout_, "costs", "annealVisitParameters/doFuncEvals");
    mxDestroyArray(NM);
    mxDestroyArray(e);
    mxDestroyArray(varargin);
    mxDestroyArray(func);
    mxDestroyArray(models);
    mxDestroyArray(canUseMatrix);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return costs;
    /*
     * 
     * 
     */
}

/*
 * The function "MannealVisitParameters_sampleFromPdf" is the implementation
 * version of the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function s = sampleFromPdf(temp, costs)
 */
static mxArray * MannealVisitParameters_sampleFromPdf(int nargout_,
                                                      mxArray * temp,
                                                      mxArray * costs) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * s = mclGetUninitializedArray();
    mxArray * cutoff = mclGetUninitializedArray();
    mxArray * dist = mclGetUninitializedArray();
    mclCopyArray(&temp);
    mclCopyArray(&costs);
    /*
     * dist = makePDF(temp, costs);
     */
    mlfAssign(
      &dist,
      mlfAnnealVisitParameters_makePDF(
        mclVa(temp, "temp"), mclVa(costs, "costs")));
    /*
     * cutoff = rand;
     */
    mlfAssign(&cutoff, mlfNRand(1, NULL));
    /*
     * s = find(cumsum(dist) >= cutoff);
     */
    mlfAssign(
      &s,
      mlfFind(
        NULL,
        NULL,
        mclGe(
          mclVe(mlfCumsum(mclVv(dist, "dist"), NULL)),
          mclVv(cutoff, "cutoff"))));
    /*
     * s = s(1);
     */
    mlfAssign(&s, mclIntArrayRef1(mclVsv(s, "s"), 1));
    mclValidateOutput(
      s, 1, nargout_, "s", "annealVisitParameters/sampleFromPdf");
    mxDestroyArray(dist);
    mxDestroyArray(cutoff);
    mxDestroyArray(costs);
    mxDestroyArray(temp);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return s;
    /*
     * 
     * 
     */
}

/*
 * The function "MannealVisitParameters_makePDF" is the implementation version
 * of the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = makePDF(temp, costs)
 */
static mxArray * MannealVisitParameters_makePDF(int nargout_,
                                                mxArray * temp,
                                                mxArray * costs) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * pdf = mclGetUninitializedArray();
    mxArray * ans = mclGetUninitializedArray();
    mxArray * w = mclGetUninitializedArray();
    mxArray * good = mclGetUninitializedArray();
    mxArray * bad = mclGetUninitializedArray();
    mclCopyArray(&temp);
    mclCopyArray(&costs);
    /*
     * % Forms exponential probability distribution given a temperature and vector of
     * % costs.  Internal function for simulated annealing algorithm.
     * 
     * bad = find(isnan(costs));
     */
    mlfAssign(
      &bad, mlfFind(NULL, NULL, mclVe(mlfIsnan(mclVa(costs, "costs")))));
    /*
     * if isempty(bad)
     */
    if (mlfTobool(mclVe(mlfIsempty(mclVv(bad, "bad"))))) {
        /*
         * pdf = eprob(temp, costs);
         */
        mlfAssign(
          &pdf,
          mlfAnnealVisitParameters_eprob(
            mclVa(temp, "temp"), mclVa(costs, "costs")));
    /*
     * else
     */
    } else {
        /*
         * good = find(~isnan(costs));
         */
        mlfAssign(
          &good,
          mlfFind(NULL, NULL, mclNot(mclVe(mlfIsnan(mclVa(costs, "costs"))))));
        /*
         * w = costs(good);
         */
        mlfAssign(
          &w, mclArrayRef1(mclVsa(costs, "costs"), mclVsv(good, "good")));
        /*
         * pdf = eprob(temp, w);
         */
        mlfAssign(
          &pdf,
          mlfAnnealVisitParameters_eprob(mclVa(temp, "temp"), mclVv(w, "w")));
        /*
         * pdf(good) = pdf;
         */
        mclArrayAssign1(&pdf, mclVsv(pdf, "pdf"), mclVsv(good, "good"));
        /*
         * pdf(bad) = 0;
         */
        mclArrayAssign1(&pdf, _mxarray20_, mclVsv(bad, "bad"));

        mexPrintf("Warning: the cost function generated a NaN for one or more models.");
    /*
     * end
     */
    }
    mclValidateOutput(pdf, 1, nargout_, "pdf", "annealVisitParameters/makePDF");
    mxDestroyArray(bad);
    mxDestroyArray(good);
    mxDestroyArray(w);
    mxDestroyArray(ans);
    mxDestroyArray(costs);
    mxDestroyArray(temp);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return pdf;
    /*
     * 
     */
}

/*
 * The function "MannealVisitParameters_eprob" is the implementation version of
 * the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = eprob(temp, costs)
 */
static mxArray * MannealVisitParameters_eprob(int nargout_,
                                              mxArray * temp,
                                              mxArray * costs) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * pdf = mclGetUninitializedArray();
    mxArray * scale = mclGetUninitializedArray();
    mxArray * mpdf = mclGetUninitializedArray();
    mxArray * toobig = mclGetUninitializedArray();
    mclCopyArray(&temp);
    mclCopyArray(&costs);
    /*
     * % Scales cost vector and calculates exponential probability distribution.  The
     * % scaling is necessary to permit wide ranges in temperature.  Internal
     * % function for simulated annealing algorithm.
     * 
     * toobig = 708.3964185322641;
     */
    mlfAssign(&toobig, _mxarray26_);
    /*
     * pdf = costs/temp;
     */
    mlfAssign(&pdf, mclMrdivide(mclVa(costs, "costs"), mclVa(temp, "temp")));
    /*
     * mpdf = max(pdf);
     */
    mlfAssign(&mpdf, mlfMax(NULL, mclVv(pdf, "pdf"), NULL, NULL));
    /*
     * 
     * if mpdf>toobig
     */
    if (mclGtBool(mclVv(mpdf, "mpdf"), mclVv(toobig, "toobig"))) {
        /*
         * scale = mpdf/toobig;
         */
        mlfAssign(
          &scale, mclMrdivide(mclVv(mpdf, "mpdf"), mclVv(toobig, "toobig")));
        /*
         * pdf = exp(-pdf/scale);
         */
        mlfAssign(
          &pdf,
          mlfExp(
            mclMrdivide(mclUminus(mclVv(pdf, "pdf")), mclVv(scale, "scale"))));
        /*
         * pdf = pdf/max(pdf);
         */
        mlfAssign(
          &pdf,
          mclMrdivide(
            mclVv(pdf, "pdf"),
            mclVe(mlfMax(NULL, mclVv(pdf, "pdf"), NULL, NULL))));
        /*
         * pdf = pdf.^scale;
         */
        mlfAssign(&pdf, mlfPower(mclVv(pdf, "pdf"), mclVv(scale, "scale")));
    /*
     * else
     */
    } else {
        /*
         * pdf = exp(-pdf);
         */
        mlfAssign(&pdf, mlfExp(mclUminus(mclVv(pdf, "pdf"))));
        /*
         * pdf = pdf/max(pdf);
         */
        mlfAssign(
          &pdf,
          mclMrdivide(
            mclVv(pdf, "pdf"),
            mclVe(mlfMax(NULL, mclVv(pdf, "pdf"), NULL, NULL))));
    /*
     * end
     */
    }
    /*
     * 
     * pdf = pdf/sum(pdf);
     */
    mlfAssign(
      &pdf,
      mclMrdivide(mclVv(pdf, "pdf"), mclVe(mlfSum(mclVv(pdf, "pdf"), NULL))));
    mclValidateOutput(pdf, 1, nargout_, "pdf", "annealVisitParameters/eprob");
    mxDestroyArray(toobig);
    mxDestroyArray(mpdf);
    mxDestroyArray(scale);
    mxDestroyArray(costs);
    mxDestroyArray(temp);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return pdf;
}

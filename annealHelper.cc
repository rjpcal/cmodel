/*
 * MATLAB Compiler: 2.1
 * Date: Fri Mar 23 17:17:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "-h" "annealVisitParameters" 
 */
#include "annealVisitParameters.h"

#include "rutil.h"

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

static mxArray * mlfsampleFromPdf(mxArray * temp,
											 mxArray * costs);

static void mlxsampleFromPdf(int nlhs,
									  mxArray * plhs[],
									  int nrhs,
									  mxArray * prhs[]);

static mxArray * mlfmakePDF(mxArray * temp,
									 mxArray * costs);

static void mlxmakePDF(int nlhs,
							  mxArray * plhs[],
							  int nrhs,
							  mxArray * prhs[]);

static mxArray * mlfeprob(mxArray * temp,
								  mxArray * costs);

static void mlxeprob(int nlhs,
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

static mxArray * MmakeTestModels(int nargout_,
											mxArray * x,
											mxArray * bestModel,
											mxArray * valueScalingRange,
											mxArray * deltas,
											mxArray * bounds);

static mxArray * MdoFuncEvals(int nargout_,
										mxArray * canUseMatrix,
										mxArray * models,
										mxArray * func,
										mxArray * varargin);

static mxArray * MsampleFromPdf(int nargout_,
										  mxArray * temp,
										  mxArray * costs);

static mxArray * MmakePDF(int nargout_,
								  mxArray * temp,
								  mxArray * costs);

static mxArray * Meprob(int nargout_,
								mxArray * temp,
								mxArray * costs);

static mexFunctionTableEntry local_function_table_[5]
  = { { "makeTestModels",
        NULL /*mlxmakeTestModels*/, 5, 1, NULL },
      { "doFuncEvals", NULL/*mlxdoFuncEvals*/, -4, 1, NULL },
      { "sampleFromPdf", mlxsampleFromPdf, 2, 1, NULL },
      { "makePDF", mlxmakePDF, 2, 1, NULL },
      { "eprob", mlxeprob, 2, 1, NULL } };

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
 * The function "mlfsampleFromPdf" contains the normal
 * interface for the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfsampleFromPdf(mxArray * temp,
                                                        mxArray * costs) {
    int nargout = 1;
    mxArray * s = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    s = MsampleFromPdf(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(s);
}

/*
 * The function "mlxsampleFromPdf" contains the feval
 * interface for the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). The feval function calls the implementation version of
 * annealVisitParameters/sampleFromPdf through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
static void mlxsampleFromPdf(int nlhs,
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
    mplhs[0] = MsampleFromPdf(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfmakePDF" contains the normal
 * interface for the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfmakePDF(mxArray * temp,
                                                  mxArray * costs) {
    int nargout = 1;
    mxArray * pdf = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    pdf = MmakePDF(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(pdf);
}

/*
 * The function "mlxmakePDF" contains the feval interface
 * for the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). The feval function calls the implementation version of
 * annealVisitParameters/makePDF through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
static void mlxmakePDF(int nlhs,
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
    mplhs[0] = MmakePDF(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "mlfeprob" contains the normal interface
 * for the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
static mxArray * mlfeprob(mxArray * temp,
                                                mxArray * costs) {
    int nargout = 1;
    mxArray * pdf = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, temp, costs);
    pdf = Meprob(nargout, temp, costs);
    mlfRestorePreviousContext(0, 2, temp, costs);
    return mlfReturnValue(pdf);
}

/*
 * The function "mlxeprob" contains the feval interface
 * for the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). The feval function calls the implementation version of
 * annealVisitParameters/eprob through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
static void mlxeprob(int nlhs,
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
    mplhs[0] = Meprob(nlhs, mprhs[0], mprhs[1]);
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
// XXX
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

	 validateInput(bestModel);
    mclCopyArray(&bestModel);

	 validateInput(valueScalingRange);
    mclCopyArray(&valueScalingRange);

	 validateInput(deltas);
    mclCopyArray(&deltas);

	 validateInput(bounds);
    mclCopyArray(&bounds);

	 validateInput(canUseMatrix);
    mclCopyArray(&canUseMatrix);

	 validateInput(FUN);
    mclCopyArray(&FUN);

	 validateInput(temp);
    mclCopyArray(&temp);

	 validateInput(varargin);
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
					mlfFind(
							  NULL,
							  NULL,
							  mclNe(mlfCtranspose(deltas), _mxarray20_)),
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
              MmakeTestModels(
										1,
										x,
										bestModel,
										valueScalingRange,
										deltas,
										bounds));
            /*
             * deltas, bounds);
             * 
             * costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
             */
            mlfAssign(
              &costs,
				  MdoFuncEvals(1,
                canUseMatrix,
                modelmatrix,
                FUN,
					 mlfIndexRef(
									 varargin,
									 "{?}",
									 mlfCreateColonIndex())
									));
            /*
             * S.nevals = S.nevals + length(costs);
             */
            mlfIndexAssign(
              &S,
              ".nevals",
              mclFeval(
                mclValueVarargout(),
                mlxPlus,
                mlfIndexRef(S, ".nevals"),
                mlfScalar(mclLengthInt(costs)),
                NULL));
            /*
             * 
             * % Sample from probability distribution
             * s = sampleFromPdf(temp, costs);
             */
            mlfAssign(
              &s,
              mlfsampleFromPdf(
                temp, costs));
            /*
             * 
             * bestModel(x, 1) = modelmatrix(x, s);
             */
            mclArrayAssign2(
              &bestModel,
              mclArrayRef2(
                modelmatrix,
                x,
                s),
              x,
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
    mlfIndexAssign(&S, ".newModel", bestModel);
    /*
     * 
     * S.cost = costs(s);
     */
    mlfIndexAssign(
      &S, ".cost", mclArrayRef1(costs, s));
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
static mxArray * MmakeTestModels(int nargout_,
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
        mclArrayRef1(bestModel, x),
        mclMtimes(
          valueScalingRange,
          mclArrayRef1(deltas, x))));
    /*
     * xv = xv(find( (xv<=bounds(x,2)) & (xv>=bounds(x,1)) ));
     */
    mlfAssign(
      &xv,
      mclArrayRef1(
        xv,
        mlfFind(
          NULL,
          NULL,
          mclAnd(
            mclLe(
              xv,
				  mclArrayRef2(
                  bounds, x, _mxarray22_)),
            mclGe(
              xv,
				  mclArrayRef2(
                  bounds, x, _mxarray21_))))));
    /*
     * models = bestModel*ones(1, length(xv));
     */
    mlfAssign(
      &models,
      mclMtimes(
        bestModel,
		  mlfOnes(
            _mxarray21_, mlfScalar(mclLengthInt(xv)), NULL)));
    /*
     * models(x,:) = xv;
     */
    mclArrayAssign2(
      &models, xv, x, mlfCreateColonIndex());
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
 * The function "MdoFuncEvals" is the implementation
 * version of the "annealVisitParameters/doFuncEvals" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 35-48). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function costs = doFuncEvals(canUseMatrix, models, func, varargin)
 */
static mxArray * MdoFuncEvals(int nargout_,
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
    if (mlfTobool(canUseMatrix)) {
        /*
         * costs = feval(func, models, varargin{:});
         */
        mlfAssign(
          &costs,
          mlfFeval(
            mclValueVarargout(),
            func,
            models,
				mlfIndexRef(
                varargin, "{?}", mlfCreateColonIndex()),
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
          mlfSize(mclValueVarargout(), models, _mxarray22_));
        /*
         * costs = zeros(NM, 1);
         */
        mlfAssign(&costs, mlfZeros(NM, _mxarray21_, NULL));
        /*
         * for e = 1:NM
         */
        {
            int v_ = mclForIntStart(1);
            int e_ = mclForIntEnd(NM);
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
                        func,
								mclArrayRef2(
												 models,
												 mlfCreateColonIndex(),
												 mlfScalar(v_)),
								mlfIndexRef(
												varargin,
												"{?}",
												mlfCreateColonIndex()),
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
 * The function "MsampleFromPdf" is the implementation
 * version of the "annealVisitParameters/sampleFromPdf" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 48-55). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function s = sampleFromPdf(temp, costs)
 */
static mxArray * MsampleFromPdf(int nargout_,
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
      mlfmakePDF(
        temp, costs));
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
          mlfCumsum(dist, NULL),
          cutoff)));
    /*
     * s = s(1);
     */
    mlfAssign(&s, mclIntArrayRef1(s, 1));
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
 * The function "MmakePDF" is the implementation version
 * of the "annealVisitParameters/makePDF" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 55-72). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = makePDF(temp, costs)
 */
static mxArray * MmakePDF(int nargout_,
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
      &bad, mlfFind(NULL, NULL, mlfIsnan(costs)));
    /*
     * if isempty(bad)
     */
    if (mlfTobool(mlfIsempty(bad))) {
        /*
         * pdf = eprob(temp, costs);
         */
        mlfAssign(
          &pdf,
          mlfeprob(
            temp, costs));
    /*
     * else
     */
    } else {
        /*
         * good = find(~isnan(costs));
         */
        mlfAssign(
          &good,
          mlfFind(NULL, NULL, mclNot(mlfIsnan(costs))));
        /*
         * w = costs(good);
         */
        mlfAssign(
          &w, mclArrayRef1(costs, good));
        /*
         * pdf = eprob(temp, w);
         */
        mlfAssign(
          &pdf,
          mlfeprob(temp, w));
        /*
         * pdf(good) = pdf;
         */
        mclArrayAssign1(&pdf, pdf, good);
        /*
         * pdf(bad) = 0;
         */
        mclArrayAssign1(&pdf, _mxarray20_, bad);

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
 * The function "Meprob" is the implementation version of
 * the "annealVisitParameters/eprob" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 72-92). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */
/*
 * function pdf = eprob(temp, costs)
 */
static mxArray * Meprob(int nargout_,
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
    mlfAssign(&pdf, mclMrdivide(costs, temp));
    /*
     * mpdf = max(pdf);
     */
    mlfAssign(&mpdf, mlfMax(NULL, pdf, NULL, NULL));
    /*
     * 
     * if mpdf>toobig
     */
    if (mclGtBool(mpdf, toobig)) {
        /*
         * scale = mpdf/toobig;
         */
        mlfAssign(
          &scale, mclMrdivide(mpdf, toobig));
        /*
         * pdf = exp(-pdf/scale);
         */
        mlfAssign(
          &pdf,
          mlfExp(
            mclMrdivide(mclUminus(pdf), scale)));
        /*
         * pdf = pdf/max(pdf);
         */
        mlfAssign(
          &pdf,
          mclMrdivide(
            pdf,
            mlfMax(NULL, pdf, NULL, NULL)));
        /*
         * pdf = pdf.^scale;
         */
        mlfAssign(&pdf, mlfPower(pdf, scale));
    /*
     * else
     */
    } else {
        /*
         * pdf = exp(-pdf);
         */
        mlfAssign(&pdf, mlfExp(mclUminus(pdf)));
        /*
         * pdf = pdf/max(pdf);
         */
        mlfAssign(
          &pdf,
          mclMrdivide(
            pdf,
            mlfMax(NULL, pdf, NULL, NULL)));
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
      mclMrdivide(pdf, mlfSum(pdf, NULL)));
    mclValidateOutput(pdf, 1, nargout_, "pdf", "annealVisitParameters/eprob");
    mxDestroyArray(toobig);
    mxDestroyArray(mpdf);
    mxDestroyArray(scale);
    mxDestroyArray(costs);
    mxDestroyArray(temp);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return pdf;
}

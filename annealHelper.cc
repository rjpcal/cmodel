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

static mxArray * MannealVisitParameters(int nargout_,
                                        mxArray * bestModel,
                                        mxArray * valueScalingRange,
                                        mxArray * deltas,
                                        mxArray * bounds,
                                        mxArray * canUseMatrix,
                                        mxArray * FUN,
                                        mxArray * temp,
                                        mxArray * varargin);

static mxArray * makeTestModels(mxArray * x,
										  mxArray * bestModel,
										  mxArray * valueScalingRange,
										  mxArray * deltas,
										  mxArray * bounds);

static mxArray * doFuncEvals(mxArray * canUseMatrix,
									  mxArray * models,
									  mxArray * func,
									  mxArray * varargin);

static mxArray * sampleFromPdf(mxArray * temp,
										 mxArray * costs);

static mxArray * makePDF(mxArray * temp,
								 mxArray * costs);

static mxArray * eprob(mxArray * temp,
							  mxArray * costs);

static mexFunctionTableEntry local_function_table_[5]
  = { { "makeTestModels",
        NULL /*mlxmakeTestModels*/, 5, 1, NULL },
      { "doFuncEvals", NULL/*mlxdoFuncEvals*/, -4, 1, NULL },
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
             * deltas, bounds);
             */
            mlfAssign(
              &modelmatrix,
              makeTestModels(
									  x,
									  bestModel,
									  valueScalingRange,
									  deltas,
									  bounds));
            /*
             * 
             * costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
             */
            mlfAssign(
              &costs,
				  doFuncEvals(
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
              sampleFromPdf(
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
static mxArray * makeTestModels(mxArray * x,
										  mxArray * bestModel,
										  mxArray * valueScalingRange,
										  mxArray * deltas,
										  mxArray * bounds) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_annealVisitParameters);
    mxArray * models = mclGetUninitializedArray();
    mxArray * xv = mclGetUninitializedArray();
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
      models, 1, 1, "models", "annealVisitParameters/makeTestModels");
    mxDestroyArray(xv);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return models;
    /*
     * 
     * 
     */
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
static mxArray * doFuncEvals(mxArray * canUseMatrix,
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
      costs, 1, 1, "costs", "annealVisitParameters/doFuncEvals");
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
static mxArray * sampleFromPdf(mxArray * temp,
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
    mlfAssign(&dist, makePDF(temp, costs));

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
      s, 1, 1, "s", "annealVisitParameters/sampleFromPdf");
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
static mxArray * makePDF(mxArray * temp,
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
        mlfAssign(&pdf, eprob(temp, costs));
    }
	 else {
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

        mlfAssign(&pdf, eprob(temp, w));

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
    mclValidateOutput(pdf, 1, 1, "pdf", "annealVisitParameters/makePDF");
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
static mxArray * eprob(mxArray * temp,
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
    mclValidateOutput(pdf, 1, 1, "pdf", "annealVisitParameters/eprob");
    mxDestroyArray(toobig);
    mxDestroyArray(mpdf);
    mxDestroyArray(scale);
    mxDestroyArray(costs);
    mxDestroyArray(temp);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return pdf;
}

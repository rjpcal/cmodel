///////////////////////////////////////////////////////////////////////
//
// doSimplex.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 06:20:45 2001
// written: Wed Apr 18 14:31:36 2001
// $Id$
//
//
// MATLAB Compiler: 2.1
// Date: Sun Apr  1 19:52:50 2001
//
// Arguments: "-B" "macro_default" "-O" "all" "-O"
// "fold_scalar_mxarrays:on" "-O" "fold_non_scalar_mxarrays:on" "-O"
// "optimize_integer_for_loops:on" "-O" "array_indexing:on" "-O"
// "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C" "-t" "-T"
// "link:mexlibrary" "libmatlbmx.mlib" "-h" "doSimplex"
//
///////////////////////////////////////////////////////////////////////

#ifndef DOSIMPLEX_CC_DEFINED
#define DOSIMPLEX_CC_DEFINED

#include "doSimplex.h"

#include "error.h"
#include "mexbuf.h"
#include "mtx.h"
#include "multivarfunction.h"
#include "strings.h"

#include <fstream.h>
#include <iomanip.h>
#include "libmatlbm.h"

#define LOCAL_PROF
#include "trace.h"

class MatlabFunction : public MultivarFunction {
  mxArray* itsFunfcn;
  mxArray* itsVarargin_ref;

  static mxArray* getref(mxArray* varargin)
	 {
		return mlfIndexRef(varargin,
								 "{?}",
								 mlfCreateColonIndex());
	 }

  double evaluate_mx(mxArray* x_mx)
  {
	 DOTRACE("evaluate_mx");

	 mxArray* mx =  mlfFeval(mclValueVarargout(),
									 itsFunfcn,
									 x_mx,
									 getref(itsVarargin_ref),
									 NULL);
	 double result = mxGetScalar(mx);
	 mxDestroyArray(mx);
	 return result;
  }

  virtual double doEvaluate(const Mtx& x)
  {
	 return evaluate_mx(x.makeMxArray());
  }

public:
  MatlabFunction(mxArray* funfcn_mx, mxArray* varargin_mx) :
	 MultivarFunction(),
	 itsFunfcn(funfcn_mx),
	 itsVarargin_ref(varargin_mx)
  {}

  virtual ~MatlabFunction() {}
};

int extractMaxIters(const mxArray* arr, int numModelParams)
{
  if (mxIsChar(arr))
	 {
		if (Mtx::extractString(arr) == "200*numberofvariables")
		  {
			 return 200*numModelParams;
		  }
		else
		  {
			 mexErrMsgTxt("Option must be an integer value "
							  "if not the default.");
		  }
	 }

  return int(mxGetScalar(arr));
}

namespace {
  MexBuf* mexBuf = 0;
}

void InitializeModule_doSimplex(void) {
  mexBuf = new MexBuf;
  cout.rdbuf(mexBuf);
  cerr.rdbuf(mexBuf);
}

void TerminateModule_doSimplex(void) {
  delete mexBuf;
}

static mxArray * MdoSimplex(mxArray * * fval,
                            mxArray * * exitflag,
                            mxArray * * output,
                            int nargout_,
                            mxArray * funfcn,
                            mxArray * x_in,
                            mxArray * printtype,
                            mxArray * tolx,
                            mxArray * tolf,
                            mxArray * maxfun,
                            mxArray * maxiter,
                            mxArray * debugFlags,
                            mxArray * varargin);

_mexLocalFunctionTable _local_function_table_doSimplex
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfDoSimplex" contains the normal interface for the
 * "doSimplex" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfDoSimplex(mxArray * * fval,
                       mxArray * * exitflag,
                       mxArray * * output,
                       mxArray * funfcn,
                       mxArray * x_in,
                       mxArray * printtype,
                       mxArray * tolx,
                       mxArray * tolf,
                       mxArray * maxfun,
                       mxArray * maxiter,
                       mxArray * debugFlags,
                       ...) {
    mxArray * varargin = NULL;
    int nargout = 1;
    mxArray * x = mclGetUninitializedArray();
    mxArray * fval__ = mclGetUninitializedArray();
    mxArray * exitflag__ = mclGetUninitializedArray();
    mxArray * output__ = mclGetUninitializedArray();
    mlfVarargin(&varargin, debugFlags, 0);
    mlfEnterNewContext(
      3,
      -9,
      fval,
      exitflag,
      output,
      funfcn,
      x_in,
      printtype,
      tolx,
      tolf,
      maxfun,
      maxiter,
      debugFlags,
      varargin);
    if (fval != NULL) {
        ++nargout;
    }
    if (exitflag != NULL) {
        ++nargout;
    }
    if (output != NULL) {
        ++nargout;
    }
    x
      = MdoSimplex(
          &fval__,
          &exitflag__,
          &output__,
          nargout,
          funfcn,
          x_in,
          printtype,
          tolx,
          tolf,
          maxfun,
          maxiter,
          debugFlags,
          varargin);
    mlfRestorePreviousContext(
      3,
      8,
      fval,
      exitflag,
      output,
      funfcn,
      x_in,
      printtype,
      tolx,
      tolf,
      maxfun,
      maxiter,
      debugFlags);
    mxDestroyArray(varargin);
    if (fval != NULL) {
        mclCopyOutputArg(fval, fval__);
    } else {
        mxDestroyArray(fval__);
    }
    if (exitflag != NULL) {
        mclCopyOutputArg(exitflag, exitflag__);
    } else {
        mxDestroyArray(exitflag__);
    }
    if (output != NULL) {
        mclCopyOutputArg(output, output__);
    } else {
        mxDestroyArray(output__);
    }
    return mlfReturnValue(x);
}

/*
 * The function "mlxDoSimplex" contains the feval interface for the "doSimplex"
 * M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * The feval function calls the implementation version of doSimplex through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxDoSimplex(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[9];
    mxArray * mplhs[4];
    int i;
    if (nlhs > 4) {
		mexErrMsgTxt("Run-time Error: File: doSimplex Line: 1 Column: 1 "
						 "The function \"doSimplex\" was called with more "
						 "than the declared number of outputs (4).");
    }
    for (i = 0; i < 4; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 8 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 8; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0,
      8,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6],
      mprhs[7]);
    mprhs[8] = NULL;
    mlfAssign(&mprhs[8], mclCreateVararginCell(nrhs - 8, prhs + 8));
    mplhs[0]
      = MdoSimplex(
          &mplhs[1],
          &mplhs[2],
          &mplhs[3],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4],
          mprhs[5],
          mprhs[6],
          mprhs[7],
          mprhs[8]);
    mlfRestorePreviousContext(
      0,
      8,
      mprhs[0],
      mprhs[1],
      mprhs[2],
      mprhs[3],
      mprhs[4],
      mprhs[5],
      mprhs[6],
      mprhs[7]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 4 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 4; ++i) {
        mxDestroyArray(mplhs[i]);
    }
    mxDestroyArray(mprhs[8]);
}

class SimplexOptimizer {
private:
  MultivarFunction& itsObjective;

  const Mtx itsInitialParams;
  const int itsPrnt;
  const double itsTolx;
  const double itsTolf;
  const int itsNparams;
  const int itsMaxFevals;
  const int itsMaxIters;

  Mtx itsSimplex;
  Mtx itsFvals;

  int itsIterCount;

  enum IterType { INITIAL, EXPAND, REFLECT, CONTRACT_OUTSIDE, CONTRACT_INSIDE, SHRINK };

  IterType itsCurIter;

  void buildInitialSimplex()
  {
	 // Following improvement suggested by L.Pfeffer at Stanford
	 // 5 percent deltas for non-zero terms
	 const double usual_delta = 0.05;

	 // Even smaller delta for zero elements of x
	 const double zero_term_delta = 0.00025;

	 for (int j = 0; j < itsNparams; ++j)
		{
		  itsSimplex.column(j+1) = itsInitialParams;

		  if (itsSimplex.at(j,j+1) != 0.0)
			 itsSimplex.at(j,j+1) *= (1.0 + usual_delta);
		  else
			 itsSimplex.at(j,j+1) = zero_term_delta;

		  itsFvals.at(0,j+1) = itsObjective.evaluate(itsSimplex.column(j+1));
		}
  }

  struct FuncPoint
  {
	 FuncPoint(const Mtx& x_, double f_) : x(x_), f(f_) {}

	 const Mtx x;
	 const double f;

	 bool betterThan (const FuncPoint& other)
	   { return (f < other.f); }
  };

  FuncPoint evaluate(const Mtx& x)
  {
	 return FuncPoint(x, itsObjective.evaluate(x));
  }

  void putInSimplex(const FuncPoint& p, int pointNumber)
  {
	 itsSimplex.column(pointNumber) = p.x.column(0);
	 itsFvals.at(0, pointNumber) = p.f;
  }

  void putInSimplex(const Mtx& params, int simplexPoint)
  {
	 putInSimplex(evaluate(params), simplexPoint);
  }

  FuncPoint simplexAt(int simplexPoint)
  {
	 return FuncPoint(itsSimplex.column(simplexPoint),
							itsFvals.at(0, simplexPoint));
  }

#if 0
  void fullSort()
  {
	 DOTRACE("fullSort");
	 // sort so itsSimplex.column(0) has the lowest function value 
	 Mtx index = itsFvals.row(0).getSortOrder();

	 itsFvals.row(0).reorder(index);
  	 itsSimplex.reorderColumns(index);
  }
#endif

  void minimalSort()
  {
	 DOTRACE("minimalSort");

	 Mtx index = itsFvals.row(0).getSortOrder();

	 const int smallest = int(index.at(0));
	 const int largest = int(index.at(itsNparams));
	 const int largest2 = int(index.at(itsNparams-1));

	 // These swaps are smart enough to check if the column numbers
	 // are the same before doing the swap

	 itsFvals.swapColumns(0, smallest);
	 itsFvals.swapColumns(largest, itsNparams);
	 itsFvals.swapColumns(largest2, itsNparams-1);

	 itsSimplex.swapColumns(0, smallest);
	 itsSimplex.swapColumns(largest, itsNparams);
	 itsSimplex.swapColumns(largest2, itsNparams-1);
  }

  bool withinTolf()
  {
	 MtxConstIter fvals = itsFvals.rowIter(0);
	 const double f0 = *fvals;
	 ++fvals;

	 for (; fvals.hasMore(); ++fvals)
		{
		  if (fabs(*fvals - f0) > itsTolf)
			 return false;
		}

	 return true;
  }

  bool withinTolx()
  {
	 const MtxConstIter col0_ = itsSimplex.columnIter(0);

	 for (int col = 1; col < itsSimplex.ncols(); ++col)
		{
		  MtxConstIter col0(col0_);
		  MtxConstIter coln = itsSimplex.columnIter(col);

		  for (; col0.hasMore(); ++col0, ++coln)
			 if ( fabs(*col0 - *coln) > itsTolx ) return false;
		}

	 return true;
  }

  bool tooManyFevals()
  {
	 if (funcCount() < itsMaxFevals)
		return false;

	 if (itsPrnt != NONE) {
		cerr << "\nExiting: Maximum number of function evaluations "
			  << "has been exceeded\n"
			  << "         - increase MaxFunEvals option.\n"
			  << "         Current function value: "
			  << bestFval()
			  << "\n\n";
	 }

	 return true;
  }

  bool tooManyIters()
  {
	 if (itsIterCount < itsMaxIters)
		return false;

	 if (itsPrnt != NONE) {
		cerr << "\nExiting: Maximum number of iterations "
			  << "has been exceeded\n"
			  << "         - increase MaxIter option.\n"
			  << "         Current function value: "
			  << bestFval()
			  << "\n\n";
	 }	 

	 return true;
  }

  void doOneIter();

  void printHeader()
  {
	 if (itsPrnt == ITER)
		cerr << "\n Iteration   Func-count     min f(x)         Procedure\n";
  }

  static const char* iterTypeString(IterType how)
  {
	 switch (how)
		{
		case INITIAL:
		  return "initial";
		  break;
		case EXPAND:
		  return "expand";
		  break;
		case REFLECT:
		  return "reflect";
		  break;
		case CONTRACT_OUTSIDE:
		  return "contract outside";
		  break;
		case CONTRACT_INSIDE:
		  return "contract inside";
		  break;
		case SHRINK:
		  return "shrink";
		  break;
		default:
		  return "unknown";
		  break;
		}
  }

  void printIter()
  {
	 if (itsPrnt == ITER)
		{
		  cerr << setw(6) << itsIterCount
				 << setw(13) << funcCount()
				 << setw(17) << setprecision(6) << itsFvals.at(0,0)
				 << "         " << iterTypeString(itsCurIter) << '\n';
		}
  }

  void printSimplex()
  {
	 if (itsPrnt == SIMPLEX)
		{
		  cerr << '\n' << iterTypeString(itsCurIter) << '\n';
		  itsSimplex.print("simplex");
		  itsFvals.print("fvals");
		  cerr << "funcEvals: " << funcCount() << '\n';
		}
  }

public:
  enum PrintType { NONE, NOTIFY, FINAL, ITER, SIMPLEX };

private:
  static PrintType extractPrinttype(const fixed_string& printtype)
  {
	 if      (printtype == "notify")   return NOTIFY;
	 else if (printtype == "none")     return NONE;
	 else if (printtype == "off")      return NONE;
	 else if (printtype == "iter")     return ITER;
	 else if (printtype == "final")    return FINAL;
	 else if (printtype == "simplex")  return SIMPLEX;

	 return NOTIFY;
  }

public:
  SimplexOptimizer(MultivarFunction& objective,
						 const Mtx& x_in,
						 const fixed_string& printtype,
						 const double tolx,
						 const double tolf,
						 const int nparams,
						 const int maxfun,
						 const int maxiter) :
	 itsObjective(objective),

	 itsInitialParams(x_in.asColumn()), // Set up simplex near the initial guess
	 itsPrnt(extractPrinttype(printtype)),
	 itsTolx(tolx),
	 itsTolf(tolf),
	 itsNparams(nparams),
	 itsMaxFevals(maxfun),
	 itsMaxIters(maxiter),

	 itsSimplex(nparams, nparams+1),
	 itsFvals(1, nparams+1),

	 itsIterCount(1)
  {
	 DOTRACE("SimplexOptimizer::SimplexOptimizer");
	 // Place input guess in the simplex! (credit L.Pfeffer at Stanford)
	 putInSimplex(itsInitialParams, 0);

	 buildInitialSimplex();

	 minimalSort();
  }

  virtual ~SimplexOptimizer() {}

  int optimize();

  virtual Mtx bestParams()
  { return itsSimplex.column(0); }

  virtual double bestFval()
  { return itsFvals.at(0,0); }

  virtual int iterCount()
  { return itsIterCount; }

  virtual int funcCount()
  { return itsObjective.evalCount(); }

  virtual const char* algorithm()
  { return "Nelder-Mead simplex direct search"; }
};


int SimplexOptimizer::optimize()
{
  printHeader();
  itsCurIter = INITIAL;
  printIter();
  printSimplex();

  // Iterate until the diameter of the simplex is less than itsTolx AND
  // the function values differ from the min by less than itsTolf, or
  // the max function evaluations are exceeded. (Cannot use OR
  // instead of AND.)

  while (!withinTolf() || !withinTolx())
	 {
		if (tooManyIters()) return 0;
		if (tooManyFevals()) return 0;

		doOneIter();
	 } // end main algorithm loop

  cerr << "\nOptimization terminated successfully:\n"
		 << " the current x satisfies the termination criteria using "
		 << "OPTIONS.TolX of " << itsTolx << '\n'
		 << " and F(X) satisfies the convergence criteria using "
		 << "OPTIONS.TolFun of " << itsTolf << '\n';

  return 1;
}


void SimplexOptimizer::doOneIter()
{
DOTRACE("SimplexOptimizer::doOneIter");

  // compute average of the itsNparams (NOT itsNparams+1) best points
  Mtx xbar(itsSimplex.columns(0,itsNparams).meanColumn());

  // Direction along which to reflect/expand/contract
  Mtx direction(xbar - itsSimplex.columns(itsNparams,1));

  // Calculate the reflection point
  FuncPoint rflPt = evaluate(xbar + direction);

  if (rflPt.betterThan(simplexAt(0)))
 	 {
 		DOTRACE("rflPt.betterThan(simplexAt(0))");

 		// Calculate the expansion point
 		FuncPoint expPt = evaluate(xbar + direction*2.0);

 		if (expPt.betterThan(rflPt))
 		  { putInSimplex(expPt, itsNparams); itsCurIter = EXPAND; }
 		else
 		  { putInSimplex(rflPt, itsNparams); itsCurIter = REFLECT; }
 	 }
  else
	 {
		DOTRACE("!addReflection");

		if (rflPt.betterThan(simplexAt(itsNparams-1)))
		  { putInSimplex(rflPt, itsNparams); itsCurIter = REFLECT; }
		else
		  {
			 // Perform contraction

			 if (rflPt.betterThan(simplexAt(itsNparams)))
				{
				  // Perform an outside contraction
				  FuncPoint ictPt = evaluate(xbar + direction*0.5);

				  if (ictPt.betterThan(rflPt))
					 { putInSimplex(ictPt, itsNparams); itsCurIter = CONTRACT_OUTSIDE; }
				  else
					 { itsCurIter = SHRINK; }
				}
			 else
				{
				  // Perform an inside contraction
				  FuncPoint octPt = evaluate(xbar - direction*0.5);

				  if (octPt.betterThan(simplexAt(itsNparams)))
					 { putInSimplex(octPt, itsNparams); itsCurIter = CONTRACT_INSIDE; }
				  else
					 { itsCurIter = SHRINK; }
				}

			 if (SHRINK == itsCurIter)
				{
				  for (int j = 1; j < itsNparams+1; ++j)
					 {
						putInSimplex(itsSimplex.columns(0,1) +
										 (itsSimplex.columns(j,1) -
										  itsSimplex.columns(0,1)) * 0.5,
										 j);
					 }
				}
		  }
	 }

  minimalSort();

  ++itsIterCount;

  printIter();
  printSimplex();
}


/*
 * The function "MdoSimplex" is the implementation version of the "doSimplex"
 * M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * It contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */

static mxArray * MdoSimplex(mxArray * * fval,
                            mxArray * * exitflag_mx,
                            mxArray * * output,
                            int nargout_,
                            mxArray * funfcn_mx,
                            mxArray * x_in,
                            mxArray * printtype_mx,
                            mxArray * tolx_mx,
                            mxArray * tolf_mx,
                            mxArray * maxfun_mx,
                            mxArray * maxiter_mx,
                            mxArray * debugFlags_mx,
                            mxArray * varargin) {

DOTRACE("MdoSimplex");
 
  mexLocalFunctionTable save_local_function_table_ =
	 mclSetCurrentLocalFunctionTable(&_local_function_table_doSimplex);

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
  if (debugFlags_mx && mxIsStruct(debugFlags_mx))
	 {
		mxArray* debugFlag = mxGetField(debugFlags_mx, 0, "debugFlag");

		if (debugFlag)
		  {
			 if (mxGetScalar(debugFlag) == -1) {
				ofstream ofs("profdata.out");
				Util::Prof::printAllProfData(ofs);
				return mxCreateScalarDouble(-1.0);
			 }

			 if (mxGetScalar(debugFlag) == -2) {
				Util::Prof::resetAllProfData();
				return mxCreateScalarDouble(-2.0);
			 }
		  }
	 }
#endif

  try {

	 // numModelParams = prod(size(x));
	 const int numModelParams = mxGetM(x_in) * mxGetN(x_in);

	 // Convert to inline function as needed.
	 // XXX Since this requires "object-oriented" programming, we can't keep this
	 // and still use the MATLAB compiler
	 // %funfcn = fcnchk(funfcn,length(varargin));

	 MatlabFunction objective(funfcn_mx, varargin);

	 SimplexOptimizer opt(objective,
								 Mtx(x_in),
								 Mtx::extractString(printtype_mx),
								 mxGetScalar(tolx_mx),
								 mxGetScalar(tolf_mx),
								 numModelParams,
								 extractMaxIters(maxfun_mx, numModelParams),
								 extractMaxIters(maxiter_mx, numModelParams)
								 );

	 int exitFlag = opt.optimize();

	 *fval = mxCreateScalarDouble(opt.bestFval());
	 *exitflag_mx = mxCreateScalarDouble(exitFlag);

	 mlfIndexAssign(output, ".iterations",
						 mxCreateScalarDouble(opt.iterCount()));

	 mlfIndexAssign(output, ".funcCount",
						 mxCreateScalarDouble(opt.funcCount()));

	 mlfIndexAssign(output, ".algorithm", mxCreateString(opt.algorithm()));

	 mclSetCurrentLocalFunctionTable(save_local_function_table_);

	 return opt.bestParams().makeMxArray();
  }
  catch (ErrorWithMsg& err) {
	 mexErrMsgTxt(err.msg_cstr());
  }
  catch (...) {
	 mexErrMsgTxt("an unknown C++ exception occurred.");
  }

  return (mxArray*) 0; // can't happen, but placate compiler
}

static const char vcid_doSimplex_cc[] = "$Header$";
#endif // !DOSIMPLEX_CC_DEFINED

/*
 * MATLAB Compiler: 2.1
 * Date: Sun Apr  1 19:52:50 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "-h" "doSimplex" 
 */
#include "doSimplex.h"

#include "mtx.h"
#include "strings.h"

#include <fstream.h>
#include "libmatlbm.h"

#include "trace.h"

namespace {
  template <class T>
  inline void localswap(T& t1, T& t2)
  {
	 T t1_(t1);
	 t1 = t2;
	 t2 = t1_;
  }
}

class FuncEvaluator {
  mxArray* itsFunfcn;
  mxArray* itsVarargin_ref;

  static mxArray* getref(mxArray* varargin)
	 {
		return mlfIndexRef(varargin,
								 "{?}",
								 mlfCreateColonIndex());
	 }

public:
  FuncEvaluator(mxArray* funfcn_mx, mxArray* varargin_mx) :
	 itsFunfcn(funfcn_mx),
	 itsVarargin_ref(varargin_mx)
  {
  }

  ~FuncEvaluator()
  {
  }

  mxArray* evaluate(mxArray* x_mx)
  {
	 DOTRACE("evaluate");
	 return mlfFeval(mclValueVarargout(),
						  itsFunfcn,
						  x_mx,
						  getref(itsVarargin_ref),
						  NULL);
  }
};

// ? max(abs(funcVals(1)-funcVals(two2np1))) <= tolf
bool withinTolf(mxArray* funcVals_mx, const double tolf)
{
  const Mtx funcVals_ref(funcVals_mx, Mtx::BORROW);

  MtxConstIter fvals = funcVals_ref.rowIter(0);
  const double f0 = *fvals;
  ++fvals;

  for (; fvals.hasMore(); ++fvals)
	 {
		if (fabs(*fvals - f0) > tolf)
		  return false;
	 }

  return true;
}

// ? max(max(abs(theSimplex(:,two2np1)-theSimplex(:,onesn)))) <= tolx
bool withinTolx(mxArray* simplex_mx, const double tolx)
{
  const Mtx simplex(simplex_mx, Mtx::BORROW);

  const MtxConstIter col0_ = simplex.columnIter(0);

  for (int col = 1; col < simplex.ncols(); ++col)
	 {
		MtxConstIter col0(col0_);
		MtxConstIter coln = simplex.columnIter(col);

		for (; col0.hasMore(); ++col0, ++coln)
		  if ( fabs(*col0 - *coln) > tolx ) return false;
	 }

  return true;
}


static mxChar _array1_[136] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'd', 'o', 'S', 'i', 'm',
                                'p', 'l', 'e', 'x', ' ', 'L', 'i', 'n', 'e',
                                ':', ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm',
                                'n', ':', ' ', '1', ' ', 'T', 'h', 'e', ' ',
                                'f', 'u', 'n', 'c', 't', 'i', 'o', 'n', ' ',
                                '"', 'd', 'o', 'S', 'i', 'm', 'p', 'l', 'e',
                                'x', '"', ' ', 'w', 'a', 's', ' ', 'c', 'a',
                                'l', 'l', 'e', 'd', ' ', 'w', 'i', 't', 'h',
                                ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h', 'a',
                                'n', ' ', 't', 'h', 'e', ' ', 'd', 'e', 'c',
                                'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u', 'm',
                                'b', 'e', 'r', ' ', 'o', 'f', ' ', 'o', 'u',
                                't', 'p', 'u', 't', 's', ' ', '(', '4', ')',
                                '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[21] = { '2', '0', '0', '*', 'n', 'u', 'm',
                               'b', 'e', 'r', 'o', 'f', 'v', 'a',
                               'r', 'i', 'a', 'b', 'l', 'e', 's' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;

static mxChar _array6_[65] = { 'O', 'p', 't', 'i', 'o', 'n', ' ', 0x0027, 'M',
                               'a', 'x', 'F', 'u', 'n', 'E', 'v', 'a', 'l', 's',
                               0x0027, ' ', 'm', 'u', 's', 't', ' ', 'b', 'e',
                               ' ', 'a', 'n', ' ', 'i', 'n', 't', 'e', 'g', 'e',
                               'r', ' ', 'v', 'a', 'l', 'u', 'e', ' ', 'i', 'f',
                               ' ', 'n', 'o', 't', ' ', 't', 'h', 'e', ' ', 'd',
                               'e', 'f', 'a', 'u', 'l', 't', '.' };
static mxArray * _mxarray5_;

static mxChar _array8_[61] = { 'O', 'p', 't', 'i', 'o', 'n', ' ', 0x0027,
                               'M', 'a', 'x', 'I', 't', 'e', 'r', 0x0027,
                               ' ', 'm', 'u', 's', 't', ' ', 'b', 'e', ' ',
                               'a', 'n', ' ', 'i', 'n', 't', 'e', 'g', 'e',
                               'r', ' ', 'v', 'a', 'l', 'u', 'e', ' ', 'i',
                               'f', ' ', 'n', 'o', 't', ' ', 't', 'h', 'e',
                               ' ', 'd', 'e', 'f', 'a', 'u', 'l', 't', '.' };
static mxArray * _mxarray7_;

static mxChar _array10_[6] = { 'n', 'o', 't', 'i', 'f', 'y' };
static mxArray * _mxarray9_;
static mxArray * _mxarray11_;

static mxChar _array15_[4] = { 'n', 'o', 'n', 'e' };
static mxArray * _mxarray14_;

static mxChar _array17_[3] = { 'o', 'f', 'f' };
static mxArray * _mxarray16_;

static mxArray * _array13_[2] = { NULL /*_mxarray14_*/, NULL /*_mxarray16_*/ };
static mxArray * _mxarray12_;
static mxArray * _mxarray18_;

static mxChar _array20_[4] = { 'i', 't', 'e', 'r' };
static mxArray * _mxarray19_;
static mxArray * _mxarray21_;

static mxChar _array23_[5] = { 'f', 'i', 'n', 'a', 'l' };
static mxArray * _mxarray22_;
static mxArray * _mxarray24_;

static mxChar _array26_[7] = { 's', 'i', 'm', 'p', 'l', 'e', 'x' };
static mxArray * _mxarray25_;
static mxArray * _mxarray27_;

static mxChar _array29_[54] = { ' ', 'I', 't', 'e', 'r', 'a', 't', 'i', 'o',
                                'n', ' ', ' ', ' ', 'F', 'u', 'n', 'c', '-',
                                'c', 'o', 'u', 'n', 't', ' ', ' ', ' ', ' ',
                                ' ', 'm', 'i', 'n', ' ', 'f', '(', 'x', ')',
                                ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                'P', 'r', 'o', 'c', 'e', 'd', 'u', 'r', 'e' };
static mxArray * _mxarray28_;
static mxArray * _mxarray30_;
static mxArray * _mxarray31_;
static mxArray * _mxarray32_;
static mxArray * _mxarray33_;

static mxChar _array35_[7] = { 'i', 'n', 'i', 't', 'i', 'a', 'l' };
static mxArray * _mxarray34_;

static mxChar _array37_[1] = { ' ' };
static mxArray * _mxarray36_;

static mxChar _array39_[39] = { ' ', '%', '5', '.', '0', 'f', ' ', ' ',
                                ' ', ' ', ' ', ' ', ' ', ' ', '%', '5',
                                '.', '0', 'f', ' ', ' ', ' ', ' ', ' ',
                                '%', '1', '2', '.', '6', 'g', ' ', ' ',
                                ' ', ' ', ' ', ' ', ' ', ' ', ' ' };
static mxArray * _mxarray38_;

static mxChar _array43_[6] = { 'f', 'o', 'r', 'm', 'a', 't' };
static mxArray * _mxarray42_;

static mxChar _array45_[13] = { 'f', 'o', 'r', 'm', 'a', 't', 's',
                                'p', 'a', 'c', 'i', 'n', 'g' };
static mxArray * _mxarray44_;

static mxArray * _array41_[2] = { NULL /*_mxarray42_*/, NULL /*_mxarray44_*/ };
static mxArray * _mxarray40_;

static mxChar _array47_[7] = { 'c', 'o', 'm', 'p', 'a', 'c', 't' };
static mxArray * _mxarray46_;

static mxChar _array49_[5] = { 's', 'h', 'o', 'r', 't' };
static mxArray * _mxarray48_;

static mxChar _array51_[1] = { 'e' };
static mxArray * _mxarray50_;
static mxArray * _mxarray52_;

static mxChar _array54_[6] = { 'e', 'x', 'p', 'a', 'n', 'd' };
static mxArray * _mxarray53_;

static mxChar _array56_[7] = { 'r', 'e', 'f', 'l', 'e', 'c', 't' };
static mxArray * _mxarray55_;

static mxChar _array58_[16] = { 'c', 'o', 'n', 't', 'r', 'a', 'c', 't',
                                ' ', 'o', 'u', 't', 's', 'i', 'd', 'e' };
static mxArray * _mxarray57_;

static mxChar _array60_[6] = { 's', 'h', 'r', 'i', 'n', 'k' };
static mxArray * _mxarray59_;

static mxChar _array62_[15] = { 'c', 'o', 'n', 't', 'r', 'a', 'c', 't',
                                ' ', 'i', 'n', 's', 'i', 'd', 'e' };
static mxArray * _mxarray61_;

static mxChar _array64_[33] = { 'N', 'e', 'l', 'd', 'e', 'r', '-', 'M', 'e',
                                'a', 'd', ' ', 's', 'i', 'm', 'p', 'l', 'e',
                                'x', ' ', 'd', 'i', 'r', 'e', 'c', 't', ' ',
                                's', 'e', 'a', 'r', 'c', 'h' };
static mxArray * _mxarray63_;

static mxChar _array66_[65] = { 'E', 'x', 'i', 't', 'i', 'n', 'g', ':', ' ',
                                'M', 'a', 'x', 'i', 'm', 'u', 'm', ' ', 'n',
                                'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                'f', 'u', 'n', 'c', 't', 'i', 'o', 'n', ' ',
                                'e', 'v', 'a', 'l', 'u', 'a', 't', 'i', 'o',
                                'n', 's', ' ', 'h', 'a', 's', ' ', 'b', 'e',
                                'e', 'n', ' ', 'e', 'x', 'c', 'e', 'e', 'd',
                                'e', 'd' };
static mxArray * _mxarray65_;

static mxChar _array68_[39] = { ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                ' ', '-', ' ', 'i', 'n', 'c', 'r', 'e',
                                'a', 's', 'e', ' ', 'M', 'a', 'x', 'F',
                                'u', 'n', 'E', 'v', 'a', 'l', 's', ' ',
                                'o', 'p', 't', 'i', 'o', 'n', '.' };
static mxArray * _mxarray67_;

static mxChar _array70_[38] = { ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                ' ', 'C', 'u', 'r', 'r', 'e', 'n', 't',
                                ' ', 'f', 'u', 'n', 'c', 't', 'i', 'o',
                                'n', ' ', 'v', 'a', 'l', 'u', 'e', ':',
                                ' ', '%', 'f', ' ', 0x005c, 'n' };
static mxArray * _mxarray69_;

static mxChar _array72_[55] = { 'E', 'x', 'i', 't', 'i', 'n', 'g', ':',
                                ' ', 'M', 'a', 'x', 'i', 'm', 'u', 'm',
                                ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ',
                                'o', 'f', ' ', 'i', 't', 'e', 'r', 'a',
                                't', 'i', 'o', 'n', 's', ' ', 'h', 'a',
                                's', ' ', 'b', 'e', 'e', 'n', ' ', 'e',
                                'x', 'c', 'e', 'e', 'd', 'e', 'd' };
static mxArray * _mxarray71_;

static mxChar _array74_[35] = { ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                '-', ' ', 'i', 'n', 'c', 'r', 'e', 'a', 's',
                                'e', ' ', 'M', 'a', 'x', 'I', 't', 'e', 'r',
                                ' ', 'o', 'p', 't', 'i', 'o', 'n', '.' };
static mxArray * _mxarray73_;

static mxChar _array76_[192] = { 0x005c, 'n', 'O', 'p', 't', 'i', 'm', 'i',
                                 'z', 'a', 't', 'i', 'o', 'n', ' ', 't', 'e',
                                 'r', 'm', 'i', 'n', 'a', 't', 'e', 'd', ' ',
                                 's', 'u', 'c', 'c', 'e', 's', 's', 'f', 'u',
                                 'l', 'l', 'y', ':', 0x005c, 'n', ' ', 't',
                                 'h', 'e', ' ', 'c', 'u', 'r', 'r', 'e', 'n',
                                 't', ' ', 'x', ' ', 's', 'a', 't', 'i', 's',
                                 'f', 'i', 'e', 's', ' ', 't', 'h', 'e', ' ',
                                 't', 'e', 'r', 'm', 'i', 'n', 'a', 't', 'i',
                                 'o', 'n', ' ', 'c', 'r', 'i', 't', 'e', 'r',
                                 'i', 'a', ' ', 'u', 's', 'i', 'n', 'g', ' ',
                                 'O', 'P', 'T', 'I', 'O', 'N', 'S', '.', 'T',
                                 'o', 'l', 'X', ' ', 'o', 'f', ' ', '%', 'e',
                                 ' ', 0x005c, 'n', ' ', 'a', 'n', 'd', ' ',
                                 'F', '(', 'X', ')', ' ', 's', 'a', 't', 'i',
                                 's', 'f', 'i', 'e', 's', ' ', 't', 'h', 'e',
                                 ' ', 'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e',
                                 'n', 'c', 'e', ' ', 'c', 'r', 'i', 't', 'e',
                                 'r', 'i', 'a', ' ', 'u', 's', 'i', 'n', 'g',
                                 ' ', 'O', 'P', 'T', 'I', 'O', 'N', 'S', '.',
                                 'T', 'o', 'l', 'F', 'u', 'n', ' ', 'o', 'f',
                                 ' ', '%', 'e', ' ', 0x005c, 'n' };
static mxArray * _mxarray75_;

void InitializeModule_doSimplex(void) {
    _mxarray0_ = mclInitializeString(136, _array1_);
    _mxarray2_ = mclInitializeString(21, _array3_);
    _mxarray4_ = mclInitializeDouble(200.0);
    _mxarray5_ = mclInitializeString(65, _array6_);
    _mxarray7_ = mclInitializeString(61, _array8_);
    _mxarray9_ = mclInitializeString(6, _array10_);
    _mxarray11_ = mclInitializeDouble(1.0);
    _mxarray14_ = mclInitializeString(4, _array15_);
    _array13_[0] = _mxarray14_;
    _mxarray16_ = mclInitializeString(3, _array17_);
    _array13_[1] = _mxarray16_;
    _mxarray12_ = mclInitializeCellVector(1, 2, _array13_);
    _mxarray18_ = mclInitializeDouble(0.0);
    _mxarray19_ = mclInitializeString(4, _array20_);
    _mxarray21_ = mclInitializeDouble(3.0);
    _mxarray22_ = mclInitializeString(5, _array23_);
    _mxarray24_ = mclInitializeDouble(2.0);
    _mxarray25_ = mclInitializeString(7, _array26_);
    _mxarray27_ = mclInitializeDouble(4.0);
    _mxarray28_ = mclInitializeString(54, _array29_);
    _mxarray30_ = mclInitializeDouble(.5);
    _mxarray31_ = mclInitializeDouble(.05);
    _mxarray32_ = mclInitializeDouble(.00025);
    _mxarray33_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray34_ = mclInitializeString(7, _array35_);
    _mxarray36_ = mclInitializeString(1, _array37_);
    _mxarray38_ = mclInitializeString(39, _array39_);
    _mxarray42_ = mclInitializeString(6, _array43_);
    _array41_[0] = _mxarray42_;
    _mxarray44_ = mclInitializeString(13, _array45_);
    _array41_[1] = _mxarray44_;
    _mxarray40_ = mclInitializeCellVector(1, 2, _array41_);
    _mxarray46_ = mclInitializeString(7, _array47_);
    _mxarray48_ = mclInitializeString(5, _array49_);
    _mxarray50_ = mclInitializeString(1, _array51_);
    _mxarray52_ = mclInitializeCharVector(0, 0, (mxChar *)NULL);
    _mxarray53_ = mclInitializeString(6, _array54_);
    _mxarray55_ = mclInitializeString(7, _array56_);
    _mxarray57_ = mclInitializeString(16, _array58_);
    _mxarray59_ = mclInitializeString(6, _array60_);
    _mxarray61_ = mclInitializeString(15, _array62_);
    _mxarray63_ = mclInitializeString(33, _array64_);
    _mxarray65_ = mclInitializeString(65, _array66_);
    _mxarray67_ = mclInitializeString(39, _array68_);
    _mxarray69_ = mclInitializeString(38, _array70_);
    _mxarray71_ = mclInitializeString(55, _array72_);
    _mxarray73_ = mclInitializeString(35, _array74_);
    _mxarray75_ = mclInitializeString(192, _array76_);
}

void TerminateModule_doSimplex(void) {
    mxDestroyArray(_mxarray75_);
    mxDestroyArray(_mxarray73_);
    mxDestroyArray(_mxarray71_);
    mxDestroyArray(_mxarray69_);
    mxDestroyArray(_mxarray67_);
    mxDestroyArray(_mxarray65_);
    mxDestroyArray(_mxarray63_);
    mxDestroyArray(_mxarray61_);
    mxDestroyArray(_mxarray59_);
    mxDestroyArray(_mxarray57_);
    mxDestroyArray(_mxarray55_);
    mxDestroyArray(_mxarray53_);
    mxDestroyArray(_mxarray52_);
    mxDestroyArray(_mxarray50_);
    mxDestroyArray(_mxarray48_);
    mxDestroyArray(_mxarray46_);
    mxDestroyArray(_mxarray40_);
    mxDestroyArray(_mxarray44_);
    mxDestroyArray(_mxarray42_);
    mxDestroyArray(_mxarray38_);
    mxDestroyArray(_mxarray36_);
    mxDestroyArray(_mxarray34_);
    mxDestroyArray(_mxarray33_);
    mxDestroyArray(_mxarray32_);
    mxDestroyArray(_mxarray31_);
    mxDestroyArray(_mxarray30_);
    mxDestroyArray(_mxarray28_);
    mxDestroyArray(_mxarray27_);
    mxDestroyArray(_mxarray25_);
    mxDestroyArray(_mxarray24_);
    mxDestroyArray(_mxarray22_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
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
        mlfError(_mxarray0_);
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

static mxArray * doSimplexImpl(mxArray * * fval,
										 mxArray * * exitflag,
										 mxArray * * output,
										 int nargout_,
										 mxArray * funfcn_mx,
										 mxArray * x_in,
										 mxArray * printtype,
										 mxArray * tolx_mx,
										 mxArray * tolf_mx,
										 mxArray * maxfun_mx,
										 mxArray * maxiter_mx,
										 mxArray * debugFlags_mx,
										 mxArray * varargin)
{

    mxArray * x = mclGetUninitializedArray();
    mxArray * convmsg1 = mclGetUninitializedArray();
    mxArray * msg = mclGetUninitializedArray();
    mxArray * fxcc = mclGetUninitializedArray();
    mxArray * xcc = mclGetUninitializedArray();
    mxArray * fxc = mclGetUninitializedArray();
    mxArray * xc = mclGetUninitializedArray();
    mxArray * fxe = mclGetUninitializedArray();
    mxArray * xe = mclGetUninitializedArray();
    mxArray * fxr = mclGetUninitializedArray();
    mxArray * xr = mclGetUninitializedArray();
    mxArray * xbar = mclGetUninitializedArray();
    mxArray * formatsave = mclGetUninitializedArray();
    mxArray * itercount = mclGetUninitializedArray();
    mxArray * how = mclGetUninitializedArray();
    mxArray * f = mclGetUninitializedArray();
    mxArray * y = mclGetUninitializedArray();
    mxArray * j = mclGetUninitializedArray();
    mxArray * zero_term_delta = mclGetUninitializedArray();
    mxArray * usual_delta = mclGetUninitializedArray();
    mxArray * funcVals = mclGetUninitializedArray();
    mxArray * theSimplex = mclGetUninitializedArray();
    mxArray * initialParams = mclGetUninitializedArray();
    mxArray * two2np1 = mclGetUninitializedArray();
    mxArray * onesn = mclGetUninitializedArray();
    mxArray * sigma = mclGetUninitializedArray();
    mxArray * psi = mclGetUninitializedArray();
    mxArray * chi = mclGetUninitializedArray();
    mxArray * rho = mclGetUninitializedArray();
    mxArray * header = mclGetUninitializedArray();
    mxArray * prnt = mclGetUninitializedArray();
    mxArray * ans = mclGetUninitializedArray();

    mclCopyArray(&funfcn_mx);
    mclCopyInputArg(&x, x_in);
    mclCopyArray(&printtype);
    mclCopyArray(&tolx_mx);
    mclCopyArray(&tolf_mx);
    mclCopyArray(&maxfun_mx);
    mclCopyArray(&maxiter_mx);
    mclCopyArray(&debugFlags_mx);
    mclCopyArray(&varargin);


	 const double tolx = mxGetScalar(tolx_mx);
	 const double tolf = mxGetScalar(tolf_mx);

	 FuncEvaluator fevaluator(funfcn_mx, varargin);

    // numModelParams = prod(size(x));
	 const int numModelParams = mxGetM(x) * mxGetN(x);

    mxArray* const numModelParams_mx =
		mclInitialize(mxCreateScalarDouble(numModelParams));

	 int maxfun = 0;

    // In case the defaults were gathered from calling: optimset('simplex'):
    if (mxIsChar(maxfun_mx))
		{
		  if (Mtx::extractString(maxfun_mx) == "200*numberofvariables")
			 {
				maxfun = 200*numModelParams;
			 }
		  else
			 {
				mexErrMsgTxt("Option 'MaxFunEvals' must be an integer value "
								 "if not the default.");
			 }
		}
	 else
		{
		  maxfun = int(mxGetScalar(maxfun_mx));
		}

	 int maxiter = 0;

    // if ischar(maxiter)
    if (mxIsChar(maxiter_mx))
		{
        if (Mtx::extractString(maxiter_mx) == "200*numberofvariables")
			 {
				maxiter = 200*numModelParams;
			 }
		  else
			 {
            // error('')
				mexErrMsgTxt("Option 'MaxIter' must be an integer value "
								 "if not the default.");
        }
    }
	 else
		{
		  maxiter = int(mxGetScalar(maxiter_mx));
		}

    // 
    // switch printtype
    {
        mxArray * v_ = mclInitialize(mclVa(printtype, "printtype"));
        if (mclSwitchCompare(v_, _mxarray9_)) {

            // case 'notify'
            // prnt = 1;
            mlfAssign(&prnt, _mxarray11_);

        // case {'none','off'}
        } else if (mclSwitchCompare(v_, _mxarray12_)) {

            // prnt = 0;
            mlfAssign(&prnt, _mxarray18_);

        // case 'iter'
        } else if (mclSwitchCompare(v_, _mxarray19_)) {

            // prnt = 3;
            mlfAssign(&prnt, _mxarray21_);

        // case 'final'
        } else if (mclSwitchCompare(v_, _mxarray22_)) {

            // prnt = 2;
            mlfAssign(&prnt, _mxarray24_);

        // case 'simplex'
        } else if (mclSwitchCompare(v_, _mxarray25_)) {

            // prnt = 4;
            mlfAssign(&prnt, _mxarray27_);

        // otherwise
        } else {

            // prnt = 1;
            mlfAssign(&prnt, _mxarray11_);

        // end
        }
        mxDestroyArray(v_);
    }

    // 
    // header = ' Iteration   Func-count     min f(x)         Procedure';
    mlfAssign(&header, _mxarray28_);

    // 
    // % Convert to inline function as needed.
    // % XXX Since this requires "object-oriented" programming, we can't keep this
    // % and still use the MATLAB compiler
    // %funfcn = fcnchk(funfcn,length(varargin));
    // 
    // % Initialize parameters
    // rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    mlfAssign(&rho, _mxarray11_);
    mlfAssign(&chi, _mxarray24_);
    mlfAssign(&psi, _mxarray30_);
    mlfAssign(&sigma, _mxarray30_);

    // onesn = ones(1,numModelParams);
    mlfAssign(
      &onesn,
      mlfOnes(_mxarray11_, numModelParams_mx, NULL));

    // two2np1 = 2:numModelParams+1;
    mlfAssign(
      &two2np1,
      mlfColon(
        _mxarray24_,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // 
    // % Set up a simplex near the initial guess.
    // initialParams = x(:); % Force initialParams to be a column vector
    mlfAssign(
      &initialParams, mclArrayRef1(mclVsa(x, "x"), mlfCreateColonIndex()));

    // theSimplex = zeros(numModelParams,numModelParams+1);
    mlfAssign(
      &theSimplex,
      mlfZeros(
        numModelParams_mx,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // funcVals = zeros(1,numModelParams+1);
    mlfAssign(
      &funcVals,
      mlfZeros(
        _mxarray11_,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // theSimplex(:,1) = initialParams;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
    mclArrayAssign2(
      &theSimplex,
      mclVsv(initialParams, "initialParams"),
      mlfCreateColonIndex(),
      _mxarray11_);

    // funcVals(:,1) = feval(funfcn,initialParams,varargin{:}); 
    mclArrayAssign2(
      &funcVals,
      mlfFeval(
        mclValueVarargout(),
        mclVa(funfcn_mx, "funfcn"),
        mclVv(initialParams, "initialParams"),
        mclVe(
          mlfIndexRef(
            mclVsa(varargin, "varargin"), "{?}", mlfCreateColonIndex())),
        NULL),
      mlfCreateColonIndex(),
      _mxarray11_);

    // 
    // % Following improvement suggested by L.Pfeffer at Stanford
    // usual_delta = 0.05;             % 5 percent deltas for non-zero terms
    mlfAssign(&usual_delta, _mxarray31_);

    // zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
    mlfAssign(&zero_term_delta, _mxarray32_);

    // for j = 1:numModelParams
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(numModelParams_mx);
        if (v_ > e_) {
            mlfAssign(&j, _mxarray33_);
        } else {

            // y = initialParams;
            // if y(j) ~= 0
            // y(j) = (1 + usual_delta)*y(j);
            // else 
            // y(j) = zero_term_delta;
            // end  
            // theSimplex(:,j+1) = y;
            // f = feval(funfcn,y,varargin{:});
            // funcVals(1,j+1) = f;
            // end     
            for (; ; ) {
                mlfAssign(&y, mclVsv(initialParams, "initialParams"));
                if (mclNeBool(
                      mclVe(mclIntArrayRef1(mclVsv(y, "y"), v_)),
                      _mxarray18_)) {
                    mclIntArrayAssign1(
                      &y,
                      mclMtimes(
                        mclPlus(_mxarray11_, mclVv(usual_delta, "usual_delta")),
                        mclVe(mclIntArrayRef1(mclVsv(y, "y"), v_))),
                      v_);
                } else {
                    mclIntArrayAssign1(
                      &y, mclVsv(zero_term_delta, "zero_term_delta"), v_);
                }
                mclArrayAssign2(
                  &theSimplex,
                  mclVsv(y, "y"),
                  mlfCreateColonIndex(),
                  mlfScalar(v_ + 1));
                mlfAssign(
                  &f,
                  mlfFeval(
                    mclValueVarargout(),
                    mclVa(funfcn_mx, "funfcn"),
                    mclVv(y, "y"),
                    mclVe(
                      mlfIndexRef(
                        mclVsa(varargin, "varargin"),
                        "{?}",
                        mlfCreateColonIndex())),
                    NULL));
                mclIntArrayAssign2(&funcVals, mclVsv(f, "f"), 1, v_ + 1);
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&j, mlfScalar(v_));
        }
    }

    // 
    // % sort so theSimplex(1,:) has the lowest function value 
    // [funcVals,j] = sort(funcVals);
    mlfAssign(&funcVals, mlfSort(&j, mclVv(funcVals, "funcVals"), NULL));

    // theSimplex = theSimplex(:,j);
    mlfAssign(
      &theSimplex,
      mclArrayRef2(
        mclVsv(theSimplex, "theSimplex"),
        mlfCreateColonIndex(),
        mclVsv(j, "j")));

    // 
    // how = 'initial';
    mlfAssign(&how, _mxarray34_);

    // itercount = 1;
    mlfAssign(&itercount, _mxarray11_);

	 int func_evals = numModelParams+1;

    // if prnt == 3
    if (mclEqBool(mclVv(prnt, "prnt"), _mxarray21_)) {

        // disp(' ')
        mlfDisp(_mxarray36_);

        // disp(header)
        mlfDisp(mclVv(header, "header"));

        // disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount, func_evals, funcVals(1)), how]) 
        mlfDisp(
          mlfHorzcat(
            mclVe(
              mlfSprintf(
                NULL,
                _mxarray38_,
                mclVv(itercount, "itercount"),
                mxCreateScalarDouble(func_evals),
                mclVe(mclIntArrayRef1(mclVsv(funcVals, "funcVals"), 1)),
                NULL)),
            mclVv(how, "how"),
            NULL));

    // elseif prnt == 4
    } else if (mclEqBool(mclVv(prnt, "prnt"), _mxarray27_)) {

        // clc
        mlfClc();

        // formatsave = get(0,{'format','formatspacing'});
        mlfAssign(&formatsave, mlfNGet(1, _mxarray18_, _mxarray40_, NULL));

        // format compact
        mlfFormat(_mxarray46_, NULL);

        // format short e
        mlfFormat(_mxarray48_, _mxarray50_);

        // disp(' ')
        mlfDisp(_mxarray36_);

        // disp(how)
        mlfDisp(mclVv(how, "how"));

        // theSimplex
        mclPrintArray(mclVsv(theSimplex, "theSimplex"), "theSimplex");

        // funcVals
        mclPrintArray(mclVsv(funcVals, "funcVals"), "funcVals");

        // func_evals
        mclPrintArray(mxCreateScalarDouble(func_evals), "func_evals");

    // end
    }

    // exitflag = 1;
    mlfAssign(exitflag, _mxarray11_);


	 //---------------------------------------------------------------------
	 //
    // Iterate until the diameter of the simplex is less than tolx AND
    // the function values differ from the min by less than tolf, or
    // the max function evaluations are exceeded. (Cannot use OR
    // instead of AND.)
	 //
	 //---------------------------------------------------------------------


    // while func_evals < maxfun & itercount < maxiter
    for (;;) { // Main algorithm
		DOTRACE("Main algorithm");

		{DOTRACE("Main loop condition"); // 0.5%
//  		mxArray * a_ = mclInitialize(
//  											  mclLt(
//  													  mxCreateScalarDouble(func_evals),
//  													  mclVa(maxfun_mx, "maxfun")));

		bool okEvals = (func_evals < maxfun);

		if (okEvals
			 && mxGetScalar(itercount) < mxGetScalar(maxiter_mx)) {
//  		  mxDestroyArray(a_);
		} else {
//  		  mxDestroyArray(a_);
		  break;
		}
		}

		{DOTRACE("check if done");

		  if (withinTolf(funcVals, tolf) && withinTolx(theSimplex, tolx))
			 break; // main loop
		}

      // how = '';
      mlfAssign(&how, _mxarray52_);

		{DOTRACE("compute reflection point");

      // 
      // Compute the reflection point
      // 

      // xbar = average of the numModelParams (NOT numModelParams+1) best points
      // xbar = sum(theSimplex(:,one2n), 2)/numModelParams;

		{DOTRACE("compute xbar");
      mlfAssign(
        &xbar,
		  mxCreateDoubleMatrix(numModelParams,1,mxREAL));

		Mtx xbar_ref(xbar, Mtx::REFER);
		const Mtx simplex(Mtx(theSimplex, Mtx::BORROW)
								.columns(0,numModelParams));

		const double numparams_inv = 1.0/numModelParams;
		MtxIter xbar_itr_ = xbar_ref.columnIter(0);

		for (int r = 0; r < numModelParams; ++r, ++xbar_itr_)
		  *xbar_itr_ = simplex.row(r).sum() * numparams_inv;

		}

		{DOTRACE("compute xr");
      // xr = (1 + rho)*xbar - rho*theSimplex(:,end);
      mlfAssign(
        &xr,
        mclMinus(
          mclMtimes(
            mclPlus(_mxarray11_, mclVv(rho, "rho")), mclVv(xbar, "xbar")),
          mclMtimes(
            mclVv(rho, "rho"),
            mclVe(
              mclArrayRef2(
                theSimplex,
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(theSimplex, "theSimplex"),
                  _mxarray24_,
                  _mxarray24_))))));
		}

		mlfAssign(&fxr, fevaluator.evaluate(xr));

		++func_evals;
		}

      // 
      // if fxr < funcVals(:,1)
      if (mclLtBool(
            mclVv(fxr, "fxr"),
            mclVe(
              mclArrayRef2(
                mclVsv(funcVals, "funcVals"),
                mlfCreateColonIndex(),
                _mxarray11_)))) {
			 DOTRACE("if fxr < funcVals(:,1)");

          // % Calculate the expansion point
          // xe = (1 + rho*chi)*xbar - rho*chi*theSimplex(:,end);
          mlfAssign(
            &xe,
            mclMinus(
              mclMtimes(
                mclPlus(
                  _mxarray11_,
                  mclMtimes(mclVv(rho, "rho"), mclVv(chi, "chi"))),
                mclVv(xbar, "xbar")),
              mclMtimes(
                mclMtimes(mclVv(rho, "rho"), mclVv(chi, "chi")),
                mclVe(
                  mclArrayRef2(
                    theSimplex,
                    mlfCreateColonIndex(),
                    mlfEnd(
                      mclVv(theSimplex, "theSimplex"),
                      _mxarray24_,
                      _mxarray24_))))));


          mlfAssign(&fxe, fevaluator.evaluate(xe));

			 ++func_evals;

          // if fxe < fxr
          if (mclLtBool(mclVv(fxe, "fxe"), mclVv(fxr, "fxr"))) {

              // theSimplex(:,end) = xe;
              mclArrayAssign2(
                &theSimplex,
                mclVsv(xe, "xe"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(theSimplex, "theSimplex"), _mxarray24_, _mxarray24_));

              // funcVals(:,end) = fxe;
              mclArrayAssign2(
                &funcVals,
                mclVsv(fxe, "fxe"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals, "funcVals"), _mxarray24_, _mxarray24_));

              // how = 'expand';
              mlfAssign(&how, _mxarray53_);

          // else
          } else {

              // theSimplex(:,end) = xr; 
              mclArrayAssign2(
                &theSimplex,
                mclVsv(xr, "xr"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(theSimplex, "theSimplex"), _mxarray24_, _mxarray24_));

              // funcVals(:,end) = fxr;
              mclArrayAssign2(
                &funcVals,
                mclVsv(fxr, "fxr"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals, "funcVals"), _mxarray24_, _mxarray24_));

              // how = 'reflect';
              mlfAssign(&how, _mxarray55_);

          // end
          }

      // else % funcVals(:,1) <= fxr
      } else {
			 DOTRACE("else funcVals(:,1) <= fxr");

          // if fxr < funcVals(:,numModelParams)
          if (mclLtBool(
                mclVv(fxr, "fxr"),
                mclVe(
                  mclArrayRef2(
                    mclVsv(funcVals, "funcVals"),
                    mlfCreateColonIndex(),
                    numModelParams_mx)))) {

              // theSimplex(:,end) = xr; 
              mclArrayAssign2(
                &theSimplex,
                mclVsv(xr, "xr"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(theSimplex, "theSimplex"), _mxarray24_, _mxarray24_));

              // funcVals(:,end) = fxr;
              mclArrayAssign2(
                &funcVals,
                mclVsv(fxr, "fxr"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals, "funcVals"), _mxarray24_, _mxarray24_));

              // how = 'reflect';
              mlfAssign(&how, _mxarray55_);

          // else % fxr >= funcVals(:,numModelParams) 
          } else {

              // % Perform contraction
              // if fxr < funcVals(:,end)
              if (mclLtBool(
                    mclVv(fxr, "fxr"),
                    mclVe(
                      mclArrayRef2(
                        mclVsv(funcVals, "funcVals"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals, "funcVals"),
                          _mxarray24_,
                          _mxarray24_))))) {

                  // % Perform an outside contraction
                  // xc = (1 + psi*rho)*xbar - psi*rho*theSimplex(:,end);
                  mlfAssign(
                    &xc,
                    mclMinus(
                      mclMtimes(
                        mclPlus(
                          _mxarray11_,
                          mclMtimes(mclVv(psi, "psi"), mclVv(rho, "rho"))),
                        mclVv(xbar, "xbar")),
                      mclMtimes(
                        mclMtimes(mclVv(psi, "psi"), mclVv(rho, "rho")),
                        mclVe(
                          mclArrayRef2(
                            theSimplex,
                            mlfCreateColonIndex(),
                            mlfEnd(
                              mclVv(theSimplex, "theSimplex"),
                              _mxarray24_,
                              _mxarray24_))))));


                  mlfAssign(&fxc, fevaluator.evaluate(xc));

						++func_evals;

                  // 
                  // if fxc <= fxr
                  if (mclLeBool(mclVv(fxc, "fxc"), mclVv(fxr, "fxr"))) {

                      // theSimplex(:,end) = xc; 
                      mclArrayAssign2(
                        &theSimplex,
                        mclVsv(xc, "xc"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(theSimplex, "theSimplex"),
                          _mxarray24_,
                          _mxarray24_));

                      // funcVals(:,end) = fxc;
                      mclArrayAssign2(
                        &funcVals,
                        mclVsv(fxc, "fxc"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals, "funcVals"),
                          _mxarray24_,
                          _mxarray24_));

                      // how = 'contract outside';
                      mlfAssign(&how, _mxarray57_);

                  // else
                  } else {

                      // % perform a shrink
                      // how = 'shrink'; 
                      mlfAssign(&how, _mxarray59_);

                  // end
                  }

              // else
              } else {

                  // % Perform an inside contraction
                  // xcc = (1-psi)*xbar + psi*theSimplex(:,end);
                  mlfAssign(
                    &xcc,
                    mclPlus(
                      mclMtimes(
                        mclMinus(_mxarray11_, mclVv(psi, "psi")),
                        mclVv(xbar, "xbar")),
                      mclMtimes(
                        mclVv(psi, "psi"),
                        mclVe(
                          mclArrayRef2(
                            theSimplex,
                            mlfCreateColonIndex(),
                            mlfEnd(
                              mclVv(theSimplex, "theSimplex"),
                              _mxarray24_,
                              _mxarray24_))))));


                  mlfAssign(&fxcc, fevaluator.evaluate(xcc));

						++func_evals;

                  // 
                  // if fxcc < funcVals(:,end)
                  if (mclLtBool(
                        mclVv(fxcc, "fxcc"),
                        mclVe(
                          mclArrayRef2(
                            mclVsv(funcVals, "funcVals"),
                            mlfCreateColonIndex(),
                            mlfEnd(
                              mclVv(funcVals, "funcVals"),
                              _mxarray24_,
                              _mxarray24_))))) {

                      // theSimplex(:,end) = xcc;
                      mclArrayAssign2(
                        &theSimplex,
                        mclVsv(xcc, "xcc"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(theSimplex, "theSimplex"),
                          _mxarray24_,
                          _mxarray24_));

                      // funcVals(:,end) = fxcc;
                      mclArrayAssign2(
                        &funcVals,
                        mclVsv(fxcc, "fxcc"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals, "funcVals"),
                          _mxarray24_,
                          _mxarray24_));

                      // how = 'contract inside';
                      mlfAssign(&how, _mxarray61_);

                  // else
                  } else {

                      // % perform a shrink
                      // how = 'shrink';
                      mlfAssign(&how, _mxarray59_);

                  // end
                  }

              // end
              }

              // if strcmp(how,'shrink')
              if (mlfTobool(
                    mclVe(mlfStrcmp(mclVv(how, "how"), _mxarray59_)))) {
                  mclForLoopIterator viter__;

                  // for j=two2np1
                  for (mclForStart(
                         &viter__, mclVv(two2np1, "two2np1"), NULL, NULL);
                       mclForNext(&viter__, &j);
                       ) {

                      // theSimplex(:,j)=theSimplex(:,1)+sigma*(theSimplex(:,j) - theSimplex(:,1));
                      mclArrayAssign2(
                        &theSimplex,
                        mclPlus(
                          mclVe(
                            mclArrayRef2(
                              theSimplex,
                              mlfCreateColonIndex(),
                              _mxarray11_)),
                          mclMtimes(
                            mclVv(sigma, "sigma"),
                            mclMinus(
                              mclVe(
                                mclArrayRef2(
                                  theSimplex,
                                  mlfCreateColonIndex(),
                                  mclVsv(j, "j"))),
                              mclVe(
                                mclArrayRef2(
                                  theSimplex,
                                  mlfCreateColonIndex(),
                                  _mxarray11_))))),
                        mlfCreateColonIndex(),
                        mclVsv(j, "j"));

                      // funcVals(:,j) = feval(funfcn,theSimplex(:,j),varargin{:});
                      mclArrayAssign2(
                        &funcVals,
								fevaluator.evaluate(mclArrayRef2(theSimplex,
																mlfCreateColonIndex(),
																mclVsv(j, "j"))),
                        mlfCreateColonIndex(),
                        mclVsv(j, "j"));

                  }
                  mclDestroyForLoopIterator(viter__);

						func_evals += numModelParams;

              }

          }

      }

		{DOTRACE("sort funcVals");
      // [funcVals,j] = sort(funcVals);

		// Throw away the actual sorted result; just keep the indices
		mxDestroyArray(mlfSort(&j, funcVals, NULL));
		}

		{DOTRACE("reorder simplex");
		Mtx jref(j, Mtx::BORROW);
		const int smallest = int(jref.at(0)) - 1;
		const int largest = int(jref.at(numModelParams)) - 1;
		const int largest2 = int(jref.at(numModelParams-1)) - 1;

		Mtx simref(theSimplex, Mtx::REFER);
		Mtx fvref(funcVals, Mtx::REFER);

		// These swaps are smart enough to check if the column numbers
		// are the same before doing the swap

		fvref.swapColumns(0, smallest);
		fvref.swapColumns(largest, numModelParams);
		fvref.swapColumns(largest2, numModelParams-1);

		simref.swapColumns(0, smallest);
		simref.swapColumns(largest, numModelParams);
		simref.swapColumns(largest2, numModelParams-1);
		}

      // itercount = itercount + 1;
      mlfAssign(
        &itercount, mclPlus(mclVv(itercount, "itercount"), _mxarray11_));

      // if prnt == 3
      if (mclEqBool(mclVv(prnt, "prnt"), _mxarray21_)) {

          // disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount, func_evals, funcVals(1)), how]) 
          mlfDisp(
            mlfHorzcat(
              mclVe(
                mlfSprintf(
                  NULL,
                  _mxarray38_,
                  mclVv(itercount, "itercount"),
                  mxCreateScalarDouble(func_evals),
                  mclVe(mclIntArrayRef1(mclVsv(funcVals, "funcVals"), 1)),
                  NULL)),
              mclVv(how, "how"),
              NULL));

      // elseif prnt == 4
      } else if (mclEqBool(mclVv(prnt, "prnt"), _mxarray27_)) {

          // disp(' ')
          mlfDisp(_mxarray36_);

          // disp(how)
          mlfDisp(mclVv(how, "how"));

          // theSimplex
          mclPrintArray(theSimplex, "theSimplex");

          // funcVals
          mclPrintArray(mclVsv(funcVals, "funcVals"), "funcVals");

          // func_evals
          mclPrintArray(mxCreateScalarDouble(func_evals), "func_evals");

      // end  
      }

    // end   % while
    }

    // 
    // 
    // x(:) = theSimplex(:,1);
    mclArrayAssign1(
      &x,
      mclArrayRef2(
        theSimplex, mlfCreateColonIndex(), _mxarray11_),
      mlfCreateColonIndex());

    // if prnt == 4,
    if (mclEqBool(mclVv(prnt, "prnt"), _mxarray27_)) {

        // % reset format
        // set(0,{'format','formatspacing'},formatsave);
        mclAssignAns(
          &ans,
          mlfNSet(
            0,
            _mxarray18_,
            _mxarray40_,
            mclVv(formatsave, "formatsave"),
            NULL));

    } // end Main algorithm

    // output.iterations = itercount;
    mlfIndexAssign(output, ".iterations", mclVsv(itercount, "itercount"));

    // output.funcCount = func_evals;
    mlfIndexAssign(output, ".funcCount", mxCreateScalarDouble(func_evals));

    // output.algorithm = 'Nelder-Mead simplex direct search';
    mlfIndexAssign(output, ".algorithm", _mxarray63_);

    // 
    // fval = min(funcVals);
    mlfAssign(fval, mlfMin(NULL, mclVv(funcVals, "funcVals"), NULL, NULL));

    // if func_evals >= maxfun 
    if (mclGeBool(mxCreateScalarDouble(func_evals), mclVa(maxfun_mx, "maxfun"))) {

        // if prnt > 0
        if (mclGtBool(mclVv(prnt, "prnt"), _mxarray18_)) {

            // disp(' ')
            mlfDisp(_mxarray36_);

            // disp('Exiting: Maximum number of function evaluations has been exceeded')
            mlfDisp(_mxarray65_);

            // disp('         - increase MaxFunEvals option.')
            mlfDisp(_mxarray67_);

            // msg = sprintf('         Current function value: %f \n', fval);
            mlfAssign(
              &msg, mlfSprintf(NULL, _mxarray69_, mclVv(*fval, "fval"), NULL));

            // disp(msg)
            mlfDisp(mclVv(msg, "msg"));

        // end
        }

        // exitflag = 0;
        mlfAssign(exitflag, _mxarray18_);

    // elseif itercount >= maxiter 
    } else if (mclGeBool(
                 mclVv(itercount, "itercount"), mclVa(maxiter_mx, "maxiter"))) {

        // if prnt > 0
        if (mclGtBool(mclVv(prnt, "prnt"), _mxarray18_)) {

            // disp(' ')
            mlfDisp(_mxarray36_);

            // disp('Exiting: Maximum number of iterations has been exceeded')
            mlfDisp(_mxarray71_);

            // disp('         - increase MaxIter option.')
            mlfDisp(_mxarray73_);

            // msg = sprintf('         Current function value: %f \n', fval);
            mlfAssign(
              &msg, mlfSprintf(NULL, _mxarray69_, mclVv(*fval, "fval"), NULL));

            // disp(msg)
            mlfDisp(mclVv(msg, "msg"));

        // end
        }

        // exitflag = 0; 
        mlfAssign(exitflag, _mxarray18_);

    // else
    } else {

        // if prnt > 1
        if (mclGtBool(mclVv(prnt, "prnt"), _mxarray11_)) {

            // convmsg1 = sprintf([ ...
            mlfAssign(
              &convmsg1,
              mlfSprintf(
                NULL,
                _mxarray75_,
                mclVa(tolx_mx, "tolx"),
                mclVa(tolf_mx, "tolf"),
                NULL));


            // '\nOptimization terminated successfully:\n',...
            // ' the current x satisfies the termination criteria using OPTIONS.TolX of %e \n',...
            // ' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of %e \n'
            // ],tolx, tolf);
            // disp(convmsg1)
            mlfDisp(mclVv(convmsg1, "convmsg1"));


        // end
        }


        // exitflag = 1;
        mlfAssign(exitflag, _mxarray11_);


    // end
    }

    mclVo(&x);
    mclValidateOutput(x, 1, nargout_, "x", "doSimplex");
    mclValidateOutput(*fval, 2, nargout_, "fval", "doSimplex");
    mclValidateOutput(*exitflag, 3, nargout_, "exitflag", "doSimplex");
    mclValidateOutput(*output, 4, nargout_, "output", "doSimplex");

    mxDestroyArray(numModelParams_mx);
    mxDestroyArray(ans);
    mxDestroyArray(prnt);
    mxDestroyArray(header);
    mxDestroyArray(rho);
    mxDestroyArray(chi);
    mxDestroyArray(psi);
    mxDestroyArray(sigma);
    mxDestroyArray(onesn);
    mxDestroyArray(two2np1);
    mxDestroyArray(initialParams);
    mxDestroyArray(theSimplex);
    mxDestroyArray(funcVals);
    mxDestroyArray(usual_delta);
    mxDestroyArray(zero_term_delta);
    mxDestroyArray(j);
    mxDestroyArray(y);
    mxDestroyArray(f);
    mxDestroyArray(how);
    mxDestroyArray(itercount);
    mxDestroyArray(formatsave);
    mxDestroyArray(xbar);
    mxDestroyArray(xr);
    mxDestroyArray(fxr);
    mxDestroyArray(xe);
    mxDestroyArray(fxe);
    mxDestroyArray(xc);
    mxDestroyArray(fxc);
    mxDestroyArray(xcc);
    mxDestroyArray(fxcc);
    mxDestroyArray(msg);
    mxDestroyArray(convmsg1);
    mxDestroyArray(varargin);
    mxDestroyArray(debugFlags_mx);
    mxDestroyArray(maxiter_mx);
    mxDestroyArray(maxfun_mx);
    mxDestroyArray(tolf_mx);
    mxDestroyArray(tolx_mx);
    mxDestroyArray(printtype);
    mxDestroyArray(funfcn_mx);

    return x;
}


/*
 * The function "MdoSimplex" is the implementation version of the "doSimplex"
 * M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * It contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function [x,fval,exitflag,output] = ...
 */

static mxArray * MdoSimplex(mxArray * * fval,
                            mxArray * * exitflag,
                            mxArray * * output,
                            int nargout_,
                            mxArray * funfcn_mx,
                            mxArray * x_in,
                            mxArray * printtype,
                            mxArray * tolx_mx,
                            mxArray * tolf_mx,
                            mxArray * maxfun_mx,
                            mxArray * maxiter_mx,
                            mxArray * debugFlags_mx,
                            mxArray * varargin) {

DOTRACE("MdoSimplex");

// HOME

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

	 mxArray* result = doSimplexImpl(fval, exitflag, output, nargout_,
												funfcn_mx,
												x_in,
												printtype,
												tolx_mx,
												tolf_mx,
												maxfun_mx,
												maxiter_mx,
												debugFlags_mx,
												varargin);

    mclSetCurrentLocalFunctionTable(save_local_function_table_);

	 return result;
}

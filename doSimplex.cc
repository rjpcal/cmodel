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

#define LOCAL_PROF
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

class DualRepMtx {
private:
  mutable Mtx itsMtx;
  mutable bool itsMtxIsValid;

  mutable mxArray* itsArray;
  mutable bool itsArrayIsValid;

  void updateMtx() const
  {
	 itsMtx = Mtx(itsArray);
	 itsMtxIsValid = true;
  }

  void updateArray() const
  {
	 mlfAssign(&itsArray, itsMtx.makeMxArray());
	 itsArrayIsValid = true;
  }

public:
  DualRepMtx() :
	 itsMtx(0,0),
	 itsMtxIsValid(true),
	 itsArray(mclGetUninitializedArray()),
	 itsArrayIsValid(false)
  {}

  DualRepMtx(const Mtx& other) :
	 itsMtx(other),
	 itsMtxIsValid(true),
	 itsArray(mclGetUninitializedArray()),
	 itsArrayIsValid(false)
  {}

  ~DualRepMtx()
  {
	 mxDestroyArray(itsArray);
  }

  void assignArray(mxArray* rhs)
  {
	 mlfAssign(&itsArray, rhs);
	 itsArrayIsValid = true;
	 itsMtxIsValid = false;
  }

  void assignMtx(const Mtx& rhs)
  {
	 itsMtx = rhs;
	 itsMtxIsValid = true;
	 itsArrayIsValid = false;
  }

  mxArray* asArray() const
  {
	 if (!itsArrayIsValid) updateArray();
	 return itsArray;
  }

  const Mtx& asMtx() const
  {
	 if (!itsMtxIsValid) updateMtx();
	 return itsMtx;
  }

  Mtx& ncMtx()
  {
	 if (!itsMtxIsValid) updateMtx();
	 itsArrayIsValid = false;
	 return itsMtx;
  }
};

class FuncEvaluator {
  int itsEvalCount;
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
	 itsEvalCount(0),
	 itsFunfcn(funfcn_mx),
	 itsVarargin_ref(varargin_mx)
  {
  }

  ~FuncEvaluator()
  {
  }

  int evalCount() const { return itsEvalCount; }

  mxArray* evaluate_mx(mxArray* x_mx)
  {
	 DOTRACE("evaluate_mx");

	 ++itsEvalCount;

	 return mlfFeval(mclValueVarargout(),
						  itsFunfcn,
						  x_mx,
						  getref(itsVarargin_ref),
						  NULL);
  }

  double evaluate(mxArray* x_mx)
  {
	 DOTRACE("evaluate");

	 ++itsEvalCount;

	 mxArray* mx =  mlfFeval(mclValueVarargout(),
									 itsFunfcn,
									 x_mx,
									 getref(itsVarargin_ref),
									 NULL);
	 double result = mxGetScalar(mx);
	 mxDestroyArray(mx);
	 return result;
  }

  double evaluate(const Mtx& x)
  {
	 return evaluate(x.makeMxArray());
  }
};

// ? max(abs(funcVals(1)-funcVals(two2np1))) <= tolf
bool withinTolf(const Mtx& funcVals, const double tolf)
{
  MtxConstIter fvals = funcVals.rowIter(0);
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
bool withinTolx(const Mtx& simplex, const double tolx)
{
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

int extractPrinttype(const mxArray* printtype_mx)
{
  fixed_string printtype = Mtx::extractString(printtype_mx);

  if (printtype == "notify")
	 {
		return 1;
	 }
  else if (printtype == "none" || printtype == "off")
	 {
		return 0;
	 }
  else if (printtype == "iter")
	 {
		return 3;
	 }
  else if (printtype == "final")
	 {
		return 2;
	 }
  else if (printtype == "simplex")
	 {
		return 4;
	 }

  return 1;
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
										 FuncEvaluator& objective,
										 const Mtx& x_in,
										 const int prnt,
										 const double tolx,
										 const double tolf,
										 const int numModelParams,
										 const int maxfun,
										 const int maxiter)
{

// HOME

  double fxcc = 0.0;
  double fxc = 0.0;
  double fxe = 0.0;
  double fxr = 0.0;

  mxArray* j_mx = mclGetUninitializedArray();

  fixed_string how;

  DualRepMtx xr_dr;
  DualRepMtx xe_dr;
  DualRepMtx xc_dr;
  DualRepMtx xcc_dr;
  DualRepMtx xbar_dr;

  Mtx xx(x_in);


  const double rho = 1.0;
  const double chi = 2.0;
  const double psi = 0.5;
  const double sigma = 0.5;

  // Set up a simplex near the initial guess.
  DualRepMtx initialParams_dr(xx.asColumn());

  // theSimplex = zeros(numModelParams,numModelParams+1);
  DualRepMtx theSimplex_dr(Mtx(numModelParams,numModelParams+1));

  // funcVals = zeros(1,numModelParams+1);
  DualRepMtx funcVals_dr(Mtx(1,numModelParams+1));

  // theSimplex(:,1) = initialParams;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
  theSimplex_dr.ncMtx().column(0) = initialParams_dr.asMtx().columns(0,1);

  // funcVals(:,1) = feval(funfcn,initialParams,varargin{:}); 
  funcVals_dr.ncMtx().at(0,0) = objective.evaluate(initialParams_dr.asMtx());

  {
	 // Following improvement suggested by L.Pfeffer at Stanford
	 // 5 percent deltas for non-zero terms
	 const double usual_delta = 0.05;

	 // Even smaller delta for zero elements of x
	 const double zero_term_delta = 0.00025;

	 for (int j_zero_based = 0; j_zero_based < numModelParams;
			++j_zero_based)
		{
		  Mtx y = initialParams_dr.asMtx();

		  if (y.at(j_zero_based) != 0.0)
			 {
				y.at(j_zero_based) *= (1.0 + usual_delta);
			 }
		  else
			 {
				y.at(j_zero_based) = zero_term_delta;
			 }

		  theSimplex_dr.ncMtx().column(j_zero_based+1) = y;

		  funcVals_dr.ncMtx().at(0,j_zero_based+1) = objective.evaluate(y);
		}
  }

  // 
  // % sort so theSimplex(1,:) has the lowest function value 
  // [funcVals,j_mx] = sort(funcVals);
  funcVals_dr.assignArray(mlfSort(&j_mx, funcVals_dr.asArray(), NULL));

  // theSimplex = theSimplex(:,j_mx);
  theSimplex_dr.assignArray(mclArrayRef2(theSimplex_dr.asArray(),
													  mlfCreateColonIndex(),
													  j_mx));

  how = "initial";

  int itercount = 1;

  if (prnt == 3)
	 {
		mexPrintf("\n Iteration   Func-count     min f(x)         Procedure\n");
		mexPrintf(" %5d        %5d     %12.6g         ",
					 itercount, objective.evalCount(), double(funcVals_dr.asMtx().at(0,0)));
		mexPrintf(how.c_str());
		mexPrintf("\n");
	 }
  else if (prnt == 4)
	 {
		mexPrintf("\n%s\n", how.c_str());
		theSimplex_dr.asMtx().print("theSimplex");
		funcVals_dr.asMtx().print("funcVals");
		mexPrintf("func_evals: %d\n", objective.evalCount());
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


  // while objective.evalCount() < maxfun & itercount_mx < maxiter
  for (;;) { // Main algorithm
	 DOTRACE("Main algorithm");

	 {DOTRACE("check feval and iter counts");
	 if ((objective.evalCount() >= maxfun) || (itercount >= maxiter))
		break; // out of main loop
	 }

	 {DOTRACE("check if done");
	 if (withinTolf(funcVals_dr.asMtx(), tolf) &&
		  withinTolx(theSimplex_dr.asMtx(), tolx))
		break; // out of main loop
	 }

	 how = "";

	 {DOTRACE("compute reflection point");

	 // 
	 // Compute the reflection point
	 // 

	 // xbar = average of the numModelParams (NOT numModelParams+1) best points
	 // xbar = sum(theSimplex(:,one2n), 2)/numModelParams;

	 {DOTRACE("compute xbar");

	 Mtx xbar_new(numModelParams,1);
	 const Mtx simplex_1_n(theSimplex_dr.asMtx().columns(0,numModelParams));

	 const double numparams_inv = 1.0/numModelParams;
	 MtxIter xbar_itr = xbar_new.columnIter(0);

	 for (int r = 0; r < numModelParams; ++r, ++xbar_itr)
		*xbar_itr = simplex_1_n.row(r).sum() * numparams_inv;


	 xbar_dr.assignMtx(xbar_new);
	 }

	 {DOTRACE("compute xr");
	 // xr_mx = (1 + rho)*xbar - rho*theSimplex(:,end);
	 xr_dr.assignMtx(xbar_dr.asMtx()*(1.0+rho) - 
						  theSimplex_dr.asMtx().columns(numModelParams,1) * rho);
	 }

	 fxr = objective.evaluate(xr_dr.asMtx());
	 }

	 // 
	 // if fxr < funcVals(:,1)
	 if (fxr < funcVals_dr.asMtx().at(0,0))
		{
		  DOTRACE("if fxr < funcVals(:,1)");

		  // % Calculate the expansion point
		  // xe = (1 + rho*chi)*xbar - rho*chi*theSimplex(:,end);
		  xe_dr.assignMtx((xbar_dr.asMtx() * (1 + rho*chi))
								-
								(theSimplex_dr.asMtx().columns(numModelParams,1)
								 *(rho*chi)));

		  fxe = objective.evaluate(xe_dr.asMtx());

		  // if fxe < fxr
		  if (fxe < fxr) {

			 // theSimplex(:,end) = xe;
			 theSimplex_dr.ncMtx().column(numModelParams) = xe_dr.asMtx();

			 // funcVals(:,end) = fxe;
			 funcVals_dr.ncMtx().at(0,numModelParams) = fxe;

			 how = "expand";

          // else
		  } else {

			 // theSimplex(:,end) = xr; 
			 theSimplex_dr.ncMtx().column(numModelParams) = xr_dr.asMtx();

			 // funcVals(:,end) = fxr;
			 funcVals_dr.ncMtx().at(0,numModelParams) = fxr;

			 how = "reflect";

          // end
		  }

		  // else % funcVals(:,1) <= fxr
      } else {
		  DOTRACE("else funcVals(:,1) <= fxr");

		  // if fxr < funcVals(:,numModelParams)
		  if (fxr < funcVals_dr.asMtx().at(0,numModelParams-1)) {

			 // theSimplex(:,end) = xr; 
			 theSimplex_dr.ncMtx().column(numModelParams) = xr_dr.asMtx();

			 // funcVals(:,end) = fxr;
			 funcVals_dr.ncMtx().at(0,numModelParams) = fxr;

			 how = "reflect";

          // else % fxr >= funcVals(:,numModelParams) 
		  } else {

			 // % Perform contraction
			 // if fxr < funcVals(:,end)
			 if (fxr < funcVals_dr.asMtx().at(0,numModelParams)) {

				// % Perform an outside contraction
				// xc = (1 + psi*rho)*xbar -
				//            psi*rho*theSimplex(:,end);
				xc_dr.assignMtx(xbar_dr.asMtx()*(1.0 + psi*
															rho)
									 - (theSimplex_dr.asMtx().columns(numModelParams,1)
										 * (psi * rho)));


				fxc = objective.evaluate(xc_dr.asMtx());

				// if fxc <= fxr
				if (fxc <= fxr) {

				  // theSimplex(:,end) = xc; 
				  theSimplex_dr.ncMtx().column(numModelParams) = xc_dr.asMtx();

				  // funcVals(:,end) = fxc;
				  funcVals_dr.ncMtx().at(0,numModelParams) = fxc;

				  how = "contract outside";

				} else {

				  how = "shrink";

				}

			 } else {

				// % Perform an inside contraction
				// xcc = (1-psi)*xbar + psi*theSimplex(:,end);
				xcc_dr.assignMtx(xbar_dr.asMtx()*(1.0-psi)
									  + 
									  (theSimplex_dr.asMtx().columns(numModelParams,1)
										* psi));

				fxcc = objective.evaluate(xcc_dr.asMtx());

				// if fxcc < funcVals(:,end)
				if (fxcc < funcVals_dr.asMtx().at(0,numModelParams)) {

				  // theSimplex(:,end) = xcc;
				  theSimplex_dr.ncMtx().column(numModelParams) = xcc_dr.asMtx();

				  // funcVals(:,end) = fxcc;
				  funcVals_dr.ncMtx().at(0,numModelParams) = fxcc;

				  how = "contract inside";

				}
				else
				  {
					 how = "shrink";
				  }

			 }

			 if (how == "shrink")
				{
				  // zero based index
				  for (int j = 1; j < numModelParams+1; ++j)
					 {
						// theSimplex(:,j_mx)=theSimplex(:,1)+sigma*(theSimplex(:,j_mx) - theSimplex(:,1));
						theSimplex_dr.ncMtx().column(j) =
						  theSimplex_dr.asMtx().columns(0,1) +
						  (theSimplex_dr.asMtx().columns(j,1) -
							theSimplex_dr.asMtx().columns(0,1)) * sigma;

						// funcVals(:,j_mx) = feval(funfcn,theSimplex(:,j_mx),varargin{:});
						funcVals_dr.ncMtx().at(0,j) =
						  objective.evaluate(theSimplex_dr.asMtx()
													.columns(j,1)
													.makeMxArray()
													);
					 }
				}
		  }
      }

	 {DOTRACE("sort funcVals");
	 // [funcVals,j_mx] = sort(funcVals);

	 // Throw away the actual sorted result; just keep the indices
	 mxDestroyArray(mlfSort(&j_mx, funcVals_dr.asArray(), NULL));
	 }

	 {DOTRACE("reorder simplex");
	 Mtx jref(j_mx, Mtx::BORROW);
	 const int smallest = int(jref.at(0)) - 1;
	 const int largest = int(jref.at(numModelParams)) - 1;
	 const int largest2 = int(jref.at(numModelParams-1)) - 1;

	 // These swaps are smart enough to check if the column numbers
	 // are the same before doing the swap

	 funcVals_dr.ncMtx().swapColumns(0, smallest);
	 funcVals_dr.ncMtx().swapColumns(largest, numModelParams);
	 funcVals_dr.ncMtx().swapColumns(largest2, numModelParams-1);

	 theSimplex_dr.ncMtx().swapColumns(0, smallest);
	 theSimplex_dr.ncMtx().swapColumns(largest, numModelParams);
	 theSimplex_dr.ncMtx().swapColumns(largest2, numModelParams-1);
	 }


	 ++itercount;

	 if (prnt == 3)
		{
		  mexPrintf(" %5d        %5d     %12.6g         ",
						itercount, objective.evalCount(), double(funcVals_dr.asMtx().at(0,0)));
		  mexPrintf(how.c_str());
		  mexPrintf("\n");
		}
	 else if (prnt == 4)
		{
		  mexPrintf("\n%s\n", how.c_str());
		  theSimplex_dr.asMtx().print("theSimplex");
		  funcVals_dr.asMtx().print("funcVals");
		  mexPrintf("func_evals: %d\n", objective.evalCount());
		}
  }

  xx = theSimplex_dr.asMtx().columns(0,1);

  // end Main algorithm

  mlfIndexAssign(output, ".iterations", mxCreateScalarDouble(itercount));

  mlfIndexAssign(output, ".funcCount", mxCreateScalarDouble(objective.evalCount()));

  // output.algorithm = 'Nelder-Mead simplex direct search';
  mlfIndexAssign(output, ".algorithm", _mxarray63_);

  mlfAssign(fval, mlfMin(NULL, funcVals_dr.asArray(), NULL, NULL));

  if (objective.evalCount() >= maxfun)
	 {
		if (prnt > 0) {
		  mexPrintf("\nExiting: Maximum number of function evaluations "
						"has been exceeded\n");
		  mexPrintf("         - increase MaxFunEvals option.\n");
		  mexPrintf("         Current function value: %f \n\n",
						mxGetScalar(*fval));
		}

		mlfAssign(exitflag, mxCreateScalarDouble(0));
	 }
  else if (itercount >= maxiter)
	 {
		if (prnt > 0) {
		  mexPrintf("\nExiting: Maximum number of iterations "
						"has been exceeded\n");
		  mexPrintf("         - increase MaxIter option.\n");
		  mexPrintf("         Current function value: %f \n\n",
						mxGetScalar(*fval));
		}

		mlfAssign(exitflag, mxCreateScalarDouble(0));
	 }
  else
	 {
		if (prnt > 1)
		  {
			 const char* format = 
				"\nOptimization terminated successfully:\n"
				" the current x satisfies the termination criteria using "
				"OPTIONS.TolX of %e \n"
				" and F(X) satisfies the convergence criteria using "
				"OPTIONS.TolFun of %e \n";

			 mexPrintf(format, tolx, tolf);
		  }

		mlfAssign(exitflag, mxCreateScalarDouble(1.0));
  }

  mclValidateOutput(*fval, 2, nargout_, "fval", "doSimplex");
  mclValidateOutput(*exitflag, 3, nargout_, "exitflag", "doSimplex");
  mclValidateOutput(*output, 4, nargout_, "output", "doSimplex");

  mxDestroyArray(j_mx);

  return xx.makeMxArray();
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

  // numModelParams = prod(size(x));
  const int numModelParams = mxGetM(x_in) * mxGetN(x_in);

  // Convert to inline function as needed.
  // XXX Since this requires "object-oriented" programming, we can't keep this
  // and still use the MATLAB compiler
  // %funfcn = fcnchk(funfcn,length(varargin));

  FuncEvaluator objective(funfcn_mx, varargin);

  mxArray* result = doSimplexImpl(fval, exitflag, output, nargout_,
											 objective,
											 Mtx(x_in),
											 extractPrinttype(printtype_mx),
											 mxGetScalar(tolx_mx),
											 mxGetScalar(tolf_mx),
											 numModelParams,
											 extractMaxIters(maxfun_mx, numModelParams),
											 extractMaxIters(maxiter_mx, numModelParams)
											 );

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  return result;
}

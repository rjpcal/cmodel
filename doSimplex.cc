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

  mxArray* evaluate_mx(mxArray* x_mx)
  {
	 DOTRACE("evaluate_mx");
	 return mlfFeval(mclValueVarargout(),
						  itsFunfcn,
						  x_mx,
						  getref(itsVarargin_ref),
						  NULL);
  }

  double evaluate(mxArray* x_mx)
  {
	 DOTRACE("evaluate");
	 mxArray* mx =  mlfFeval(mclValueVarargout(),
									 itsFunfcn,
									 x_mx,
									 getref(itsVarargin_ref),
									 NULL);
	 double result = mxGetScalar(mx);
	 mxDestroyArray(mx);
	 return result;
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
										 mxArray * funfcn_mx,
										 mxArray * x_in,
										 const int prnt,
										 const double tolx,
										 const double tolf,
										 const int numModelParams,
										 const int maxfun,
										 const int maxiter,
										 mxArray * debugFlags_mx,
										 mxArray * varargin)
{

// HOME

  double fxcc = 0.0;
  mxArray* xcc_mx = mclGetUninitializedArray();
  mxArray* fxc_mx = mclGetUninitializedArray();
  mxArray* xc_mx = mclGetUninitializedArray();
  mxArray* fxe_mx = mclGetUninitializedArray();
  mxArray* xe_mx = mclGetUninitializedArray();
  mxArray* fxr_mx = mclGetUninitializedArray();
  mxArray* xr_mx = mclGetUninitializedArray();
  mxArray* formatsave_mx = mclGetUninitializedArray();
  mxArray* how_mx = mclGetUninitializedArray();
  mxArray* f_mx = mclGetUninitializedArray();
  mxArray* y_mx = mclGetUninitializedArray();
  mxArray* j_mx = mclGetUninitializedArray();
  mxArray* zero_term_delta_mx = mclGetUninitializedArray();
  mxArray* usual_delta_mx = mclGetUninitializedArray();
  mxArray* funcVals_mx = mclGetUninitializedArray();
//    mxArray* theSimplex_mx = mclGetUninitializedArray();
  mxArray* two2np1_mx = mclGetUninitializedArray();
  mxArray* header_mx = mclGetUninitializedArray();
  mxArray* ans_mx = mclGetUninitializedArray();

  DualRepMtx xbar_dr;
  DualRepMtx theSimplex_dr;

  mclCopyArray(&funfcn_mx);
  Mtx xx(x_in);
  mclCopyArray(&debugFlags_mx);
  mclCopyArray(&varargin);


	 FuncEvaluator fevaluator(funfcn_mx, varargin);

    mxArray* const numModelParams_mx =
		mclInitialize(mxCreateScalarDouble(numModelParams));

    // 
    // header_mx = ' Iteration   Func-count     min f(x)         Procedure';
    mlfAssign(&header_mx, _mxarray28_);

    // 
    // % Convert to inline function as needed.
    // % XXX Since this requires "object-oriented" programming, we can't keep this
    // % and still use the MATLAB compiler
    // %funfcn = fcnchk(funfcn,length(varargin));
    // 
    // % Initialize parameters
    // rho_mx = 1; chi_mx = 2; psi_mx = 0.5; sigma_mx = 0.5;
    mxArray* rho_mx = mclInitialize(_mxarray11_);
    mxArray* chi_mx = mclInitialize(_mxarray24_);
    mxArray* psi_mx = mclInitialize(_mxarray30_);
    mxArray* sigma_mx = mclInitialize(_mxarray30_);

    // two2np1_mx = 2:numModelParams+1;
    mlfAssign(
      &two2np1_mx,
      mlfColon(
        _mxarray24_,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // 
    // % Set up a simplex near the initial guess.
    // initialParams = x(:); % Force initialParams to be a column vector
	 DualRepMtx initialParams_dr;
	 {
		Mtx xx_copy(xx);
		xx_copy.reshape(xx.nelems(),1);
		initialParams_dr.assignMtx(xx_copy);
	 }

    // theSimplex_mx = zeros(numModelParams,numModelParams+1);
//      mlfAssign(
//        &theSimplex_mx,
//        mlfZeros(
//          numModelParams_mx,
//  		  mxCreateScalarDouble(numModelParams+1),
//          NULL));
    theSimplex_dr.assignArray(
      mlfZeros(
        numModelParams_mx,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // funcVals_mx = zeros(1,numModelParams+1);
    mlfAssign(
      &funcVals_mx,
      mlfZeros(
        _mxarray11_,
		  mxCreateScalarDouble(numModelParams+1),
        NULL));

    // theSimplex_mx(:,1) = initialParams;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
	 theSimplex_dr.ncMtx().column(0) = initialParams_dr.asMtx().columns(0,1);
//      mclArrayAssign2(
//        &theSimplex_mx,
//        initialParams_dr.asArray(),
//        mlfCreateColonIndex(),
//        _mxarray11_);

    // funcVals_mx(:,1) = feval(funfcn,initialParams,varargin{:}); 
    mclArrayAssign2(
      &funcVals_mx,
      mlfFeval(
        mclValueVarargout(),
        funfcn_mx,
        initialParams_dr.asArray(),
		  mlfIndexRef(mclVsa(varargin, "varargin"), "{?}", mlfCreateColonIndex()),
        NULL),
      mlfCreateColonIndex(),
      _mxarray11_);

    // 
    // % Following improvement suggested by L.Pfeffer at Stanford
    // usual_delta_mx = 0.05;             % 5 percent deltas for non-zero terms
    mlfAssign(&usual_delta_mx, _mxarray31_);

    // zero_term_delta_mx = 0.00025;      % Even smaller delta for zero elements of x
    mlfAssign(&zero_term_delta_mx, _mxarray32_);

    // for j_mx = 1:numModelParams
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(numModelParams_mx);
        if (v_ > e_) {
            mlfAssign(&j_mx, _mxarray33_);
        } else {

            for (int j_one_based = 1; j_one_based <= numModelParams;
					  ++j_one_based) {
				  // y = initialParams;
				  mlfAssign(&y_mx, initialParams_dr.asArray());
				  // if y(j_mx) ~= 0
				  if (mclNeBool(
									 mclVe(mclIntArrayRef1(y_mx, v_)),
									 _mxarray18_)) {
					 //   y(j_mx) = (1 + usual_delta_mx)*y(j_mx);
					 mclIntArrayAssign1(
											  &y_mx,
											  mclMtimes(
															mclPlus(_mxarray11_, mclVv(usual_delta_mx, "usual_delta_mx")),
															mclVe(mclIntArrayRef1(y_mx, v_))),
											  v_);
					 // else 
				  } else {
					 //   y(j_mx) = zero_term_delta_mx;
					 mclIntArrayAssign1(
											  &y_mx, mclVsv(zero_term_delta_mx, "zero_term_delta_mx"), v_);
					 // end  
				  }
				  // theSimplex_mx(:,j_mx+1) = y;
				  // end     
//  				  mclArrayAssign2(
//  										&theSimplex_mx,
//  										y_mx,
//  										mlfCreateColonIndex(),
//  										mlfScalar(v_ + 1));
				  theSimplex_dr.ncMtx().column(j_one_based-1+1) = Mtx(y_mx);

				  // f = feval(funfcn,y,varargin{:});
				  mlfAssign(
								&f_mx,
								mlfFeval(
											mclValueVarargout(),
											mclVa(funfcn_mx, "funfcn"),
											mclVv(y_mx, "y_mx"),
											mclVe(
													mlfIndexRef(
																	mclVsa(varargin, "varargin"),
																	"{?}",
																	mlfCreateColonIndex())),
											NULL));
				  // funcVals_mx(1,j_mx+1) = f;
				  mclIntArrayAssign2(&funcVals_mx, mclVsv(f_mx, "f_mx"), 1, v_ + 1);
				  if (v_ == e_) {
					 break;
				  }
				  ++v_;
            }
            mlfAssign(&j_mx, mlfScalar(v_));
        }
    }

    // 
    // % sort so theSimplex_mx(1,:) has the lowest function value 
    // [funcVals_mx,j_mx] = sort(funcVals_mx);
    mlfAssign(&funcVals_mx, mlfSort(&j_mx, mclVv(funcVals_mx, "funcVals_mx"), NULL));

    // theSimplex_mx = theSimplex_mx(:,j_mx);
//      mlfAssign(
//        &theSimplex_mx,
//        mclArrayRef2(
//          mclVsv(theSimplex_mx, "theSimplex_mx"),
//          mlfCreateColonIndex(),
//          mclVsv(j_mx, "j_mx")));
    theSimplex_dr.assignArray(
      mclArrayRef2(
//          mclVsv(theSimplex_mx, "theSimplex_mx"),
						 theSimplex_dr.asArray(),
						 mlfCreateColonIndex(),
						 mclVsv(j_mx, "j_mx")));

    // 
    // how_mx = 'initial';
    mlfAssign(&how_mx, _mxarray34_);

	 int itercount = 1;

	 int func_evals = numModelParams+1;

    if (prnt == 3) {

        // disp(' ')
        mlfDisp(_mxarray36_);

        // disp(header_mx)
        mlfDisp(mclVv(header_mx, "header_mx"));

        // disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount_mx, func_evals, funcVals_mx(1)), how]) 
        mlfDisp(
          mlfHorzcat(
            mclVe(
              mlfSprintf(
                NULL,
                _mxarray38_,
                mxCreateScalarDouble(itercount),
                mxCreateScalarDouble(func_evals),
                mclVe(mclIntArrayRef1(mclVsv(funcVals_mx, "funcVals_mx"), 1)),
                NULL)),
            mclVv(how_mx, "how_mx"),
            NULL));

    } else if (prnt == 4) {

        // clc
        mlfClc();

        // formatsave_mx = get(0,{'format','formatspacing'});
        mlfAssign(&formatsave_mx, mlfNGet(1, _mxarray18_, _mxarray40_, NULL));

        // format compact
        mlfFormat(_mxarray46_, NULL);

        // format short e
        mlfFormat(_mxarray48_, _mxarray50_);

        // disp(' ')
        mlfDisp(_mxarray36_);

        // disp(how_mx)
        mlfDisp(mclVv(how_mx, "how_mx"));

        // theSimplex_mx
        mclPrintArray(theSimplex_dr.asArray(), "theSimplex");

        // funcVals_mx
        mclPrintArray(mclVsv(funcVals_mx, "funcVals_mx"), "funcVals_mx");

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


    // while func_evals < maxfun & itercount_mx < maxiter
    for (;;) { // Main algorithm
		DOTRACE("Main algorithm");

		{DOTRACE("check feval and iter counts");
		if ((func_evals >= maxfun) || (itercount >= maxiter))
		  break; // out of main loop
		}

		{DOTRACE("check if done");
		if (withinTolf(funcVals_mx, tolf) && withinTolx(theSimplex_dr.asMtx(), tolx))
		  break; // out of main loop
		}

      // how_mx = '';
      mlfAssign(&how_mx, _mxarray52_);

		{DOTRACE("compute reflection point");

      // 
      // Compute the reflection point
      // 

      // xbar = average of the numModelParams (NOT numModelParams+1) best points
      // xbar = sum(theSimplex_mx(:,one2n), 2)/numModelParams;

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
      // xr_mx = (1 + rho_mx)*xbar - rho_mx*theSimplex_mx(:,end);
		const Mtx xbar(xbar_dr.asMtx());
//  		const Mtx theSimplex(theSimplex_mx, Mtx::BORROW);

		Mtx xr = xbar*(1.0+mxGetScalar(rho_mx)) - 
		  theSimplex_dr.asMtx().columns(numModelParams,1) * mxGetScalar(rho_mx);

		mlfAssign(&xr_mx, xr.makeMxArray());
		}

		mlfAssign(&fxr_mx, fevaluator.evaluate_mx(xr_mx));

		++func_evals;
		}

      // 
      // if fxr_mx < funcVals_mx(:,1)
      if (mclLtBool(
            mclVv(fxr_mx, "fxr_mx"),
            mclVe(
              mclArrayRef2(
                mclVsv(funcVals_mx, "funcVals_mx"),
                mlfCreateColonIndex(),
                _mxarray11_)))) {
			 DOTRACE("if fxr_mx < funcVals_mx(:,1)");

          // % Calculate the expansion point
          // xe_mx = (1 + rho_mx*chi_mx)*xbar - rho_mx*chi_mx*theSimplex_mx(:,end);
          mlfAssign(
            &xe_mx,
            mclMinus(
              mclMtimes(
                mclPlus(
                  _mxarray11_,
                  mclMtimes(mclVv(rho_mx, "rho_mx"), mclVv(chi_mx, "chi_mx"))),
                xbar_dr.asArray()),
              mclMtimes(
                mclMtimes(mclVv(rho_mx, "rho_mx"), mclVv(chi_mx, "chi_mx")),
                mclVe(
//                    mclArrayRef2(
//                      theSimplex_mx,
//                      mlfCreateColonIndex(),
//                      mlfEnd(
//                        mclVv(theSimplex_mx, "theSimplex_mx"),
//                        _mxarray24_,
//                        _mxarray24_))
							 theSimplex_dr.asMtx()
							 .columns(numModelParams,1).makeMxArray()
							 )
					 )));


          mlfAssign(&fxe_mx, fevaluator.evaluate_mx(xe_mx));

			 ++func_evals;

          // if fxe_mx < fxr_mx
          if (mclLtBool(mclVv(fxe_mx, "fxe_mx"), mclVv(fxr_mx, "fxr_mx"))) {

              // theSimplex_mx(:,end) = xe_mx;
//                mclArrayAssign2(
//                  &theSimplex_mx,
//                  mclVsv(xe_mx, "xe_mx"),
//                  mlfCreateColonIndex(),
//                  mlfEnd(
//                    mclVv(theSimplex_mx, "theSimplex_mx"), _mxarray24_, _mxarray24_));
				theSimplex_dr.ncMtx().column(numModelParams) = Mtx(xe_mx);

              // funcVals_mx(:,end) = fxe_mx;
              mclArrayAssign2(
                &funcVals_mx,
                mclVsv(fxe_mx, "fxe_mx"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals_mx, "funcVals_mx"), _mxarray24_, _mxarray24_));

              // how_mx = 'expand';
              mlfAssign(&how_mx, _mxarray53_);

          // else
          } else {

              // theSimplex_mx(:,end) = xr_mx; 
//                mclArrayAssign2(
//                  &theSimplex_mx,
//                  mclVsv(xr_mx, "xr_mx"),
//                  mlfCreateColonIndex(),
//                  mlfEnd(
//                    mclVv(theSimplex_mx, "theSimplex_mx"), _mxarray24_, _mxarray24_));
				theSimplex_dr.ncMtx().column(numModelParams) = Mtx(xr_mx);

              // funcVals_mx(:,end) = fxr_mx;
              mclArrayAssign2(
                &funcVals_mx,
                mclVsv(fxr_mx, "fxr_mx"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals_mx, "funcVals_mx"), _mxarray24_, _mxarray24_));

              // how_mx = 'reflect';
              mlfAssign(&how_mx, _mxarray55_);

          // end
          }

      // else % funcVals_mx(:,1) <= fxr_mx
      } else {
			 DOTRACE("else funcVals_mx(:,1) <= fxr_mx");

          // if fxr_mx < funcVals_mx(:,numModelParams)
          if (mclLtBool(
                mclVv(fxr_mx, "fxr_mx"),
                mclVe(
                  mclArrayRef2(
                    mclVsv(funcVals_mx, "funcVals_mx"),
                    mlfCreateColonIndex(),
                    numModelParams_mx)))) {

              // theSimplex_mx(:,end) = xr_mx; 
//                mclArrayAssign2(
//                  &theSimplex_mx,
//                  mclVsv(xr_mx, "xr_mx"),
//                  mlfCreateColonIndex(),
//                  mlfEnd(
//                    mclVv(theSimplex_mx, "theSimplex_mx"), _mxarray24_, _mxarray24_));
				theSimplex_dr.ncMtx().column(numModelParams) = Mtx(xr_mx);

              // funcVals_mx(:,end) = fxr_mx;
              mclArrayAssign2(
                &funcVals_mx,
                mclVsv(fxr_mx, "fxr_mx"),
                mlfCreateColonIndex(),
                mlfEnd(
                  mclVv(funcVals_mx, "funcVals_mx"), _mxarray24_, _mxarray24_));

              // how_mx = 'reflect';
              mlfAssign(&how_mx, _mxarray55_);

          // else % fxr_mx >= funcVals_mx(:,numModelParams) 
          } else {

              // % Perform contraction
              // if fxr_mx < funcVals_mx(:,end)
              if (mclLtBool(
                    mclVv(fxr_mx, "fxr_mx"),
                    mclVe(
                      mclArrayRef2(
                        mclVsv(funcVals_mx, "funcVals_mx"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals_mx, "funcVals_mx"),
                          _mxarray24_,
                          _mxarray24_))))) {

                  // % Perform an outside contraction
                  // xc_mx = (1 + psi_mx*rho_mx)*xbar -
					   //            psi_mx*rho_mx*theSimplex_mx(:,end);
                  mlfAssign(
                    &xc_mx,
                    mclMinus(
                      mclMtimes(
                        mclPlus(
                          _mxarray11_,
                          mclMtimes(mclVv(psi_mx, "psi_mx"), mclVv(rho_mx, "rho_mx"))),
                        xbar_dr.asArray()),
                      mclMtimes(
                        mclMtimes(mclVv(psi_mx, "psi_mx"), mclVv(rho_mx, "rho_mx")),
                        mclVe(
//                            mclArrayRef2(
//                              theSimplex_mx,
//                              mlfCreateColonIndex(),
//                              mlfEnd(
//                                mclVv(theSimplex_mx, "theSimplex_mx"),
//                                _mxarray24_,
//                                _mxarray24_))
										theSimplex_dr.asMtx().columns(numModelParams,1)
										.makeMxArray()
										)
								)));


                  mlfAssign(&fxc_mx, fevaluator.evaluate_mx(xc_mx));

						++func_evals;

                  // 
                  // if fxc_mx <= fxr_mx
                  if (mclLeBool(mclVv(fxc_mx, "fxc_mx"), mclVv(fxr_mx, "fxr_mx"))) {

                      // theSimplex_mx(:,end) = xc_mx; 
//                        mclArrayAssign2(
//                          &theSimplex_mx,
//                          mclVsv(xc_mx, "xc_mx"),
//                          mlfCreateColonIndex(),
//                          mlfEnd(
//                            mclVv(theSimplex_mx, "theSimplex_mx"),
//                            _mxarray24_,
//                            _mxarray24_));
						  theSimplex_dr.ncMtx().column(numModelParams) = Mtx(xc_mx);

                      // funcVals_mx(:,end) = fxc_mx;
                      mclArrayAssign2(
                        &funcVals_mx,
                        mclVsv(fxc_mx, "fxc_mx"),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals_mx, "funcVals_mx"),
                          _mxarray24_,
                          _mxarray24_));

                      // how_mx = 'contract outside';
                      mlfAssign(&how_mx, _mxarray57_);

                  // else
                  } else {

                      // % perform a shrink
                      // how_mx = 'shrink'; 
                      mlfAssign(&how_mx, _mxarray59_);

                  // end
                  }

              // else
              } else {

                  // % Perform an inside contraction
                  // xcc_mx = (1-psi_mx)*xbar + psi_mx*theSimplex_mx(:,end);
                  mlfAssign(
                    &xcc_mx,
                    mclPlus(
                      mclMtimes(
                        mclMinus(_mxarray11_, mclVv(psi_mx, "psi_mx")),
                        xbar_dr.asArray()),
                      mclMtimes(
                        mclVv(psi_mx, "psi_mx"),
                        mclVe(
//                            mclArrayRef2(
//                              theSimplex_mx,
//                              mlfCreateColonIndex(),
//                              mlfEnd(
//                                mclVv(theSimplex_mx, "theSimplex_mx"),
//                                _mxarray24_,
//                                _mxarray24_))
										theSimplex_dr.asMtx().columns(numModelParams,1)
										.makeMxArray()
										)
								)));

						fxcc = fevaluator.evaluate(xcc_mx);

						++func_evals;

                  // 
                  // if fxcc < funcVals_mx(:,end)
                  if (mclLtBool(
                        mxCreateScalarDouble(fxcc),
                        mclVe(
                          mclArrayRef2(
                            mclVsv(funcVals_mx, "funcVals_mx"),
                            mlfCreateColonIndex(),
                            mlfEnd(
                              mclVv(funcVals_mx, "funcVals_mx"),
                              _mxarray24_,
                              _mxarray24_))))) {

                      // theSimplex_mx(:,end) = xcc_mx;
//                        mclArrayAssign2(
//                          &theSimplex_mx,
//                          mclVsv(xcc_mx, "xcc_mx"),
//                          mlfCreateColonIndex(),
//                          mlfEnd(
//                            mclVv(theSimplex_mx, "theSimplex_mx"),
//                            _mxarray24_,
//                            _mxarray24_));
						  theSimplex_dr.ncMtx().column(numModelParams) = Mtx(xcc_mx);

                      // funcVals_mx(:,end) = fxcc;
                      mclArrayAssign2(
                        &funcVals_mx,
                        mxCreateScalarDouble(fxcc),
                        mlfCreateColonIndex(),
                        mlfEnd(
                          mclVv(funcVals_mx, "funcVals_mx"),
                          _mxarray24_,
                          _mxarray24_));

                      // how_mx = 'contract inside';
                      mlfAssign(&how_mx, _mxarray61_);

                  // else
                  } else {

                      // % perform a shrink
                      // how_mx = 'shrink';
                      mlfAssign(&how_mx, _mxarray59_);

                  // end
                  }

              // end
              }

              // if strcmp(how_mx,'shrink')
              if (mlfTobool(
                    mclVe(mlfStrcmp(mclVv(how_mx, "how_mx"), _mxarray59_)))) {
                  mclForLoopIterator viter__;

                  // for j_mx=two2np1_mx
                  for (mclForStart(
                         &viter__, mclVv(two2np1_mx, "two2np1_mx"), NULL, NULL);
                       mclForNext(&viter__, &j_mx);
                       ) {

						  int j_zero_based = int(mxGetScalar(j_mx))-1;
                      // theSimplex_mx(:,j_mx)=theSimplex_mx(:,1)+sigma_mx*(theSimplex_mx(:,j_mx) - theSimplex_mx(:,1));
						  theSimplex_dr.ncMtx().column(j_zero_based) =
							 theSimplex_dr.asMtx().columns(0,1) +
							 (theSimplex_dr.asMtx().columns(j_zero_based,1) -
							  theSimplex_dr.asMtx().columns(0,1)) * mxGetScalar(sigma_mx);
//                        mclArrayAssign2(
//                          &theSimplex_mx,
//                          mclPlus(
//                            mclVe(
//                              mclArrayRef2(
//                                theSimplex_mx,
//                                mlfCreateColonIndex(),
//                                _mxarray11_)),
//                            mclMtimes(
//                              mclVv(sigma_mx, "sigma_mx"),
//                              mclMinus(
//                                mclVe(
//                                  mclArrayRef2(
//                                    theSimplex_mx,
//                                    mlfCreateColonIndex(),
//                                    mclVsv(j_mx, "j_mx"))),
//                                mclVe(
//                                  mclArrayRef2(
//                                    theSimplex_mx,
//                                    mlfCreateColonIndex(),
//                                    _mxarray11_))))),
//                          mlfCreateColonIndex(),
//                          mclVsv(j_mx, "j_mx"));

                      // funcVals_mx(:,j_mx) = feval(funfcn,theSimplex_mx(:,j_mx),varargin{:});
                      mclArrayAssign2(
                        &funcVals_mx,
								fevaluator.evaluate_mx(
//  															  mclArrayRef2(theSimplex_mx,
//  																mlfCreateColonIndex(),
//  																mclVsv(j_mx, "j_mx"))
															  theSimplex_dr.asMtx()
															  .columns(j_zero_based,1)
															  .makeMxArray()
															  ),
                        mlfCreateColonIndex(),
                        mclVsv(j_mx, "j_mx"));

                  }
                  mclDestroyForLoopIterator(viter__);

						func_evals += numModelParams;

              }

          }

      }

		{DOTRACE("sort funcVals_mx");
      // [funcVals_mx,j_mx] = sort(funcVals_mx);

		// Throw away the actual sorted result; just keep the indices
		mxDestroyArray(mlfSort(&j_mx, funcVals_mx, NULL));
		}

		{DOTRACE("reorder simplex");
		Mtx jref(j_mx, Mtx::BORROW);
		const int smallest = int(jref.at(0)) - 1;
		const int largest = int(jref.at(numModelParams)) - 1;
		const int largest2 = int(jref.at(numModelParams-1)) - 1;

//  		Mtx simref(theSimplex_mx, Mtx::REFER);
		Mtx fvref(funcVals_mx, Mtx::REFER);

		// These swaps are smart enough to check if the column numbers
		// are the same before doing the swap

		fvref.swapColumns(0, smallest);
		fvref.swapColumns(largest, numModelParams);
		fvref.swapColumns(largest2, numModelParams-1);

		theSimplex_dr.ncMtx().swapColumns(0, smallest);
		theSimplex_dr.ncMtx().swapColumns(largest, numModelParams);
		theSimplex_dr.ncMtx().swapColumns(largest2, numModelParams-1);
		}


		++itercount;

      if (prnt == 3) {

          // disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount_mx, func_evals, funcVals_mx(1)), how_mx]) 
          mlfDisp(
            mlfHorzcat(
              mclVe(
                mlfSprintf(
                  NULL,
                  _mxarray38_,
                  mxCreateScalarDouble(itercount),
                  mxCreateScalarDouble(func_evals),
                  mclVe(mclIntArrayRef1(mclVsv(funcVals_mx, "funcVals_mx"), 1)),
                  NULL)),
              mclVv(how_mx, "how_mx"),
              NULL));

      } else if (prnt == 4) {

          // disp(' ')
          mlfDisp(_mxarray36_);

          // disp(how_mx)
          mlfDisp(mclVv(how_mx, "how_mx"));

          // theSimplex_mx
          mclPrintArray(theSimplex_dr.asArray(), "theSimplex_mx");

          // funcVals_mx
          mclPrintArray(mclVsv(funcVals_mx, "funcVals_mx"), "funcVals_mx");

          // func_evals
          mclPrintArray(mxCreateScalarDouble(func_evals), "func_evals");

      // end  
      }

    // end   % while
    }

    // 
    // 
    // x(:) = theSimplex_mx(:,1);
	 {
//  		Mtx theSimplex_ref(theSimplex_mx, Mtx::BORROW);
		xx = theSimplex_dr.asMtx().columns(0,1);
	 }

    if (prnt == 4) {

        // % reset format
        // set(0,{'format','formatspacing'},formatsave_mx);
        mclAssignAns(
          &ans_mx,
          mlfNSet(
            0,
            _mxarray18_,
            _mxarray40_,
            mclVv(formatsave_mx, "formatsave_mx"),
            NULL));

    } // end Main algorithm

    // output.iterations = itercount_mx;
    mlfIndexAssign(output, ".iterations", mxCreateScalarDouble(itercount));

    // output.funcCount = func_evals;
    mlfIndexAssign(output, ".funcCount", mxCreateScalarDouble(func_evals));

    // output.algorithm = 'Nelder-Mead simplex direct search';
    mlfIndexAssign(output, ".algorithm", _mxarray63_);

    // 
    // fval = min(funcVals_mx);
    mlfAssign(fval, mlfMin(NULL, mclVv(funcVals_mx, "funcVals_mx"), NULL, NULL));

    // if func_evals >= maxfun 
    if (func_evals >= maxfun) {

        if (prnt > 0) {
			 mexPrintf("\nExiting: Maximum number of function evaluations "
						  "has been exceeded\n");
			 mexPrintf("         - increase MaxFunEvals option.\n");
			 mexPrintf("         Current function value: %f \n\n",
						  mxGetScalar(*fval));
        }

        // exitflag = 0;
        mlfAssign(exitflag, _mxarray18_);

    }
	 else if (itercount >= maxiter) {

        if (prnt > 0) {
			 mexPrintf("\nExiting: Maximum number of iterations "
						  "has been exceeded\n");
			 mexPrintf("         - increase MaxIter option.\n");
			 mexPrintf("         Current function value: %f \n\n",
						  mxGetScalar(*fval));
        }

        // exitflag = 0; 
        mlfAssign(exitflag, _mxarray18_);

    // else
    } else {

        if (prnt > 1) {

			 const char* format = 
				"\nOptimization terminated successfully:\n"
				" the current x satisfies the termination criteria using "
				"OPTIONS.TolX of %e \n"
				" and F(X) satisfies the convergence criteria using "
				"OPTIONS.TolFun of %e \n";

			 mexPrintf(format, tolx, tolf);
        }


        // exitflag = 1;
        mlfAssign(exitflag, _mxarray11_);


    // end
    }

    mclValidateOutput(*fval, 2, nargout_, "fval", "doSimplex");
    mclValidateOutput(*exitflag, 3, nargout_, "exitflag", "doSimplex");
    mclValidateOutput(*output, 4, nargout_, "output", "doSimplex");

    mxDestroyArray(numModelParams_mx);
    mxDestroyArray(ans_mx);
    mxDestroyArray(header_mx);
    mxDestroyArray(rho_mx);
    mxDestroyArray(chi_mx);
    mxDestroyArray(psi_mx);
    mxDestroyArray(sigma_mx);
    mxDestroyArray(two2np1_mx);
//      mxDestroyArray(theSimplex_mx);
    mxDestroyArray(funcVals_mx);
    mxDestroyArray(usual_delta_mx);
    mxDestroyArray(zero_term_delta_mx);
    mxDestroyArray(j_mx);
    mxDestroyArray(y_mx);
    mxDestroyArray(f_mx);
    mxDestroyArray(how_mx);
    mxDestroyArray(formatsave_mx);
    mxDestroyArray(xr_mx);
    mxDestroyArray(fxr_mx);
    mxDestroyArray(xe_mx);
    mxDestroyArray(fxe_mx);
    mxDestroyArray(xc_mx);
    mxDestroyArray(fxc_mx);
    mxDestroyArray(xcc_mx);
    mxDestroyArray(varargin);
    mxDestroyArray(debugFlags_mx);
    mxDestroyArray(funfcn_mx);

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

	 mxArray* result = doSimplexImpl(fval, exitflag, output, nargout_,
												funfcn_mx,
												x_in,
												extractPrinttype(printtype_mx),
												mxGetScalar(tolx_mx),
												mxGetScalar(tolf_mx),
												numModelParams,
												extractMaxIters(maxfun_mx, numModelParams),
												extractMaxIters(maxiter_mx, numModelParams),
												debugFlags_mx,
												varargin);

    mclSetCurrentLocalFunctionTable(save_local_function_table_);

	 return result;
}

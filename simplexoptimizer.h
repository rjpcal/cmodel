///////////////////////////////////////////////////////////////////////
//
// simplexoptimizer.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 14:52:30 2001
// written: Thu Feb 14 11:57:09 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef SIMPLEXOPTIMIZER_H_DEFINED
#define SIMPLEXOPTIMIZER_H_DEFINED

#include <iostream>

#include "mtx.h"
#include "multivarfunction.h"

class fstring;

class SimplexOptimizer
{
public:
  enum PrintType { NONE, NOTIFY, FINAL, ITER, SIMPLEX };

  SimplexOptimizer(MultivarFunction& objective,
                   const Mtx& x_in,
                   const fstring& printtype,
                   const int nparams,
                   const int maxfun = 0, // default is 200*nparams
                   const int maxiter = 0, // default is 200*nparams
                   const double tolx = 1e-4,
                   const double tolf = 1e-4);

  virtual ~SimplexOptimizer();

  virtual int optimize();

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

private:
  MultivarFunction& itsObjective;

  const Mtx itsInitialParams;
  const int itsPrnt;
  const int itsNparams;
  const int itsMaxFevals;
  const int itsMaxIters;
  const double itsTolx;
  const double itsTolf;

  Mtx itsSimplex;
  Mtx itsFvals;

  int itsIterCount;

  enum IterType
    {
      INITIAL,
      EXPAND,
      REFLECT,
      CONTRACT_OUTSIDE,
      CONTRACT_INSIDE,
      SHRINK
    };

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

  void putInSimplex(const FuncPoint& p, int pointNumber);

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
    // sort so itsSimplex.column(0) has the lowest function value
    Mtx index = itsFvals.row(0).getSortOrder();

    itsFvals.row(0).reorder(index);
    itsSimplex.reorderColumns(index);
  }
#endif

  void minimalSort();

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

    if (itsPrnt != NONE)
      {
        std::cerr << "\nExiting: Maximum number of function evaluations "
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

    if (itsPrnt != NONE)
      {
        std::cerr << "\nExiting: Maximum number of iterations "
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
      std::cerr << "\n Iteration   Func-count     min f(x)         Procedure\n";
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

  void printIter();

  void printSimplex();
};

static const char vcid_simplexoptimizer_h[] = "$Header$";
#endif // !SIMPLEXOPTIMIZER_H_DEFINED

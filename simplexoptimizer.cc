///////////////////////////////////////////////////////////////////////
//
// simplexoptimizer.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 14:52:57 2001
// written: Wed Jul 31 15:11:11 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef SIMPLEXOPTIMIZER_CC_DEFINED
#define SIMPLEXOPTIMIZER_CC_DEFINED

#include "cmodel/simplexoptimizer.h"

#include "util/strings.h"

#include <iostream>
#include <iomanip>

#include "util/trace.h"

namespace
{
  SimplexOptimizer::PrintType extractPrinttype(const fstring& printtype)
  {
    if      (printtype == "notify")   return SimplexOptimizer::NOTIFY;
    else if (printtype == "none")     return SimplexOptimizer::NONE;
    else if (printtype == "off")      return SimplexOptimizer::NONE;
    else if (printtype == "iter")     return SimplexOptimizer::ITER;
    else if (printtype == "final")    return SimplexOptimizer::FINAL;
    else if (printtype == "simplex")  return SimplexOptimizer::SIMPLEX;

    return SimplexOptimizer::NOTIFY;
  }
}

SimplexOptimizer::SimplexOptimizer(MultivarFunction& objective,
                                   const Mtx& x_in,
                                   const fstring& printtype,
                                   const int nparams,
                                   const int maxfun,
                                   const int maxiter,
                                   const double tolx,
                                   const double tolf) :
  itsObjective(objective),

  itsInitialParams(x_in.asColumn()), // Set up simplex near the initial guess
  itsPrnt(extractPrinttype(printtype)),
  itsNparams(nparams),
  itsInitialFevals(objective.evalCount()),
  itsMaxFevals(maxfun > 0 ? maxfun : 200*nparams),
  itsMaxIters(maxiter > 0 ? maxiter : 200*nparams),
  itsTolx(tolx),
  itsTolf(tolf),

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

SimplexOptimizer::~SimplexOptimizer() {}

void SimplexOptimizer::putInSimplex(const FuncPoint& p, int pointNumber)
{
  itsSimplex.column(pointNumber) = p.x.column(0);
  itsFvals.at(0, pointNumber) = p.f;
}

void SimplexOptimizer::printIter()
{
  if (itsPrnt == ITER)
    {
      std::cerr << std::setw(6) << itsIterCount
                << std::setw(13) << funcCount()
                << std::setw(17) << std::setprecision(6)
                << double(itsFvals.at(0,0))
                << "         " << iterTypeString(itsCurIter) << '\n';
    }
}

void SimplexOptimizer::printSimplex()
{
  if (itsPrnt == SIMPLEX)
    {
      std::cerr << '\n' << iterTypeString(itsCurIter) << '\n';
      itsSimplex.print("simplex");
      itsFvals.print("fvals");
      std::cerr << "funcEvals: " << funcCount() << '\n';
    }
}

void SimplexOptimizer::minimalSort()
{
DOTRACE("SimplexOptimizer::minimalSort");

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

int SimplexOptimizer::optimize()
{
DOTRACE("SimplexOptimizer::optimize");

  itsInitialFevals = itsObjective.evalCount();

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

  std::cerr << "\nOptimization terminated successfully:\n"
            << " the current x satisfies the termination criteria using "
            << "OPTIONS.TolX of " << itsTolx << '\n'
            << " and F(X) satisfies the convergence criteria using "
            << "OPTIONS.TolFun of " << itsTolf << '\n';

  return 1;
}

bool SimplexOptimizer::withinTolf() const
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

bool SimplexOptimizer::withinTolx() const
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

bool SimplexOptimizer::tooManyFevals() const
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

bool SimplexOptimizer::tooManyIters() const
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

void SimplexOptimizer::doOneIter()
{
DOTRACE("SimplexOptimizer::doOneIter");

  // compute average of the itsNparams (NOT itsNparams+1) best points
  Mtx xbar(itsSimplex(col_range_n(0,itsNparams)).meanColumn());

  // Direction along which to reflect/expand/contract
  Mtx direction(xbar - itsSimplex(col_range_n(itsNparams,1)));

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
                  putInSimplex(itsSimplex(col_range_n(0,1)) +
                               (itsSimplex(col_range_n(j,1)) -
                                itsSimplex(col_range_n(0,1))) * 0.5,
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

void SimplexOptimizer::printHeader()
{
  if (itsPrnt == ITER)
    std::cerr << "\n Iteration   Func-count     min f(x)         Procedure\n";
}

static const char vcid_simplexoptimizer_cc[] = "$Header$";
#endif // !SIMPLEXOPTIMIZER_CC_DEFINED

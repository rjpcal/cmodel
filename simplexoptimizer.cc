///////////////////////////////////////////////////////////////////////
//
// simplexoptimizer.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 14:52:57 2001
// written: Wed Apr 18 15:07:44 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef SIMPLEXOPTIMIZER_CC_DEFINED
#define SIMPLEXOPTIMIZER_CC_DEFINED

#include "simplexoptimizer.h"

#include "strings.h"

#include <iostream.h>
#include <iomanip.h>

#include "trace.h"

namespace
{
  SimplexOptimizer::PrintType extractPrinttype(const fixed_string& printtype)
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
											  const fixed_string& printtype,
											  const int nparams,
											  const int maxfun,
											  const int maxiter,
											  const double tolx,
											  const double tolf) :
  itsObjective(objective),

  itsInitialParams(x_in.asColumn()), // Set up simplex near the initial guess
  itsPrnt(extractPrinttype(printtype)),
  itsNparams(nparams),
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

void SimplexOptimizer::printIter()
{
  if (itsPrnt == ITER)
	 {
		cerr << setw(6) << itsIterCount
			  << setw(13) << funcCount()
			  << setw(17) << setprecision(6) << itsFvals.at(0,0)
			  << "         " << iterTypeString(itsCurIter) << '\n';
	 }
}

void SimplexOptimizer::printSimplex()
{
  if (itsPrnt == SIMPLEX)
	 {
		cerr << '\n' << iterTypeString(itsCurIter) << '\n';
		itsSimplex.print("simplex");
		itsFvals.print("fvals");
		cerr << "funcEvals: " << funcCount() << '\n';
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


static const char vcid_simplexoptimizer_cc[] = "$Header$";
#endif // !SIMPLEXOPTIMIZER_CC_DEFINED

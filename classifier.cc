///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Thu Mar  8 09:54:37 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

/*
  forwardProbit should be a member function, should:
  return a double*, OR
  return its result into a data member, OR
  return its result into a reference parameter?

  the pure virtual function should be:
  virtual void void computeMu(Rat& modelParams, double* muOut) = 0;

  then functions in the base class could be:
  double loglikelihoodFor(Rat& modelParams);
  double devianceFor(Rat& modelParams);


  can call from MATLAB like this:
  classifier(modelParams, 'cssm', 'll', objParams, observedIncidence, ...);
  OR
  classifier(modelParams, 'cssm', 'dev', objParams, observedIncidence, ...);
 */

#include "classifier.h"
#include "cssm_scaleWeights1.h"
#include "forwardProbit1.h"
#include "loglikelihood1.h"

#include "rutil.h"

#include "trace.h"

#include <cmath>
#include <fstream.h>
#include "libmatlbm.h"

void InitializeModule_classifier(void) {
}

void TerminateModule_classifier(void) {
}

#if 0
// A quick and dirty approximation (WHICH IS CURRENTLY GIVING
// INCORRECT RESULTS!) to the exponential function using the first
// eight terms of the Taylor series
double fastexp(double x) {
  return (1.0 +
			 x * (1.0 +
					x * (0.5 +
						  x * (0.1666666667 +
								 x * (0.0416666667 +
										x * (0.0083333333 +
											  x * (0.0013888889 +
													 x * 0.0001984127)))))));
}
#endif

double minkDist(const double* w, int nelems,
                const double* x1, int stride1,
                const double* x2, int stride2,
                double r, double r_inv)
{
  double wt_sum = 0.0;
  for (int k = 0; k < nelems; ++k)
	 {
		wt_sum +=
		  w[k] *
		  pow( abs( x1[k*stride1] - x2[k*stride2]), r);
	 }
  return pow(wt_sum, r_inv);
}

//
// This is just a functor that binds 4 of the 5 arguments to
// minkDist2, so that we don't have to keep passing 5 arguments to
// functions all the time in the time-critical inner loop.
//

class MinkDist2Binder {
public:
  MinkDist2Binder(const double* attWeights, int nelems,
						int stride1,
						const double* x2) :
	 itsAttWeights(attWeights),
	 itsNelems(nelems),
	 itsStride1(stride1),
	 itsX2(x2)
  {}

  // Specialized Minkowski distance for r==2
  double minkDist2(const double* x1) const
  {
	 double wt_sum = 0.0;
	 const double* x2 = itsX2;
	 const double* w = itsAttWeights;
	 for (int k = 0; k < itsNelems; ++k, x1 += itsStride1, ++x2, ++w)
		{
		  wt_sum +=
			 (*w) * 
			 ((*x1) - (*x2)) * ((*x1) - (*x2));
		}
	 return sqrt(wt_sum);	 
  }

private:
  const double* const itsAttWeights;
  int itsNelems;
  int itsStride1;
  const double* const itsX2;
};

void linearCombo(int nelems, const double* w, int w_stride,
                 const double* const* elems, int elems_stride, int dim,
                 double* result)
{
DOTRACE("linearCombo");
  // e.g nelems == 3   dim == 4
  //
  //               | e11  e12  e13  e14 |
  // [w1 w2 w3] *  | e21  e22  e23  e24 | = 
  //               | e31  e32  e33  e34 |
  //
  //
  // [ w1*e11+w2*e21+w3*e31  w1*e12+w2*e22+w3*e32  ... ]
  for (int d = 0; d < dim; ++d)
    {
      result[d] = 0.0;
      for (int elem = 0; elem < nelems; ++elem)
        result[d] += w[elem*w_stride] * elems[elem][d*elems_stride];
    }
}

class ModelCssm {
private:
  const Rat& objParams;
  const Rat& observedIncidence;
  const int numStoredExemplars;
  const int num1;
  const int num2;
  const int numTrainingExemplars;
  const int numAllExemplars;

  const int dimObjParams;

  double* const mu;

  typedef const double* constDblPtr;

  constDblPtr* const cat1;
  constDblPtr* const cat2;

  // Count the category training exemplars
  int countCategory(const Rat& params, int category)
  {
	 int n = 0;
    for (int i = 0; i < params.mrows(); ++i)
      {
        if (int(params.at(i,0)) == category)
          ++n;
      }
	 return n;
  }

  void resetMu()
  {
    for (int i = 0; i < numAllExemplars; ++i)
      mu[i] = 0.0;
  }

  void computeSimilarity(const double* attWeights,
								 const double* storedExemplar1,
								 const double* storedExemplar2,
								 double minkPower,
								 double minkPowerInv)
  {
	 DOTRACE("computeSimilarity");

	 for (int y = 0; y < numAllExemplars; ++y) {

		// compute similarity of ex-y to stored-1-x
		double sim1 =
		  minkDist(attWeights, dimObjParams,
					  objParams.data()+y+numAllExemplars, numAllExemplars,
					  storedExemplar1, 1,
					  minkPower, minkPowerInv);

		mu[y] += exp(-sim1);

		// compute similarity of ex-y to stored-2-x
		double sim2 =
		  minkDist(attWeights, dimObjParams,
					  objParams.data()+y+numAllExemplars, numAllExemplars,
					  storedExemplar2, 1,
					  minkPower, minkPowerInv);

		mu[y] -= exp(-sim2);
	 }
  }

  void computeSimilarity2(const double* attWeights,
								  const double* storedExemplar1,
								  const double* storedExemplar2)
  {
	 DOTRACE("computeSimilarity2");

	 MinkDist2Binder binder1(attWeights, dimObjParams,
									 numAllExemplars,
									 storedExemplar1);

	 MinkDist2Binder binder2(attWeights, dimObjParams,
									 numAllExemplars,
									 storedExemplar2);

	 // This finds the first testExemplar data point (skipping the
	 // first column of objParams, which contains category labels)
	 const double* testExemplar = objParams.data()+numAllExemplars;

	 for (int y = 0; y < numAllExemplars; ++y, ++testExemplar) {

		// compute similarity of ex-y to stored-1-x
		const double sim1 = binder1.minkDist2(testExemplar);

		// compute similarity of ex-y to stored-2-x
		const double sim2 = binder2.minkDist2(testExemplar);

		mu[y] += exp(-sim1) - exp(-sim2);
	 }
  }

public:
  ModelCssm(const Rat& objParams_,
				const Rat& observedIncidence_,
				int numStoredExemplars_) :
	 objParams(objParams_),
	 observedIncidence(observedIncidence_),
	 numStoredExemplars(numStoredExemplars_),
	 num1(countCategory(objParams_, 0)),
	 num2(countCategory(objParams_, 1)),
	 numTrainingExemplars(num1),
	 numAllExemplars(objParams_.mrows()),
	 dimObjParams(4),
	 mu(new double[numAllExemplars]),
	 cat1(new constDblPtr[num1]),
	 cat2(new constDblPtr[num2])
  {
	 if (num1 != num2) {
		mexErrMsgTxt("the two categories must have the "
						 "same number of training exemplars");
	 }

	 // Find the category 1 and category 2 training exemplars
    int c1=0,c2=0;
    for (int i = 0; i < objParams.mrows(); ++i)
      {
        if (int(objParams.at(i,0)) == 0)
          cat1[c1++] = objParams.address(i,1);
        else if (int(objParams.at(i,0)) == 1)
          cat2[c2++] = objParams.address(i,1);
      }
  }

  ~ModelCssm()
  {
	 delete [] cat2;
	 delete [] cat1;
	 delete [] mu;
  }

  double loglikelihoodFor(Rat& modelParams);
};


_mexLocalFunctionTable _local_function_table_classifier
  = { 0, (mexFunctionTableEntry *)NULL };

///////////////////////////////////////////////////////////////////////
//
// compute the minus loglikelihood for the constrained summed
// similarity model (cssm)
//
// based on the MATLAB function llcssm 
//
//     function ll = llcssm(modelParams, objParams, ...
//                          observedIncidence, numStoredExemplars)
//
///////////////////////////////////////////////////////////////////////

double ModelCssm::loglikelihoodFor(Rat& modelParams) {
DOTRACE("ModelCssm::loglikelihoodFor");

  if (modelParams.length() !=
		(2*numTrainingExemplars*numStoredExemplars + 6)) {
    mexErrMsgTxt("wrong number of model parameters");
  }

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights,
  //
  //---------------------------------------------------------------------

  double* attWeights = modelParams.data();

  {
    for (int i = 0; i < dimObjParams; ++i)
      attWeights[i] = abs(attWeights[i]);
  }

  //---------------------------------------------------------------------
  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //
  //---------------------------------------------------------------------

  double* rawWeights = modelParams.data() + 6;
  const int numRawWeights = modelParams.nelems() - 6;

  cssm_scaleWeights(rawWeights, numRawWeights,
                    numStoredExemplars, numTrainingExemplars);

  const double* scaledWeights = rawWeights;

  //---------------------------------------------------------------------
  //
  // Compute mu, a matrix of differences of summed similarities.
  //
  //---------------------------------------------------------------------

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  resetMu();

  double stored1[dimObjParams];
  double stored2[dimObjParams];

  for (int x = 0; x < numStoredExemplars; ++x) {

    linearCombo(numTrainingExemplars, scaledWeights+x,
                2*numStoredExemplars,
                cat1, numAllExemplars, dimObjParams,
                &stored1[0]);

    linearCombo(numTrainingExemplars, scaledWeights+x+numStoredExemplars,
                2*numStoredExemplars,
                cat2, numAllExemplars, dimObjParams,
                &stored2[0]);

	 {
		if (minkPower == 2.0) {
		  computeSimilarity2(attWeights, &stored1[0], &stored2[0]);
		}
		else {
		  computeSimilarity(attWeights, &stored1[0], &stored2[0],
								  minkPower, minkPowerInv);
		}
	 }
  }


  //---------------------------------------------------------------------
  //
  // Compute the predicted probabilities based on the similarity
  // differences in mu.
  //
  //---------------------------------------------------------------------

  // thresh = modelParams(5);
  const double thresh = modelParams.data()[4];

  // sigmaNoise = sqrt(n1+n2)*modelParams(6);
  const double sigmaNoise =
    sqrt(numStoredExemplars*2)*modelParams.data()[5];

  // predictedProbability = forwardProbit(mu, thresh, sigmaNoise);
  mxArray* predictedProbability_mx =
	 forwardProbit(mu, numAllExemplars, thresh, sigmaNoise);

  validateVariable(predictedProbability_mx);

  //---------------------------------------------------------------------
  //
  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.
  //
  //---------------------------------------------------------------------

  double ll = -loglikelihood(mxGetPr(predictedProbability_mx),
									  numAllExemplars,
									  observedIncidence.data(),
									  observedIncidence.data()+numAllExemplars);

  mxDestroyArray(predictedProbability_mx);

  return ll;
}



static mxArray * Mclassifier(int /* nargout_ */,
									  mxArray * modelParams_mx,
									  mxArray * objParams_mx,
									  mxArray * observedIncidence_mx,
									  mxArray * numStoredExemplars_mx)
{
DOTRACE("Mclassifier");

#ifdef LOCAL_PROF
  if (mxGetScalar(numStoredExemplars_mx) == -1) {
	 ofstream ofs("profdata.out");
	 Util::Prof::printAllProfData(ofs);
	 return mxCreateScalarDouble(-1.0);
  }

  if (mxGetScalar(numStoredExemplars_mx) == -2) {
	 Util::Prof::resetAllProfData();
	 return mxCreateScalarDouble(-2.0);
  }
#endif

  //---------------------------------------------------------------------
  //
  // Set up
  //
  //---------------------------------------------------------------------

  mexLocalFunctionTable save_local_function_table_ =
    mclSetCurrentLocalFunctionTable(&_local_function_table_classifier);

  mclCopyArray(&modelParams_mx); // THIS IS REQUIRED
  validateInput(modelParams_mx);

  validateInput(objParams_mx);
  const Rat objParams(objParams_mx);

  validateInput(observedIncidence_mx);
  const Rat observedIncidence(observedIncidence_mx);

  int numStoredExemplars = int(mxGetScalar(numStoredExemplars_mx));

  //---------------------------------------------------------------------
  //
  // Call the real computational function for each set of model params
  //
  //---------------------------------------------------------------------

  Rat allModelParams(modelParams_mx);

  mxArray* ll_mx = mxCreateDoubleMatrix(allModelParams.ncols(), 1, mxREAL);
  Rat ll(ll_mx);

  ModelCssm model(objParams,
						observedIncidence,
						numStoredExemplars);

  for (int i = 0; i < allModelParams.ncols(); ++i)
	 {
		Rat modelParams(allModelParams.columnSlice(i));
		ll.at(i) = model.loglikelihoodFor(modelParams);
	 }

  //---------------------------------------------------------------------
  //
  // Clean up
  //
  //---------------------------------------------------------------------

  mxDestroyArray(modelParams_mx);

  mclSetCurrentLocalFunctionTable(save_local_function_table_);

  return ll_mx;
}

/*
 * The function "mlfClassifier" contains the normal interface for the
 * "classifier" function. This function processes any input arguments
 * and passes them to the implementation version of the function,
 * appearing above.
*/
mxArray* mlfClassifier(mxArray* modelParams_mx,
								mxArray* objParams_mx,
								mxArray* observedIncidence_mx,
								mxArray* numStoredExemplars_mx)
{
DOTRACE("mlfClassifier");

  int nargout = 1;

  mlfEnterNewContext(0, 4,
							modelParams_mx, objParams_mx,
							observedIncidence_mx, numStoredExemplars_mx);

  mxArray * ll = Mclassifier(nargout, modelParams_mx, objParams_mx,
									  observedIncidence_mx, numStoredExemplars_mx);

  mlfRestorePreviousContext(0, 4,
									 modelParams_mx, objParams_mx,
									 observedIncidence_mx, numStoredExemplars_mx);

  return mlfReturnValue(ll);
}

/*
 * The function "mlxClassifier" contains the feval interface for the
 * "classifier" function. The feval function calls the implementation
 * version of classifier through this function. This function
 * processes any input arguments and passes them to the implementation
 * version of the function, appearing above.  */
void mlxClassifier(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
DOTRACE("mlxClassifier");

  mxArray * mprhs[4];
  mxArray * mplhs[1];
  int i;
  if (nlhs > 1) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of outputs (1).");
  }
  if (nrhs > 4) {
	 mexErrMsgTxt("Run-time Error: "
					  "The function \"classifier\" was called with more "
					  "than the declared number of inputs (4).");
  }
  for (i = 0; i < 1; ++i) {
	 mplhs[i] = mclGetUninitializedArray();
  }
  for (i = 0; i < 4 && i < nrhs; ++i) {
	 mprhs[i] = prhs[i];
  }
  for (; i < 4; ++i) {
	 mprhs[i] = NULL;
  }

  mlfEnterNewContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  mplhs[0] = Mclassifier(nlhs, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  mlfRestorePreviousContext(0, 4, mprhs[0], mprhs[1], mprhs[2], mprhs[3]);

  plhs[0] = mplhs[0];
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED

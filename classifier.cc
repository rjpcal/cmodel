///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Thu Mar  8 11:10:05 2001
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

#include "rutil.h"

#include "trace.h"

#include <cmath>
#include <fstream.h>

class Num {
public:
  static double erfc(double x) {
	 double z = fabs(x);

	 double t = 1.0/(1.0+0.5*fabs(x));

	 double ans =
      t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
      t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
      t*(-0.82215223+t*0.17087277)))))))))
		;

	 return x >= 0.0 ? ans : 2.0-ans;
  }

  static double gammaln(double xx)
  {
	 double tol = 1e-50;

	 int ival = int(xx);

	 if ( (ival < TABLE_SIZE) && (xx - double(ival) < tol) )
		{
		  if (!filled) fillTable();
		  return lookup[ival];
		}
	 else 
		{
		  // mexPrintf("couldn't use lookup\n");
		  return gammalnEngine(xx);
		}
  }

  static const double SQRT_2 = 1.41421356237;

private:
  static bool filled;
  static const int TABLE_SIZE = 101;
  static double lookup[TABLE_SIZE];

  static void fillTable()
  {
	 for (int i = 0; i < TABLE_SIZE; ++i)
		lookup[i] = gammalnEngine(double(i));

	 filled = true;
  }

  static double gammalnEngine(double xx)
  {
  DOTRACE("gammalnEngine");
	 static double cof[6] = {76.18009172947146,     -86.50532032941677,
									 24.01409824083091,     -1.231738572450155,
									 0.1208650973866179e-2, -0.5395239384953e-5};
	 double x,y,tmp,ser;
	 int j;

	 y = x = xx;
	 tmp = x+5.5;
	 tmp -= (x+0.5)*log(tmp);
	 ser=1.000000000190015;
	 for (j=0;j<=5;j++) ser += cof[j]/++y;
	 return -tmp+log(2.5066282746310005*ser/x);
  }
};

bool Num::filled = false;
double Num::lookup[TABLE_SIZE] = { 0.0 };

void Classifier::forwardProbit(double* diffEvidence,
										 int num_points,
										 double thresh_,
										 double sigmaNoise_,
										 double* ppOut)
{
DOTRACE("Classifier::forwardProbit");

  for (int i = 0; i < num_points; ++i) {
	 double alpha_val = (thresh_-diffEvidence[i]) / sigmaNoise_;

	 ppOut[i] = 0.5*Num::erfc(alpha_val / Num::SQRT_2);
  }
}

///////////////////////////////////////////////////////////////////////
//
// ll = sum(gammaln(1+sum(observedIncidence,2))) - ...
//      sum(sum(gammaln(observedIncidence+1)))  + ...
//      sum(sum(observedIncidence.*log(predictedProbability)));
//
///////////////////////////////////////////////////////////////////////

double Classifier::loglikelihood(const double* predictedProbability1,
											int numPoints,
											const double* observedIncidence1,
											const double* observedIncidence2)
{
DOTRACE("Classifier::loglikelihood");
  double ll = 0.0;

  const double LOG_10_MINUS_50 = -115.1293;

  for(int k = 0; k < numPoints; ++k) {
	 double oi1 = observedIncidence1[k];
	 double oi2 = observedIncidence2[k];

	 // term 1
	 ll += Num::gammaln(1.0 + oi1 + oi2);

	 // term 2
	 ll -= Num::gammaln(1.0+oi1);
	 ll -= Num::gammaln(1.0+oi2);

	 // term3
	 double pp_val = predictedProbability1[k];

	 if (pp_val < 1e-50) ll += oi1 * LOG_10_MINUS_50;
	 else                ll += oi1 * log(pp_val);

	 pp_val = 1.0 - pp_val;

	 if (pp_val < 1e-50) ll += oi2 * LOG_10_MINUS_50;
	 else                ll += oi2 * log(pp_val);
  }

  return ll;
}

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



///////////////////////////////////////////////////////////////////////
//
// ModelCssm member definitions
//
///////////////////////////////////////////////////////////////////////

ModelCssm::ModelCssm(const Rat& objParams_,
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

ModelCssm::~ModelCssm()
{
  delete [] cat2;
  delete [] cat1;
  delete [] mu;
}

// Count the category training exemplars
int ModelCssm::countCategory(const Rat& params, int category) {
  int n = 0;
  for (int i = 0; i < params.mrows(); ++i)
	 {
		if (int(params.at(i,0)) == category)
		  ++n;
	 }
  return n;
}


void ModelCssm::scaleWeights(double* weights,
									  int numRawWeights,
									  int numStoredExemplars,
									  int numTrainingExemplars)
{
DOTRACE("ModelCssm::scaleWeights");

  int mrows = numStoredExemplars*2;
  int ncols = numTrainingExemplars;

  if ( numRawWeights != (mrows*ncols) )
	 mexErrMsgTxt("weights must have "
					  "2*numStoredExemplars*numTrainingExemplars elements");

  for (int i = 0; i < mrows; ++i)
	 {
		double sum_wt = 0.0;
		{
		  for (int ix = i; ix < numRawWeights; ix+=mrows)
			 sum_wt += abs(weights[ix]);
		}
		{
		  for (int ix = i; ix < numRawWeights; ix+=mrows)
			 weights[ix] = abs(weights[ix]) / sum_wt;
		}
	 }
}


void ModelCssm::computeSimilarity(const double* attWeights,
											 const double* storedExemplar1,
											 const double* storedExemplar2,
											 double minkPower,
											 double minkPowerInv)
{
DOTRACE("ModelCssm::computeSimilarity");

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


void ModelCssm::computeSimilarity2(const double* attWeights,
											  const double* storedExemplar1,
											  const double* storedExemplar2)
{
DOTRACE("ModelCssm::computeSimilarity2");

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

  scaleWeights(rawWeights, numRawWeights,
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
  double* predictedProbability = new double[numAllExemplars];
  forwardProbit(mu, numAllExemplars, thresh, sigmaNoise,
					 predictedProbability);

  //---------------------------------------------------------------------
  //
  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.
  //
  //---------------------------------------------------------------------

  double ll = -loglikelihood(predictedProbability,
									  numAllExemplars,
									  observedIncidence.data(),
									  observedIncidence.data()+numAllExemplars);

  delete [] predictedProbability;

  return ll;
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED

///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Thu Mar  8 09:58:28 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"

extern void InitializeModule_classifier(void);
extern void TerminateModule_classifier(void);
extern _mexLocalFunctionTable _local_function_table_classifier;

extern mxArray* mlfClassifier(mxArray* modelParams,
										mxArray* objParams,
										mxArray* observedIncidence,
										mxArray* numExemplars);

extern void mlxClassifier(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]);

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED

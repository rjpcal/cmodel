///////////////////////////////////////////////////////////////////////
//
// rutil.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 19 16:20:06 2001
// written: Thu Feb 14 11:57:01 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef RUTIL_H_DEFINED
#define RUTIL_H_DEFINED

#define validateInput(x) (x) = mlfAssign(&(x), mclVa((x), #x))

#define validateVariable(x) (x) = mlfAssign(&(x), mclVv((x), #x))

static const char vcid_rutil_h[] = "$Header$";
#endif // !RUTIL_H_DEFINED

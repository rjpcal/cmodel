///////////////////////////////////////////////////////////////////////
//
// nrutil.h
//
// created: Thu Apr 19 17:31:40 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef NRUTIL_H_DEFINED
#define NRUTIL_H_DEFINED


#define NR_END 1

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

inline double FMAX(double a, double b)
{ return (a > b) ? a : b; }

inline int IMIN(int a, int b)
{ return (a < b) ? a : b; }

inline double SQR(double a) { return a*a; }

inline double pythag(double a, double b)
{
  double absa,absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

inline double* vector(long nl, long nh)
{
  double* v;
  v = (double*) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  return v-nl+NR_END;
}

inline void free_vector(double* v, long nl, long nh)
{
  free(v+nl-NR_END);
}

inline double** matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  m = (double**) malloc((size_t)((nrow+NR_END)*sizeof(double*)));

  m += NR_END;
  m -= nrl;

  m[nrl]=(double*) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));

  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;++i)
         m[i]=m[i-1]+ncol;

  return m;
}

inline void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free(m[nrl]+ncl-NR_END);
  free(m+nrl-NR_END);
}

static const char vcid_nrutil_h[] = "$Id$ $URL$";
#endif // !NRUTIL_H_DEFINED

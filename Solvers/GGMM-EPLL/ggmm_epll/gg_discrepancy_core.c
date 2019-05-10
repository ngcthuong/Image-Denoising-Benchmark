/*  Code used by gg_discrepancy.m
**
**  Citation:
**  If you use this code please cite:
**  C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
**  restoration with generalized Gaussian mixture model patch
**  priors", arXiv.
**
**  License details as in license.txt
**  ________________________________________
*/

#include <mex.h>
#include <stdio.h>
#include <math.h>

static void usage()
{
  char str[1024];
  sprintf(str, "usage: gg_distance(x, nu, a1, a2, b1, b2, h)\n");
  mexErrMsgTxt(str);
}

// fprintf('%f, ', log1p(exp(-7:7)))
static const double log1pexm_lookup[] =
  { 0.000911, 0.002476, 0.006715, 0.018150, 0.048587, 0.126928, 0.313262, 0.693147, 1.313262, 2.126928, 3.048587, 4.018150, 5.006715, 6.002476, 7.000911 };

static double log1pexp(double x)
{
  int xf;
  if (x <= -7)
    return 0;
  if (x >= 7)
    return x;
  xf = floorf(x);
  return
    (1 - (x - xf)) * log1pexm_lookup[xf + 7] +
    ((x - xf)) * log1pexm_lookup[xf + 8];
}


void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int k, d, N;
  const double* x;
  const double* nu;
  const double* a1;
  const double* a2;
  const double* b1;
  const double* b2;
  const double* h;
  double* ih;
  const mwSize* sx;
  double* f;

  int i;
  double lx, f1, f2, tmp;

  if (nrhs < 7 || nlhs > 1)
    {
      usage();
      return;
    }
  for (k = 0; k < 7; ++k)
    if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]))
      {
	usage();
	return;
      }
  x  = mxGetData(prhs[0]);
  nu = mxGetData(prhs[1]);
  a1 = mxGetData(prhs[2]);
  a2 = mxGetData(prhs[3]);
  b1 = mxGetData(prhs[4]);
  b2 = mxGetData(prhs[5]);
  h  = mxGetData(prhs[6]);

  sx = mxGetDimensions(prhs[0]);

  d  = sx[0];
  N  = sx[1];

  plhs[0] = mxCreateDoubleMatrix(d, N, mxREAL);
  f = (double*) mxGetPr(plhs[0]);

  ih = malloc(d * sizeof(double));
  for (k = 0; k < d; ++k)
    ih[k] = 1 / h[k];
  for (i = 0; i < N; ++i)
    for (k = 0; k < d; ++k)
      {
	lx = x[i * d + k];
	if (isinf(lx) || isnan(lx))
	  {
	    f[i * d + k] = -INFINITY;
	    continue;
	  }
	else
	  {
	    f1 = a1[k] * lx + b1[k];
	    f2 = a2[k] * lx + b2[k];
	    tmp = (f2 - f1) * ih[k];
	    if (nu[k] < 2)
	      tmp = f2 - h[k] * log1pexp(tmp);
	    else
	      if (nu[k] > 2)
		tmp = f1 + h[k] * log1pexp(tmp);
	      else
		tmp = f1;
	    f[i * d + k] = tmp;
	  }
      }
  free(ih);
}

/*  Code to extract patches from an image
**  For usage/input/output details, refer getpatches.m
**
**  Citation:
**  If you use this code please cite:
**  S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
**  GMM-based patch priors for image restoration: Three ingredients for a
**  100x speed-up", arXiv.
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
  sprintf(str, "usage: patches = getpatches(img, P, mask)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int M, N, P;
  int i, j, ip, jp;
  int k, l;
  long Q, q;
  double* patches;
  const char* mask;
  const mwSize* sp;
  const double* img;

  if (nrhs < 3|| nlhs > 1)
    {
      usage();
      return;
    }
  for (k = 0; k < 2; ++k)
    if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]))
      {
	usage();
	return;
      }
  if (!mxIsLogical(prhs[2]))
    {
      usage();
      return;
    }
  img = mxGetData(prhs[0]);
  P       = (int) *((double*) mxGetData(prhs[1]));
  mask    = (char*) mxGetLogicals(prhs[2]);


  sp      = mxGetDimensions(prhs[0]);
  M       = (int) sp[0];
  N       = (int) sp[1];

  Q = 0;
  for (i = 0; i < N; ++i)
    for (j = 0; j < M; ++j)
      if (mask[i * M + j])
	++Q;

  plhs[0] = mxCreateDoubleMatrix(P*P, Q, mxREAL);
  patches = (double*) mxGetPr(plhs[0]);

  q = 0;
  for (i = 0; i < N; ++i)
    for (j = 0; j < M; ++j)
      if (mask[i * M + j])
	{
	  for (k = 0; k < P; ++k)
	    for (l = 0; l < P; ++l)
	      {
		ip = (i + l) % N;
		jp = (j + k) % M;
		patches[(q * P + k) * P + l] = img[ip * M + jp];
	      }
	  ++q;
	}
}

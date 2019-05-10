/*  Code to project patches to make an image
**  For usage/input/output details, refer projpatches.m
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
  sprintf(str, "usage: [img, norm] = projpatches(patches, M, N, mask)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  int M, N, P;
  int i, j, ip, jp;
  int k, l;
  long q;
  const double* patches;
  const char* mask;
  const mwSize* sp;
  double* img;
  double* norm;

  if (nrhs < 4|| nlhs > 2)
    {
      usage();
      return;
    }
  for (k = 0; k < 3; ++k)
    if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]))
      {
	usage();
	return;
      }
  if (!mxIsLogical(prhs[3]))
    {
      usage();
      return;
    }
  patches = mxGetData(prhs[0]);
  M       = (int) *((double*) mxGetData(prhs[1]));
  N       = (int) *((double*) mxGetData(prhs[2]));
  mask    = (char*) mxGetLogicals(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
  img     = (double*) mxGetPr(plhs[0]);
  norm    = (double*) mxGetPr(plhs[1]);

  sp      = mxGetDimensions(prhs[0]);
  P       = (int) sqrt(sp[0]);

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
		img[ip * M + jp] += patches[(q * P + k) * P + l];
		++norm[ip * M + jp];
	      }
	  ++q;
	}
  for (i = 0; i < N; ++i)
    for (j = 0; j < M; ++j)
      if (norm[i * M + j] > 0)
	img[i * M + j] /= norm[i * M + j];
}

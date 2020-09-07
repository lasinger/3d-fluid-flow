/*
Copyright (c) 2020 ETH Zurich

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Author: Katrin Lasinger
*/

#include <vector>
#include <algorithm>
#include <omp.h>
#include "mex.h"
#include <math.h>

double phi_lin(int idx, double x){
	if(idx == 0)
		return 1.0-x;
	return x;
}

/// matlab calling
void interpolateFlowToPart ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *u_grid = (Scalar*)   mxGetPr(prhs[0]); //flow volume: 3x(NxMxL)
  Scalar *part = (Scalar*)   mxGetPr(prhs[1]); //3d particles
  size_t N       = (size_t)  (*mxGetPr(prhs[2]));
  size_t M       = (size_t)  (*mxGetPr(prhs[3]));
  size_t L       = (size_t)  (*mxGetPr(prhs[4]));

  const mwSize* dims = mxGetDimensions(prhs[1]);
  size_t numpart = dims[1]; // list of points
  size_t partDim = dims[0]; // featureDim
 
  plhs[0] = mxCreateDoubleMatrix( 3, numpart, mxREAL); 

  Scalar* u_part = (Scalar*)mxGetPr( plhs[0] );

#pragma omp parallel for //schedule (static)
for (long long i=0; i<numpart; i++ )
{
  //2d part
  Scalar px_i = part[i*partDim];
  Scalar py_i = part[i*partDim+1];
  Scalar pz_i = part[i*partDim+2];

  long xInt = (long)floor(px_i);
  long yInt = (long)floor(py_i);
  long zInt = (long)floor(pz_i);
    
  for (long ix=0; ix<=1; ix++)
	  for (long iy=0; iy<=1; iy++)
		  for (long iz=0; iz<=1; iz++){
			  long xgrid = xInt + ix;
              long ygrid = yInt + iy;
              long zgrid = zInt + iz;
			  
			  if (xgrid<0 || xgrid>N-1 || ygrid<0 || ygrid>M-1 || zgrid<0 || zgrid>L-1)
                    continue;
              
			  Scalar w = phi_lin(ix,px_i-(Scalar)xInt)*phi_lin(iy,py_i-(Scalar)yInt)*phi_lin(iz,pz_i-(Scalar)zInt);
			  size_t idx = zgrid + ygrid*L + xgrid*M*L;

			  u_part[i*3] = u_part[i*3] + u_grid[idx*3]*w;
			  u_part[i*3+1] = u_part[i*3+1] + u_grid[idx*3+1]*w;
			  u_part[i*3+2] = u_part[i*3+2] + u_grid[idx*3+2]*w;
		  }

}

}
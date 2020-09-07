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

double phi0(double x) {
	return 2.0*(x - 0.5)*(x - 1.0);
}
double phi1(double x) {
	return 4.0*x*(1.0 - x);
}
double phi2(double x) {
	return 2.0*x*(x - 0.5);
}

double phi(int idx, double x) {
	if (idx == 0)
		return phi0(x);
	if (idx == 1)
		return phi1(x);
	return phi2(x);
}

double phi0Prime(double x) {
	return 4.0*x - 3.0;
}
double phi1Prime(double x) {
	return 4.0 - 8.0*x;
}
double phi2Prime(double x) {
	return 4.0*x - 1.0;
}

double phiPrime(int idx, double x) {
	if (idx == 0)
		return phi0Prime(x);
	if (idx == 1)
		return phi1Prime(x);
	return phi2Prime(x);
}



/// matlab calling
void derivativeFlowPartP2P1 ( int nlhs, mxArray *plhs[],
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
 
  //input: u_grid, part, N, M, L
  //output: u_part

  plhs[0] = mxCreateDoubleMatrix( 3, numpart, mxREAL);
  plhs[1] = mxCreateDoubleMatrix( 3, numpart, mxREAL);
  plhs[2] = mxCreateDoubleMatrix( 3, numpart, mxREAL);

  Scalar* dudx_part = (Scalar*)mxGetPr( plhs[0] );
  Scalar* dudy_part = (Scalar*)mxGetPr( plhs[1] );
  Scalar* dudz_part = (Scalar*)mxGetPr( plhs[2] );

#pragma omp parallel for //schedule (static)
for (long long i=0; i<numpart; i++ )
{
  //2d part
  Scalar px_i = part[i*partDim]/2.0;
  Scalar py_i = part[i*partDim+1]/2.0;
  Scalar pz_i = part[i*partDim+2]/2.0;

  long xInt = (long)floor(px_i);
  long yInt = (long)floor(py_i);
  long zInt = (long)floor(pz_i);
    
  for (long ix=0; ix<=2; ix++)
	  for (long iy=0; iy<=2; iy++)
		  for (long iz=0; iz<=2; iz++){
			  long xgrid = xInt * 2 + ix;
			  long ygrid = yInt * 2 + iy;
			  long zgrid = zInt * 2 + iz;
			  
			  if (xgrid<0 || xgrid>N-1 || ygrid<0 || ygrid>M-1 || zgrid<0 || zgrid>L-1)
                    continue;
              
			  Scalar dx = phiPrime(ix,px_i-(Scalar)xInt)*phi(iy,py_i-(Scalar)yInt)*phi(iz,pz_i-(Scalar)zInt);
			  Scalar dy = phi(ix,px_i-(Scalar)xInt)*phiPrime(iy,py_i-(Scalar)yInt)*phi(iz,pz_i-(Scalar)zInt);
			  Scalar dz = phi(ix,px_i-(Scalar)xInt)*phi(iy,py_i-(Scalar)yInt)*phiPrime(iz,pz_i-(Scalar)zInt);

			  size_t idx = zgrid + ygrid*L + xgrid*M*L;

			  dudx_part[i*3] = dudx_part[i*3] + u_grid[idx*3]*dx;
			  dudx_part[i*3+1] = dudx_part[i*3+1] + u_grid[idx*3+1]*dx;
			  dudx_part[i*3+2] = dudx_part[i*3+2] + u_grid[idx*3+2]*dx;

			  dudy_part[i*3] = dudy_part[i*3] + u_grid[idx*3]*dy;
			  dudy_part[i*3+1] = dudy_part[i*3+1] + u_grid[idx*3+1]*dy;
			  dudy_part[i*3+2] = dudy_part[i*3+2] + u_grid[idx*3+2]*dy;

			  dudz_part[i*3] = dudz_part[i*3] + u_grid[idx*3]*dz;
			  dudz_part[i*3+1] = dudz_part[i*3+1] + u_grid[idx*3+1]*dz;
			  dudz_part[i*3+2] = dudz_part[i*3+2] + u_grid[idx*3+2]*dz;
		  }

}


}
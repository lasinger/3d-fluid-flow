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

#define _USE_MATH_DEFINES

#include <numeric>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "mex.h"
#include <math.h>
#include "FastExp64.h"

/// matlab calling
void computeGradient_Polycam ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *output1;

  Scalar *part_I_x = (Scalar*)   mxGetPr(prhs[0]); //x coord of 2d part pos
  Scalar *part_I_y = (Scalar*)   mxGetPr(prhs[1]); //y coord of 2d part pos
  Scalar *part_X = (Scalar*)   mxGetPr(prhs[2]); //x coord of 3d part pos
  Scalar *part_Y = (Scalar*)   mxGetPr(prhs[3]); //y coord of 3d part pos
  Scalar *part_Z = (Scalar*)   mxGetPr(prhs[4]); //z coord of 3d part pos
  Scalar *part_int = (Scalar*)   mxGetPr(prhs[5]); //particle intensity
  size_t numpart       = (size_t)  (*mxGetPr(prhs[6]));
  Scalar sigma       = (Scalar)  (*mxGetPr(prhs[7]));
  Scalar *residualImg = (Scalar*)   mxGetPr(prhs[8]); //residual img
  size_t N       = (size_t)  (*mxGetPr(prhs[9]));
  size_t M       = (size_t)  (*mxGetPr(prhs[10]));
  Scalar *a_x = (Scalar*)   mxGetPr(prhs[11]); //parameters of polynomial camera for x
  Scalar *a_y = (Scalar*)   mxGetPr(prhs[12]); //parameters of polynomial camera for y
 

  //input: particle 2d pos, 3d pos, sigma, numpart
  //output: gradient per particle (for one camera)

  Scalar radius = 3*sigma;
  Scalar radiusSq = radius*radius;
  Scalar sigmaSq = sigma*sigma;
  Scalar invSigmaSq = 1.0/sigmaSq;
  
  Scalar integralWeight = 1.0/sqrt(2.0*M_PI*sigmaSq);

  plhs[0] = mxCreateDoubleMatrix( numpart, 1, mxREAL);
  Scalar* part_grad_x = (Scalar*)mxGetPr( plhs[0] );


  plhs[1] = mxCreateDoubleMatrix( numpart, 1, mxREAL);
  Scalar* part_grad_y = (Scalar*)mxGetPr( plhs[1] );

  plhs[2] = mxCreateDoubleMatrix( numpart, 1, mxREAL);
  Scalar* part_grad_z = (Scalar*)mxGetPr( plhs[2] );

#pragma omp parallel for //schedule (static)
for (long long i=0; i<numpart; i++ )
{
  //2d part
  Scalar px_i = part_I_x[i];
  Scalar py_i = part_I_y[i];

  //3d part
  Scalar pX = part_X[i];
  Scalar pY = part_Y[i];
  Scalar pZ = part_Z[i];
  
  Scalar A_dX[] = {0, 1, 0, 0, 2*pX, pY, 0, pZ, 0, 0, 3*pX*pX, 2*pX*pY, pY*pY, 0, 2*pX*pZ, pY*pZ, 0, pZ*pZ, 0};
  Scalar A_dY[] = {0, 0, 1, 0, 0, pX, 2*pY, 0, pZ, 0, 0, pX*pX, pX*2*pY, 3*pY*pY, 0, pX*pZ, 2*pY*pZ, 0, pZ*pZ};
  Scalar A_dZ[] = {0, 0, 0, 1, 0, 0, 0, pX, pY, 2*pZ, 0, 0, 0, 0, pX*pX, pX*pY, pY*pY, pX*2*pZ, pY*2*pZ};

  // particle intensity
  Scalar partintPart = -2.0f/sigmaSq*part_int[i];

  for (long qx = (long)floor(px_i-radius); qx<=(long)ceil(px_i+radius); qx++ )
    for (long qy = (long)floor(py_i-radius); qy<=(long)ceil(py_i+radius); qy++ )  
    {
	  if(qx<0 || qx>=N || qy<0 || qy>=M)
		  continue;
		
	  Scalar xDiff = qx-px_i;
	  Scalar yDiff = qy-py_i;
		
	  Scalar distSq = xDiff*xDiff + yDiff*yDiff;
      if(distSq > radiusSq)
		  continue;

	  Scalar weight_pi_q = integralWeight* fast_exp_64(-distSq*invSigmaSq) * residualImg[qy + qx*M];
	
	  part_grad_x[i] += (-std::inner_product(a_x,a_x+19,A_dX,0.0)*xDiff-std::inner_product(a_y,a_y+19,A_dX,0.0)*yDiff) * weight_pi_q * partintPart; 
	  part_grad_y[i] += (-std::inner_product(a_x,a_x+19,A_dY,0.0)*xDiff-std::inner_product(a_y,a_y+19,A_dY,0.0)*yDiff) * weight_pi_q * partintPart; 
	  part_grad_z[i] += (-std::inner_product(a_x,a_x+19,A_dZ,0.0)*xDiff-std::inner_product(a_y,a_y+19,A_dZ,0.0)*yDiff) * weight_pi_q * partintPart; 

	}
}

}
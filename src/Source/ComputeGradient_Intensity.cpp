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

#include <vector>
#include <algorithm>
#include <omp.h>
#include "mex.h"
#include <math.h>
#include "FastExp64.h"

/// matlab calling
void computeGradient_Intensity ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *output1;

  Scalar *part_I_x = (Scalar*)   mxGetPr(prhs[0]); //x coord of 2d part pos
  Scalar *part_I_y = (Scalar*)   mxGetPr(prhs[1]); //y coord of 2d part pos
  size_t numpart       = (size_t)  (*mxGetPr(prhs[2]));
  Scalar sigma       = (Scalar)  (*mxGetPr(prhs[3]));
  Scalar *residualImg = (Scalar*)   mxGetPr(prhs[4]); //residual img
  size_t N       = (size_t)  (*mxGetPr(prhs[5]));
  size_t M       = (size_t)  (*mxGetPr(prhs[6]));
 
  //input: particle 2d pos, 3d pos, sigma, numpart
  //output: gradient per particle (for one camera)

  Scalar radius = 3*sigma;
  Scalar radiusSq = radius*radius;
  Scalar sigmaSq = sigma*sigma;
  Scalar invSigmaSq = 1.0/sigmaSq;
  
  Scalar integralWeight = 1.0/sqrt(2.0*M_PI*sigmaSq);

  plhs[0] = mxCreateDoubleMatrix( numpart, 1, mxREAL);
  Scalar* part_grad_x = (Scalar*)mxGetPr( plhs[0] );

#pragma omp parallel for //schedule (static)
for (long long i=0; i<numpart; i++ )
{
  //2d part
  Scalar px_i = part_I_x[i];
  Scalar py_i = part_I_y[i];


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


      part_grad_x[i] += integralWeight* fast_exp_64(-distSq*invSigmaSq) * residualImg[qy + qx*M];
	
	}
}


}
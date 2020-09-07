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
void computeGradient ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *output1;

  Scalar *part_I_x = (Scalar*)   mxGetPr(prhs[0]); //x coord of 2d part pos
  Scalar *part_I_y = (Scalar*)   mxGetPr(prhs[1]); //y coord of 2d part pos
  Scalar *part_int = (Scalar*)   mxGetPr(prhs[2]); //particle intensity
  size_t numpart       = (size_t)  (*mxGetPr(prhs[3]));
  Scalar sigma       = (Scalar)  (*mxGetPr(prhs[4]));
  Scalar *ProjectionPart1_I_x = (Scalar*)   mxGetPr(prhs[5]);
  Scalar *ProjectionPart1_I_y = (Scalar*)   mxGetPr(prhs[6]);
  Scalar *ProjectionPart2_I_x = (Scalar*)   mxGetPr(prhs[7]);
  Scalar *ProjectionPart2_I_y = (Scalar*)   mxGetPr(prhs[8]);
  Scalar *ProjectionPart3_I_x = (Scalar*)   mxGetPr(prhs[9]);
  Scalar *ProjectionPart3_I_y = (Scalar*)   mxGetPr(prhs[10]);
  Scalar *residualImg = (Scalar*)   mxGetPr(prhs[11]); //residual img
  size_t N       = (size_t)  (*mxGetPr(prhs[12]));
  size_t M       = (size_t)  (*mxGetPr(prhs[13]));

  //input: particle 2d pos, 3d pos, sigma, numpart
  //output: gradient per particle (for one camera)

  Scalar radius = 3*sigma;
  Scalar radiusSq = radius*radius;
  Scalar sigmaSq = sigma*sigma;
  Scalar invSigmaSq = 1.0/sigmaSq;
  
  Scalar integralWeight = 1.0/sqrt(2.0*M_PI*sigmaSq);
  plhs[0] = mxCreateDoubleMatrix(3, numpart, mxREAL);
  Scalar* part_grad = (Scalar*)mxGetPr(plhs[0]);


#pragma omp parallel for //schedule (static)
for (long long i=0; i<numpart; i++ )
{
  //2d part
  Scalar px_i = part_I_x[i];
  Scalar py_i = part_I_y[i];

  Scalar P1xi = ProjectionPart1_I_x[i];
  Scalar P1yi = ProjectionPart1_I_y[i];
  Scalar P2xi = ProjectionPart2_I_x[i];
  Scalar P2yi = ProjectionPart2_I_y[i];
  Scalar P3xi = ProjectionPart3_I_x[i];
  Scalar P3yi = ProjectionPart3_I_y[i];

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

	  part_grad[i*3] += (P1xi*xDiff+P1yi*yDiff) * weight_pi_q * partintPart;
	  part_grad[i*3+1] += (P2xi*xDiff+P2yi*yDiff) * weight_pi_q * partintPart;
	  part_grad[i*3+2] += (P3xi*xDiff+P3yi*yDiff) * weight_pi_q * partintPart;
	
	}
}


}
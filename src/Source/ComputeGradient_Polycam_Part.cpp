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
void computeGradient_Polycam_Part ( int nlhs, mxArray *plhs[],
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
  Scalar *dudpx = (Scalar*)mxGetPr(prhs[13]); //
  Scalar *dudpy = (Scalar*)mxGetPr(prhs[14]); //
  Scalar *dudpz = (Scalar*)mxGetPr(prhs[15]); //

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
  Scalar duxdpx = dudpx[i * 3];
  Scalar duydpx = dudpx[i * 3 + 1];
  Scalar duzdpx = dudpx[i * 3 + 2];
  Scalar duxdpy = dudpy[i * 3 ];
  Scalar duydpy = dudpy[i * 3 + 1];
  Scalar duzdpy = dudpy[i * 3 + 2];
  Scalar duxdpz = dudpz[i * 3];
  Scalar duydpz = dudpz[i * 3 + 1];
  Scalar duzdpz = dudpz[i * 3 + 2];
  
  Scalar A_dX[] = { 0, 1+duxdpx, duydpx,   duzdpx,   2*pX*(1+duxdpx), (1+duxdpx)*pY+duydpx*pX, 2*pY*duydpx,     (1+duxdpx)*pZ+duzdpx*pX, duydpx*pZ+duzdpx*pY,     2*pZ*duzdpx,     3*pX*pX*(1+duxdpx), 2*pX*pY*(1+duxdpx)+pX*pX*duydpx, pY*pY*(1+duxdpx)+2*pX*pY*duydpx, 3*pY*pY*duydpx,     2*pX*pZ*(1+duxdpx)+pX*pX*duzdpx, pY*pZ*(1+duxdpx)+pX*pZ*duydpx+pX*pY*duzdpx, 2*pY*pZ*duydpx+pY*pY*duzdpx,     pZ*pZ*(1+duxdpx)+2*pX*pZ*duzdpx, pZ*pZ*duydpx+2*pY*pZ*duzdpx     };
  Scalar A_dY[] = { 0, duxdpy,   1+duydpy, duzdpy,   2*pX*duxdpy,     (1+duydpy)*pX+duxdpy*pY, 2*pY*(1+duydpy), duxdpy*pZ+duzdpy*pX,     (1+duydpy)*pZ+duzdpy*pY, 2*pZ*duzdpy,     3*pX*pX*duxdpy,     2*pX*pY*duxdpy+pX*pX*(1+duydpy), pY*pY*duxdpy+2*pX*pY*(1+duydpy), 3*pY*pY*(1+duydpy), 2*pX*pZ*duxdpy+pX*pX*duzdpy,     pY*pZ*duxdpy+pX*pZ*(1+duydpy)+pX*pY*duzdpy, 2*pY*pZ*(1+duydpy)+pY*pY*duzdpy, pZ*pZ*duxdpy+2*pX*pZ*duzdpy,     pZ*pZ*(1+duydpy)+2*pY*pZ*duzdpy };
  Scalar A_dZ[] = { 0, duxdpz,   duydpz,   1+duzdpz, 2*pX*duxdpz,     duydpz*pX+duxdpz*pY,     2*pY*duydpz,     (1+duzdpz)*pX+duxdpz*pZ, (1+duzdpz)*pY+duydpz*pZ, 2*pZ*(1+duzdpz), 3*pX*pX*duxdpz,     2*pX*pY*duxdpz+pX*pY*duydpz,     pY*pY*duxdpz+2*pX*pY*duydpz,     3*pY*pY*duydpz,     2*pX*pZ*duxdpz+pX*pX*(1+duzdpz), pY*pZ*duxdpz+pX*pZ*duydpz+pX*pY*(1+duzdpz), 2*pY*pY*duydpz+pY*pY*(1+duzdpz), pZ*pZ*duxdpz+2*pX*pZ*(1+duzdpz), pZ*pZ*duydpz+2*pY*pZ*(1+duzdpz) };


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
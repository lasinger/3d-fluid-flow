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
void render2dPart_allcam ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  const mxArray *cell_array_ptr_in; // per cam 2d array of xy coordinates of 2d part pos
  const mwSize *dims;
  std::vector< Scalar*> part;
  const mxArray *cell_element_ptr;
  mwIndex jcell;
  cell_array_ptr_in = prhs[0];

  size_t numcam = mxGetNumberOfElements(cell_array_ptr_in);

  
  part.resize(numcam);

  for (jcell=0; jcell<numcam; jcell++) {
	  cell_element_ptr = mxGetCell(cell_array_ptr_in,jcell);

	  part[jcell] = (Scalar*) mxGetPr(cell_element_ptr);
  }

  Scalar *part_int = (Scalar*)   mxGetPr(prhs[1]); //particle intensity
  size_t numpart       = (size_t)  (*mxGetPr(prhs[2]));
  Scalar sigma       = (Scalar)  (*mxGetPr(prhs[3]));
  size_t N       = (size_t)  (*mxGetPr(prhs[4]));
  size_t M       = (size_t)  (*mxGetPr(prhs[5]));
  Scalar renderRadFactor = 4.0f;
  if (nrhs>=7) {
	  renderRadFactor = (Scalar)(*mxGetPr(prhs[6]));
  }
  Scalar *sigma_list;
  if(sigma < 0){
	  if(nrhs<8){
		  sigma = 1.0;
		  mexPrintf("not enough input args, so sigma was set to 1");
	  }
	  sigma_list = (Scalar*) mxGetPr(prhs[7]); // no check here if really exists... could check and set sigma fixed in case not exist
  }


  //input: particle 2d pos, 3d pos, sigma, numpart
  //output: gradient per particle (for one camera)
  //sigma = 1.0;

  Scalar radius = renderRadFactor*sigma;
  Scalar radiusSq = radius*radius;
  Scalar sigmaSq = sigma*sigma;
  Scalar invSigmaSq = 1.0/sigmaSq;
  
  Scalar integralWeight = 1.0/sqrt(2.0*M_PI*sigmaSq);

  std::vector< mxArray*> renderedImg(numcam, NULL);
  for (size_t i = 0; i < (mwIndex)numcam; i++)
	renderedImg[i] = mxCreateDoubleMatrix( M, N, mxREAL); 

  // for imagewatch use: @mem(((((renderedImg)._Myfirst)[0])._Myfirst),  FLOAT64, 1, 700, 400, 700*8)

#pragma omp parallel for schedule(static)
for (int cam = 0; cam < numcam;cam++){
	Scalar* renderedImg_cam = (Scalar*)mxGetPr( renderedImg[cam] );
	for (size_t i = 0; i < numpart; i++)
	{
		//no need to render particles with 0 intensity
		if(part_int[i]<=1e-5)
			continue;

		//2d part
		Scalar px_i = part[cam][i*2];
		Scalar py_i = part[cam][i*2+1];

		if (sigma == -1){
			radius = renderRadFactor*sigma_list[i];
			radiusSq = radius*radius;
			sigmaSq = sigma_list[i]*sigma_list[i];
			invSigmaSq = 1.0/sigmaSq;
		}


		long floorX = (long)floor(px_i - radius);
		long ceilX = (long)ceil(px_i + radius);
		long floorY = (long)floor(py_i - radius);
		long ceilY = (long)ceil(py_i + radius);
		
		for (long qx = floorX; qx <= ceilX; qx++)
			for (long qy = floorY; qy <= ceilY; qy++)
			{
				if (qx < 0 || qx >= N || qy < 0 || qy >= M)
					continue;

				Scalar xDiff = (Scalar)qx - px_i;
				Scalar yDiff = (Scalar)qy - py_i;

				Scalar distSq = xDiff*xDiff + yDiff*yDiff;
				if (distSq > radiusSq)
					continue;

				size_t idx1d = qy + qx*M;

				Scalar intensity = integralWeight * part_int[i] * fast_exp_64(-distSq * invSigmaSq);
				
				renderedImg_cam[idx1d] = renderedImg_cam[idx1d] + intensity;
			}
	}
}



  //------------------------------------------------------------------------------
  //write back to matlab

   mwSize dimsOutCell[1];
  dimsOutCell[0] = numcam;
  // indices
  mxArray *cell_array_ptr;
  cell_array_ptr = mxCreateCellArray(1, dimsOutCell);
  
  for (size_t i = 0; i < (mwIndex)numcam; i++) {
	  mxSetCell(cell_array_ptr, i, mxDuplicateArray(renderedImg[i]));
	  mxDestroyArray(renderedImg[i]);
  }

  plhs[0] = cell_array_ptr;

}
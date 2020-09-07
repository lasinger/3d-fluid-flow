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

#include <Eigen/Dense>
using namespace Eigen;

/// matlab calling
void proj2d ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *part = (Scalar*)mxGetPr(prhs[0]);
 
  Scalar *P_vec = (Scalar*)mxGetPr(prhs[1]);

  Matrix<Scalar, 3, 4> P;

  for(int j=0;j<4;j++)
    for (int i=0;i<3;i++)
		P(i,j) = P_vec[i+j*3];

  const mwSize* dims = mxGetDimensions(prhs[0]);
  size_t numpart = dims[1]; // list of points
  size_t partDim = dims[0]; // featureDim

  plhs[0] = mxCreateDoubleMatrix(2,numpart, mxREAL);
  Scalar* part2d = (Scalar*)mxGetPr(plhs[0]);


  plhs[1] = mxCreateDoubleMatrix(3,numpart, mxREAL);
  Scalar* xhom = (Scalar*)mxGetPr(plhs[1]);
 
#pragma omp parallel for //schedule (static)
  for (long long i = 0; i < numpart; i++)
  {

	Vector4d pt;
	pt << part[i*3],part[i*3+1], part[i*3+2], 1;


	Vector3d pt2d_hom;
	pt2d_hom = P * pt;

	xhom[i*3] =pt2d_hom(0);
	xhom[i*3+1] = pt2d_hom(1);
	xhom[i*3+2] = pt2d_hom(2);
	part2d[i*2] = pt2d_hom(0) / pt2d_hom(2);
	part2d[i*2+1] = pt2d_hom(1) / pt2d_hom(2);


  }

}
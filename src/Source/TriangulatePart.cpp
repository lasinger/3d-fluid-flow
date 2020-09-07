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

#include <Eigen/Dense>
using namespace Eigen;


// Code adapted from Johannes L. Schoeneberger
// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
//https://github.com/colmap/colmap/blob/master/src/base/triangulation.cc
Eigen::Vector3d TriangulatePoint(const Eigen::Matrix<double, 3, 4>& proj_matrix1,
                                 const Eigen::Matrix<double, 3, 4>& proj_matrix2,
                                 const double point1x, const double point1y,
                                 const double point2x, const double point2y) {
  Eigen::Matrix4d A;

  A.row(0) = point1x * proj_matrix1.row(2) - proj_matrix1.row(0);
  A.row(1) = point1y * proj_matrix1.row(2) - proj_matrix1.row(1);
  A.row(2) = point2x * proj_matrix2.row(2) - proj_matrix2.row(0);
  A.row(3) = point2y * proj_matrix2.row(2) - proj_matrix2.row(1);

  Eigen::JacobiSVD<Eigen::Matrix4d> svd(A, Eigen::ComputeFullV);

  return svd.matrixV().col(3).hnormalized();
}


 Vector2d polynomialCameraForward( Matrix<double, 19, 1> a_x, Matrix<double, 19, 1> a_y, Vector3d pt3d ){
	double pX = pt3d[0];
	double pY = pt3d[1];
	double pZ = pt3d[2];

	Matrix<double, 19, 1> A;
	A << 1, pX, pY, pZ, pX*pX, pX*pY, pY*pY, pX*pZ, pY*pZ, pZ*pZ, pX*pX*pX, pX*pX*pY, pX*pY*pY, pY*pY*pY, pX*pX*pZ, pX*pY*pZ, pY*pY*pZ, pX*pZ*pZ, pY*pZ*pZ;

	Vector2d pt2d;
	pt2d[0] = A.dot(a_x);
	pt2d[1] = A.dot(a_y);

	return pt2d;

 }

template <class Vec19_>
Vector3d TriangulatePoint_poly(Vec19_ a_x0, Vec19_ a_y0, Vec19_ a_x1, Vec19_ a_y1, Vector2d p0, Vector2d p1)
{
	double thresh = 0.001;
	Vector3d pt3d(0,0,0); //seems to be sufficient enough to just initialize with 0

	for(int iterations = 0; iterations<1000; iterations++ )
	{

        // check new projection
        Vector2d pt2dEst_1 = polynomialCameraForward(a_x0,a_y0,pt3d);
        Vector2d pt2dEst_2 = polynomialCameraForward(a_x1,a_y1,pt3d);

        Vector2d diff1 = p0-pt2dEst_1;
        Vector2d diff2 = p1-pt2dEst_2;

        if (diff1.norm()<thresh && diff2.norm()<thresh){
            break;
		}

        double pX = pt3d[0];
		double pY = pt3d[1];
		double pZ = pt3d[2];

        Vec19_ A_dX, A_dY, A_dZ;
	    A_dX << 0, 1, 0, 0, 2*pX, pY, 0, pZ, 0, 0, 3*pX*pX, 2*pX*pY, pY*pY, 0, 2*pX*pZ, pY*pZ, 0, pZ*pZ, 0;
	    A_dY << 0, 0, 1, 0, 0, pX, 2*pY, 0, pZ, 0, 0, pX*pX, pX*2*pY, 3*pY*pY, 0, pX*pZ, 2*pY*pZ, 0, pZ*pZ;
	    A_dZ << 0, 0, 0, 1, 0, 0, 0, pX, pY, 2*pZ, 0, 0, 0, 0, pX*pX, pX*pY, pY*pY, pX*2*pZ, pY*2*pZ;

		MatrixXd A(4,3);

        A(0,0) = A_dX.dot(a_x0);
        A(1,0) = A_dX.dot(a_y0);
		A(2,0) = A_dX.dot(a_x1);
        A(3,0) = A_dX.dot(a_y1);
		A(0,1) = A_dY.dot(a_x0);
        A(1,1) = A_dY.dot(a_y0);
		A(2,1) = A_dY.dot(a_x1);
        A(3,1) = A_dY.dot(a_y1);
        A(0,2) = A_dZ.dot(a_x0);
        A(1,2) = A_dZ.dot(a_y0);
		A(2,2) = A_dZ.dot(a_x1);
        A(3,2) = A_dZ.dot(a_y1);

        Vector4d delta_x;
		delta_x << diff1(0), diff1(1), diff2(0), diff2(1);

        //delta_X = A\delta_x;
		Vector3d delta_X = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(delta_x);

        pt3d = pt3d+delta_X;

	}

	return pt3d;

}

Vector3d polynomialCameraImg2World(Matrix<double, 19, 1> a_x, Matrix<double, 19, 1> a_y, Vector2d pt2d, double depth)
{
	double thresh = 0.001;
	Vector3d pt3d(0,0,depth); //seems to be sufficient enough to just initialize with 0

	for(int iterations = 0; iterations<1000; iterations++ )
	{

        // check new projection
        Vector2d pt2dEst = polynomialCameraForward(a_x,a_y,pt3d);

        Vector2d diff = pt2d-pt2dEst;

        if (diff.norm()<thresh){
            break;
		}

        double pX = pt3d[0];
		double pY = pt3d[1];
		double pZ = pt3d[2];

        Matrix<double, 19, 1> A_dX, A_dY, A_dZ;
	    A_dX << 0, 1, 0, 0, 2*pX, pY, 0, pZ, 0, 0, 3*pX*pX, 2*pX*pY, pY*pY, 0, 2*pX*pZ, pY*pZ, 0, pZ*pZ, 0;
	    A_dY << 0, 0, 1, 0, 0, pX, 2*pY, 0, pZ, 0, 0, pX*pX, pX*2*pY, 3*pY*pY, 0, pX*pZ, 2*pY*pZ, 0, pZ*pZ;
	    //A_dZ << 0, 0, 0, 1, 0, 0, 0, pX, pY, 2*pZ, 0, 0, 0, 0, pX*pX, pX*pY, pY*pY, pX*2*pZ, pY*2*pZ;

		MatrixXd A(2,2);

        A(0,0) = A_dX.dot(a_x);
        A(1,0) = A_dX.dot(a_y);
		A(0,1) = A_dY.dot(a_x);
        A(1,1) = A_dY.dot(a_y);

        Vector2d delta_x;
		delta_x << diff(0), diff(1);

		Vector2d delta_X = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(delta_x);

        pt3d(0) = pt3d(0)+delta_X(0);
        pt3d(1) = pt3d(1)+delta_X(1);

	}

	return pt3d;
}


/// matlab calling - polynomial camera model
void triangulatePartPoly ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;
  typedef Matrix<double, 19, 1> Vector19d;
  

  Scalar *output1;

  const mxArray *cell_array_ptr_in = prhs[0]; // per cam 2d array of xy coordinates+intensity of 2d part pos
  std::vector< Scalar*> part;
  const mxArray *cell_element_ptr;
  const mxArray *cell_element_ptr2;
  mwIndex jcell;

  int numcam = mxGetNumberOfElements(cell_array_ptr_in);

  
  part.resize(numcam);
  std::vector<int> numPart(numcam,0);

  for (jcell=0; jcell<numcam; jcell++) {
	  cell_element_ptr = mxGetCell(cell_array_ptr_in,jcell);

	  const mwSize* dims = mxGetDimensions(cell_element_ptr);

	  numPart[jcell] = dims[1];

	  part[jcell] = (Scalar*) mxGetPr(cell_element_ptr);
  }

  const mxArray *cell_array_ptr_in_a = prhs[1]; // per cam: 19x2 coefficients of camera model
  std::vector< Vector19d > a_x;
  std::vector< Vector19d > a_y;
  a_x.resize(numcam);
  a_y.resize(numcam);

  for (jcell=0; jcell<numcam; jcell++) {
	  cell_element_ptr = mxGetCell(cell_array_ptr_in_a,jcell);

	  Scalar* a_vec = (Scalar*) mxGetPr(cell_element_ptr);
	  for (int i=0;i<19;i++){
		a_x[jcell](i) = a_vec[i];
		a_y[jcell](i) = a_vec[i+19];
	  }

  }
  int N       = (int)  (*mxGetPr(prhs[2]));
  int M       = (int)  (*mxGetPr(prhs[3]));
  int L       = (int)  (*mxGetPr(prhs[4]));
  Scalar triangError = (Scalar)  (*mxGetPr(prhs[5]));
  Scalar triangErrorSq = triangError*triangError;

  Scalar integralWeight = 1.0/sqrt(2.0*M_PI);

  std::vector< Scalar > part3d; //list of 3d pts incl intensity (so 4d)

  std::vector<std::vector<Scalar> > local(omp_get_max_threads());
  std::vector<int> numTriangulatedPts(omp_get_max_threads(),0);
#pragma omp parallel //num_threads(1)
{
	int np = omp_get_num_threads();
	
	int currThreadNum = omp_get_thread_num();
	
  //iterate over 2d peaks of cam0
  //#pragma omp parallel for
#pragma omp for //schedule(static)
  for(int i=0;i<numPart[0];i++){

	Vector2d pt_ref;
	pt_ref << part[0][i*3],part[0][i*3+1];
	Scalar distSq = pow(floor(pt_ref(0)+0.5)-pt_ref(0),2)+pow(floor(pt_ref(1)+0.5)-pt_ref(1),2);
	Scalar int_ref = part[0][i*3+2]/integralWeight/ exp(-distSq); //sigma=1

    //first slice (at z=0)
	Vector3d pt3d_front = polynomialCameraImg2World( a_x[0], a_y[0], pt_ref, 0 );
    //last slice (at z=L-1)
    Vector3d pt3d_back = polynomialCameraImg2World( a_x[0], a_y[0], pt_ref, L-1 );
	
    //project 3d start and end pt to image
    Vector2d front_epi = polynomialCameraForward( a_x[1], a_y[1], pt3d_front );
	Vector2d back_epi = polynomialCameraForward( a_x[1], a_y[1], pt3d_back );

    Scalar maxBB_x = std::max(front_epi(0),back_epi(0)) + triangError;
    Scalar maxBB_y = std::max(front_epi(1),back_epi(1)) + triangError;
    Scalar minBB_x = std::min(front_epi(0),back_epi(0)) - triangError;
    Scalar minBB_y = std::min(front_epi(1),back_epi(1)) - triangError;

	Scalar lowerTerm = std::sqrt(std::pow(back_epi(0)-front_epi(0),2)+std::pow(back_epi(1)-front_epi(1),2));

    int numPtsOnLine=0;
    int numFound3dpts=0;
    int idxCurrPt = numTriangulatedPts[currThreadNum];

    //get distance of pts to epipolar line
    for (int j=1; j<numPart[1];j++){


		if(part[1][j*3]<minBB_x || part[1][j*3]>maxBB_x || part[1][j*3+1]<minBB_y || part[1][j*3+1]>maxBB_y)
			continue;

       // compute distance
       Scalar upperTerm = abs((back_epi(0)-front_epi(0))*(front_epi(1)-part[1][j*3+1])-(front_epi(0)-part[1][j*3])*(back_epi(1)-front_epi(1)));
       Scalar d = upperTerm/lowerTerm;

       if (d > triangError)
           continue;

       numPtsOnLine++;

	   Vector3d pt3d = TriangulatePoint_poly(a_x[0],a_y[0],a_x[1],a_y[1],pt_ref,Vector2d(part[1][j*3],part[1][j*3+1]));
	   

       // can add a check if pt is actually within volume (since my
       // epipolar line is alway from z=0 to z=L-1 it could be that part
       // of the ray is outside of the volume to the sides)
       // possible check for x and y should be sufficient
	   // assuming that N is already size including padding (so minus also for upper bound)
       if( pt3d(0) < 0 || pt3d(0) > N-1 || pt3d(1) < 0 || pt3d(1) > M-1 || pt3d(2) < 0 || pt3d(2) > L-1)
           continue;

		bool match=true;

	   // loop over remaining cams
	   for (int c=2;c<numcam;c++){
		   if(match == false) //if no match in prev other cam
			   continue;
		   match = false;
		   Vector2d pt2d_proj = polynomialCameraForward( a_x[c], a_y[c], pt3d );

		   Scalar maxBB_x_oc = pt2d_proj(0) + triangError;
           Scalar maxBB_y_oc = pt2d_proj(1) + triangError;
           Scalar minBB_x_oc = pt2d_proj(0) - triangError;
           Scalar minBB_y_oc = pt2d_proj(1) - triangError;

		   for (int k=1; k<numPart[c];k++){

			   if(part[c][k*3]<minBB_x_oc || part[c][k*3]>maxBB_x_oc || part[c][k*3+1]<minBB_y_oc || part[c][k*3+1]>maxBB_y_oc)
				   continue;

			   if (pow(part[c][k*3]-pt2d_proj(0),2) + (pow(part[c][k*3+1]-pt2d_proj(1),2)) >triangErrorSq)
					continue;

			   match = true;
		   }
	   }

	   if(match){
		   numTriangulatedPts[currThreadNum]++;
		   numFound3dpts++; 
		   local[currThreadNum].push_back(pt3d(0));
		   local[currThreadNum].push_back(pt3d(1));
		   local[currThreadNum].push_back(pt3d(2));
		   local[currThreadNum].push_back(int_ref);

	   }

	}
    
    //distribute intensity over found ref pts
    if(numFound3dpts>1){
		for(int ii=0;ii<numFound3dpts;ii++)
			local[currThreadNum][(idxCurrPt+ii)*4+3] = int_ref*(Scalar)4/(Scalar)(3+numFound3dpts);
	}
  }
}

for (int p = 0; p < omp_get_max_threads(); ++p){
	if (numTriangulatedPts[p] >0)
		part3d.insert(part3d.end(),local[p].begin(),local[p].end());
}

int numTriangulatedPtsAll = part3d.size()/4;

  //------------------------------------------------------------------------------
  //write back to matlab

  plhs[0] = mxCreateDoubleMatrix( 4, numTriangulatedPtsAll, mxREAL);
  output1  = mxGetPr(plhs[0]);
  for (int i=0;i<part3d.size();i++)
      output1[i] = part3d[i];

 

}

/// matlab calling
void triangulatePart ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *output1;

  

  const mxArray *cell_array_ptr_in = prhs[0]; // per cam 2d array of xy coordinates+intensity of 2d part pos
  std::vector< Scalar*> part;
  const mxArray *cell_element_ptr;
  const mxArray *cell_element_ptr2;
  mwIndex jcell;

  int numcam = mxGetNumberOfElements(cell_array_ptr_in);

  
  part.resize(numcam);
  std::vector<int> numPart(numcam,0);

  for (jcell=0; jcell<numcam; jcell++) {
	  cell_element_ptr = mxGetCell(cell_array_ptr_in,jcell);

	  const mwSize* dims = mxGetDimensions(cell_element_ptr);

	  numPart[jcell] = dims[1];

	  part[jcell] = (Scalar*) mxGetPr(cell_element_ptr);
  }

  const mxArray *cell_array_ptr_in_P = prhs[1]; // per cam 2d array of xy coordinates+intensity of 2d part pos
  const mxArray *cell_array_ptr_in_C = prhs[2];
  std::vector< Matrix<Scalar, 3, 4> > P;
  std::vector< Vector3d > C;
  P.resize(numcam);
  C.resize(numcam);

  for (jcell=0; jcell<numcam; jcell++) {
	  cell_element_ptr = mxGetCell(cell_array_ptr_in_P,jcell);
	  cell_element_ptr2 = mxGetCell(cell_array_ptr_in_C,jcell);

	  Scalar* p_vec = (Scalar*) mxGetPr(cell_element_ptr);
	  for(int j=0;j<4;j++)
		for (int i=0;i<3;i++)
			  P[jcell](i,j) = p_vec[i+j*3];

	  Scalar* c_vec = (Scalar*) mxGetPr(cell_element_ptr2);
		for (int i=0;i<3;i++)
			  C[jcell](i) = c_vec[i];
  }
  int N       = (int)  (*mxGetPr(prhs[3]));
  int M       = (int)  (*mxGetPr(prhs[4]));
  int L       = (int)  (*mxGetPr(prhs[5]));
  Scalar triangError = (Scalar)  (*mxGetPr(prhs[6]));
  Scalar triangErrorSq = triangError*triangError;

	  

 
  Scalar test = P[0](1,2);
  Scalar integralWeight = 1.0/sqrt(2.0*M_PI);

  std::vector< Scalar > part3d; //list of 3d pts incl intensity (so 4d)
  
  Matrix3d M_inv = P[0].block<3,3>(0,0).inverse();
  Matrix3d M2_inv = P[1].block<3,3>(0,0).inverse();
  Vector3d Pcol = P[0].col(3);

  std::vector<std::vector<Scalar> > local(omp_get_max_threads());
  std::vector<int> numTriangulatedPts(omp_get_max_threads(),0);
#pragma omp parallel //num_threads(1)
{

	int np = omp_get_num_threads();
	
	int currThreadNum = omp_get_thread_num();
	
  //iterate over 2d peaks of cam0
  //#pragma omp parallel for
#pragma omp for //schedule(static)
  for(int i=0;i<numPart[0];i++){

	Vector3d pt_ref;
	pt_ref << part[0][i*3],part[0][i*3+1],1;
	Scalar distSq = pow(floor(pt_ref(0)+0.5)-pt_ref(0),2)+pow(floor(pt_ref(1)+0.5)-pt_ref(1),2);
	Scalar int_ref = part[0][i*3+2]/integralWeight/ exp(-distSq); //sigma=1

	Vector3d pt;
	pt = M_inv * (pt_ref - Pcol);

    Vector3d x1=C[0];
    Vector3d x2=pt;
    Vector3d x21 = x2-x1;
    Scalar dvN = -x21(2);

    //first slice (at z=0)
    Scalar dx1Z = x1(2);
    Scalar t = dx1Z/dvN;
    Vector4d pt3d_front;
	pt3d_front << x1 + t*x21,1;
    //last slice (at z=L-1)
    dx1Z = x1(2)-(L-1);
    t = dx1Z/dvN;
    Vector4d pt3d_back;
	pt3d_back << x1 + t*x21,1;



    //project 3d start and end pt to image
    Vector3d front_epi = P[1]*pt3d_front;
    front_epi = front_epi/front_epi(2);
    Vector3d back_epi = P[1]*pt3d_back;
    back_epi = back_epi/back_epi(2);

    Scalar maxBB_x = std::max(front_epi(0),back_epi(0)) + triangError;
    Scalar maxBB_y = std::max(front_epi(1),back_epi(1)) + triangError;
    Scalar minBB_x = std::min(front_epi(0),back_epi(0)) - triangError;
    Scalar minBB_y = std::min(front_epi(1),back_epi(1)) - triangError;

	Scalar lowerTerm = std::sqrt(std::pow(back_epi(0)-front_epi(0),2)+std::pow(back_epi(1)-front_epi(1),2));

    int numPtsOnLine=0;
    int numFound3dpts=0;
    int idxCurrPt = numTriangulatedPts[currThreadNum];

    //get distance of pts to epipolar line
    for (int j=1; j<numPart[1];j++){

		
		// if outside bounding box ignore
		if(part[1][j*3]<minBB_x || part[1][j*3]>maxBB_x || part[1][j*3+1]<minBB_y || part[1][j*3+1]>maxBB_y)
			continue;

       // compute distance
       Scalar upperTerm = abs((back_epi(0)-front_epi(0))*(front_epi(1)-part[1][j*3+1])-(front_epi(0)-part[1][j*3])*(back_epi(1)-front_epi(1)));
       Scalar d = upperTerm/lowerTerm;

       if (d > triangError)
           continue;

       numPtsOnLine++;

	   Vector3d pt3d = TriangulatePoint(P[0],P[1],pt_ref(0),pt_ref(1),part[1][j*3],part[1][j*3+1]);
	   

       // can add a check if pt is actually within volume (since my
       // epipolar line is alway from z=0 to z=L-1 it could be that part
       // of the ray is outside of the volume to the sides)
       // possible check for x and y should be sufficient
	   // assuming that N is already size including padding (so minus also for upper bound)
       if( pt3d(0) < 0 || pt3d(0) > N-1 || pt3d(1) < 0 || pt3d(1) > M-1 || pt3d(2) < 0 || pt3d(2) > L-1)
           continue;

	   
		Vector4d pt3d_hom;
		pt3d_hom << pt3d,1;

		bool match=true;

	   // loop over remaining cams
	   for (int c=2;c<numcam;c++){
		   if(match == false) //if no match in prev other cam
			   continue;
		   match = false;
		   Vector3d pt2d_proj = P[c]*pt3d_hom;
		   pt2d_proj = pt2d_proj/pt2d_proj(2);

		   Scalar maxBB_x_oc = pt2d_proj(0) + triangError;
           Scalar maxBB_y_oc = pt2d_proj(1) + triangError;
           Scalar minBB_x_oc = pt2d_proj(0) - triangError;
           Scalar minBB_y_oc = pt2d_proj(1) - triangError;

		   for (int k=1; k<numPart[c];k++){
			   //Vector3d pt_cam_other;
			   //pt_cam_other << part[c][k*3],part[c][k*3+1],1;

			   //if(pt_cam_other(0)<minBB_x_oc || pt_cam_other(0)>maxBB_x_oc || pt_cam_other(1)<minBB_y_oc || pt_cam_other(1)>maxBB_y_oc)
			//	   continue;

			  // if (pow(pt_cam_other(0)-pt2d_proj(0),2) + (pow(pt_cam_other(1)-pt2d_proj(1),2)) >triangErrorSq)
				//	continue;

			   if(part[c][k*3]<minBB_x_oc || part[c][k*3]>maxBB_x_oc || part[c][k*3+1]<minBB_y_oc || part[c][k*3+1]>maxBB_y_oc)
				   continue;

			   if (pow(part[c][k*3]-pt2d_proj(0),2) + (pow(part[c][k*3+1]-pt2d_proj(1),2)) >triangErrorSq)
					continue;

			   match = true;
		   }
	   }

	   if(match){
		   numTriangulatedPts[currThreadNum]++;
		   numFound3dpts++;
		   local[currThreadNum].push_back(pt3d(0));
		   local[currThreadNum].push_back(pt3d(1));
		   local[currThreadNum].push_back(pt3d(2));
		   local[currThreadNum].push_back(int_ref);
	   }

	}
    
    //distribute intensity over found ref pts
    if(numFound3dpts>1){
		for(int ii=0;ii<numFound3dpts;ii++)
			local[currThreadNum][(idxCurrPt+ii)*4+3] = int_ref*(Scalar)4/(Scalar)(3+numFound3dpts);
	}


  }


}

for (int p = 0; p < omp_get_max_threads(); ++p){
	if (numTriangulatedPts[p] >0)
		part3d.insert(part3d.end(),local[p].begin(),local[p].end());
}



int numTriangulatedPtsAll = part3d.size()/4;

  //input: cell array of 2d pts per cam incl intensity, P (and others?) 
  //output: triangulated 3d pts (get 0:n part per 2d pt)

  //------------------------------------------------------------------------------
  //write back to matlab

  plhs[0] = mxCreateDoubleMatrix( 4, numTriangulatedPtsAll, mxREAL);
  output1  = mxGetPr(plhs[0]);
  for (int i=0;i<part3d.size();i++)
      output1[i] = part3d[i];

}
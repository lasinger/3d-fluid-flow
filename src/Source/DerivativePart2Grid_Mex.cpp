#include <iostream>
#include <string.h>

#include "mex.h"

#include "DerivativePart2Grid_lin.cpp"
#include "DerivativePart2Grid_quad.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
	int doP2P1 = 0;
	if (nrhs >= 6) {
		doP2P1 = (int)(*mxGetPr(prhs[5]));
	}
	if (doP2P1 <= 0) {
		derivativePart2Grid_lin(nlhs, plhs, nrhs, prhs);
	}
	else {
		derivativePart2Grid_quad(nlhs, plhs, nrhs, prhs);
	}
}

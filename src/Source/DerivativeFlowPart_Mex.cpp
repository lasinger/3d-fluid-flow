#include <iostream>
#include <string.h>

#include "mex.h"

#include "DerivativeFlowPart.cpp"
#include "DerivativeFlowPartP2P1.cpp"


void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
	int doP2P1 = 0;
	if (nrhs >= 6) {
		doP2P1 = (int)(*mxGetPr(prhs[5]));
	}
	if (doP2P1 <= 0) {
		derivativeFlowPart(nlhs, plhs, nrhs, prhs);
	}
	else {
		derivativeFlowPartP2P1(nlhs, plhs, nrhs, prhs);
	}
}

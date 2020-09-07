#include <iostream>
#include <string.h>

#include "mex.h"
#include "InterpolateFlowToPart.cpp"
#include "InterpolateFlowToPartP2P1.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
	int doP2P1 = 0;
	if(nrhs >= 6){
		doP2P1       = (int)  (*mxGetPr(prhs[5]));
	}
	if(doP2P1 <= 0){
		interpolateFlowToPart( nlhs, plhs, nrhs, prhs );
	}else{
		interpolateFlowToPartP2P1( nlhs, plhs, nrhs, prhs );
	}
}

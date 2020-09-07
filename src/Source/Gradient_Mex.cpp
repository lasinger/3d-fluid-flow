#include <iostream>
#include <string.h>

#include "mex.h"
#include "ComputeGradient.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    computeGradient( nlhs, plhs, nrhs, prhs );
}

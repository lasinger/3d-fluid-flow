#include <iostream>
#include <string.h>

#include "mex.h"
#include "ComputeGradient_Polycam.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    computeGradient_Polycam( nlhs, plhs, nrhs, prhs );
}

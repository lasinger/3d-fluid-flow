#include <iostream>
#include <string.h>

#include "mex.h"
#include "ComputeGradient_Intensity.cpp"


void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    computeGradient_Intensity( nlhs, plhs, nrhs, prhs );
}

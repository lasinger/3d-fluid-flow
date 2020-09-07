#include <iostream>
#include <string.h>

#include "mex.h"
#include "ComputeGradient_Polycam_Part.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    computeGradient_Polycam_Part( nlhs, plhs, nrhs, prhs );
}

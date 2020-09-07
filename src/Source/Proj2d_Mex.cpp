#include <iostream>
#include <string.h>

#include "mex.h"
#include "Proj2d.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    proj2d( nlhs, plhs, nrhs, prhs );
}

#include <iostream>
#include <string.h>

#include "mex.h"
#include "TriangulatePart.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    triangulatePart( nlhs, plhs, nrhs, prhs );
}

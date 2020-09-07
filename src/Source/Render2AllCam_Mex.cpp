#include <iostream>
#include <string.h>

#include "mex.h"
#include "Render2dPart_allcam.cpp"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    render2dPart_allcam( nlhs, plhs, nrhs, prhs );
}

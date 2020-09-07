% Copyright (c) 2020 ETH Zurich
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
% Author: Katrin Lasinger

function [estU_grid, estU_prev, estU_part, reproj1, part2d1, xhom1 ] = inertial_flow(estU_grid, estU_prev, estPart, tau_k, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize,numpart, femSpace)

if ~exist('femSpace','var')
	femSpace = 0;
end

estU_i = estU_grid + tau_k.*(estU_grid-estU_prev);
estU_prev = estU_grid;
estU_grid = estU_i;

part2d1 = cell(1,numCam);
xhom1 = cell(1,numCam);

estU_part  = interpFlowToPartmex( estU_grid, estPart(1:3,:)./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace );


for cam=1:numCam
    [part2d1{cam}, xhom1{cam}] = proj2d(estPart(1:3,:)+estU_part,P{cam});
end
numpart = size(estPart,2);
reproj1 = render2dallcammex(part2d1,estPart(4,:),numpart,sigma,N2d,M2d );

end


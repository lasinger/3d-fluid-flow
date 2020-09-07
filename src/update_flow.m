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

function [ rightTerm, leftTerm, reprojTmp1, part2dTmp1, xhomTmp1, estTmp_u_grid, estTmp_pressure, estTmp_u_part ] = update_flow( estPart, Lu, flow_k_DataDeriv, flow_k, Lpl_flow, estPressure, costOrig, reproj0, img0, img1, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, lambda, kinetic, iterations, Nabla_incomp, Div_incomp, Lpl_incomp, Nabla,divMethod,femSpace, alpha )

    if ~exist('femSpace','var')
        femSpace = 0;
    end

    numpart = size(estPart,2);
    
    tol = 1e-4;
    

    flow_k_bar = flow_k - (1/Lu) * ( flow_k_DataDeriv );

    if divMethod == 0 % no proxmap
        flow_new = flow_k_bar;
        estTmp_pressure = estPressure;
    else
        [estTmp_pressure,flow_new] = solveDiv0Constraint( iterations, Lds,Mds,Nds, flow_k_bar, estPressure, Lu, kinetic, tol, Nabla_incomp, Div_incomp, Lpl_incomp );
    end
    estTmp_u_grid = reshape(flow_new,Lds*Mds*Nds,3)';

    diff =  1/16 * lambda * (sum(([Nabla * flow_k(1:Nds*Mds*Lds); Nabla * flow_k(1+Nds*Mds*Lds : 2*Nds*Mds*Lds); Nabla * flow_k(1 + 2*Nds*Mds*Lds : end) ]).^2) - sum(([Nabla * flow_new(1:Nds*Mds*Lds); Nabla * flow_new(1+Nds*Mds*Lds : 2*Nds*Mds*Lds); Nabla * flow_new(1 + 2*Nds*Mds*Lds : end) ]).^2));
    
    if divMethod == 0 && alpha > 0
        diff = diff + 1/16 * alpha * lambda*((sum(Div_incomp*flow_k).^2 - sum(Div_incomp*flow_new).^2));
    end
    
    part2dTmp1 = cell(1,numCam);
    xhomTmp1 = cell(1,numCam);

    estTmp_u_part  = interpFlowToPartmex( estTmp_u_grid, estPart(1:3,:)./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace );

    for cam=1:numCam
        [part2dTmp1{cam}, xhomTmp1{cam}] = proj2d(estPart(1:3,:)+estTmp_u_part,P{cam});
    end
    reprojTmp1 = render2dallcammex(part2dTmp1,estPart(4,:),numpart,sigma,N2d,M2d );

    rightTerm = costOrig + sum(flow_k_DataDeriv.*(flow_new-flow_k)) + Lu/2 * sum((flow_new-flow_k).^2) + diff; 
    leftTerm = cost_dataterm(img0,reproj0,img1,reprojTmp1);

end


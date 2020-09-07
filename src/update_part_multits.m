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

function [ rightTerm, leftTerm, estTmp_p, reproj, part2d, xhom, estU_part ] = update_part_multits ( estPart, estPartInt,ref_idx, Lp, Dc, estU_grid, costOrig, imgs, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, femSpace )

if ~exist('femSpace','var')
	femSpace = 0;
end

estTmp_p = estPart - (1/Lp).*Dc;

numpart = size(estTmp_p,2);

% when particle moves also flow at particle location changes
% (while flow on grid stays constant)
num_ts = length(imgs);
estU_part = cell(1,num_ts-1);
for i=1:num_ts-1
    estU_part{i} = interpFlowToPartmex( estU_grid{i}, estTmp_p(1:3,:)./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace );
end

part2d = cell(1,num_ts);
xhom = cell(1,num_ts);
reproj = cell(1,num_ts);

for i=1:num_ts
    part2d{i} = cell(1,numCam);
    xhom{i} = cell(1,numCam);
    reproj{i} = cell(1,numCam);

    for cam=1:numCam
        if i<ref_idx
            [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estTmp_p(1:3,:)+estU_part{i},P{cam});
        elseif i>ref_idx
            [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estTmp_p(1:3,:)+estU_part{i-1},P{cam});
        else
            [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estTmp_p(1:3,:),P{cam});
        end
    end
    reproj{i} = render2dallcammex(part2d{i},estPartInt,numpart,sigma,N2d,M2d );
end
rightTerm = costOrig + sum(sum(Dc.*(estTmp_p-estPart))) + Lp/2 * sum((estTmp_p(:)-estPart(:)).^2); %svds((estTmp_p-estPart),1); %
leftTerm = cost_dataterm_multits(imgs,reproj);

end
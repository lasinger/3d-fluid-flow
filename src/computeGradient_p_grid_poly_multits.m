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

function [ Dc ] = computeGradient_p_grid_poly_multits( imgs, reproj, part2d, estPart, estU_grid, estU_part, numCam, N2d, M2d, sigma, P, numpart, ref_idx, Nds, Mds, Lds, stepsize, femSpace )

Dc = zeros(3,numpart);

for i=1:length(imgs)
    for cam=1:numCam

        % compute gradient
        diff = imgs{i}{cam} - reproj{i}{cam};
        
        dudpx = zeros(3,numpart);
        dudpy = zeros(3,numpart);
        dudpz = zeros(3,numpart);
        if i~=ref_idx
            idx=i;
            if i>ref_idx
                idx=i-1;
            end
            [dudpx, dudpy, dudpz] = derivativeFlowPartmex( estU_grid{idx}./stepsize, (estPart(1:3,:)+estU_part{idx}(1:3,:))./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace );
        end
        
        if i<ref_idx
            partTmp = estPart(1:3,:)+estU_part{i}(1:3,:);
        elseif i>ref_idx
            partTmp = estPart(1:3,:)+estU_part{i-1}(1:3,:);
        else
            partTmp = estPart(1:3,:);
        end
        

        [part_grad_x0, part_grad_y0, part_grad_z0] = gradientpolyPartmex(part2d{i}{cam}(1,:),part2d{i}{cam}(2,:),partTmp(1,:),partTmp(2,:),partTmp(3,:),estPart(4,:),numpart,sigma,diff(:),N2d,M2d,P{cam}(:,1),P{cam}(:,2),dudpx, dudpy, dudpz);
        %[part_grad_x0, part_grad_y0, part_grad_z0] = gradientpartpolymex(part2d{i}{cam}(1,:),part2d{i}{cam}(2,:),partTmp(1,:),partTmp(2,:),partTmp(3,:),estPart(4,:),numpart,sigma,diff{i}(:),N2d,M2d,P{cam}(:,1),P{cam}(:,2));
        
        Dc = Dc + [part_grad_x0';part_grad_y0';part_grad_z0'];

    end
end

Dc = Dc/numCam/length(imgs);
Dc = -Dc;

end


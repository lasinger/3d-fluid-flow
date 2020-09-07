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

function [ gradUgrid ] = computeGradient_gridU_poly( img1, reproj1, part2d1, estPart, estU_part, numCam, N2d, M2d, sigma, P, numpart, part, Nds, Mds, Lds, femSpace )
 
if ~exist('femSpace','var')
	femSpace = 0;
end

gradUpart = zeros(3,numpart); % p, c, u

for cam=1:numCam
            
    % compute gradient
    diff1 = img1{cam} - reproj1{cam};

    [part_grad_x1, part_grad_y1, part_grad_z1] = gradientpolymex(part2d1{cam}(1,:),part2d1{cam}(2,:),estPart(1,:)+estU_part(1,:),estPart(2,:)+estU_part(2,:),estPart(3,:)+estU_part(3,:),estPart(4,:),numpart,sigma,diff1(:),N2d,M2d,P{cam}(:,1),P{cam}(:,2));
    
    gradUpart = gradUpart + [part_grad_x1';part_grad_y1';part_grad_z1'];
end

gradUpart = gradUpart/numCam;
gradUpart = -gradUpart;


if femSpace == 0
    gradUgrid = interp3c_partToGrid_trilinear(gradUpart, numpart, part, Nds, Mds, Lds);
else
    gradUgrid = interp3c_partToGrid_triquadratic(gradUpart, numpart, part, Nds, Mds, Lds );
end


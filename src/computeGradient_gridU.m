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

function [ gradUgrid ] = computeGradient_gridU( img1, reproj1, part2d1, partInt, xhom1, numCam, N2d, M2d, sigma, P, numpart, part, Nds, Mds, Lds, femSpace )
 
if ~exist('femSpace','var')
	femSpace = 0;
end

gradUpart = zeros(3,numpart); % p, c, u

for cam=1:numCam
            
    % compute gradient
    diff1 = img1{cam} - reproj1{cam};
    
    xhom11=xhom1{cam}(1,:);
    xhom12=xhom1{cam}(2,:);
    xhom13=xhom1{cam}(3,:);
    xhom13_2=xhom13.^2;
    Pcam = P{cam};
    ProjPart1_I_x = ((Pcam(3,1).*xhom11)./xhom13_2 - Pcam(1,1)./xhom13);
    ProjPart1_I_y = ((Pcam(3,1).*xhom12)./xhom13_2 - Pcam(2,1)./xhom13);
    ProjPart2_I_x = ((Pcam(3,2).*xhom11)./xhom13_2 - Pcam(1,2)./xhom13);
    ProjPart2_I_y = ((Pcam(3,2).*xhom12)./xhom13_2 - Pcam(2,2)./xhom13);
    ProjPart3_I_x = ((Pcam(3,3).*xhom11)./xhom13_2 - Pcam(1,3)./xhom13);
    ProjPart3_I_y = ((Pcam(3,3).*xhom12)./xhom13_2 - Pcam(2,3)./xhom13);
    part_grad = gradientmex(part2d1{cam}(1,:),part2d1{cam}(2,:),partInt,numpart,sigma,ProjPart1_I_x,ProjPart1_I_y,ProjPart2_I_x,ProjPart2_I_y,ProjPart3_I_x,ProjPart3_I_y,diff1(:),N2d,M2d);
    
    gradUpart = gradUpart + part_grad;

end

gradUpart = gradUpart/numCam;

gradUpart = -gradUpart;

gradUgrid = derivativePart2Gridmex(gradUpart, part, Nds, Mds, Lds, femSpace);




end

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

function [ Dc ] = computeGradient_intensity_multits( imgs, reproj, part2d, numCam, N2d, M2d, sigma, numpart )

Dc = zeros(1,numpart); % p, c, u

for i=1:length(imgs)
    for cam=1:numCam

        % compute gradient
        diff = imgs{i}{cam} - reproj{i}{cam};

       [part_grad0] = gradientintmex(part2d{i}{cam}(1,:),part2d{i}{cam}(2,:),numpart,sigma,diff(:),N2d,M2d);

        Dc = Dc + part_grad0';


    end
end
Dc = Dc/numCam/length(imgs);
Dc = -Dc;

end


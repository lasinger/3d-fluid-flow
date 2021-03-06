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

function [estPart, estPart_prev, reproj ] = inertial_intensity_multits(estPart, estPart_prev, part2d, tau_k, numCam, sigma, N2d, M2d)

    idx=4;
    estPart_i = estPart(idx,:) + tau_k.*(estPart(idx,:)-estPart_prev(idx,:));
    estPart_prev(idx,:) = estPart(idx,:);
    estPart(idx,:) = estPart_i;
    estPart(idx,estPart(idx,:)<0) = 0;
    
    numpart = size(estPart,2);
    num_ts=length(part2d);
    reproj = cell(1,num_ts);
    for i=1:num_ts
        reproj{i} = render2dallcammex(part2d{i},estPart(idx,:),numpart,sigma,N2d,M2d );
    end

end


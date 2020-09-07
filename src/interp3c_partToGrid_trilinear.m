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

function [grid3c] = interp3c_partToGrid_trilinear(part_val3c, numpart, part, Nds, Mds, Lds )
% bilinear interpolation to grid

grid3c = zeros(3,Nds*Mds*Lds);

for i=1:numpart
        
    % get relevant grid nbs and their weights

    xInt = floor(part(1,i));
    yInt = floor(part(2,i));
    zInt = floor(part(3,i));

    for ix=0:1
        for iy=0:1
            for iz=0:1
                
                xgrid = xInt + ix;
                ygrid = yInt + iy;
                zgrid = zInt + iz;
                if xgrid<0 || xgrid>Nds-1 || ygrid<0 || ygrid>Mds-1 || zgrid<0 || zgrid>Lds-1
                    continue;
                end
                
                w = (1-abs(xgrid-part(1,i))) * (1-abs(ygrid-part(2,i))) * (1-abs(zgrid-part(3,i)));
                
                % we have matlab indices starting from 1 (and particle pos starting from 0)
                idx = zgrid + ygrid*Lds + (xgrid)*Mds*Lds + 1;    
                grid3c(:,idx) = grid3c(:,idx) + part_val3c(:,i)*w;
                
            end
        end

    end




end
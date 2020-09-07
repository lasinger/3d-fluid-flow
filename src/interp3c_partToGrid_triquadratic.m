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

function [grid3c] = interp3c_partToGrid_triquadratic(part_val3c, numpart, part, Nds, Mds, Lds )
% triquadratic interpolation to grid

grid3c = zeros(3,Nds*Mds*Lds);

for i=1:numpart
        
    % get relevant grid nbs and their weights

    % divide by 2 to have pressure grid as reference
    xInt = floor(part(1,i)/2);
    yInt = floor(part(2,i)/2);
    zInt = floor(part(3,i)/2);

    for ix=0:2
        for iy=0:2
            for iz=0:2
                
                xgrid = xInt*2 + ix;
                ygrid = yInt*2 + iy;
                zgrid = zInt*2 + iz;
                if xgrid<0 || xgrid>Nds-1 || ygrid<0 || ygrid>Mds-1 || zgrid<0 || zgrid>Lds-1
                    continue;
                end
                
                w = phi(ix,part(1,i)/2-xInt)*phi(iy,part(2,i)/2-yInt)*phi(iz,part(3,i)/2-zInt);
                
                % we have matlab indices starting from 1 (and particle pos starting from 0)
                idx = zgrid + ygrid*Lds + (xgrid)*Mds*Lds + 1;   
                grid3c(:,idx) = grid3c(:,idx) + part_val3c(:,i)*w;
                
            end
        end

    end
end

end

function [val] = phi0(x)
	val = 2.0*(x-0.5)*(x-1.0);
end
function [val] = phi1(x)
	val = 4.0*x*(1.0-x);
end
function [val] = phi2(x)
	val = 2.0*x*(x-0.5);
end

function [val] = phi(idx, x)
    if idx == 0
		val = phi0(x);
    elseif idx == 1
		val = phi1(x);
    else
		val = phi2(x);
    end
end

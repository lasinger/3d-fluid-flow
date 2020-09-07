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

function [ rightTerm, leftTerm, reproj, estTmp_c ] = update_intensity_multits( estPartInt, part2d, Lc, Dc, costOrig, imgs, numCam, sigma, N2d, M2d, sparsityNorm, my)

    estTmp_c = estPartInt - (1/Lc).*Dc;

    if sparsityNorm == 0
        % 0-norm
        zeroPartMask = Lc/2 * estTmp_c.^2 < my | estTmp_c<0;
        estTmp_c(zeroPartMask) = 0;
    elseif sparsityNorm == 1
        % 1-norm
        zeroPartMask = estTmp_c<my/Lc;
        estTmp_c(zeroPartMask) = 0;
        estTmp_c(~zeroPartMask) = estTmp_c(~zeroPartMask) - my/Lc;
    else
        %no sparsity
        estTmp_c = max(estTmp_c,0);
    end

    numpart = size(estTmp_c,2);
    num_ts=length(part2d);
    reproj = cell(1,num_ts);
    for i=1:num_ts
        reproj{i} = render2dallcammex(part2d{i},estTmp_c,numpart,sigma,N2d,M2d );
    end
    
    %add diff of old sparsity minus new sparsity (since it's part of the
    %energy)
    diff = 0;
    if sparsityNorm == 0
        diff = my * (nnz(estPartInt) - nnz(estTmp_c));
    elseif sparsityNorm == 1
        diff = my * (sum(estPartInt) - sum(estTmp_c));
    end
           
    rightTerm = costOrig + diff + sum(Dc.*(estTmp_c-estPartInt)) + Lc/2 * sum((estTmp_c-estPartInt).^2);
    leftTerm = cost_dataterm_multits(imgs,reproj);
end


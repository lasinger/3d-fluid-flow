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

function [ pts ] = get2Dpeaks( img, epsilon )

    img=img/max(img(:));
    maxImg = ordfilt2(img,9,ones(3,3));

    %use only voxel that are bright enough and suppress non-maxima
    [ppY,ppX] = ind2sub(size(img),find(maxImg==img & img>epsilon));
    vals1 = zeros(length(ppY),1);
    for i=1:length(ppY)
        vals1(i) = img(ppY(i),ppX(i));
    end

    %subpixel
    [ppY,ppX] = subpixFit2d(img,ppY,ppX,1);
    pts=[ppX,ppY,vals1]';

end


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

function [reproj, residual, pts3d] = triangulate_part(img, P, triangError, N, M, L, sigma, intThresh,displayFigures)
% input should be residual images (orig img - already reconstructed
% particles)
% then detect peaks and finally triangulate

numCam = size(P,1);
[M2d,N2d] = size(img{1});

numTriangulatedPts = 0;
pts3d=[];
triangErrorSq = triangError*triangError;

% peak detection in all 2D images - only for one timestep
for c=1:numCam
    % pts in pixel coordinates (starting with 1)
    pts1_2d{c} = get2Dpeaks(img{c},intThresh);
    % get pts starting from 0
    pts1_2d{c}(1:2,:) = pts1_2d{c}(1:2,:) - 1;
end

if length(P{1}) == 4 % pinhole camera

    for i=1:numCam
        [~,~,V] = svd(P{i});
        C{i}  = V(1:end-1,end)/V(end,end);
        vv{i} = getViewVector(P{i},C{i},N2d/2,M2d/2);

        MM = P{i}(:,1:3);
        M_inv{i} = inv(MM);
        [Rinv,Kinv] = qr(M_inv{i});
        R{i} = Rinv'; K{i} = inv(Kinv);

        [r1 r2 r3] = dcm2angle(R{i});
        view_angle = abs(r2*180/pi);
    end
    
     tic
     pts3d = triangulatePartmex(pts1_2d,P,C,N,M,L,triangError);
     toc
 
else % polynomial camera
    tic
    pts3d = triangulatePartPolymex(pts1_2d,P,N,M,L,triangError);
    toc
end

% generate rendered particle images and residual images
if ~isempty(pts3d)
    for cam=1:numCam
        [part2d{cam}, xhom{cam}] = proj2d(pts3d(1:3,:),P{cam});
    end
    reproj = render2dallcammex(part2d,pts3d(4,:),size(pts3d,2),sigma,N2d,M2d );
else
    for cam=1:numCam
        reproj{cam} = zeros(size(img{cam}));
    end
end

for cam=1:numCam
   residual{cam} = img{cam} - reproj{cam};
end

if displayFigures
    figure;
    imshow(img{1},[]);
    figure;
    imshow(reproj{1},[]);
    figure;
    imshow(residual{1},[]);
end
end
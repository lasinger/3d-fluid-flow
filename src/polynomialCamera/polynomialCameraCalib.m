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

function [a] = polynomialCameraCalib(pts3d,pts2d)
% see "The Cubic Rational Polynomial Camera Model" by Hartley and Saxena
% esp for details on data normalization

numpts = size(pts2d,1);

% get translation & scaling (Tx and TX matrix) for normalization of 2d and
% 3d data
m2d = mean(pts2d);
s = sum(sum((pts2d-m2d).^2))^0.5;
sxh = (sqrt(2)*numpts)/s;

m3d = mean(pts3d);
S = sum(sum((pts3d-m3d).^2))^0.5;
SXh = (sqrt(3)*numpts)/S;

Tx = sxh*[1 0 -m2d(1); 0 1 -m2d(2); 0 0 1/sxh];
TX = SXh*[1 0 0 -m3d(1); 0 1 0 -m3d(2); 0 0 1 -m3d(3); 0 0 0 1/SXh];

pts3d(:,4) = ones(numpts,1);
pts2d(:,3) = ones(numpts,1);

pts2dOrig = pts2d;
pts3dOrig = pts3d;
pts2d = (Tx * pts2d')';
pts3d = (TX * pts3d')';

A = [ones(numpts,1),pts3d(:,1), pts3d(:,2),pts3d(:,3),pts3d(:,1).^2,...
    pts3d(:,1).*pts3d(:,2),pts3d(:,2).^2,pts3d(:,1).*pts3d(:,3),pts3d(:,2).*pts3d(:,3),pts3d(:,3).^2,...
    pts3d(:,1).^3,pts3d(:,1).^2.*pts3d(:,2),pts3d(:,1).*pts3d(:,2).^2,pts3d(:,2).^3,pts3d(:,1).^2.*pts3d(:,3),...
    pts3d(:,1).*pts3d(:,2).*pts3d(:,3),pts3d(:,2).^2.*pts3d(:,3),pts3d(:,1).*pts3d(:,3).^2,pts3d(:,2).*pts3d(:,3).^2];

b1 = pts2d(:,1);
b2 = pts2d(:,2);
% get coefficients for normalized input and output
a2(:,1) = A\b1;
a2(:,2) = A\b2;

% iterative Levenberg-Marquardt method
opts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', ...
                                   'Diagnostics','on', ...
                                   'Display','iter', ...
                                   'FinDiffType','central', ...
                                   'FunValCheck','on', ...
                                   'TolFun',1e-10, ...
                                   'TolX',1e-10);

rp_u = @(aOptim) A*aOptim-b1;
a2(:,1) = lsqnonlin(rp_u,a2(:,1),[],[],opts);
rp_v = @(aOptim) A*aOptim-b2; 
a2(:,2) = lsqnonlin(rp_v,a2(:,2),[],[],opts);

if(length(a2) < 20)
    a2(20,:) = 0;
end

a_ = a2(:,1);

aMat1 = zeros(4,4,4);
aMat1(:,:,1) = [a_(11) a_(12) a_(15) a_(5);...
               a_(12) a_(13) a_(16) a_(6);...
               a_(15) a_(16) a_(18) a_(8);...
               a_(5)  a_(6)  a_(8)  a_(2)]; 
aMat1(:,:,2) = [a_(12) a_(13) a_(16) a_(6);...
               a_(13) a_(14) a_(17) a_(7);...
               a_(16) a_(17) a_(19) a_(9);...
               a_(6)  a_(7)  a_(9)  a_(3)];       
aMat1(:,:,3) = [a_(15) a_(16) a_(18) a_(8);...
               a_(16) a_(17) a_(19) a_(9);...
               a_(18) a_(19) a_(20) a_(10);...
               a_(8)  a_(9)  a_(10) a_(4)];
aMat1(:,:,4) = [a_(5) a_(6) a_(8) a_(2);...
               a_(6) a_(7) a_(9) a_(3);...
               a_(8) a_(9) a_(10) a_(4);...
               a_(2)  a_(3)  a_(4)  a_(1)];
           

a_ = a2(:,2);
aMat2 = zeros(4,4,4);
aMat2(:,:,1) = [a_(11) a_(12) a_(15) a_(5);...
               a_(12) a_(13) a_(16) a_(6);...
               a_(15) a_(16) a_(18) a_(8);...
               a_(5)  a_(6)  a_(8)  a_(2)]; 
aMat2(:,:,2) = [a_(12) a_(13) a_(16) a_(6);...
               a_(13) a_(14) a_(17) a_(7);...
               a_(16) a_(17) a_(19) a_(9);...
               a_(6)  a_(7)  a_(9)  a_(3)];       
aMat2(:,:,3) = [a_(15) a_(16) a_(18) a_(8);...
               a_(16) a_(17) a_(19) a_(9);...
               a_(18) a_(19) a_(20) a_(10);...
               a_(8)  a_(9)  a_(10) a_(4)];
aMat2(:,:,4) = [a_(5) a_(6) a_(8) a_(2);...
               a_(6) a_(7) a_(9) a_(3);...
               a_(8) a_(9) a_(10) a_(4);...
               a_(2)  a_(3)  a_(4)  a_(1)];
           
weightMat = ones(4,4,4);
weightMat(:,:,1) = [1   1/3 1/3 1/3;...
                    1/3 1/3 1/6 1/6;...
                    1/3 1/6 1/3 1/6;...
                    1/3 1/6 1/6 1/3];
weightMat(:,:,2) = [1/3 1/3 1/6 1/6;...
                    1/3   1 1/3 1/3;...
                    1/6 1/3 1/3 1/6;...
                    1/6 1/3 1/6 1/3]; 
weightMat(:,:,3) = [1/3 1/6 1/3 1/6;...
                    1/6 1/3 1/3 1/6;...
                    1/3 1/3   1 1/3;...
                    1/6 1/6 1/3 1/3]; 
weightMat(:,:,4) = [1/3 1/6 1/6 1/3;...
                    1/6 1/3 1/6 1/3;...
                    1/6 1/6 1/3 1/3;...
                    1/3 1/3 1/3   1]; 

% tensor for normalized data
N1_scale = aMat1.*weightMat;
N2_scale = aMat2.*weightMat;

N1 = zeros(4,4,4);
N2 = zeros(4,4,4);

% compute P'*TX
for i=1:4
    for j=1:4
        for k=1:4            
            
            % compute P'*TX
            for p=1:4
                for q=1:4
                    for r=1:4   
                        N1(q,p,r) = N1(q,p,r) + N1_scale(j,i,k) * TX(i,q) * TX(j,p) * TX(k,r);
                        N2(q,p,r) = N2(q,p,r) + N2_scale(j,i,k) * TX(i,q) * TX(j,p) * TX(k,r);
                    end
                end
            end
        end
    end
end

% get coefficients a back
aAlt = zeros(19,2);
aAlt(1,1) = N1(4,4,4);
aAlt(2,1) = N1(4,1,4)+N1(1,4,4)+N1(4,4,1);
aAlt(3,1) = N1(4,2,4)+N1(2,4,4)+N1(4,4,2);
aAlt(4,1) = N1(4,3,4)+N1(3,4,4)+N1(4,4,3);
aAlt(5,1) = N1(4,1,1)+N1(1,4,1)+N1(1,1,4);
aAlt(6,1) = N1(4,2,1)+N1(2,4,1)+N1(4,1,2)+N1(1,4,2)+N1(1,2,4)+N1(2,1,4);
aAlt(7,1) = N1(4,2,2)+N1(2,4,2)+N1(2,2,4);
aAlt(8,1) = N1(4,3,1)+N1(3,4,1)+N1(4,1,3)+N1(1,4,3)+N1(1,3,4)+N1(3,1,4);
aAlt(9,1) = N1(4,3,2)+N1(3,4,2)+N1(4,2,3)+N1(2,4,3)+N1(2,3,4)+N1(3,2,4);
aAlt(10,1) = N1(4,3,3)+N1(3,4,3)+N1(3,3,4);
aAlt(11,1) = N1(1,1,1);
aAlt(12,1) = N1(1,2,1)+N1(2,1,1)+N1(1,1,2);
aAlt(13,1) = N1(2,2,1)+N1(2,1,2)+N1(1,2,2);
aAlt(14,1) = N1(2,2,2);
aAlt(15,1) = N1(1,3,1)+N1(3,1,1)+N1(1,1,3);
aAlt(16,1) = N1(3,2,1)+N1(2,3,1)+N1(3,1,2)+N1(1,3,2)+N1(1,2,3)+N1(2,1,3);
aAlt(17,1) = N1(2,3,2)+N1(3,2,2)+N1(2,2,3);
aAlt(18,1) = N1(3,3,1)+N1(1,3,3)+N1(3,1,3);
aAlt(19,1) = N1(3,3,2)+N1(2,3,3)+N1(3,2,3);
aAlt(1,2) = N2(4,4,4);
aAlt(2,2) = N2(4,1,4)+N2(1,4,4)+N2(4,4,1);
aAlt(3,2) = N2(4,2,4)+N2(2,4,4)+N2(4,4,2);
aAlt(4,2) = N2(4,3,4)+N2(3,4,4)+N2(4,4,3);
aAlt(5,2) = N2(4,1,1)+N2(1,4,1)+N2(1,1,4);
aAlt(6,2) = N2(4,2,1)+N2(2,4,1)+N2(4,1,2)+N2(1,4,2)+N2(1,2,4)+N2(2,1,4);
aAlt(7,2) = N2(4,2,2)+N2(2,4,2)+N2(2,2,4);
aAlt(8,2) = N2(4,3,1)+N2(3,4,1)+N2(4,1,3)+N2(1,4,3)+N2(1,3,4)+N2(3,1,4);
aAlt(9,2) = N2(4,3,2)+N2(3,4,2)+N2(4,2,3)+N2(2,4,3)+N2(2,3,4)+N2(3,2,4);
aAlt(10,2) = N2(4,3,3)+N2(3,4,3)+N2(3,3,4);
aAlt(11,2) = N2(1,1,1);
aAlt(12,2) = N2(1,2,1)+N2(2,1,1)+N2(1,1,2);
aAlt(13,2) = N2(2,2,1)+N2(2,1,2)+N2(1,2,2);
aAlt(14,2) = N2(2,2,2);
aAlt(15,2) = N2(1,3,1)+N2(3,1,1)+N2(1,1,3);
aAlt(16,2) = N2(3,2,1)+N2(2,3,1)+N2(3,1,2)+N2(1,3,2)+N2(1,2,3)+N2(2,1,3);
aAlt(17,2) = N2(2,3,2)+N2(3,2,2)+N2(2,2,3);
aAlt(18,2) = N2(3,3,1)+N2(1,3,3)+N2(3,1,3);
aAlt(19,2) = N2(3,3,2)+N2(2,3,3)+N2(3,2,3);

%transform back inv(Tx)*P''
a = zeros(19,2);
Tx_inv = inv(Tx);
a(:,1) = Tx_inv(1,1).*aAlt(:,1); 
a(1,1) = a(1,1)+Tx_inv(1,3);
a(:,2) = Tx_inv(2,2).*aAlt(:,2);
a(1,2) = a(1,2)+Tx_inv(2,3);

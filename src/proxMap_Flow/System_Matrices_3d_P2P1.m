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

function [Nabla, Div, Lpl] = System_Matrices_3d_P2P1( L, M, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDER IS Z,Y,X !!

%coefficients
f = 1; %3;
% phi = lin (1-x),x
% gsi = quad 2(x-0.5)(x-1),4x(1-x),2x(x-0.5)
% gsi' = (4x-3),(4-8x),(4x-1)
% integral phi_i*gsi_j'
a(1,1) = f*-5/6;
a(1,2) = f*2/3;
a(1,3) = f*1/6;
a(2,1) = f*-1/6;
a(2,2) = f*-2/3;
a(2,3) = f*5/6;

% integral phi_i*gsi_j
b(1,1) = f*1/6;
b(1,2) = f*1/3;
b(1,3) = 0;
b(2,1) = 0;
b(2,2) = f*1/3;
b(2,3) = f*1/6;

A_u = zeros(8,27);
A_v = zeros(8,27);
A_w = zeros(8,27);

%actually [0,0.5,1]
for r=1:2
    for s=1:2
        for t=1:2
            for i=1:3
                for j=1:3
                    for k=1:3
                        % convert 3d to 1d idx in correct order
                        idx1 = (r-1)*2*2+(s-1)*2+t;
                        idx2 = (i-1)*3*3+(j-1)*3+k;
                        A_u(idx1,idx2) = a(r,i)*b(s,j)*b(t,k);
                        A_v(idx1,idx2) = b(r,i)*a(s,j)*b(t,k);
                        A_w(idx1,idx2) = b(r,i)*b(s,j)*a(t,k);
                    end
                end
            end
        end
    end
end

% here I need to group 4er pairs per voxel .. 
% so it is .. id 1:N-1 and 2:N

% still x and x+1 together but with y and z via 00, 01, 10, 11. displaced
% 

Lp = floor((L+1)/2);
Mp = floor((M+1)/2);
Np = floor((N+1)/2);

idxMatrix = reshape(1:N*M*L, L,M,N); % component order 
idxMatrix_p = reshape(1:Np*Mp*Lp, Lp,Mp,Np);

% for each vertex location build row vector of 5x5x5=125
            % entries from M
            % each row in M corresponds to a 3x3x3 matrix that is put ont o
            % the current vertex (borders need to be handled accordingly)

%loop over vertices
I = zeros(1,Np*Mp*Lp*8*27);
J = zeros(1,Np*Mp*Lp*8*27);
V = zeros(1,Np*Mp*Lp*8*27);
idx=1;
for r=1:Np
    for s=1:Mp
        for t=1:Lp
            % looping over all 8x27 elements of M (2x2x2, 3x3x3) and setting
            % them at the correct indices... (+border handling)
            
            % per vertex loop over 2x2x2 locations where the 3x3x3 row
            % vector of M is put (in total 8x27, which is the matrix A)
            idxRowLpl = idxMatrix_p(t,s,r);
            
            
            for rr=0:-1:-1
                for sr=0:-1:-1
                    for tr=0:-1:-1
                        if r+rr < 1 || r+rr+1 > Np
                            continue
                        end
                        if s+sr < 1 || s+sr+1 > Mp
                            continue
                        end
                        if t+tr < 1 || t+tr+1 > Lp
                            continue
                        end
                        
                        idxRowM = -rr*2*2-sr*2-tr+1;

                        I(idx:idx+27*3-1) = idxRowLpl;
                        Jmat = idxMatrix((t+tr-1)*2+1:(t+tr-1)*2+3,(s+sr-1)*2+1:(s+sr-1)*2+3,(r+rr-1)*2+1:(r+rr-1)*2+3);
                        Jmat = Jmat(:);
                        J(idx:idx+26) = Jmat;
                        V(idx:idx+26) = A_u(idxRowM,:);
                        idx=idx+27;
                        J(idx:idx+26) = Jmat+N*M*L;
                        V(idx:idx+26) = A_v(idxRowM,:);
                        idx=idx+27;
                        J(idx:idx+26) = Jmat+2*N*M*L;
                        V(idx:idx+26) = A_w(idxRowM,:);
                        idx=idx+27;

                    end
                end
            end
        end
    end
end

I(idx:end) = [];
J(idx:end) = [];
V(idx:end) = [];

Div = sparse(I,J,V,Np*Mp*Lp,3*N*M*L);
Nabla = Div';
Lpl   = Div*Nabla;



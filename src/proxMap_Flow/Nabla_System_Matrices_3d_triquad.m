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

function [Nabla] = Nabla_System_Matrices_3d_triquad( L, M, N)

%coefficients
f = 15;
% integral phi_i*phi_j'
a(1,1) = f*-1/2; %phi0'phi0
a(1,2) = f*2/3; %phi1'phi0
a(1,3) = f*-1/6; %phi2'phi0
a(2,1) = f*-2/3;
a(2,2) = 0;
a(2,3) = f*2/3;
a(3,1) = f*1/6;
a(3,2) = f*-2/3;
a(3,3) = f*1/2; 

% integral phi_i*phi_j
b(1,1) = f*2/15;
b(1,2) = f*1/15;
b(1,3) = f*-1/30;
b(2,1) = b(1,2);
b(2,2) = f*8/15;
b(2,3) = f*1/15;
b(3,1) = b(1,3);
b(3,2) = b(2,3);
b(3,3) = f*2/15;


A_u = zeros(27,27);
A_v = zeros(27,27);
A_w = zeros(27,27);

%actually [0,0.5,1]
for r=1:3
    for s=1:3
        for t=1:3
            for i=1:3
                for j=1:3
                    for k=1:3
                        % convert 3d to 1d idx in correct order
                        idx1 = (r-1)*3*3+(s-1)*3+t;
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

idxMatrix = reshape(1:N*M*L, L,M,N); % component order 

% for each vertex location build row vector of 5x5x5=125
            % entries from M
            % each row in M corresponds to a 3x3x3 matrix that is put ont o
            % the current vertex (borders need to be handled accordingly)

Lp = floor((L+1)/2);
Mp = floor((M+1)/2);
Np = floor((N+1)/2);
            
%loop over vertices
I = zeros(1,(Np-1)*(Mp-1)*(Lp-1)*27*27);
J = zeros(1,(Np-1)*(Mp-1)*(Lp-1)*27*27);
V = zeros(1,(Np-1)*(Mp-1)*(Lp-1)*27*27);
idx=1;
for ip=0:Np-2
    for jp=0:Mp-2
        for kp=0:Lp-2

            Jmat = idxMatrix(kp*2+1:kp*2+3,jp*2+1:jp*2+3,ip*2+1:ip*2+3);
            Jmat = Jmat(:);
            
            % iterate over 3x3x3 voxel elements
            for iv=1:3
                for jv=1:3
                    for kv=1:3
                        
                        i=ip*2+iv;
                        j=jp*2+jv;
                        k=kp*2+kv;
                        
                        idxRowLpl = idxMatrix(k,j,i);
                        idxRowM = (iv-1)*3*3+(jv-1)*3+kv;

                        I(idx:idx+26) = idxRowLpl;
                        
                        J(idx:idx+26) = Jmat;
                        V(idx:idx+26) = A_u(idxRowM,:);
                        idx=idx+27;
                        I(idx:idx+26) = idxRowLpl+N*M*L;
                        J(idx:idx+26) = Jmat;
                        V(idx:idx+26) = A_v(idxRowM,:);
                        idx=idx+27;
                        I(idx:idx+26) = idxRowLpl+2*N*M*L;
                        J(idx:idx+26) = Jmat;
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
Nabla = sparse(I,J,V,3*N*M*L,N*M*L);
Nabla = Nabla/f.^3;
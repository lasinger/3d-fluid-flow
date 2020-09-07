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

function [Nabla] = Nabla_System_Matrices_3d_trilin( L, M, N)

%coefficients
f = 1;
% integral phi_i*phi_j' 
% lin: phi0=(1-x), phi1=x, phi0'=-1, phi1'=1
a(1,1) = -0.5; %phi0'phi0
a(1,2) = 0.5; %phi1'phi'0
a(2,1) = -0.5; %phi0'phi1
a(2,2) = 0.5; %phi1'phi1

%lin: phi_0 = 1-x. phi_1 = x
% integral phi_i*phi_j
b(1,1) = 0.25;
b(1,2) = 0.25;
b(2,1) = 0.25;
b(2,2) = 0.25;


%b = b';

A_u = zeros(8,8);
A_v = zeros(8,8);
A_w = zeros(8,8);

%actually [0,0.5,1]
for r=1:2
    for s=1:2
        for t=1:2
            for i=1:2
                for j=1:2
                    for k=1:2
                        % convert 3d to 1d idx in correct order
                        idx1 = (r-1)*2*2+(s-1)*2+t;
                        idx2 = (i-1)*2*2+(j-1)*2+k;
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

Lp = L; 
Mp = M; 
Np = N; 
            
%loop over vertices
I = zeros(1,(N-1)*(M-1)*(L-1)*8*8);
J = zeros(1,(N-1)*(M-1)*(L-1)*8*8);
V = zeros(1,(N-1)*(M-1)*(L-1)*8*8);
idx=1;
for ip=0:Np-2
    for jp=0:Mp-2
        for kp=0:Lp-2

            Jmat = idxMatrix(kp+1:kp+2,jp+1:jp+2,ip+1:ip+2);
            Jmat = Jmat(:);
            
            % iterate over 3x3x3 voxel elements
            for iv=1:2
                for jv=1:2
                    for kv=1:2
                        
                        i=ip+iv;
                        j=jp+jv;
                        k=kp+kv;
                        
                        idxRowLpl = idxMatrix(k,j,i);
                        idxRowM = (iv-1)*2*2+(jv-1)*2+kv;

                        I(idx:idx+7) = idxRowLpl;
                       
                        J(idx:idx+7) = Jmat;
                        V(idx:idx+7) = A_u(idxRowM,:);
                        idx=idx+8;
                        I(idx:idx+7) = idxRowLpl+N*M*L;
                        J(idx:idx+7) = Jmat;
                        V(idx:idx+7) = A_v(idxRowM,:);
                        idx=idx+8;
                        I(idx:idx+7) = idxRowLpl+2*N*M*L;
                        J(idx:idx+7) = Jmat;
                        V(idx:idx+7) = A_w(idxRowM,:);
                        idx=idx+8;
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
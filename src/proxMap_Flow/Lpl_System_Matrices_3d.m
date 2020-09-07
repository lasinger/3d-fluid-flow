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
% Author: Katrin Lasinger, Christoph Vogel

function [Lpl] = Lpl_System_Matrices_3d( L, M, N, useX, useY, useZ) % flipped N L order here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M=10;
%L=10;
%N=10; % z direction

if ~exist('useY','var')
    useY = 1;
end
if ~exist('useZ','var')
    useZ = 1;
end
if ~exist('useX','var')
    useX = 1;
end

idxMatrix = reshape(1:N*M*L, L,M,N); % component order: 

%dx 
I=speye( N*M*L, N*M*L );

if useX
  ijX=idxMatrix([2:L,L], :, :);
  E=sparse( 1:N*M*L, ijX(:), ones(N*M*L,1), N*M*L, N*M*L );
  dX = I-E;
  clear('E');
  clear('ijX');  
else
  dX = sparse( N*M*L, N*M*L );
end
Lpl   = dX'*dX;
clear('dX');

if useY
  ijY=idxMatrix(:, [2:M,M], :);
  S=sparse( 1:N*M*L, ijY(:), ones(N*M*L,1), N*M*L, N*M*L );
  dY = I-S;
  clear('S');  
  clear('ijY');
else
  dY = sparse( N*M*L, N*M*L );
end
Lpl   = Lpl + dY'*dY;
clear('dY');

if useZ
  ijZ=idxMatrix(:, :, [2:N,N]);
  D=sparse( 1:N*M*L, ijZ(:), ones(N*M*L,1), N*M*L, N*M*L );    
  dZ = I-D;
  clear('D');
  clear('ijZ');
else
  dZ = sparse( N*M*L, N*M*L );
end
clear('I');
clear('idxMatrix');
Lpl   = Lpl + dZ'*dZ;
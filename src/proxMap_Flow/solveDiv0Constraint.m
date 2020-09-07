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

% input: 
% iterations are the number of steps we take here.  
% since we are doing this all the time and kick start from the old solution
% this could be just 10 iterations; with a large tolerance .. 
% size of volume in x:L, y:M, z:N 
% solution for the flow of last iteration + gradient : v_k + L * \nabla_H|v_k
% tol : up to which tolerance the solution is computed .. 1e-4 or 1e-3
% might be sufficient (not sure .. )
% old solution for the pressure term (start with 0)
%
% Lu: Lipshitz constant used (only need to compute v_new from new pressure solution)
%
% Beware: the Lipshitz test has still to be done afterwards.
% one might want to split the problem 
% additional kinetic penalty as in Capture to Simulation (Gregson) does
% this
function [pressure_new,v_new] = solveDiv0Constraint( iterations, L,M,N, flow_k_bar, pressure_k, Lu, kinetic, tol, Nabla, Div, Lpl )

if L > 180
  truncPCG = 1e-0;iterations = int32(iterations * 2.0); % less mem need more iterations but also faster .. 
else
  truncPCG = 1e-1;iterations = int32(iterations * 1.3);
end
%truncPCG = 1e-2;  % slower but more accurate .. and more mem -- yet still
%should be prefered as long as memory is there!

if ~exist('tol', 'var')
  tol = 1e-3;
end

if ~exist('pressure_k', 'var')
  pressure_k = zeros(L*M*N,1);
end

if ~exist('kinetic', 'var')
  kinetic = 0;
end

Divflow_k_bar = Lu*Div*flow_k_bar;

Li = ichol(Lpl+(1e-7)*speye(size(Lpl)), struct('type','ict','droptol', truncPCG));% faster .. needs to multiply with this matrix per iteration

[pressure_new,flag,relres,iter,resvec] = pcg( Lpl+(1e-7)*speye(size(Lpl)), Divflow_k_bar, tol, iterations, Li, Li', pressure_k );

clear('Li');

v_new   = (Lu/(kinetic+Lu)) * flow_k_bar - (1./(kinetic+Lu)) * (Nabla * pressure_new); % no need to store Nabla now

return

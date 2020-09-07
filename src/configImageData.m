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

function [params] = configImageData(varargin)

%algorithm params
params.sigma = 1; % sigma rendering Gaussin blobs
params.intThresh = 0.02; % intensity threshold for when a particle vanishes
params.my = 0.0001; % i think needs to be smaller for 0norm than for 1norm (0.0001 seems to be good for 0norm - but might depend)
params.lambda = 0.000004;
params.divMethod = 1; %0..soft constraint (depending on alpha), 1..hard constraint
params.femSpace = -1; % 0...P1P0, 1...P2P1, 2...P2P0, -1...P1P0 (different Nabla, best results for flow, artifacts on pressure, default)
params.sparsityNorm = 0; %1,0 norm or no sparsity (just making sure that all are positive)
params.its = 40;
params.kinetic = 0;
params.pcgIts = 20;
params.lipThresh = 1e-6;
params.triangErrorList = [0.8,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.5,1.5,1.5,2.0];
params.stepsizeInt = 10;
params.pyramid_factor = 0.94;
params.pyramid_levels = 10;
params.stopThreshFine = 0.01;
params.alpha = 40;
params.lambdaPressure = 0;

% îo
params.calibPtsPath = '../data/calibpoints.txt';
params.imgPath = '../data/ppp0175/';
params.parPathResults = '../results/';
params.imageFormat = 'partImg_%04d_%d.tif'; %indices: 1st: time (dependent on params.timesteps), 2nd: cam (starting with 0)
params.timesteps = [1,2]; % ids of used time steps (currently 2 or 3 time steps supported, refence idx = ceil(num_ts/2)

% data
params.methodName = 'testData';
params.factor_vox_mm = 20;
params.width_mm = 51.2;
params.height_mm = 25.6;
params.depth_mm = 17.6;
params.box_offset_x = params.width_mm/2;
params.box_offset_y = params.height_mm/2;
params.box_offset_z = params.depth_mm/2;
params.isPolynomialCamModel = 0;

% others
params.displayFigures = 1;

% Loading optional arguments
if length(varargin) == 1
    varargin = varargin{1};
end
for i=1:2:length(varargin)-1
    params.(varargin{i}) = varargin{i+1};
end

end
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

function[] = sceneFlow_imageInput(varargin)

params = configImageData(varargin);

addpath('proxMap_Flow/');
addpath('polynomialCamera/');

% load calibration data
formatSpec = '%f %f %f %f %f %f %f %f %f'; % x y z u0 v0 u1 v1 u2 v2 u3 v3
A = readmatrix(params.calibPtsPath,'NumHeaderLines',1);
worldPts = [A(:,1:3)];

numCam = (size(A,2)-3)/2;
for cam=1:numCam
    imagePts{cam} = A(:,4+(cam-1)*2:5+(cam-1)*2)';
end

ts = params.timesteps;
num_ts = length(ts);    

ref_idx = ceil(num_ts/2);

% conversion mm to voxel (voxel size approx. pixel size, starting at 1)
boxDimensions = struct('factor_vox_mm',params.factor_vox_mm,'width_mm',params.width_mm,'height_mm',params.height_mm,'depth_mm',params.depth_mm);
boxOffset =[params.box_offset_x,params.box_offset_y,params.box_offset_z];
worldPtsVox = worldPts + repmat(boxOffset,length(worldPts),1);% mm to vox
worldPtsVox = worldPtsVox * boxDimensions.factor_vox_mm+1;

N=boxDimensions.factor_vox_mm*boxDimensions.width_mm;
M=boxDimensions.factor_vox_mm*boxDimensions.height_mm;
L=boxDimensions.factor_vox_mm*boxDimensions.depth_mm;

numPart = length(worldPtsVox);
imagePtsMat = zeros(numPart,2,numCam);

% get camera model
P = cell(numCam,1);
for cam=1:numCam
    if params.isPolynomialCamModel
        P{cam} = polynomialCameraCalib(worldPtsVox,imagePts{cam}');
    else
        P{cam} = dlt(imagePts{cam},worldPtsVox);
    end
end

% read image data and normalize
maxIntensity = 0;
imgs = cell(1,num_ts);
for cam=1:numCam
    
    for t=1:num_ts
        imgFilename = append(params.imgPath,sprintf(params.imageFormat,ts(t),cam-1));
        imgs{t}{cam} = double(imread(imgFilename));
        maxIntensity = max(maxIntensity,max(imgs{t}{cam}(:)));
    end
end
for cam=1:numCam
    for t=1:num_ts
        imgs{t}{cam} = min(imgs{t}{cam}/maxIntensity+0.001,1);
    end
end

N2d = size(imgs{1}{1},2);
M2d = size(imgs{1}{1},1);
estPart = [];

% write config data to txt file
formatOut = 'yymmdd_HHMM_SSFFF';
nowStr = datestr(now,formatOut);

fileID = fopen([params.parPathResults 'flowResult_',params.methodName,'_',nowStr,'.txt'],'w');
fprintf(fileID, 'Parameters:\r\n');
fprintf(fileID, 'Img Filename: %s\r\n',imgFilename);
fprintf(fileID, 'FEM space: %.6f\r\n',params.femSpace);
fprintf(fileID, 'Sigma: %.6f\r\n',params.sigma);
fprintf(fileID, 'Intensity thresh: %.6f\r\n',params.intThresh);
fprintf(fileID, 'My (sparsity): %.6f\r\n',params.my);
fprintf(fileID, 'Lambda (quad. reg): %.10f\r\n',params.lambda);
fprintf(fileID, 'Alpha: %.6f\r\n',params.alpha);
fprintf(fileID, 'Div Method: %d\r\n',params.divMethod);
fprintf(fileID, 'Sparsity norm: %d\r\n',params.sparsityNorm);
fprintf(fileID, 'Iterations: %d\r\n',params.its);
fprintf(fileID, 'PCG Iterations: %d\r\n',params.pcgIts);
fprintf(fileID, 'Kinetic: %.2f\r\n',params.kinetic);
fprintf(fileID, 'Lip thresh: %f\r\n',params.lipThresh);
triangErrorListOutput=sprintf('%.2f ', params.triangErrorList);
fprintf(fileID, 'Triang error list: [%s]\r\n',triangErrorListOutput);
fprintf(fileID, 'stop thresh fine: %f\r\n',params.stopThreshFine);
fprintf(fileID, 'Stepsize: %.2f\r\n',params.stepsizeInt);
fprintf(fileID, 'Pyr factor: %.2f\r\n',params.pyramid_factor);
fprintf(fileID, 'Pyr levels: %d\r\n',params.pyramid_levels);
fprintf(fileID, 'N: %d\r\n',N);
fprintf(fileID, 'M: %d\r\n',M);
fprintf(fileID, 'L: %d\r\n',L);
fprintf(fileID, 'polynomialModel: %f\r\n',params.isPolynomialCamModel);
fprintf(fileID, 'ts: ');
for t=1:num_ts
    fprintf(fileID, '%d,',ts(t));
end
fprintf(fileID, '\r\n');

% flow estimation and particle reconstruction
totalTime = tic;
[estPart, estU_grid, estU_part, estPressure, Nds, Mds, Lds, cost] = sceneParticleFlow(imgs, N, M, L, num_ts, ref_idx, N2d, M2d, P, params, fileID, estPart);
elapsedTime = toc(totalTime)


% convert particles & flow from voxel to mm

estPart_ref(1:3,:) = (estPart(1:3,:)-1)/boxDimensions.factor_vox_mm-repmat(boxOffset',1,length(estPart));
estU_part_ref = estU_part{ref_idx}/boxDimensions.factor_vox_mm;

% stepsize compared to voxel dim
stepsize = [ (N-1)/(Nds-1); (M-1)/(Mds-1); (L-1)/(Lds-1) ];

[z_mat, y_mat, x_mat] = ndgrid((0:stepsize(3):(L-1))/(L-1), (0:stepsize(2):(M-1))/(M-1), (0:stepsize(1):(N-1))/(N-1));
    
est_flow = struct;
est_flow.U = reshape(estU_grid{ref_idx}(1,:),Lds,Mds,Nds)/boxDimensions.factor_vox_mm;
est_flow.V = reshape(estU_grid{ref_idx}(2,:),Lds,Mds,Nds)/boxDimensions.factor_vox_mm;
est_flow.W = reshape(estU_grid{ref_idx}(3,:),Lds,Mds,Nds)/boxDimensions.factor_vox_mm;
est_flow.X = x_mat*boxDimensions.width_mm-boxOffset(1);
est_flow.Y = x_mat*boxDimensions.height_mm-boxOffset(2);
est_flow.Z = x_mat*boxDimensions.depth_mm-boxOffset(3);

% write output to mat file
save([params.parPathResults 'flowResult_',params.methodName,'_',nowStr,'.mat'], '-mat', '-v7.3','est_flow','estU_part_ref','estPart_ref','Lds','Mds','Nds','L','M','N','stepsize','ts','ref_idx','P','boxDimensions','boxOffset','params');


% output visualizations
zSlice = 78;
aepeF = cost;

% single slice
fig_h=figure;ax = gca;imshow(imresize(squeeze(est_flow.U(round(zSlice/stepsize(1)),:,:)),4,'nearest'),[]); colormap(ax,jet(255)); colorbar;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(fig_h,[params.parPathResults 'flowResult_',params.methodName,'_',nowStr,'_flowDS.png']);
close(fig_h);

% residual image
for cam=1:numCam
    part2d{cam} = proj2d(estPart(1:3,:),P{cam});
end
estImgIt = render2dallcammex(part2d,estPart(4,:),size(estPart,2),params.sigma,N2d,M2d );
residual01 = imgs{ref_idx}{1}-estImgIt{1};
fig_h=figure;ax = gca;imshow(residual01,[]);
saveas(fig_h,[params.parPathResults 'flowResult_',params.methodName,'_',nowStr,'_residual.png']);
close(fig_h);


[~, Div_incomp, ~] = System_Matrices_3d_P1P0( Lds, Mds, Nds);
flow_k = reshape(estU_grid{ref_idx}',3*Nds*Mds*Lds,1);
div=Div_incomp*flow_k;
aad_ds = mean(abs(div));
clear('Div_incomp');
clear('div');
clear('flow_k');

fprintf(fileID, '\r\n\r\nElapsed time: %f\r\n', elapsedTime);
fprintf(fileID, '\r\nFinal results:\n\r');
fprintf(fileID, 'Num. reconstructed particles:%d\r\n',size(estPart,2));
fprintf(fileID, 'Cost:%f\r\n',cost);
fprintf(fileID, 'Divergence:%f\r\n',aad_ds);

fclose(fileID);

end
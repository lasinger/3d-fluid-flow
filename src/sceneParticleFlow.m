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

function [estPart, estU_grid, estU_part, estPressure, Nds, Mds, Lds, cost] = sceneParticleFlow(imgs, N, M, L, num_ts, ref_idx, N2d, M2d, P, params, fileID, estPart, gtData, part)

writeToFile = false;
if exist('fileID','var')
    writeToFile = true;
end
if ~exist('estPart','var')
    estPart = [];
end
if ~exist('gtData','var')
    gtData = false;
end

numCam = length(P);

%flow grid
Nds = ceil(N/params.stepsizeInt);
Mds = ceil(M/params.stepsizeInt);
Lds = ceil(L/params.stepsizeInt);

%pressure grid
Nds_p = Nds-1;
Mds_p = Mds-1;
Lds_p = Lds-1;

if params.femSpace == 1 %P2P1
    % need to be odd numbers
    Nds = roundOdd(N/params.stepsizeInt);
    Mds = roundOdd(M/params.stepsizeInt);
    Lds = roundOdd(L/params.stepsizeInt);

    Nds_p = (Nds+1)/2;
    Mds_p = (Mds+1)/2;
    Lds_p = (Lds+1)/2;
elseif params.femSpace == 2 %P2P0
    % need to be odd numbers
    Nds = roundOdd(N/params.stepsizeInt);
    Mds = roundOdd(M/params.stepsizeInt);
    Lds = roundOdd(L/params.stepsizeInt);

    Nds_p = (Nds-1)/2;
    Mds_p = (Mds-1)/2;
    Lds_p = (Lds-1)/2;
end

if params.femSpace <= 0
    params.lambdaPressure = 0;
end

stepsize = [ (N-1)/(Nds-1); (M-1)/(Mds-1); (L-1)/(Lds-1) ];

interpMethod = 'cubic';

stepsizeLip = 2;
Lc = 1;
Lp = 1;
Lu = 1;

width_Pyramid{1} = Nds;
height_Pyramid{1} = Mds;
depth_Pyramid{1} = Lds;
sigma_Pyramid{1} = params.sigma;


for i=2:params.pyramid_levels
    if params.femSpace <= 0
        width_Pyramid{i}  = round(params.pyramid_factor*width_Pyramid{i-1});
        height_Pyramid{i} = round(params.pyramid_factor*height_Pyramid{i-1});
        depth_Pyramid{i} = round(params.pyramid_factor*depth_Pyramid{i-1});
    else
        width_Pyramid{i}  = roundOdd(params.pyramid_factor*width_Pyramid{i-1});
        height_Pyramid{i} = roundOdd(params.pyramid_factor*height_Pyramid{i-1});
        depth_Pyramid{i} = roundOdd(params.pyramid_factor*depth_Pyramid{i-1});
    end
    sigma_Pyramid{i} = sigma_Pyramid{i-1}/params.pyramid_factor;
end

% init on highest pyr lvl
estU_grid = cell(1,num_ts-1);
estPressure = cell(1,num_ts-1);
estU_part = cell(1,num_ts-1);
pressure = [];
for i=1:num_ts-1
    estU_grid{i} = zeros(3,Nds*Mds*Lds);
    
    if params.divMethod == 0 % soft constraint
        estPressure{i} = 0;
    else
        estPressure{i} = zeros(Nds_p*Mds_p*Lds_p,1); 
    end
end

triangIt = 0;
numTriangIts = length(params.triangErrorList);

% triangulation loop
for triangError=params.triangErrorList
    tic
    triangIt=triangIt+1;
    stopThresh = params.stopThreshFine*(numTriangIts+2-triangIt)/2; %0.5*(1/triangIt);
    
    disp(['triang error: ' num2str(triangError) ', triang it: ' num2str(triangIt)]);
    if ~isempty(estPart)
        for cam=1:numCam
            part2d{cam} = proj2d(estPart(1:3,:),P{cam});
        end
        estImgIt0 = render2dallcammex(part2d,estPart(4,:),size(estPart,2),params.sigma,N2d,M2d );
        for cam=1:numCam
            residual{cam} = max(imgs{ref_idx}{cam} - estImgIt0{cam},0);
        end
    else
        residual = imgs{ref_idx};
    end

    [reproj, residual, estPartNew] = triangulate_part(residual,P,triangError,N,M,L,1,params.intThresh,params.displayFigures);

    disp(['num part prev:' num2str(size(estPart,2))]);
    disp(['added particles:' num2str(size(estPartNew,2))]);
    estPart = [estPart,estPartNew];
    
    
    numpart = size(estPart,2);
    disp(['num part:' num2str(size(estPart,2))]);
    if params.displayFigures
        figure;
    end
    flowSigma = sigma_Pyramid{params.pyramid_levels};
    cost = 0;
    
    for level = params.pyramid_levels:-1:1
        sigmaPyr = sigma_Pyramid{level};

        Nds = width_Pyramid{level};
        Mds = height_Pyramid{level};
        Lds = depth_Pyramid{level};

        stepsizePyr = [ (N-1)/(Nds-1); (M-1)/(Mds-1); (L-1)/(Lds-1) ];

        if level == params.pyramid_levels

            %downsample to lowest pyr lvl
            [z_mat_orig, y_mat_orig, x_mat_orig] = ndgrid((0:stepsize(3):(L-1))/(L-1), (0:stepsize(2):(M-1))/(M-1), (0:stepsize(1):(N-1))/(N-1));
            [z_mat_interp, y_mat_interp, x_mat_interp] = ndgrid((0:stepsizePyr(3):(L-1))/(L-1), (0:stepsizePyr(2):(M-1))/(M-1), (0:stepsizePyr(1):(N-1))/(N-1));
            
            for i=1:num_ts-1
                u = reshape(estU_grid{i}(1,:),depth_Pyramid{1},height_Pyramid{1},width_Pyramid{1});
                v = reshape(estU_grid{i}(2,:),depth_Pyramid{1},height_Pyramid{1},width_Pyramid{1});
                w = reshape(estU_grid{i}(3,:),depth_Pyramid{1},height_Pyramid{1},width_Pyramid{1});
                u = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,u,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_u);
                v = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,v,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_v);
                w = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,w,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_w);
                estU_grid{i} = [u(:)';v(:)';w(:)'];
            end
            if params.divMethod
                if params.femSpace <= 0
                    z_mat_orig = z_mat_orig(1:end-1,1:end-1,1:end-1);
                    y_mat_orig = y_mat_orig(1:end-1,1:end-1,1:end-1);
                    x_mat_orig = x_mat_orig(1:end-1,1:end-1,1:end-1);
                    z_mat_interp = z_mat_interp(1:end-1,1:end-1,1:end-1);
                    y_mat_interp = y_mat_interp(1:end-1,1:end-1,1:end-1);
                    x_mat_interp = x_mat_interp(1:end-1,1:end-1,1:end-1);
                elseif params.femSpace == 2
                    z_mat_orig = z_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    y_mat_orig = y_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    x_mat_orig = x_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    z_mat_interp = z_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                    y_mat_interp = y_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                    x_mat_interp = x_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                else
                    z_mat_orig = z_mat_orig(1:2:end,1:2:end,1:2:end);
                    y_mat_orig = y_mat_orig(1:2:end,1:2:end,1:2:end);
                    x_mat_orig = x_mat_orig(1:2:end,1:2:end,1:2:end);
                    z_mat_interp = z_mat_interp(1:2:end,1:2:end,1:2:end);
                    y_mat_interp = y_mat_interp(1:2:end,1:2:end,1:2:end);
                    x_mat_interp = x_mat_interp(1:2:end,1:2:end,1:2:end);
                end
                for i=1:num_ts-1
                    if params.femSpace <= 0
                        pressure = reshape(estPressure{i},depth_Pyramid{1}-1,height_Pyramid{1}-1,width_Pyramid{1}-1);
                    elseif params.femSpace == 2
                        pressure = reshape(estPressure{i},(depth_Pyramid{1}-1)/2,(height_Pyramid{1}-1)/2,(width_Pyramid{1}-1)/2);
                    else
                        pressure = reshape(estPressure{i},(depth_Pyramid{1}+1)/2,(height_Pyramid{1}+1)/2,(width_Pyramid{1}+1)/2);
                    end
                    pressure = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,pressure,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_w);
                    estPressure{i} = pressure(:);
                end
            end
        else

            [z_mat_orig, y_mat_orig, x_mat_orig] = ndgrid((0:stepsizePyr_prev(3):(L-1))/(L-1), (0:stepsizePyr_prev(2):(M-1))/(M-1), (0:stepsizePyr_prev(1):(N-1))/(N-1));
            [z_mat_interp, y_mat_interp, x_mat_interp] = ndgrid((0:stepsizePyr(3):(L-1))/(L-1), (0:stepsizePyr(2):(M-1))/(M-1), (0:stepsizePyr(1):(N-1))/(N-1));

            for i=1:num_ts-1
                u = reshape(estU_grid{i}(1,:),Lprev,Mprev,Nprev);
                v = reshape(estU_grid{i}(2,:),Lprev,Mprev,Nprev);
                w = reshape(estU_grid{i}(3,:),Lprev,Mprev,Nprev);
                u = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,u,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_u);
                v = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,v,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_v);
                w = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,w,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_w);
                estU_grid{i} = [u(:)';v(:)';w(:)'];
            end
            if params.divMethod
                if params.femSpace <= 0
                    z_mat_orig = z_mat_orig(1:end-1,1:end-1,1:end-1);
                    y_mat_orig = y_mat_orig(1:end-1,1:end-1,1:end-1);
                    x_mat_orig = x_mat_orig(1:end-1,1:end-1,1:end-1);
                    z_mat_interp = z_mat_interp(1:end-1,1:end-1,1:end-1);
                    y_mat_interp = y_mat_interp(1:end-1,1:end-1,1:end-1);
                    x_mat_interp = x_mat_interp(1:end-1,1:end-1,1:end-1);
                elseif params.femSpace == 2
                    z_mat_orig = z_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    y_mat_orig = y_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    x_mat_orig = x_mat_orig(1:2:end-1,1:2:end-1,1:2:end-1);
                    z_mat_interp = z_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                    y_mat_interp = y_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                    x_mat_interp = x_mat_interp(1:2:end-1,1:2:end-1,1:2:end-1);
                else
                    z_mat_orig = z_mat_orig(1:2:end,1:2:end,1:2:end);
                    y_mat_orig = y_mat_orig(1:2:end,1:2:end,1:2:end);
                    x_mat_orig = x_mat_orig(1:2:end,1:2:end,1:2:end);
                    z_mat_interp = z_mat_interp(1:2:end,1:2:end,1:2:end);
                    y_mat_interp = y_mat_interp(1:2:end,1:2:end,1:2:end);
                    x_mat_interp = x_mat_interp(1:2:end,1:2:end,1:2:end);
                end
                for i=1:num_ts-1
                    if params.femSpace <= 0
                        pressure = reshape(estPressure{i},Lprev-1,Mprev-1,Nprev-1);
                    elseif params.femSpace == 2
                        pressure = reshape(estPressure{i},(Lprev-1)/2,(Mprev-1)/2,(Nprev-1)/2);
                    else
                        pressure = reshape(estPressure{i},(Lprev+1)/2,(Mprev+1)/2,(Nprev+1)/2);
                    end
                    pressure = (interp3(y_mat_orig,z_mat_orig,x_mat_orig,pressure,y_mat_interp, z_mat_interp, x_mat_interp,interpMethod,0)); %/rescale_factor_w);
                    estPressure{i} = pressure(:);

                end
            end
        end

        for i=1:num_ts-1
            estU_part{i} = interpFlowToPartmex( estU_grid{i}, (estPart(1:3,:))./repmat(stepsizePyr,1,numpart), Nds, Mds, Lds, params.femSpace );
        end
        disp(['pyr lvl: ' num2str(level)]);
        disp(['mean int before: ' num2str(mean(estPart(4,:)))]);
        
        myPyr = params.my; % * 1/(level);
        itsPyr = params.its; %/sigmaPyr
        lambdaPyr = params.lambda*stepsizePyr(1)^3;
        %lambdaPyr = params.lambda*((N-1)/(width_Pyramid{params.pyramid_levels}-1))^3;
        %lambdaPyr = params.lambda*params.stepsizeInt;
        %lambdaPyr = params.lambda*sigmaPyr;
        [ estPart, estU_grid, estU_part, estPressure, Lp, Lu, Lc, cost ] = partFlowOptimization_multits( imgs, ref_idx, N2d, M2d, sigmaPyr, P, estPart, estU_part, estU_grid, estPressure, N, M, L, Nds, Mds, Lds, itsPyr, Lp, Lu, Lc, stepsizeLip, myPyr, params.sparsityNorm, flowSigma,stepsizePyr,cost,params.divMethod,params.femSpace,params.alpha,stopThresh,lambdaPyr,params.kinetic,params.pcgIts,params.lipThresh, params.lambdaPressure,params.displayFigures); 
        disp(['mean int after: ' num2str(mean(estPart(4,estPart(4,:)>0)))]);

        Lprev = Lds;
        Mprev = Mds;
        Nprev = Nds;
        stepsizePyr_prev = stepsizePyr;

        numpart = size(estPart,2);

    end
    disp(['num vanished: ' num2str(sum(estPart(4,:)<=0))]);
    for i=1:num_ts-1
        estU_part{i}(:,estPart(4,:)<=0) = [];
    end
    estPart(:,estPart(4,:)<=0) = [];
    numpart = size(estPart,2);
    %remove particles that are outside of the measurement volume
    maskOut = estPart(1,:)<1 | estPart(1,:)>N | estPart(2,:)<1 | estPart(2,:)>M | estPart(3,:)<1 | estPart(3,:)>L;
    disp(['num outside: ' num2str(sum(maskOut))]);
    for i=1:num_ts-1
        estU_part{i}(:,maskOut) = [];
    end
    estPart(:,maskOut) = [];
    toc
    
    if writeToFile
        fprintf(fileID, 'Iteration: %d, triang error: %.2f\r\n',triangIt,triangError);
        fprintf(fileID, 'Cost:%f\r\n',cost);
    end
    if gtData
	
		distThresh = 1.0;
		[ numTruePart, numGhostPart, numVanished ] = evaluatePositionError(part,estPart,params.intThresh,distThresh);
		fprintf(fileID, 'Position Error: num true part: %d, num found: %d, num ghost: %d, num vanished: %d\r\n',size(part,2),numTruePart,numGhostPart,numVanished);
		[ numTruePart, numGhostPart, numVanished ] = evaluatePositionError(part,estPart,params.my,distThresh);
		fprintf(fileID, 'Position Error: num true part: %d, num found: %d, num ghost: %d, num vanished: %d\r\n',size(part,2),numTruePart,numGhostPart,numVanished);
    end
end
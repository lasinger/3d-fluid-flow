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

function [ estPart, estU_grid, estU_part, estPressure, Lp, Lu, Lc, costOrig ] = partFlowOptimization_multits( imgs, ref_idx, N2d, M2d, sigma, P, estPart, estU_part, estU_grid, estPressure, N, M, L, Nds, Mds, Lds, its, Lp, Lu, Lc, stepsizeLip, my, sparsityNorm, flowSigma,stepsize, costPrev, divMethod,femSpace,alpha,stopThresh,lambda,kinetic,iterations,lipThresh,lambdaPressure,displayFigures) %, f_residual, f_flow )

num_ts = length(imgs);

numCam = length(P);
numpart = size(estPart,2);

part2d = cell(1,num_ts);
reproj = cell(1,num_ts);
xhom = cell(1,num_ts);
for i=1:num_ts
	part2d{i} = cell(1,numCam);
	reproj{i} = cell(1,numCam);
	xhom{i} = cell(1,numCam);
end

for i=1:num_ts
	for cam=1:numCam
        if i<ref_idx
            [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estPart(1:3,:)+estU_part{i},P{cam});
        elseif i>ref_idx
            [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estPart(1:3,:)+estU_part{i-1},P{cam});
        else
           [part2d{i}{cam}, xhom{i}{cam}] = proj2d(estPart(1:3,:),P{cam});
        end
	end
	reproj{i} = render2dallcammex(part2d{i},estPart(4,:),numpart,sigma,N2d,M2d );
end

costOrig = cost_dataterm_multits(imgs,reproj);
estPart_prev = estPart;
estU_prev = estU_grid;

if divMethod == 0 && alpha == 0 % no div constraint
    Nabla_incomp=0;
    Div_incomp=0;
    Lpl_incomp=0;
elseif femSpace <= 0 %P1P0
    [Nabla_incomp, Div_incomp, Lpl_incomp] = System_Matrices_3d_P1P0( Lds, Mds, Nds); 
elseif femSpace == 2 %P2P0
    [Nabla_incomp, Div_incomp, Lpl_incomp] = System_Matrices_3d_P2P0( Lds, Mds, Nds);
else
    [Nabla_incomp, Div_incomp, Lpl_incomp] = System_Matrices_3d_P2P1( Lds, Mds, Nds);
    if lambdaPressure > 0
        LplLinPressure = lambdaPressure * Lpl_System_Matrices_3d( floor((Lds+1)/2), floor((Mds+1)/2), floor((Nds+1)/2) );
        Lpl_incomp = Lpl_incomp+LplLinPressure;
    end
end

if femSpace == 0
    Nabla = Nabla_System_Matrices_3d_trilin( Lds, Mds, Nds);
    Lpl = Nabla'*Nabla;
elseif femSpace == -1
    Nabla = Nabla_System_Matrices_3d( Lds, Mds, Nds);
    Lpl = Lpl_System_Matrices_3d( Lds, Mds, Nds );
else
    Nabla = Nabla_System_Matrices_3d_triquad(Lds, Mds, Nds);
    Lpl = Lpl_System_Matrices_3d_triquad(Lds, Mds, Nds);
    
end

for warp=1:its

    for cam=1:numCam

        %visualize change every now and then
        if cam==1 && mod(warp,4) == 1 && displayFigures
            subaxis(2,2,1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
            imshow(imresize(imgs{num_ts}{cam}-reproj{num_ts}{cam},2,'nearest'),[]);
            title(['sigma: ' num2str(sigma) ', it: ' num2str(warp) ] );
            xlabel(['Cost: ' num2str(costOrig) ', Diff: ' num2str(costPrev-costOrig) ' , Lp: ' num2str(Lp) ', Lc: ' num2str(Lc) ', Lu: ' num2str(Lu), ', Lambda: ' num2str(lambda)]);
            subaxis(2,2,2, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
            flowU_est = reshape(estU_grid{num_ts-1}(1,:),Lds,Mds,Nds);
            zSlice = 10;
            if L>100
                zSlice = 78;
            end
            ax = gca;
            imshow(imresize(squeeze(flowU_est(round(zSlice/stepsize(3)),:,:)),8),[]);
            colormap(ax,jet);colorbar;
            xlabel('flow U');
            if divMethod
                flow_k = reshape(estU_grid{num_ts-1}',3*Nds*Mds*Lds,1);
                divergence=Div_incomp*flow_k;
                
                if femSpace <= 0
                    flowU_est = reshape(estPressure{num_ts-1},Lds-1,Mds-1,Nds-1);
                    divergence = reshape(divergence,Lds-1,Mds-1,Nds-1);
                elseif femSpace == 2
                    flowU_est = reshape(estPressure{num_ts-1},(Lds-1)/2,(Mds-1)/2,(Nds-1)/2);
                    divergence = reshape(divergence,(Lds-1)/2,(Mds-1)/2,(Nds-1)/2);
                    zSlice = floor(zSlice/2);
                else
                    flowU_est = reshape(estPressure{num_ts-1},(Lds+1)/2,(Mds+1)/2,(Nds+1)/2);
                    divergence = reshape(divergence,(Lds+1)/2,(Mds+1)/2,(Nds+1)/2);
                    zSlice = ceil(zSlice/2);
                end
                
                
                subaxis(2,2,3, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
                ax = gca;
                imshow(imresize(squeeze(divergence(round(zSlice/stepsize(3)),:,:)),8),[]);
                colormap(ax,jet);colorbar;
                xlabel('divergence');
                
                subaxis(2,2,4, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
                ax = gca;
                imshow(imresize(squeeze(flowU_est(round(zSlice/stepsize(3)),:,:)),8),[]);
                colormap(ax,jet);colorbar;
                xlabel('pressure');
            end
            
            
            drawnow
        end
    end
    
    costPrev = costOrig;

    tau_k = 1/sqrt(2);
	
	% particles

    % inertial step
    [estPart, estPart_prev, estU_part, reproj, part2d, xhom ] = inertial_part_multits(estPart, estPart_prev, ref_idx, estU_grid, tau_k, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize,femSpace);
    
    %if pinhole (3x4 matrix)
    if length(P{1}) == 4
        Dc = computeGradient_p_grid_multits( imgs, reproj, part2d, estPart(4,:), xhom, numCam, N2d, M2d, sigma, P, numpart, estPart, estU_grid, estU_part, ref_idx, Nds, Mds, Lds, stepsize, femSpace );
    else % polynomial camera model
        Dc = computeGradient_p_grid_poly_multits( imgs, reproj, part2d, estPart, estU_grid, estU_part, numCam, N2d, M2d, sigma, P, numpart, ref_idx, Nds, Mds, Lds, stepsize, femSpace );
    end
    costOrig = cost_dataterm_multits(imgs,reproj); 
   
    %Lp
    idx = 1:3;
    Lp_old = Lp;
    update_p = true;

    % update step
    [ rightTerm, leftTerm, estTmp_p, reprojTmp, part2dTmp, xhomTmp, estTmp_u_part ] = update_part_multits(estPart(idx,:), estPart(4,:), ref_idx, Lp, Dc(idx,:), estU_grid, costOrig, imgs, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, femSpace);

    if warp == 1 
        while (leftTerm-rightTerm) <= 0
            Lp = Lp/stepsizeLip;
            [ rightTerm, leftTerm, estTmp_p, reprojTmp, part2dTmp, xhomTmp, estTmp_u_part ] = update_part_multits(estPart(idx,:), estPart(4,:), ref_idx, Lp, Dc(idx,:), estU_grid, costOrig, imgs, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, femSpace);
        end
    end
    while (leftTerm-rightTerm) > 0 
        Lp = Lp*stepsizeLip;
        
        if max(sqrt(sum(Dc(idx,:).^2,1)))/Lp < lipThresh
            Lp=Lp_old;
            update_p=false;
            break;
        end
        
        [ rightTerm, leftTerm, estTmp_p, reprojTmp, part2dTmp, xhomTmp, estTmp_u_part ] = update_part_multits(estPart(idx,:), estPart(4,:), ref_idx, Lp, Dc(idx,:), estU_grid, costOrig, imgs, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, femSpace);
    end
    if update_p
        estPart(1:3,:) = estTmp_p;
        estU_part = estTmp_u_part;
        reproj = reprojTmp;
        part2d = part2dTmp;
        xhom = xhomTmp;
        costOrig = leftTerm; % only dataterm
    end

    %inertial step
    [estPart, estPart_prev, reproj ] = inertial_intensity_multits(estPart, estPart_prev, part2d, tau_k, numCam, sigma, N2d, M2d);
    
    %Lc
    Dc =computeGradient_intensity_multits( imgs, reproj, part2d, numCam, N2d, M2d, sigma, numpart );
    costOrig = cost_dataterm_multits(imgs,reproj);
    idx = 4;
    Lc_old = Lc;
    update_c = true;
    
    [ rightTerm, leftTerm, reprojTmp, estTmp_c ] = update_intensity_multits( estPart(idx,:), part2d, Lc, Dc, costOrig, imgs, numCam, sigma, N2d, M2d, sparsityNorm, my); %, lambda, Nabla, estU_grid );
    
    if warp == 1 
        while (leftTerm-rightTerm) <= 0
            Lc = Lc/stepsizeLip;
            [ rightTerm, leftTerm, reprojTmp, estTmp_c ] = update_intensity_multits( estPart(idx,:), part2d, Lc, Dc, costOrig, imgs, numCam, sigma, N2d, M2d, sparsityNorm, my); %, lambda, Nabla, estU_grid );
        end
    end
    while (leftTerm-rightTerm) > 0
        Lc = Lc*stepsizeLip;
        
        % if maximum update is smaller than threshold then don't update
        if max(abs(Dc))/Lc < lipThresh 
            Lc=Lc_old;
            update_c=false;
            break;
        end
        
        [ rightTerm, leftTerm, reprojTmp, estTmp_c ] = update_intensity_multits( estPart(idx,:), part2d, Lc, Dc, costOrig, imgs, numCam, sigma, N2d, M2d, sparsityNorm, my); %, lambda, Nabla, estU_grid );
    end
    if update_c
        estPart(4,:) = estTmp_c;
        reproj = reprojTmp;
        costOrig = leftTerm;
    end


    if sigma <= flowSigma
        for i=1:num_ts
			flow_idx = i;
			if i==ref_idx
				continue
			end
			if i>ref_idx
				flow_idx = i-1;
			end
                
			%inertial step
			[estU_grid{flow_idx}, estU_prev{flow_idx}, estU_part{flow_idx}, reproj{i}, part2d{i}, xhom{i} ] = inertial_flow(estU_grid{flow_idx}, estU_prev{flow_idx}, estPart, tau_k, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize,numpart, femSpace);
			costOrig = cost_dataterm(imgs{ref_idx}, reproj{ref_idx},imgs{i},reproj{i});
            
			% get gradient at grid positions
			%if pinhole (3x4 matrix)
            if length(P{1}) == 4
				gradUgrid = computeGradient_gridU(imgs{i}, reproj{i}, part2d{i}, estPart(4,:), xhom{i}, numCam, N2d, M2d, sigma, P, numpart, estPart(1:3,:)./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace);
            else % polynomial camera model
				gradUgrid = computeGradient_gridU_poly(imgs{i}, reproj{i}, part2d{i}, estPart, estU_part{flow_idx}, numCam, N2d, M2d, sigma, P, numpart, estPart(1:3,:)./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace);
            end
            gradUgrid = gradUgrid/num_ts;

			flow_k = reshape(estU_grid{flow_idx}',3*Nds*Mds*Lds,1);
			Lpl_flow = 0.125*lambda*[Lpl * flow_k(1:Nds*Mds*Lds); Lpl * flow_k(1+Nds*Mds*Lds : 2*Nds*Mds*Lds); Lpl * flow_k(1 + 2*Nds*Mds*Lds : end) ];
			alphaDiv = 0;
			if divMethod == 0 && alpha ~= 0
				%alphaDiv = 0.125*alpha*[Lpl_incomp*flow_k(1:Nds*Mds*Lds); Lpl_incomp * flow_k(1+Nds*Mds*Lds : 2*Nds*Mds*Lds); Lpl_incomp * flow_k(1 + 2*Nds*Mds*Lds : end)];
				alphaDiv = 0.125*alpha*lambda*Div_incomp'*(Div_incomp*flow_k);
			end
			flow_k_DataDeriv = reshape(gradUgrid',3*Nds*Mds*Lds,1) + Lpl_flow + alphaDiv;
			
			Lu_old = Lu;
			update_u = true;
			
			[ rightTerm, leftTerm, reprojTmp1, part2dTmp1, xhomTmp1, estTmp_u_grid, estTmp_pressure, estTmp_u_part ] = update_flow( estPart, Lu, flow_k_DataDeriv, flow_k, Lpl_flow, estPressure{flow_idx}, costOrig, reproj{ref_idx}, imgs{ref_idx}, imgs{i}, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, lambda,kinetic, iterations,  Nabla_incomp, Div_incomp, Lpl_incomp, Nabla,divMethod,femSpace,alpha);
			if warp == 1 
				while (leftTerm-rightTerm) <= 0 && Lu>0.25 % hard block of Lu getting to small cause sometimes that seems to happen (either bcs of to small radius for rendering or bcs of limited pcg iterations - but if it happens behaviour is weird)
					Lu = Lu/stepsizeLip;
					[ rightTerm, leftTerm, reprojTmp1, part2dTmp1, xhomTmp1, estTmp_u_grid, estTmp_pressure, estTmp_u_part ] = update_flow( estPart, Lu, flow_k_DataDeriv, flow_k, Lpl_flow, estPressure{flow_idx}, costOrig, reproj{ref_idx}, imgs{ref_idx}, imgs{i}, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, lambda,kinetic, iterations, Nabla_incomp, Div_incomp, Lpl_incomp, Nabla,divMethod,femSpace,alpha);
				end
			end
			while (leftTerm-rightTerm) > 0
				Lu = Lu*stepsizeLip;
				
				% if maximum update is smaller than threshold then don't update
				if max(sqrt(flow_k_DataDeriv(1:Nds*Mds*Lds).^2+flow_k_DataDeriv(Nds*Mds*Lds+1:2*Nds*Mds*Lds).^2+flow_k_DataDeriv(2*Nds*Mds*Lds+1:end).^2))/Lu < lipThresh
					Lu=Lu_old;
					update_u=false;
					break;
                end
				
				[ rightTerm, leftTerm, reprojTmp1, part2dTmp1, xhomTmp1, estTmp_u_grid, estTmp_pressure, estTmp_u_part ] = update_flow( estPart, Lu, flow_k_DataDeriv, flow_k, Lpl_flow, estPressure{flow_idx}, costOrig, reproj{ref_idx}, imgs{ref_idx}, imgs{i}, P, numCam, sigma, N2d, M2d, Nds, Mds, Lds, stepsize, lambda,kinetic, iterations, Nabla_incomp, Div_incomp, Lpl_incomp, Nabla, divMethod, femSpace, alpha );
				
			end        
			if update_u
				estU_grid{flow_idx} = estTmp_u_grid;
				estU_part{flow_idx} = estTmp_u_part;
				estPressure{flow_idx} = estTmp_pressure;
				reproj{i} = reprojTmp1;
				xhom{i} = xhomTmp1;
				part2d{i} = part2dTmp1;
				costOrig = leftTerm;
			end
		
        end
        
        costOrig = cost_dataterm_multits(imgs,reproj);

    end
   
    if (~update_p && ~update_c && ~update_u) || (costPrev-costOrig) < stopThresh
        disp(['stop at it ' num2str(warp)]);
        break;
    end
    

end

disp(['sigma: ' num2str(sigma)]);
disp(['Lp: ' num2str(Lp) ', Lc: ' num2str(Lc) ', Lu: ' num2str(Lu)]);
disp(['cost: ' num2str(leftTerm) ', diff: ' num2str(costOrig-costPrev)]);

disp(['avg part intensity:' num2str(mean(estPart(4,:)))]);
disp(['avg intensity of part>0:' num2str(mean(estPart(4,estPart(4,:)>0)))]);
disp(['num part 0: ' num2str(sum(estPart(4,:)<=0))]);


end


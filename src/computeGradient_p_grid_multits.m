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

function [ Dc ] = computeGradient_p_grid_multits( imgs, reproj, part2d, partInt, xhom, numCam, N2d, M2d, sigma, P, numpart, estPart, estU_grid, estU_part, ref_idx, Nds, Mds, Lds, stepsize, femSpace )

Dc = zeros(3,numpart);

for i=1:length(imgs)
    for cam=1:numCam

        % compute gradient
        diff = imgs{i}{cam} - reproj{i}{cam};

        % derive u,v,w parts in x,y and z direction
        dudpx = zeros(3,numpart);
        dudpy = zeros(3,numpart);
        dudpz = zeros(3,numpart);
        if i~=ref_idx
            idx=i;
            if i>ref_idx
                idx=i-1;
            end
            % bring all to the same coordinate system - lower res grid -->
            % also displacement and thus gradient has to be rescaled... -
            [dudpx, dudpy, dudpz] = derivativeFlowPartmex( estU_grid{idx}./stepsize, (estPart(1:3,:)+estU_part{idx}(1:3,:))./repmat(stepsize,1,numpart), Nds, Mds, Lds, femSpace );
            
        end
        
        xhom1=xhom{i}{cam}(1,:);
        xhom2=xhom{i}{cam}(2,:);
        xhom3=xhom{i}{cam}(3,:);
        xhom3_2=xhom3.^2;
        Pcam = P{cam};
        
        dfxdpx = Pcam(1,1).*(1+dudpx(1,:))+Pcam(1,2).*(dudpx(2,:))+Pcam(1,3).*(dudpx(3,:));
        dfxdpy = Pcam(1,1).*(dudpy(1,:))+Pcam(1,2).*(1+dudpy(2,:))+Pcam(1,3).*(dudpy(3,:));
        dfxdpz = Pcam(1,1).*(dudpz(1,:))+Pcam(1,2).*(dudpz(2,:))+Pcam(1,3).*(1+dudpz(3,:));
        
        dfydpx = Pcam(2,1).*(1+dudpx(1,:))+Pcam(2,2).*(dudpx(2,:))+Pcam(2,3).*(dudpx(3,:));
        dfydpy = Pcam(2,1).*(dudpy(1,:))+Pcam(2,2).*(1+dudpy(2,:))+Pcam(2,3).*(dudpy(3,:));
        dfydpz = Pcam(2,1).*(dudpz(1,:))+Pcam(2,2).*(dudpz(2,:))+Pcam(2,3).*(1+dudpz(3,:));
        
        dgdpx = Pcam(3,1).*(1+dudpx(1,:))+Pcam(3,2).*(dudpx(2,:))+Pcam(3,3).*(dudpx(3,:));
        dgdpy = Pcam(3,1).*(dudpy(1,:))+Pcam(3,2).*(1+dudpy(2,:))+Pcam(3,3).*(dudpy(3,:));
        dgdpz = Pcam(3,1).*(dudpz(1,:))+Pcam(3,2).*(dudpz(2,:))+Pcam(3,3).*(1+dudpz(3,:));
        
        ProjPart1_I_x = ((dgdpx.*xhom1)./xhom3_2 - dfxdpx./xhom3);
        ProjPart1_I_y = ((dgdpx.*xhom2)./xhom3_2 - dfydpx./xhom3);
        ProjPart2_I_x = ((dgdpy.*xhom1)./xhom3_2 - dfxdpy./xhom3);
        ProjPart2_I_y = ((dgdpy.*xhom2)./xhom3_2 - dfydpy./xhom3);
        ProjPart3_I_x = ((dgdpz.*xhom1)./xhom3_2 - dfxdpz./xhom3);
        ProjPart3_I_y = ((dgdpz.*xhom2)./xhom3_2 - dfydpz./xhom3);

        
        part_grad = gradientmex(part2d{i}{cam}(1,:),part2d{i}{cam}(2,:),partInt,numpart,sigma,ProjPart1_I_x,ProjPart1_I_y,ProjPart2_I_x,ProjPart2_I_y,ProjPart3_I_x,ProjPart3_I_y,diff(:),N2d,M2d);
       
        Dc = Dc + part_grad;
    end
end

Dc = Dc/numCam/length(imgs);

Dc = -Dc;

end


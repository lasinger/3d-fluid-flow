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

function [ppY,ppX] = subpixFit2d( map, ppY,ppX,perDim )

% map: 3x3 nbh around peak

A=zeros(3*3,5);

row =1;
for dyy=-1:1
   for dxx=-1:1
       A(row,:) = [dxx*dxx dyy*dyy dxx dyy dxx*dyy ];
       row=row+1;
   end
end


%per dimension
[ymax,xmax]= size(map);

for i=1:length(ppY)

    x0=ppX(i);
    y0=ppY(i);
    
    if(perDim)
        %perDim
       

        % quadratic interpolation in y-direction
        %yy= y0 - round((ymax+1)/2);
        if (y0>1) && (y0<ymax)
            pp= mypolyfit([-1:1:1],log(squeeze(map(y0-1:1:y0+1,x0))+eps),2);
            %yy= yy - pp(2)/(2*pp(1));
            ppY(i)=y0 - pp(2)/(2*pp(1));
        end

        % quadratic interpolation in x-direction
        %xx= x0 - round((xmax+1)/2);
        if (x0>1) && (x0<xmax)
            pp= mypolyfit([-1:1:1]',log(squeeze(map(y0,x0-1:1:x0+1))+eps),2);
            %xx= xx - pp(2)/(2*pp(1));
            ppX(i) = x0 - pp(2)/(2*pp(1));
        end
    else
        %all dimensions
        b=zeros(3*3,1);
        row =1;
       for dyy=-1:1
           for dxx=-1:1
               xNb = x0+dxx;
               yNb = y0+dyy;
               if(xNb<1 || xNb>xmax)
                   xNb= x0-dxx;
               end
               if(yNb<1 || yNb>ymax)
                   yNb= y0-dyy;
               end
               b(row) = 1-map(yNb,xNb);
               row=row+1;
           end
       end
        fp = A\b;

        %set up function and get minimum
        fun = @(x) ( fp(1)*x(1)*x(1) + fp(2)*x(2)*x(2) + fp(3)*x(1) + fp(4)*x(2) + fp(5)*x(1)*x(2) );
        x00 = [0,0];
        [x,mc] = fminsearch(fun,x00);

        ppX(i)=ppX(i)+x(1);
        ppY(i)=ppY(i)+x(2);
    end

end



%end

function [ p ] = mypolyfit(x,y,n)
% basically copy of polyfit without all that extra stuff
x = x(:);
y = y(:);
% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);

p = R\(Q'*y);    % Same as p = V\y;
p = p.'; 

%end
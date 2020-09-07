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

function P = dlt(x,X)

    if(size(X,1)>size(X,2))
        X =X';
    end
    if(size(x,1)>size(x,2))
        x = x';
    end
    
    my_x = mean(x')';
    my_X = mean(X')';
    
    s_x = sqrt(2)*length(x)/(sum(sqrt(sum((x-repmat(my_x,1,length(x))).^2))));
    s_X = sqrt(2)*length(X)/(sum(sqrt(sum((X-repmat(my_X,1,length(X))).^2))));
    
    Tx = s_x * [eye(2),-my_x; 0, 0, 1/s_x];
    TX = s_X * [eye(3),-my_X; 0, 0,0, 1/s_X];
    
    x = [ x; ones(1,size(x,2))];
    X = [ X; ones(1,size(X,2))];
    
    y = Tx*x;
    Y = TX*X;
    
    A = zeros(size(X,2)*2,12);
    for i = 1:size(X,2)
        A(2*(i-1)+1,:) = [  zeros(1,4)    , -Y(:,i)'       ,  y(2,i)*Y(:,i)'];
        A(2*(i-1)+2,:) = [  Y(:,i)'       ,  zeros(1,4)    , -y(1,i)*Y(:,i)'];
    end
    
    [~,~,V] = svd(A);
    Q = reshape(V(:,end),[4,3]);
    Q = Q';
    P = Tx\Q*TX;
    
    if det(P(1:3,1:3))<0
        P = -P;
    end
    
end






















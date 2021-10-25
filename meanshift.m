% MIT License
% 
% Copyright (c) 2021 Fraunhofer-Chalmers Centre
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
% Implementation of the mean shift algorithm
% Cheng, Yizong (1995). Mean Shift, Mode Seeking, and Clustering. 
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 17, 790-799.
%
% [mode,SE] = meanshift(X,x0,K,n)
% 
% inputs:
% X     - data set
% x0    - starting value for mode search
% K     - smoothing function
% n     - number of iterations
%
% outputs: 
% mode  - mode of smoothed distribution
% SE    - standard deviation normalized by mode

function [mode,SE] = meanshift(X,x0,K,n)

X_std = std(X);
X_min = min(X);

X = X - X_min;

w = x0-X_min;

for i = 1 : n
    Kn = K(w);
    Kn = Kn-logsum(Kn);
    for j = 1 : size(X,2)
        w(j) =  exp(logsum( log(X(:,j)) + Kn )); 
    end
end

mode = w + X_min;
SE = abs(std(X)./mode);

function S = logsum(L)

% uses the identity log(a + b) = log(a * (1 + b/a)) = log a + log(1 + b/a)
% to get log(sum([p1 p2 p3 p4 ...])) given log([p1 p2 p3 p4 ...])
% without losing information through underflow in exp(log([p1 p2 p3 p4 ...]))
%
% S = logsum(L)
% -------------
% L - logarithms of values  
% S - logarithm of sum of values

L(L==-inf) = [];
L = sort(L); % important for numerical precision
if isempty(L)
    S = -inf;
    return
end
S = L(1);
for i = 2 : length(L)
    S = L(i) + log(1+exp(S-L(i)));
end


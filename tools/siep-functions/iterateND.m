function f = iterateND(X,q,Q,theta)
% SYNTAX:
%
%     f = iterateND(X,q,Q) : (theta) optional
%
% DESCRIPTION:
%
%    Implementation of iteration map for SVD power algorithm applied to 
%    ND tensor
% 
% INPUTS:
%   X:  tensor of order N and dimension (n1,n2,...,nN)
%   q:  vector of dimension (sum(n)+1) used as initial conditon of
%       iteraton map
%   Q: cell containing projection matrices Q{1}, ..., Q{N} of dimensions
%   n1^2, n2^2,...,nN^2
%
% OUTPUTS:
%   F
%   
% AUTHOR:
%     Siep Weiland
%
% VERSION:
%     Created June 2019
%     Last update: June 16, 2019

% define tensor dimensions (3rd order only)
n = size(X);

% define index ranges:
N = length(n);
ind = cell(N,1);

w = cell(N,1);

for i = 1:N
    if i == 1
        nmin = 1; nmax = n(1);
    else
        nmin = sum(n(1:i-1))+1;
        nmax = sum(n(1:i));
    end
    ind{i} = nmin:nmax;

    w{i} = q(nmin:nmax);
    w{i} = Q{i}*w{i}./norm(Q{i}*w{i}); % Project and normalize the iteration vectors on im Q_i

end
sigma = abs(q(end));

% Compute the iteration result
X = tensor(X);

fcell = cell(N,1);
for i = 1:N
    ind_w = 1:N;
    ind_w(i) = [];
    wi = {w{ind_w}};
    fcell{i} = double(ttv(X, wi,-i));
end
fend = double(ttv(X,w,1:N));

f = zeros(size(q));
for i = 1:N
    f(ind{i}) = (1/sigma).*fcell{i};
end
f(end) = fend; 

f = theta*f+(1-theta)*q; % Is this needed? Also, if theta is not provided is it right to assume theta = 1?


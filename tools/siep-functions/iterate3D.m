function f = iterate3D(X,q,Q1,Q2,Q3,theta)
% SYNTAX:
%
%     f = iterate3D(X,q,Q1,Q2,Q3)
%
% DESCRIPTION:
%
%    Implementation of iteration map for SVD power algorithm applied to 
%    3D tensor
% 
% INPUTS:
%   X:  tensor of order 3 and dimension (n1,n2,n3)
%   q:  vector of dimension (n1+n2+n3+1) used as initial conditon of
%       iteraton map
%   Q1, Q2, Q3: projection matrices of dimensions n1^2, n2^2, n3^3
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

[n1,n2,n3] = size(X);

ind_1=1:n1; 
ind_2=(n1+1):(n1+n2);
ind_3=(n1+n2+1):(n1+n2+n3);

w1 = q(ind_1);
w2 = q(ind_2);
w3 = q(ind_3);
sigma = abs(q(end));

% Project and normalize the iteration vectors on im Q_i
w1 = Q1*w1./norm(Q1*w1);
w2 = Q2*w2./norm(Q2*w2);
w3 = Q3*w3./norm(Q3*w3);

% Compute the iteration result
X = tensor(X);

f1 = double(ttv(X,{w2,w3},-1));
f2 = double(ttv(X,{w1,w3},-2));
f3 = double(ttv(X,{w1,w2},-3));
f4 = double(ttv(X,{w1,w2,w3},[1 2 3]));

f = [(1/sigma).*f1; (1/sigma).*f2; (1/sigma).*f3; f4];
% f = theta*f+(1-theta)*q; % I (Bogoljub) commented this


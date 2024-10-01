function qout = iterate3D_W(X,q,W1,W2,W3,theta)
% SYNTAX:
%
%     qout = iterate3D(X,q,W1,W2,W3)
%
% DESCRIPTION:
%
%    Implementation of iteration map for SVD power algorithm applied to 
%    3D tensor with weighted inner products. 
%    Implements the iteration map q \mapsto G(q) with
%
%            G(q) = \frac{1}{T(q)} * \inv(W) * \nabla T(q) 
%
% INPUTS:
%   X:  tensor of order 3 and dimension (n1,n2,n3)
%   q:  vector of dimension (n1+n2+n3+1) used as initial conditon of
%       iteraton map to be conditoned as 
%       q=\col(x_i/\|x_i\|_{W_i} | i=1,2,3 )
%   W1, W2, W3: weighting matrices of dimensions n1^2, n2^2, n3^3
%
% OUTPUTS:
%   qout = G(q);
%   
% AUTHOR:
%     Siep Weiland
%
% VERSION:
%     Created November 2019
%     Last update: November 5, 2019
%     Note: to be adjested to allow for projections in higher order 
%     singular values.

[n1,n2,n3] = size(X);

ind_1=1:n1; 
ind_2=(n1+1):(n1+n2);
ind_3=(n1+n2+1):(n1+n2+n3);

w1 = q(ind_1);
w2 = q(ind_2);
w3 = q(ind_3);
sigma = abs(q(end));

% Normalize the iteration vector w_i w.r.t. weighted inner product
w1 = w1./sqrt(w1'*W1*w1);
w2 = w2./sqrt(w2'*W2*w2);
w3 = w3./sqrt(w3'*W3*w3);

% Compute the iteration result
X = tensor(X);

% KKT Lagrange multiplier
lambda = double(ttv(X,{w1,w2,w3},[1 2 3])); 

% KKT Gradient conditions
f1 = (1/lambda).*inv(W1)*double(ttv(X,{w2,w3},[2 3]));
f2 = (1/lambda).*inv(W2)*double(ttv(X,{w1,w3},[1 3]));
f3 = (1/lambda).*inv(W3)*double(ttv(X,{w1,w2},[1 2]));

qout = theta*[f1; f2; f3; lambda]+(1-theta)*[f1; f2; f3; lambda];


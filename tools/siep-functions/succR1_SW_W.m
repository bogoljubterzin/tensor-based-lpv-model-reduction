function [x1,x2,x3,gains,iterations,Approx] = succR1_SW_W(X,max_level,W1,W2,W3)
% SYNTAX:
%
% [x1,x2,x3,gains,iterations] = succR1_SW_W(X,r,W1,W2,W3)
%
% DESCRIPTION:
%
% Computation of r successive rank-one approximations for an order-N tensor
% with weighted inner products. At this moment N=3.
% 
% INPUTS:
%    X:        multidimensional array, double formatted
%    r:        number of successive rank one approximations that are computed 
%    W1,W2,W3: positive definite weigthing matrices on inner products of 
%              modal directions
%
% OUTPUTS:
%     x1,x2,x3:   orthonormal columns that define the singular vectors
%     gains:      vector of corresponding singular values
%     iteratons:  vector with number of iterations per level
%     Approx:     Approximate tensor consisting of sum of r rank-one
%     approximations
%
% REMARKS:
%    invokes:
%        'gs_orth' to compute a gram-schmidt orthonormalization
%                  of the basis functions computed.
%         'rank1tensor' to compute a rank-one tensor from three vectors and scalar. 
%         'iterate3D' to compute iterations in the fixed-point algorithm.
%
% AUTHOR:
%     Siep Weiland
% VERSION:
%     Created June 2019
%     Last update: November 5, 2019



%%%% Initialization
% define maximum number of iterations
maxiter = 1000;   

r = max_level; % I added this line

% define tensor dimensions
[n1,n2,n3] = size(X);

% define index ranges:
ind_1=1:n1; 
ind_2=(n1+1):(n1+n2);
ind_3=(n1+n2+1):(n1+n2+n3);

% define projection matrices: (NOT USED)
%Q1 = eye(n1);
%Q2 = eye(n2);
%Q3 = eye(n3);

% initialize subspaces of singular vectors:
x1 = inf(n1,r);
x2 = inf(n2,r);
x3 = inf(n3,r);

% initialize iteration numbers for all SVD levels:
iterations = inf(r,1);

% initialize gains:
gains = inf(r,1);

% Set defaults
if nargin < 5
     W1 =eye(n1); W2=eye(n2), W3=eye(n3)
end

%%%% Compute the optimal rank-one approximant at level 1

q = [rand(n1,1);rand(n2,1);rand(n3,1);rand(1,1)];
qnew = iterate3D_W(X,q,W1,W2,W3);  % first "weighted" iteration

k=1; % iteration number
while (norm(q-qnew)>1e-6)
    if k == maxiter
        disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
        break
    end
    q = qnew;
    qnew = iterate3D_W(X,q,W1,W2,W3);
    k = k+1;
end
disp('Iteration completed at level 1')
% Store number of iterations
iterations(1) = k;

% Update x1,x2,x3,
x1hat = qnew(ind_1);
x2hat = qnew(ind_2);
x3hat = qnew(ind_3);
x1(:,1) = x1hat./sqrt(x1hat'*W1*x1hat); % normalize w.r.t. weighted inner product
x2(:,1) = x2hat./sqrt(x2hat'*W2*x2hat); % normalize w.r.t. weighted inner product
x3(:,1) = x3hat./sqrt(x3hat'*W3*x3hat); % normalize w.r.t. weighted inner product

% Store singular value (or gain, or Lagrange multiplier)
gains(1) = qnew(end);

% Define approximate tensor
Approx =rank1tensor(gains(1),1,1,1,x1,x2,x3);
X_err = X - Approx;
disp(['relative error level 1:', num2str(frob(X_err)/frob(X))]);


%%% Proceed with next levels of singular value decompositions till
%%% reaching level max_level

for i=2:max_level
    
    % initialize q
    q = [rand(n1,1);rand(n2,1);rand(n3,1);rand(1,1)];
    qnew = iterate3D_W(X_err,q,W1,W2,W3); % first iteration on error
    
    k=1; % iteration number
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
            break
        end
        q = qnew;
        qnew = iterate3D_W(X_err,q,W1,W2,W3);
        k = k+1;
    end
    disp(['Iteration completed at level ', num2str(i)])
    iterations(i) = k;
    % update x1,x2,x3,result
    x1hat = qnew(ind_1);
    x2hat = qnew(ind_2);
    x3hat = qnew(ind_3);
    x1(:,i) = x1hat./sqrt(x1hat'*W1*x1hat); % normalize w.r.t. weighted inner product
    x2(:,i) = x2hat./sqrt(x2hat'*W2*x2hat); % normalize w.r.t. weighted inner product
    x3(:,i) = x3hat./sqrt(x3hat'*W3*x3hat); % normalize w.r.t. weighted inner product
    gains(i) = qnew(end);
    
    % Update approximation at this level
    Approx = Approx + rank1tensor(gains(i),i,i,i,x1,x2,x3);
    X_err = X-Approx;
end


%%% Orthonormalize the basis vectors found for each vector space
%%%   --> may like to skip this step
%x1 = gs_orth(x1);
%x2 = gs_orth(x2);
%x3 = gs_orth(x3);


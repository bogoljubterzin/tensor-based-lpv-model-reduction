function [x1,x2,x3,gains,iterations,Approx] = succR1_SW(X,max_level)
% SYNTAX:
%
% [x1,x2,x3,gains,iterations] = succR1_SW(X,num)
%
% DESCRIPTION:
%
% Computation of successive rank-one approximations for an order-N tensor.
% At the moment N=3.
% 
% INPUTS:
%    X:     multidimensional array, double formatted
%    num:   number of singular values that are computed
%
% OUTPUTS:
%     x1,x2,x3:   orthonormal columns that define the singular vectors
%     gains:     vector of corresponding singular values
%     iteraton_numbers: number of iterations per level
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
%     Last update: June 16, 2019



%%%% Initialization
% define maximum number of iterations
maxiter = 1000;   

% define tensor dimensions
[n1,n2,n3] = size(X);

% define index ranges:
ind_1=1:n1; 
ind_2=(n1+1):(n1+n2);
ind_3=(n1+n2+1):(n1+n2+n3);

% define projection matrices:
Q1 = eye(n1);
Q2 = eye(n2);
Q3 = eye(n3);

% initialize subspaces of singular vectors:
x1 = [];
x2 = [];
x3 = [];

% initialize iteration numbers for all SVD levels:
iterations = [];


%%%% Compute the optimal rank-one approximant at level 1

iterationconstant = rand(1,1);  %%%CHANGE THIS
q = [rand(n1,1);rand(n2,1);rand(n3,1);iterationconstant];
qnew = iterate3D(X,q,Q1,Q2,Q3);  % first iteration

k=1; % iteration number
while (norm(q-qnew)>1e-6)
    if k == maxiter
        disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
        break
    end
    q = qnew;
    qnew = iterate3D(X,q,Q1,Q2,Q3);
    k = k+1;
end
disp('Iteration completed at level 1')
iterations(1) = k;
% update x1,x2,x3,
x1hat = qnew(ind_1);
x2hat = qnew(ind_2);
x3hat = qnew(ind_3);
x1 = x1hat./norm(x1hat); % normalize
x2 = x2hat./norm(x2hat); % normalize
x3 = x3hat./norm(x3hat); % normalize
gains = qnew(end);

%%% Proceed with next levels of singular value decompositions till
%%% reaching level max_level

Approx =rank1tensor(gains(1),1,1,1,x1,x2,x3);
X_err=X-Approx;
for i=2:max_level
    % initialize q
    iterationconstant = rand(1,1);  %%%CHANGE THIS
    q = [rand(n1,1);rand(n2,1);rand(n3,1);iterationconstant];
    qnew = iterate3D(X_err,q,Q1,Q2,Q3); % first iteration on error
    
    k=1; % iteration number
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
            break
        end
        q = qnew;
        qnew = iterate3D(X_err,q,Q1,Q2,Q3);
        k = k+1;
    end
    disp(['Iteration completed at level ', num2str(i)])
    iterations(i) = k;
    % update x1,x2,x3,result
    x1hat = qnew(ind_1);
    x2hat = qnew(ind_2);
    x3hat = qnew(ind_3);
    x1(:,i) = x1hat./norm(x1hat);
    x2(:,i) = x2hat./norm(x2hat);
    x3(:,i) = x3hat./norm(x3hat);
    gains(i) = qnew(end);
    %Update the tensor
    Approx = Approx + rank1tensor(gains(i),i,i,i,x1,x2,x3); 
    X_err = X_err   - rank1tensor(gains(i),i,i,i,x1,x2,x3);  % = X-Approx
end

gains = gains';

%%% Orthonormalize the basis vectors found for each vector space
%%%   --> may like to skip this step
%x1 = gs_orth(x1);
%x2 = gs_orth(x2);
%x3 = gs_orth(x3);


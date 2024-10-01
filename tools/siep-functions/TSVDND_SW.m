function [x,gains,iterations,Xa] = TSVDND_SW(X,num_levels,theta)
% SYNTAX:
%
% [x,result,iterations,Xa] = TSVDND_SW(X,num_levels)
%
% DESCRIPTION:
%
% Computation of complete orthogonal tensor SVD for an order-N tensor till
% level num_levels. At the moment N=3.
% 
% INPUTS:
%    X:     multidimensional array, double formatted
%    num:   number of singular values that are computed
%
% OUTPUTS:
%     x:  cell array with elements x1,x2,...,xN
%               orthonormal columns that define the singular vectors
%     gains:     vector of corresponding singular values
%     iteraton_numbers: number of iterations per level
%     Xa :  approximate tensor till level num_levels
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


%%%% Initialization
% define maximum number of iterations
if nargin == 2
    theta = 1;
end

% tol = 1e-6;
tol = 1e-10;
% X = tensor(X);

maxiter = 10000;    

% define tensor dimensions 
n = size(X);

% define index ranges:
N = length(n);
ind = cell(N,1);
for i = 1:N
    if i == 1
        nmin = 1; nmax = n(1);
    else
        nmin = sum(n(1:i-1))+1;
        nmax = sum(n(1:i));
    end
    ind{i} = nmin:nmax;
end

% define initial projection matrices:
Q = cell(N,1);
for i = 1:N
    Q{i} = eye(n(i));
end


% initialize singular vectors:
x = cell(N,1);
for i = 1:N
    x{i} = [];
end

% initialize gain vector:
gains = zeros(num_levels,1);

% initialize iteration vector:
iterations = zeros(num_levels,1);

% check if num_levels < K = min modrank(X)
R = modrank(X);
K = min(R);
if num_levels > K
    num_levels = K;
    fprintf("Number of levels has to be lower or equal to the minimum modal rank of a tensor. \n");
end
%%%% Computation of TSVD
% Compute i-th level TSVD for all i <= num_levels
for i=1:num_levels
    % initialize q
    iterationconstant = rand(1,1);  
    q = zeros(ind{N}(end)+1,1);
    for j = 1:N
        q(ind{j}) = rand(n(j),1);
    end
    q(end) = iterationconstant;

    qnew = iterateND(X,q,Q,theta);  % run first iteration
    
    k=1;   % iteration number
    while (norm(q-qnew)>tol)
        if k == maxiter
            disp(['Maximum number of iterations attained at level ', num2str(i) ])
            break
        end
        q = qnew;
        qnew = iterateND(X,q,Q,theta);
        k = k+1;
    end
    disp(['Iteration at level ', num2str(i) , ' completed....'])
    
    iterations(i) = k;
    gains(i) = (qnew(end)); 

    % update x1,x2,... xN,Q1,Q2,...,QN
    for j = 1:N
        xhat = qnew(ind{j});
        x{j}(:, i) = Q{j}*xhat./norm(Q{j}*xhat);
        Q{j} = eye(n(j)) - x{j}*pinv(x{j});
       
    end
    X = double(ttm(tensor(X),Q,1:N));
end

%%% Define approximation
Xa=squeeze(zeros(size(X))); % added squeeze
for i = 1:num_levels

   Approx =  tensor(x{1}(:,i));
   for j = 1:N-1
        Approx = (ttt(Approx,tensor(x{j+1}(:,i))));
        Approx = squeeze(Approx);
   end
   Approx = gains(i).*double(Approx);
   Xa = Xa + Approx;
end




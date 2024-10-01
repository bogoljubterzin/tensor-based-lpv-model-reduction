function [x1,x2,x3,gains,iterations,Approx,W1,W2,W3,store_k_loops] = succR1_SW_dynamic_W(X,max_level,F1,F2,F3,inputs,maxiter)
% SYNTAX:
%
% [x1,x2,x3,gains,iterations] = succR1_SW(X,num,W1,W2,W3)
%
% DESCRIPTION:
%
% Computation of successive rank-one approximations for an order-N tensor.
% At the moment N=3.
%
% INPUTS:
%    X:     multidimensional array, double formatted
%    num:   number of singular values that are computed
% .  W1,W2,W3:  positive definite weigthin matrices on modal directions
%    inputs: input signals that define the operating points on which the
%    Tensor is computed. Needed to compute a sort matrix S in order to sort SDg_int
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
% vERSION:
%     Created June 2019
%     Last update: June 16, 2019



%%%% Initialization
% define maximum number of iterations
% maxiter = 1e4; % input argument

store_k_loops = zeros(maxiter,max_level);
theta = 1;
Scale = 10^-2;

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
x1 = inf(n1,max_level);
x2 = inf(n2,max_level);
x3 = inf(n3,max_level);

% initialize iteration numbers for all SvD levels:
iterations = inf(max_level,1);
gains = inf(max_level,1);

if nargin < 5
    F1 =eye(n1); F2=eye(n2); F3=eye(n3); % filter for innproducts
end

%%%% wieghing matrices %%%%
% of the form:
% W = S-1*F*S with S a sorting matrix according to the z-argument

%%%% Compute the optimal rank-one approximant at level 1 %%%%
% --------------------------------------------------------- %

q = [rand(n1,1);rand(n2,1);rand(n3,1);rand(1,1)];

% Sorting matrix %
% extract v to compute z %
v = q(ind_2);
v = v(:);
z = v.'*inputs;
diff_z = diff(sort(z));
mean_diff_z_sq = Scale*mean(diff_z)^2;

[~,sortI] = sort(z);
S = zeros(n3,n3); % sorting matrix
for s=1:n3
    S(s,sortI(s)) = 1;
end

% wieghing matrices %
W1 = F1;
W2 = F2;
D2mean = F3./mean_diff_z_sq;
W3 = S^-1*D2mean*S;
% W3 = S^-1*F3*S;
W3 = W3'*W3;


% qnew = iterate3D(X,q,Q1*chol(W1),Q2*chol(W2),Q3*chol(W3),theta);  % first iteration
qnew = iterate3D_W(X,q,W1,W2,W3,theta);  % first iteration

% Sorting matrix %
% extract v to compute z %
v = qnew(ind_2);
v = v(:);
z = v.'*inputs;
diff_z = diff(sort(z));
mean_diff_z_sq = Scale*mean(diff_z)^2;

[~,sortI] = sort(z);
S = zeros(n3,n3); % sorting matrix
for s=1:n3
    S(s,sortI(s)) = 1;
end

% wieghing matrices %
W1new = F1;
W2new = F2;
D2mean = F3./mean_diff_z_sq;
W3new = S^-1*D2mean*S;
% W3new = S^-1*F3*S;
W3new = W3new'*W3new;




k=1; % iteration number
store_k_loops(k,1) = norm(q-qnew);
while (norm(q-qnew)>1e-6)
    if k == maxiter
        disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
        break
    end
    q = qnew;
    W1 = W1new;
    W2 = W2new;
    W3 = W3new;
    
%     qnew = iterate3D(X,q,Q1*chol(W1),Q2*chol(W2),Q3*chol(W3),theta);
    qnew = iterate3D_W(X,q,W1,W2,W3,theta);
    
    % Sorting matrix %
    % extract v to compute z %
    v = qnew(ind_2);
    v = v(:);
    z = v.'*inputs;
    diff_z = diff(sort(z));
    mean_diff_z_sq = Scale*mean(diff_z)^2;
    
    [~,sortI] = sort(z);
    S = zeros(n3,n3); % sorting matrix
    for s=1:n3
        S(s,sortI(s)) = 1;
    end
    
    % weighting matrices %
    W1new = F1;
    W2new = F2;
    D2mean = F3./mean_diff_z_sq;
    W3new = S^-1*D2mean*S; 
%     W3new = S^-1*F3*S;
    W3new = W3new'*W3new;
    
    store_k_loops(k,1) = norm(q-qnew);
    k = k+1;
    if(mod(k,100) == 0)
        disp(['Singular value 1/' num2str(max_level) ' - iter ' num2str(k) '/' num2str(maxiter)])
    end
    
end
disp('Iteration completed at level 1')
iterations(1) = k;
% update x1,x2,x3,
x1hat = qnew(ind_1);
x2hat = qnew(ind_2);
x3hat = qnew(ind_3);
x1(:,1) = x1hat./sqrt(x1hat'*W1*x1hat); % normalize w.r.t. weighted inner product
x2(:,1) = x2hat./sqrt(x2hat'*W2*x2hat); % normalize w.r.t. weighted inner product
x3(:,1) = x3hat./sqrt(x3hat'*W3*x3hat); % normalize w.r.t. weighted inner product
gains(1) = qnew(end);

%%% Proceed with next levels of singular value decompositions till
%%% reaching level max_level

Approx = rank1tensor(gains(1),1,1,1,x1,x2,x3);
X_err = X-Approx;
disp(['rel error 1: ' num2str(frob(X_err)/frob(X))])

for i=2:max_level

    % initialize q
    q = [rand(n1,1);rand(n2,1);rand(n3,1);rand(1,1)];
    
    % Sorting matrix %
    % extract v to compute z %
    v = q(ind_2);
    v = v(:);
    z = v.'*inputs;
    diff_z = diff(sort(z));
    mean_diff_z_sq = Scale*mean(diff_z)^2;
    
    [~,sortI] = sort(z);
    S = zeros(n3,n3); % sorting matrix
    for s=1:n3
        S(s,sortI(s)) = 1;
    end
    
    % weighting matrices %
    W1 = F1;
    W2 = F2;
    D2mean = F3./mean_diff_z_sq;
    W3 = S^-1*D2mean*S;
%     W3 = S^-1*F3*S;
W3 = W3'*W3;
    
%     qnew = iterate3D(X_err,q,Q1*chol(W1),Q2*chol(W2),Q3*chol(W3),theta); % first iteration on error
    qnew = iterate3D_W(X_err,q,W1,W2,W3,theta);
    
    % Sorting matrix %
    % extract v to compute z %
    v = qnew(ind_2);
    v = v(:);
    z = v.'*inputs;
    diff_z = diff(sort(z));
    mean_diff_z_sq = Scale*mean(diff_z)^2;
    
    [~,sortI] = sort(z);
    S = zeros(n3,n3); % sorting matrix
    for s=1:n3
        S(s,sortI(s)) = 1;
    end
    
    % wieghing matrices %
    W1new = F1;
    W2new = F2;
    D2mean = F3./mean_diff_z_sq;
    W3new = S^-1*D2mean*S;
%     W3new = S^-1*F3*S;
W3new = W3new'*W3new;
    
   
    k=1; % iteration number
    store_k_loops(k,i) = norm(q-qnew);
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            disp(['Maximum number of iterations attained, error equals ',num2str(norm(q-qnew))])
            break
        end
        q = qnew;
        W1 = W1new;
        W2 = W2new;
        W3 = W3new;
        
%         qnew = iterate3D(X_err,q,Q1*chol(W1),Q2*chol(W2),Q3*chol(W3),theta);
        qnew = iterate3D_W(X_err,q,W1,W2,W3,theta);
        
        % Sorting matrix %
        % extract v to compute z %
        v = qnew(ind_2);
        v = v(:);
        z = v.'*inputs;
        diff_z = diff(sort(z));
        mean_diff_z_sq = Scale*mean(diff_z)^2;
  
        [~,sortI] = sort(z);
        S = zeros(n3,n3); % sorting matrix
        for s=1:n3
            S(s,sortI(s)) = 1;
        end
        
        % wieghing matrices %
        W1new = F1;
        W2new = F2;
        D2mean = F3./mean_diff_z_sq;
        W3new = S^-1*D2mean*S;  
%         W3new = S^-1*F3*S;
        W3new = W3new'*W3new;
        k = k+1;
        store_k_loops(k,i) = norm(q-qnew);
        
        if(mod(k,100) == 0)
            disp(['Singular value ' num2str(i) '/' num2str(max_level) ' - iter ' num2str(k) '/' num2str(maxiter)])
        end
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
    
    Approx=Approx + rank1tensor(gains(i),i,i,i,x1,x2,x3);
    X_err = X-Approx;
end


% relative final error %
disp(['----------------------------------------'])
disp(['Relative final error: ' num2str(frob(X_err)/frob(X))])
disp(['----------------------------------------'])



%%% Orthonormalize the basis vectors found for each vector space
%%%   --> may like to skip this step
%x1 = gs_orth(x1);
%x2 = gs_orth(x2);
%x3 = gs_orth(x3);



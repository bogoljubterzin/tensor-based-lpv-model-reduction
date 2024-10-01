function [x1,x2,x3,gains,iterations,Xa] = TSVD3D_SW(X,num_levels)
% SYNTAX:
%
% [x1,x2,x3,result,iterations,Xa] = TSVD3D_SW(X,num_levels)
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
%     x1,x2,x3:   orthonormal columns that define the singular vectors
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
%
% VERSION:
%     Created August 2019
%     Last update: August 21, 2019

%%%% Initialization
% define maximum number of iterations
maxiter = 10000;    

% define tensor dimensions (3rd order only)
[n1,n2,n3] = size(X);

% define index ranges:
ind_1=1:n1; 
ind_2=(n1+1):(n1+n2);
ind_3=(n1+n2+1):(n1+n2+n3);

% define initial projection matrices:
Q1 = eye(n1);
Q2 = eye(n2);
Q3 = eye(n3);

% initialize singular vectors:
x1 = [];
x2 = [];
x3 = [];

% initialize gain vector:
gains = zeros(num_levels,1);

% initialize iteration vector:
iterations = zeros(num_levels,1);


%%%% Computation of TSVD
% Compute i-th level TSVD for all i <= num_levels
for i=1:num_levels
    % initialize q
    iterationconstant = rand(1,1);  
    q = [rand(n1,1);rand(n2,1);rand(n3,1);iterationconstant];
    qnew = iterate3D(X,q,Q1,Q2,Q3);  % run first iteration
    
    k=1;   % iteration number
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            disp(['Maximum number of iterations attained at level ', num2str(i) ])
            break
        end
        q = qnew;
        qnew = iterate3D(X,q,Q1,Q2,Q3);
        k = k+1;
    end
    disp(['Iteration at level ', num2str(i) , ' completed....'])
    
    % update x1,x2,x3,iterations,gains,Q1,Q2,Q3
    iterations(i)=k;
    gains(i) = qnew(end);
    x1hat = qnew(ind_1);
    x2hat = qnew(ind_2);
    x3hat = qnew(ind_3);
    
    % Define singular vectors at this level
    x1(:,i) = Q1*x1hat./norm(Q1*x1hat); % normalize
    x2(:,i) = Q2*x2hat./norm(Q2*x2hat); % normalize
    x3(:,i) = Q3*x3hat./norm(Q3*x3hat); % normalize 
    
    % Define projection matrices for next level
    Q1 = eye(n1) - x1*inv(x1'*x1)*x1';   
    Q2 = eye(n2) - x2*inv(x2'*x2)*x2';
    Q3 = eye(n3) - x3*inv(x3'*x3)*x3';
    
    % Adjust X with projection matrices
    X = double(ttm(tensor(X),{Q1,Q2,Q3},[1,2,3]));
end

%%% Define approximation
Xa=zeros(n1,n2,n3);
for i = 1:num_levels
   Approx = (ttt(tensor(x1(:,i)),tensor(x2(:,i))));
   Approx = squeeze(Approx);
   Approx = ttt(Approx,tensor(x3(:,i)));
   Approx = squeeze(Approx);
   Approx = gains(i).*double(Approx);
   Xa=Xa+Approx;
end;




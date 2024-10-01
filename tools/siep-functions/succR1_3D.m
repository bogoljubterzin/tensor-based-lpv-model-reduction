function [z,t,y,result,iter_number] = succR1_3D(G,num_iter)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computation of Successive rank-one approximations for an order-3
%%% tensor. 
%%% Arguments are 'G' and 'num_iter'. 'G' is the multi-dimensional array
%%% for which the decomposition is computed, G should be a double, 
%%% and 'num_iter' is the number of singular values that are
%%% computed.
%%%
%%% This function uses 'gs_orth' to compute a gram-schmidt orthonormalization
%%% of the basis functions computed. It also uses 'rank1tensor' to compute
%%% a rank-one tensor from three vectors and scalar. Finally, 'iterate3D'
%%% is used to compute iterations in the fixed-point algorithm.
%%%
%%% The output matrices z,t,y contain orthonormal columns that form the
%%% subspaces. 'result' contains the gains found.
%%%
%%% Created by Femke van Belzen, February 11th,  2009
%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% Initialization
maxiter = 1000;    %maximum number of iterations


%%% Initialize the procedure

[n,m,p] = size(G);

Q1 = eye(n);
Q2 = eye(m);
Q3 = eye(p);


z = [];
t = [];
y = [];
iter_number = [];
%%%% Compute the optimal rank-one approximant

iterationconstant = rand(1,1);  %%%CHANGE THIS
q = [rand(n,1);rand(m,1);rand(p,1);iterationconstant];
qnew = iterate3D(G,q,Q1,Q2,Q3);

k=0;
while (norm(q-qnew)>1e-6)
    if k == maxiter
        fprintf('Max iterations attained, error equals %e \n',norm(q-qnew))
        break
    end
    q = qnew;
    qnew = iterate3D(G,q,Q1,Q2,Q3);
    k = k+1;
end
disp('Iteration 1 completed \n')
iter_number(1) = k;
% update z,t,y,result
zhat = qnew(1:n);
that = qnew((n+1):(n+m));
yhat = qnew((n+m+1):(n+m+p));
z = zhat./norm(zhat);
t = that./norm(that);
y = yhat./norm(yhat);
result = qnew(n+m+p+1);


for i=2:num_iter
    %Update the tensor
    X_from_previous_iteration = rank1tensor(result(i-1),i-1,i-1,i-1,z,t,y);
    G = G - X_from_previous_iteration;
    % initialize q
    iterationconstant = rand(1,1);  %%%CHANGE THIS
    q = [rand(n,1);rand(m,1);rand(p,1);iterationconstant];
    qnew = iterate3D(G,q,Q1,Q2,Q3);
    k=0;
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            fprintf('Max iterations attained, error equals %e \n',norm(q-qnew))
            break
        end
        q = qnew;
        qnew = iterate3D(G,q,Q1,Q2,Q3);
        k = k+1;
    end
    fprintf('Iteration %i completed \n',i)
    iter_number(i) = k;
    % update z,t,y,result
    zhat = qnew(1:n);
    that = qnew((n+1):(n+m));
    yhat = qnew((n+m+1):(n+m+p));
    z(:,i) = zhat./norm(zhat);
    t(:,i) = that./norm(that);
    y(:,i) = yhat./norm(yhat);
    result(i) = qnew(n+m+p+1);
end

result = result';

%%% Orthonormalize the basis vectors found for each vector space
z = gs_orth(z);
t = gs_orth(t);
y = gs_orth(y);


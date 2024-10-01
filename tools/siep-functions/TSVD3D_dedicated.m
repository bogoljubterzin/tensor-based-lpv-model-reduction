function [x,y,z,result,iter_number] = TSVD3D_dedicated(G,num_iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Created by Femke van Belzen, January 28th 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialization
maxiter = 10000;    %maximum number of iterations
[n,m,p] = size(G);

Q1 = eye(n);
Q2 = eye(m);
Q3 = eye(p);

x = [];
y = [];
z = [];
result = [];
iter_number = [];

%%%% Start the computations
for i=1:num_iter
    % initialize q
    iterationconstant = rand(1,1);  %%%CHANGE THIS
    q = [rand(n,1);rand(m,1);rand(p,1);iterationconstant];
    qnew = iterate3D(G,q,Q1,Q2,Q3);
    k=0;
    while (norm(q-qnew)>1e-6)
        if k == maxiter
            disp('max iterations attained')
            break
        end
        q = qnew;
        qnew = iterate3D(G,q,Q1,Q2,Q3);
        k = k+1;
    end
    fprintf('Iteration %i completed \n',i)
    iter_number(i)=k;
    % update x,y,z,t,var,result,Q1,Q2,Q3
    xhat = qnew(1:n);
    yhat = qnew((n+1):(n+m));
    zhat = qnew((n+m+1):(n+m+p));
    x(:,i) = Q1*xhat./norm(Q1*xhat);
    y(:,i) = Q2*yhat./norm(Q2*yhat);
    z(:,i) = Q3*zhat./norm(Q3*zhat);
    result(i) = qnew(n+m+p+1);
    Q1 = eye(n) - x*x';
    Q2 = eye(m) - y*y';
    G = double(ttm(tensor(G),{Q1,Q2,Q3},[1,2,3]));
end
result=result';
iter_number = iter_number';





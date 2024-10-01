%% How to get core tensor from TSVD?

T = zeros(2,2,2);
T(1,1,1) = 1; T(2,1,1) = 1; T(1,1,2) = 1;
[w1,w2,w3,sv,iter,Tapp] = succR1_SW(T,2);

t1 = tensor(w1(:,1)); t2 = tensor(w2(:,1)); t3 = tensor(w3(:,1));
A = ttt(t1,t2);
A = squeeze(A);
A = ttt(A,t3);
A = squeeze(A);
A = sv(1).*A
% A = sv(1)*A

% S = TSVcore(T,2)
%%
lmd = 0.5 + 0.5*sqrt(5);
w11 = [lmd/sqrt(lmd^2+1); 1/sqrt(lmd^2+1)];
w1 = [w11, [w11(2); -w11(1)]];
w2 = [[1;0], [0;1]];
w3 = W1_sol;

S = zeros(size(T));

for i = 1:size(w1,2)
    for j = 1:size(w2,2)
        for k = 1:size(w3,2)
            t1 = tensor(w1(:,i));
            t2 = tensor(w2(:,j));
            t3 = tensor(w3(:,k));
            Trank1 = squeeze(ttt(t1,t2));
            Trank1 = squeeze(ttt(Trank1,t3));
            S(i,j,k) = innerprod(tensor(T),Trank1);
        end
    end
end

% a1, b1, c1 matches with the first singular vectors from the solution
% sigma1 = lambda; sigma2 = 0; => R = 1;
% T = sigma_1 * w1 x w2 x w3
% but when I compute it I get a different tensor (stored in A)
% [w1,w2,w3,sv,iter,Tapp] = succR1_SW(double(A),2);

%%

%% Example IV.5

T = zeros(2,2,2);
T(1,1,1) = 2;
T(2,2,2) = 0.5*sqrt(2);
T(1,2,2) = 0.5*sqrt(2);

[w1,w2,w3,sv,iter,Tapp] = TSVD3D_SW(T,2);
S = zeros(size(T));

for i = 1:size(w1,2)
    for j = 1:size(w2,2)
        for k = 1:size(w3,2)
            t1 = tensor(w1(:,i));
            t2 = tensor(w2(:,j));
            t3 = tensor(w3(:,k));
            Trank1 = squeeze(ttt(t1,t2));
            Trank1 = squeeze(ttt(Trank1,t3));
            S(i,j,k) = innerprod(tensor(T),Trank1);
        end
    end
end

t1 = tensor(w1(:,1));
t2 = tensor(w2(:,1));
t3 = tensor(w3(:,1));
Trank1 = squeeze(ttt(t1,t2));
Trank1 = sv(1).*squeeze(ttt(Trank1,t3));

t1 = tensor(w1(:,2));
t2 = tensor(w2(:,2));
t3 = tensor(w3(:,2));
Trank2 = squeeze(ttt(t1,t2));
Trank2 = sv(2).*squeeze(ttt(Trank2,t3));

S = TSVcore(T,2);

%%
rng(1)
A = rand(2,4,3);
[w1,w2,w3,sv,iter,Tapp] = TSVD3D_SW(A,min(size(A)));

%%
T = rand(4,4,4);
S = TSVcore(T,4);

% tensSS - tensSS model of MSD2

tensSS = tensSS_dt{1}

n = 4;
Wn_cell = cell(n,1);

Wn = ReachabilityTensors(tensSS.A,tensSS.B,n);
num_levels = 10;
%%
V = [];
Z = [];

for i = 1:n
    [x,gains,iterations,Xa] = TSVDND_SW(Wn{i},num_levels);
    V = [V, x{1}];
    for j = 2:i+1
        Z = [Z, x{j}];
    end
    tsvd_info{i}.x = x;
    tsvd_info{i}.gains = gains;
    tsvd_info{i}.iterations = iterations;
end

%%
V = orth(V);
% Z = orth(Z);
% Z = null(Z');
Z = eye(tensSS_dt{1}.Np);
W = V;
constantTerm = 1;
tensSSr = tensSS.PetrovGalerkinLPV(W,V,Z,constantTerm)
[~,~,yout,~] = tensSSr.simulateSS(utrain(t),t,V'*x0{1},constantTerm) 



%%
[lpvssR,Tx] = lpvmmred(lpvss_dt{1}, 4, 'R');
tensSSrMM = tensSS.PetrovGalerkinLPV(Tx,Tx,eye(lpvssR.Np),1)
[~,~,youtMM,~] = tensSSrMM.simulateSS(utrain(t),t,Tx'*x0{1},1) 



%%


norm(yTrain{3} - youtMM)
norm(yTrain{3} - yout)

%%

figure,
plot(t,yTrain{3}, t, youtMM)
%%
theta = 1;
[x,gains,iterations,Xa] = TSVDND_SW(Wn{2},num_levels,theta)



%%

AB = ttt(tensSS.A,tensSS.B,2,1)

%%

T = rand(6,7,4,5);
[x,gains,iterations,Xa] = TSVDND_SW(T,2,theta)


function Wn = ReachabilityTensors(A,B,n)

    Wn = cell(n,1);
    
    Wn{1} = permute(B,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Wn{i} = ttt(Ap,Wn{i-1},3,1);
    end

end



% Developing tensor-based moment matching

load('benchmark_models.mat')
lpvss_dt
tensSS_dt

msd = 3; % MSD3: Nx = 100; Np = 99;

lpvss = lpvss_dt{msd};
tensSS = tensSS_dt{msd};
eta_map = eta{msd}.map;
N = 4; % moments/ time horizon
%% LPV-SS moment matching

[lpvssR,Tx] = lpvmmred(lpvss, N, 'R');
etar_map = @(x,u) eta_map(Tx*x,u);
[~, ~, yMM] = affineLpvSim(lpvssR, etar_map, utrain(t), t, Tx'*x0{msd});

%%
constantTerm = 0;
Z = eye(tensSS.Np+1);
tensSSR = tensSS.PetrovGalerkinLPV(Tx,Tx,Z,constantTerm);
[~,~,yMMLPVGalerkin,~] = tensSSR.simulateSS(utrain(t),t,Tx'*x0{msd},constantTerm);

% I performed LPVPetrovGalerkin from tensor based just to check if it is
% implemented correctly and if it gives the same results as the lpvmmred.
% The models are the same, but there is some slight difference in
% simulation 

%%
figure(1),
plot(t,yTrain{msd},t,yMM, t, yMMLPVGalerkin); 

%%
N = 3;
Wn_cell = cell(N,1);
Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
num_levels = 3; % this cant be too high (for msd3) because of memory
                % issues (max n = 3); for other models it is ok to be higher
%%
Tstates = [];
Tsched = [];
decomp = 'hosvd';

for n = 1:N 
    if lower(decomp) == "tsvd"
        % for every time horizon perform TSVD
        [x,gains,iterations,Xa] = TSVDND_SW(Wn{n},num_levels);
    
        % separate singular vectors: 
    
        % for the state projection V, we gather the
        % vectors from first dimension x{1}
        Tstates = [Tstates, x{1}];
    
        % for the scheduling projection Z, we gather vectors from dimensions
        % x{2},..., x{n}
        for j = 2:n+1
            Tsched = [Tsched, x{j}];
            % Tsched{j-1} = null(x{j}.')
        end
        % save info from decomposition from every level
        tsvd_info{n}.x = x;
        tsvd_info{n}.gains = gains;
        tsvd_info{n}.iterations = iterations;

    elseif lower(decomp) == 'hosvd'
        T = hosvd(tensor(Wn{n}),2*sqrt(eps));
        Tstates = [Tstates, T{1}];
        for j = 2:n+1
            Tsched = [Tsched, T{j}];
            % Tsched = [Tsched, null(T{j}.')];

        end
    end 
end

%% Is V ok as a state projection? 
V = orth(Tstates);
Z = eye(tensSS.Np+1);
W = V;
constantTerm = 0;
tensSSR1 = tensSS.PetrovGalerkinLPV(W,V,Z,constantTerm)
%%
[~,~,yout1,~] = tensSSR1.simulateSS(utrain(t),t,V'*x0{msd},constantTerm) 

% NaN.. lets try with another model where we can build more tensors Wn

%% Joint-Reduction

V = orth(Tstates);

Z = orth(Tsched);
% Z = null(Tsched'); % maybe a different strategy for Z?


W = V;
constantTerm = 0;
tensSSR2 = tensSS.PetrovGalerkinLPV(W,V,Z,constantTerm)

%%
[~,~,yout2,~] = tensSSR2.simulateSS(utrain(t),t,V'*x0{msd},constantTerm) 

%%
figure(3),
plot(t,yTrain{msd},t,yout2); 

%%

e1 = yTrain{msd} - yout2;
%% Local Functions
function Wn = ReachabilityTensors(A,B,n)

    Wn = cell(n,1);
    
    Wn{1} = permute(B,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Wn{i} = ttt(Ap,Wn{i-1},3,1);
    end

end

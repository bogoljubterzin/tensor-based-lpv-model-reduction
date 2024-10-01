% Developing tensor-based moment matching
l2err = @(y,ybar) norm(y-ybar)/norm(y-mean(y-ybar))*100;

load('benchmark_models.mat')
lpvss_dt
tensSS_dt

msd = 1; % MSD1: Nx = 10; Np = 9;

lpvss = lpvss_dt{msd};
tensSS = tensSS_dt{msd};
eta_map = eta{msd}.map;
%% LPV-SS moment matching
N = 3; % moments/ time horizon
tic
[lpvssR,Tx] = lpvmmred(lpvss, N, 'T');
tMM = toc;
etar_map = @(x,u) eta_map(Tx*x,u);
[~, ~, yMM] = affineLpvSim(lpvssR, etar_map, utrain(t), t, Tx'*x0{msd});
[~, ~, yMMTest] = affineLpvSim(lpvssR, etar_map, utest(t), t, Tx'*x0Test{msd});


%%
tic
decomp = 'hosvd';
N = 7;  % this cant be too high (for msd3) because of memory
        % issues (max N = 3); for other models it is ok to be higher
Wn_cell = cell(N,1);
Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);


Tstates = [];
Tsched = [];
Tstatesq = [];
Tschedq = [];
for n = 1:N

    if lower(decomp) == "tsvd"
        % for every time horizon perform TSVD
        sizeWn = size(Wn{n});
        if size(Wn{n},n+2) == 1
            Wn{n} = reshape(Wn{n}, sizeWn(1:end-1));
        end
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
        end

        Tq = hosvd(tensor(Qn{n}),2*sqrt(eps));
        Tstatesq = [Tstatesq, Tq{n+2}];
        for j = 2:n+1
            Tschedq = [Tschedq, Tq{j}];
        end
    end
end

%% method 1: Wn
% rx = 6;
% rp = 5;

V1 = orth(Tstates);
Z1 = orth(Tsched);
W1 = V1;

%%
constantTerm = 0;
tensSSR = tensSS.PetrovGalerkinLPV(W1,V1,Z1,constantTerm);
tTMMR = toc

%%
[~,~,yout2,~] = tensSSR.simulateSS(utrain(t),t,V1.'*x0{msd},constantTerm);
[~,~,yout4,~] = tensSSR.simulateSS(utest(t),t,V1.'*x0Test{msd},constantTerm);

e1 = l2err(yTrain{1},yMM)
e2 = l2err(yTrain{1}, yout2)
e4 = l2err(yTest{1},yout4)
%% method 2: Qn
V2 = orth(Tstatesq);
Z2 = orth(Tschedq);
W2 = V2;

constantTerm = 0;
tensSSRq = tensSS.PetrovGalerkinLPV(W2,V2,Z2,constantTerm);
tTMMq = toc;


%%
[~,~,yout2q,~] = tensSSRq.simulateSS(utrain(t),t,V2.'*x0{msd},constantTerm);
[~,~,yout4q,~] = tensSSRq.simulateSS(utest(t),t,V2.'*x0Test{msd},constantTerm);

e2q = l2err(yTrain{1}, yout2q)
e4q = l2err(yTest{1}, yout4q)
%%
FigMSD1_MM = figure(5)
clf(FigMSD1_MM)
plot(t, yTrain{1}, 'LineWidth',2,'LineStyle','--'); hold on; grid on;
plot(t,yMM,'LineWidth',1.5);
plot(t, yout2,'LineWidth',1.5);
xlabel("Time (s)");
ylabel("y (m)");
legend(sprintf("FOM (N_x = %d, Np = %d)",lpvss.Nx,lpvss.Np),sprintf("LPV MM (R_x = %d, R_p = %d)",lpvssR.Nx, lpvssR.Np),sprintf("Tensor MM (R_x = %d, R_p = %d)",tensSSR.Nx, tensSSR.Np));

%%
exportgraphics(FigMSD1_MM,'FigMSD1_MM.pdf')
%%

FigMSD1_MMwq = figure(7)
clf(FigMSD1_MMwq)
plot(t, yTrain{1}, 'LineWidth',2,'LineStyle','--'); hold on; grid on;
plot(t,yout2,'LineWidth',1.5);
plot(t, yout2q,'LineWidth',1.5);
xlabel("Time (s)");
ylabel("y (m)");
legend("FOM",'Wn (case 1)','Qn (case 2)');
exportgraphics(FigMSD1_MMwq,'FigMSD1_MMwq.pdf')

%%
FigMSD1_MM_test = figure(6)
clf(FigMSD1_MM_test)
plot(t, yTest{1}, 'LineWidth',2,'LineStyle','--'); hold on; grid on;
plot(t,yMMTest,'LineWidth',1);
plot(t, yout4,'LineWidth',1);
plot(t,yout4q,'LineWidth',1)
xlabel("Time (s)");
ylabel("y (m)");
legend("FOM",'LPV MM','Tensor MM (case 1)', 'Tensor MM (case 2)');
exportgraphics(FigMSD1_MM_test,'FigMSD1_MM_test.pdf')

%% Method 3: Hankel matrix


Rn = orth(Tstates);
Nn = orth(Tstatesq);
rank(Nn'*Rn,eps^2)
[N,R] = qr(Nn'*Rn);

V3 = Rn*R;
W3 = (Nn*inv(N));
% 
% N*R
% Nn'*Rn

Z3 = Tsched;

tensSSR3 = tensSS.PetrovGalerkinLPV(W3,V3,Z3,constantTerm);
%%
[~,~,yout2opt,~] = tensSSR3.simulateSS(utrain(t),t,V3'*x0{msd},constantTerm);
[~,~,yout4opt,~] = tensSSR3.simulateSS(utest(t),t,V3'*x0Test{msd},constantTerm);

%% Local Functions

function Qn = ObservabilityTensors(A,C,n)
    Qn = cell(n,1);
    Qn{1} = permute(C,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Qn{i} = ttt(Qn{i-1},Ap,i+1,1);
    end
end
function Wn = ReachabilityTensors(A,B,n)

    Wn = cell(n,1);
    
    Wn{1} = permute(B,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Wn{i} = ttt(Ap,Wn{i-1},3,1);
    end

end


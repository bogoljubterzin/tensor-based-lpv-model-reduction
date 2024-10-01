%% Analyzing tensor-based moment matching method (using benchmark MSD1)

% Questions:
% 1) What happens if we keep the number of states and schedluing the same
% as of the original model? 
% 2) How to choose Rx and Rp?
% To answer question 2) we will provide singular values, ||T-Tr||_F/||T||_F
% relative frobenious norm of the error tensor, and ||y-yr||/||y|| relative
% L2 norm of the response error singal

load('benchmark_models.mat')

msd = 1; % MSD1: Nx = 10; Np = 9;

lpvss = lpvss_dt{msd};
tensSS = tensSS_dt{msd};
eta_map = eta{msd}.map;

x01 = x0{1}; % old value saved
x02 = zeros(tensSS.Nx,1);

[~,~,y1,p1] = tensSS.simulateSS(utrain(t),t,x01,1);
yTrain{msd} = y1;

[~,x2,y2,p2] = tensSS.simulateSS(utrain(t),t,x02,1);
yTrain{msd} = y2;

[~,x2test,y2test,p2test] = tensSS.simulateSS(utest(t),t,x02,1);


%%
N = 8;  % this cant be too high (for msd3) because of memory
        % issues (max N = 3); for other models it is ok to be higher

Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);

Wn_hosvd = cell(N,1);
Qn_hosvd = cell(N,1);

Wn_tsvd = cell(N,1);
Qn_tsvd = cell(N,1);

Tstates = cell(N,1);
Tsched = cell(N,1);
Tstatesq = cell(N,1);
Tschedq = cell(N,1);
for i = 1:N
    Tstates{i} = [];
    Tstates{i} = [];
    Tsched{i} = [];
    Tschedq{i} = [];
end
%%
decomp = 'hosvd';
U1 = cell(N,1);

tol = 1e-6;
for n = 1:N
    if strcmp(decomp,"hosvd") == 1
        T = hosvd(tensor(Wn{n}),tol);
        Wn_hosvd{n} = T;
        Tstates{n} = [Tstates{n}, T{1}];
        U1{n} = T{1};
        for j = 2:n+1
            Tsched{n} = [Tsched{n}, T{j}];
        end
    
        Tq = hosvd(tensor(Qn{n}),2*sqrt(eps));
        Qn_hosvd{n} = Tq;
        Tstatesq{n} = [Tstatesq{n}, Tq{n+2}];
        for j = 2:n+1
            Tschedq{n} = [Tschedq{n}, Tq{j}];
        end
    elseif strcmp(decomp,"tsvd") == 1
        K = min(modrank(Wn{n}));
        [x,gains,iterations,Xa] = TSVDND_SW(Wn{n},K);
        Wn_tsvd{n}.x = x;
        Wn_tsvd{n}.gains = gains;
        Wn_tsvd{n}.iterations = iterations;
        Wn_tsvd{n}.Xa = Xa;

        % for the state projection V, we gather the
        % vectors from first dimension x{1}
        Tstates{n} = [Tstates{n}, x{1}];
    
        % for the scheduling projection Z, we gather vectors from dimensions
        % x{2},..., x{n}
        for j = 2:n+1
            Tsched{n} = [Tsched{n}, x{j}];
            % Tsched{j-1} = null(x{j}.')
        end

        K = min(modrank(Qn{n}));
        [x,gains,iterations,Xa] = TSVDND_SW(Qn{n},K);
        Qn_tsvd{n}.x = x;
        Qn_tsvd{n}.gains = gains;
        Qn_tsvd{n}.iterations = iterations;
        Qn_tsvd{n}.Xa = Xa;
        
        Tstatesq{n} = [Tstatesq{n}, x{n+2}];
        for j = 2:n+1
            Tschedq{n} = [Tschedq{n}, x{j}];
        end
    end
end


%%
Vr = cell(N,1);
Zr = cell(N,1);
Vq = cell(N,1);
Zq = cell(N,1);
Wr = cell(N,1);
Wq = cell(N,1);

Vr{1} = orth(Tstates{1});
Zr{1} = orth(Tsched{1});
Vq{1} = orth(Tstatesq{1});
Zq{1} = orth(Tschedq{1});

for i = 2:N
    Vr{i} = orth([Vr{i-1}, Tstates{i}]);
    Wr{i} = Vr{i};
    Zr{i} = orth([Zr{i-1}, Tsched{i}]);
    Vq{i} = orth([Vq{i-1}, Tstatesq{i}]);
    Wq{i} = Vq{i};
    Zq{i} = orth([Zq{i-1}, Tschedq{i}]);
end

% Now we have the matrices Tx and Tp for every horizon (for both
% reachability and observability case)
% Let's see if for some N we can use the optimal case
opt_case = zeros(N,1);
for i = 1:N
    Rn = Vr{i};
    Nn = Vq{i};
    rx = size(Rn,2);
    if rank(Nn'*Rn) == rx
        opt_case(i) = true;
        [N1,R] = qr(Nn'*Rn);
        Vh{i} = Rn/R;
        Wh{i} = Nn'\N1;
        Zh{i} = orth([Zr{i},Zq{i}]);
    end
end

%%
n=7;
sysrEq = tensSS.PetrovGalerkinLPV(Wh{n},Vh{n},Zh{n},0);
[~,~,yeq_train,pr2] = sysrEq.simulateSS(unew1,t,Vh{n}'*x02,constantTerm);
[~,~,yeq_test,pr2] = sysrEq.simulateSS(unew2,t,Vh{n}'*x02,constantTerm);


[lpvssR,Tx1] = lpvmmred(lpvss, n-1, 'T');
eta1_mm = @(x,u) eta_msd1.map(Tx1*x,u);
[~, x2mm, y2mm] = affineLpvSim(lpvssR, eta1_mm, unew1, t, Tx1'*x02);
[~, x2mmtest, y2mmtest] = affineLpvSim(lpvssR, eta1_mm, unew2, t, Tx1'*x02);
%%
MSDEq1 = figure(123)
clf(MSDEq1);
plot(t,y2,'LineStyle','--','LineWidth',2);hold on; grid on;
plot(t,yeq_train,'LineStyle',':','LineWidth',1.2)
plot(t, y2mm,'LineStyle','-.','LineWidth',1.2)
legend(sprintf("FOM, Nx = %d, Np = %d",tensSS.Nx,tensSS.Np), ...
    sprintf("TensMM, Nx = %d, Np = %d",sysrEq.Nx,sysrEq.Np), ...
    sprintf("LPV MM, Nx = %d, Np = %d",lpvssR.Nx,lpvssR.Np))
title("Equivalent transformation for large $N$")

MSDEq2 = figure(124)
plot(t,y2test,t,yeq_test,t,y2mmtest)

%%
constantTerm = 0;
tensSSr = cell(N,1);
yr = cell(N,1);
pr = cell(N,1);
for i = 2:N
    tensSSr{i} = tensSS.PetrovGalerkinLPV(Wr{i},Vr{i},Zr{i},constantTerm);
    [~,~,yr2,pr2] = tensSSr{i}.simulateSS(utrain(t),t,Vr{i}'*x02,constantTerm);
    yr{i} = yr2;
    pr{i} = pr2;
end


% %%
% erN_tsvd = zeros(N,1);
% for i = 2:N
%     erN_tsvd(i) = norm(y2-yr{i})./norm(y2);
% end

erN_hosvd = zeros(N,1);
for i = 2:N
    erN_hosvd(i) = norm(y2-yr{i})./norm(y2);
end
%%
Report_DiffNError = figure(21);
clf(Report_DiffNError)
plot(1:N-1, erN_tsvd(2:N),'ro','LineWidth',2); grid on;  hold on;
plot(1:N-1,erN_hosvd(2:N),'bo','LineWidth',2)
xlabel("Horizon N (-)"); ylabel("Relative L2 error e (-)")
legend("TSVD","HOSVD");
title("Relative error - Tensor MM (MSD1)")

%%
set(0, 'DefaultLineLineWidth', 1.2);

clrs = lines(4);
Report_DifferentN = figure(10),
clf(Report_DifferentN);

tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot(t,y2,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yr{2},'Color',clrs(2,:))
plot(t,yr{4},'Color',clrs(3,:))
plot(t,yr{6},'Color',clrs(4,:))
% plot(t,yrN{5},'Color',clrs(6,:))
legend("FOM", ...
    sprintf("N=1 (Rx = %d, Rp = %d)",tensSSr{2}.Nx,tensSSr{2}.Np), ...
    sprintf("N=3 (Rx = %d, Rp = %d)",tensSSr{4}.Nx,tensSSr{4}.Np), ...
    sprintf("N=5 (Rx = %d, Rp = %d)",tensSSr{6}.Nx,tensSSr{6}.Np))
sgtitle("Different horizon N - Tensor MM (case R) TSVD")
xlabel("Time (s)"); ylabel("y (m)");

nexttile;
semilogy(t,abs(y2-yr{2})./norm(y2),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2-yr{4})./norm(y2),'Color',clrs(3,:));
semilogy(t,abs(y2-yr{6})./norm(y2),'Color',clrs(4,:));
xlabel("Time (s)"); ylabel("e (m)");
ylim([1e-8,1e0])
legend(sprintf("N=1 (Rx = %d, Rp = %d)",tensSSr{2}.Nx,tensSSr{2}.Np), ...
    sprintf("N=3 (Rx = %d, Rp = %d)",tensSSr{4}.Nx,tensSSr{4}.Np), ...
    sprintf("N=5 (Rx = %d, Rp = %d)",tensSSr{6}.Nx,tensSSr{6}.Np))

%% N=2 R vs O vs H
n = 3;
tensSSR = tensSS.PetrovGalerkinLPV(Wr{n},Vr{n},Zr{n},constantTerm);
tensSSQ = tensSS.PetrovGalerkinLPV(Wq{n},Vq{n},Zq{n},constantTerm);
tensSSH = tensSS.PetrovGalerkinLPV(Wh{n},Vh{n},Zh{n},constantTerm);

[~,~,yR,pr] = tensSSR.simulateSS(utrain(t),t,Vr{n}'*x02,constantTerm);
[~,~,yQ,pq] = tensSSQ.simulateSS(utrain(t),t,Vq{n}'*x02,constantTerm);
[~,~,yH,ph] = tensSSH.simulateSS(utrain(t),t,Vh{n}'*x02,constantTerm);

%%

set(0, 'DefaultLineLineWidth', 1.2);

clrs = lines(4);
Report_DifferentCase = figure(12),
clf(Report_DifferentCase);

tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot(t,y2,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yR,'Color',clrs(2,:))
plot(t,yQ,'Color',clrs(3,:))
plot(t,yH,'Color',clrs(4,:))
% plot(t,yrN{5},'Color',clrs(6,:))
legend("FOM", ...
    sprintf("case R (Rx = %d, Rp = %d)",tensSSR.Nx,tensSSR.Np), ...
    sprintf("case Q (Rx = %d, Rp = %d)",tensSSQ.Nx,tensSSQ.Np), ...
    sprintf("case H (Rx = %d, Rp = %d)",tensSSH.Nx,tensSSH.Np))
sgtitle("Tensor MM (different cases)")
xlabel("Time (s)"); ylabel("y (m)");

nexttile;
semilogy(t,abs(y2-yR)./norm(y2),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2-yQ)./norm(y2),'Color',clrs(3,:));
semilogy(t,abs(y2-yH)./norm(y2),'Color',clrs(4,:));
xlabel("Time (s)"); ylabel("e (m)");
ylim([1e-7,1e-1])
legend(sprintf("case R (Rx = %d, Rp = %d)",tensSSR.Nx,tensSSR.Np), ...
    sprintf("case Q (Rx = %d, Rp = %d)",tensSSQ.Nx,tensSSQ.Np), ...
    sprintf("case H (Rx = %d, Rp = %d)",tensSSH.Nx,tensSSH.Np))

%% TensSS vs LPVMM (+PCA?)
tensSSMM = cell(N,1);
TxMM = cell(N,1);
yrMM = cell(N,1);
for n = 1:N-1
    [tensSSMMi, TxMMi] = tensSS.lpvTensMMStates(n, 'R');
    tensSSMM{n} = tensSSMMi;
    [~,~,yrMMi,prMMi] = tensSSMMi.simulateSS(utrain(t),t,TxMMi'*x02,1);
    TxMM{n} = TxMMi;
    yrMM{n} = yrMMi;
end
%% Applying PCA red to tensSSMM{2}
horizon = 2;
[lpvssR,Tx1] = lpvmmred(lpvss, horizon, 'R');
eta1_mm = @(x,u) eta_msd1.map(Tx1*x,u);
[~, x1mm, y1mm] = affineLpvSim(lpvssR, eta1_mm, utrain(t), t, Tx1'*x02);
pTraj = p2'; Rp = tensSSr{horizon+1}.Np
[lpvssPCA1, mapPCA1] = lpvpcared(lpvssR, Rp, pTraj,'trajectory');
etar1 = @(x,u) mapPCA1(eta1_mm(x,u).');
[~, x1pca, y1pca] = affineLpvSim(lpvssPCA1, etar1, utrain(t), t, Tx1'*x02);

% test
[~, x1pca_test, y1pcatest] = affineLpvSim(lpvssPCA1, etar1, utest(t), t, Tx1'*x02);
[~,~,yr2test,pr2test] = tensSSr{horizon+1}.simulateSS(utest(t),t,Vr{horizon+1}'*x02,constantTerm);

%%
Report_TensMMvsMMPCA = figure(14);
clf(Report_TensMMvsMMPCA)
tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot(t,y2,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yr{horizon+1},'Color',clrs(2,:))
plot(t,y1pca,'Color',clrs(3,:))
legend("FOM", ...
    sprintf("Tensor MM (Rx = %d, Rp = %d)",tensSSr{horizon+1}.Nx,tensSSr{horizon+1}.Np), ...
    sprintf("State MM + PCA (Rx = %d, Rp = %d)",lpvssPCA1.Nx,lpvssPCA1.Np))
sgtitle("Tensor MM vs State MM + PCA (on training data)")
xlabel("Time (s)"); ylabel("y (m)");

nexttile;
semilogy(t,abs(y2-yr{horizon+1})./norm(y2),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2-y1pca)./norm(y2),'Color',clrs(3,:));
xlabel("Time (s)"); ylabel("e (m)");
ylim([1e-7,1e-1])
legend(sprintf("Tensor MM (Rx = %d, Rp = %d)",tensSSr{horizon+1}.Nx,tensSSr{horizon+1}.Np), ...
    sprintf("State MM + PCA (Rx = %d, Rp = %d)",lpvssPCA1.Nx,lpvssPCA1.Np))

%%
Report_TensMMvsMMPCA_Test = figure(15);
clf(Report_TensMMvsMMPCA_Test)
tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot(t,y2test,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yr2test,'Color',clrs(2,:))
plot(t,y1pcatest,'Color',clrs(3,:))
legend("FOM", ...
    sprintf("Tensor MM (Rx = %d, Rp = %d)",tensSSr{horizon+1}.Nx,tensSSr{horizon+1}.Np), ...
    sprintf("State MM + PCA (Rx = %d, Rp = %d)",lpvssPCA1.Nx,lpvssPCA1.Np))
sgtitle("Tensor MM vs State MM + PCA (on testing data)")
xlabel("Time (s)"); ylabel("y (m)");

nexttile;
semilogy(t,abs(y2test-yr2test)./norm(y2test),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2test-y1pcatest)./norm(y2test),'Color',clrs(3,:));
xlabel("Time (s)"); ylabel("e (m)");
ylim([1e-7,1e-1])
legend(sprintf("Tensor MM (Rx = %d, Rp = %d)",tensSSr{horizon+1}.Nx,tensSSr{horizon+1}.Np), ...
    sprintf("State MM + PCA (Rx = %d, Rp = %d)",lpvssPCA1.Nx,lpvssPCA1.Np))
%%
clrs = lines(4);
Report_TensMMvsStateMM = figure(13),
clf(Report_DifferentCase);
n = 3;
tiledlayout(2,1,'TileSpacing','tight')
nexttile;
plot(t,y2,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yr{n},'Color',clrs(2,:))
plot(t,yrMM{n-1},'Color',clrs(3,:))
% plot(t,yrN{5},'Color',clrs(6,:))
legend("FOM", ...
    sprintf("Tensor MM(Rx = %d, Rp = %d)",tensSSr{3}.Nx,tensSSr{3}.Np), ...
    sprintf("State MM (Rx = %d, Rp = %d)",tensSSMM{2}.Nx,tensSSMM{2}.Np))
sgtitle("Tensor MM vs State MM (N=2)")
xlabel("Time (s)"); ylabel("y (m)");

nexttile;
semilogy(t,abs(y2-yr{n})./norm(y2),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2-yrMM{n})./norm(y2),'Color',clrs(3,:));
xlabel("Time (s)"); ylabel("e (m)");
ylim([1e-7,1e-1])
legend(sprintf("Tensor MM(Rx = %d, Rp = %d)",tensSSr{3}.Nx,tensSSr{3}.Np), ...
    sprintf("State MM (Rx = %d, Rp = %d)",tensSSMM{2}.Nx,tensSSMM{2}.Np))

%%
er = zeros(N,1);
for i = 2:N
    er(i) = norm(y2-yr{i},2)/norm(y2,2);
    erMM(i) = norm(y2-yrMM{n-1})/norm(y2,2);
end

Report_DiffNError = figure(11);
clf(Report_DiffNError)
plot(1:N-1, [er(2:N)],'ro','LineWidth',2); grid on;  hold on;
plot(1:N-1,[erMMN{:}],'bo','LineWidth',2)
xlabel("Horizon N (-)"); ylabel("Relative L2 error e (-)")
legend("Joint MM","State MM");

%% Lets answer the first question
% 1) What happens if we keep the number of states and schedluing the same
% as of the original model? 
% Taking horizon N = 6 is enough to provide us full transformation matrices
% (case 1: Txr{6}, Tpr{6}, case 2: Txq{6}, Tpq{6})

n = 5;
V1 = Vr{n};
W1 = V1;
Z1 = Zr{n};
% Z1 = eye(tensSS.Np+1);

V2 = Vq{n};
W2 = V2;
Z2 = Zq{n};

Rn = Vr{n};
Nn = Vq{n};
[N1,R] = qr(Nn'*Rn);
V3 = Rn*inv(R);
W3 = Nn'\N1;
Z3 = Z2;
%%
constantTerm = 0; 
% Note about the meaning of constant term: 
% Original model is affine: A(p) = A0 + A1p1 + ... meaning that there is a
% constant term A0 inside of it that is not scheduling dependent.
% When I perform the reduction like this, constant parts of the system are still used
% to form reachability and observability spaces Wn and Qn. The scheduling
% transformation ends up being Tp of dimension (np+1) x rp and the new
% scheduling map takes the original (1; p(t)) and transforms it into
% pnew(t) = (pnew_1(t) ... pnew_rp(t)). The resulting matrices do not have
% a constant term!

tensSSr1 = tensSS.PetrovGalerkinLPV(W1,V1,Z1,constantTerm);
tensSSr2 = tensSS.PetrovGalerkinLPV(W2,V2,Z2,constantTerm);
tensSSr3 = tensSS.PetrovGalerkinLPV(W3,V3,Z3,constantTerm);


[~,~,yr1r,pr1r] = tensSSr1.simulateSS(utrain(t),t,V1'*x01,constantTerm); 
[~,~,yr1q,pr1q] = tensSSr2.simulateSS(utrain(t),t,V2'*x01,constantTerm);
[~,~,yr1opt,pr1opt] = tensSSr3.simulateSS(utrain(t),t,V3'*x01,constantTerm);


[~,~,yr2r,pr2r] = tensSSr1.simulateSS(utrain(t),t,V1'*x02,constantTerm);
[~,~,yr2q,pr2q] = tensSSr2.simulateSS(utrain(t),t,V2'*x02,constantTerm);
[~,~,yr2opt,pr2opt] = tensSSr3.simulateSS(utrain(t),t,V3'*x02,constantTerm);


%%

e1r = norm(y1-yr1r,2)/norm(y1,2);
e1q = norm(y1-yr1q,2)/norm(y1,2);
e1opt = norm(y1-yr1opt,2)/norm(y1,2);

e2r = norm(y2-yr2r,2)/norm(y2,2);
e2q = norm(y2-yr2q,2)/norm(y2,2);
e2opt = norm(y2-yr2opt,2)/norm(y2,2);

e1 = [e1r, e1q, e1opt];
e2 = [e2r,e2q,e2opt];

[e1star, inde1] = min(e1)
[e2star,inde2] = min(e2)

e1sr = mean(e1)
e2sr = mean(e2)
fprintf("The error signals are e1 = %f, and e2 = %f \n", e1,e2);

% visualize the trajectories
clrs = lines(4);
set(0, 'DefaultLineLineWidth', 1.1);

Fig1 = figure(1);
clf(Fig1)
tiledlayout(2,2,'TileSpacing','tight');
nexttile;
plot(t,y1,'--','Color',clrs(1,:),'LineWidth',1.5); hold on; grid on;
plot(t, yr1r,'Color',clrs(2,:));
plot(t,yr1q,'Color',clrs(3,:))
plot(t,yr1opt,'Color',clrs(4,:));
title("Experiment 1 (non-zero initial cond)");
legend('Original','MM: Case 1','MM: Case 2','MM: Case 3');
nexttile;
plot(t,y2,'Color',clrs(1,:)); hold on; grid on;
plot(t,yr2r,'Color',clrs(2,:));
plot(t,yr2q,'Color',clrs(3,:));
plot(t, yr2opt,'Color',clrs(4,:))
title("Experiment 2 (zero initial cond)");
legend('Original','MM: Case 1','MM: Case 2','MM: Case 3');

nexttile;
semilogy(t,abs(y1-yr1r)./abs(y1),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y1-yr1q)./abs(y1),'Color',clrs(3,:));
semilogy(t,abs(y1-yr1opt)./abs(y1),'Color',clrs(4,:));
title(sprintf("Min L2 norm of error signal: e = %.3f, for case %d", e1star, inde1))


nexttile;
semilogy(t,abs(y2-yr2r)./abs(y2),'Color',clrs(2,:)); grid on; hold on;
semilogy(t,abs(y2-yr2q)./abs(y2),'Color',clrs(3,:));
semilogy(t,abs(y2-yr2opt)./abs(y2),'Color',clrs(4,:));

title(sprintf("Min L2 norm of error signal: e = %.3f, for case %d", e2star, inde2))

exportgraphics(Fig1,'FigTensMM_NotReduced.pdf')
%%
Etsvd = zeros(N,1);
Ehosvd = zeros(N,1);
n = 1;
Etsvd(n) = norm(double(Wn{n})-Wn_tsvd{n}.Xa,'fro');
Ehosvd(n) = norm(double(Wn{n})-double(Wn_hosvd{n}),'fro')

%%
EQtsvd = zeros(N,1);
EQhosvd = zeros(N,1);

for n=1:N
    EQtsvd(n) = norm(double(squeeze(Qn{n}))-squeeze(Qn_tsvd{n}.Xa),'fro')
    EQhosvd(n) = norm(double(Qn{n})-double(Qn_hosvd{n}),'fro')
end

%%
squeeze(Qn{2})

clear sv
Q = Qn{5};
for j = 1:numel(size(Q))
    Qm = tens2mat(Q,j);
    sv{j} = svd(Qm)
end

%%


constantTerm = 0;
for n = 1:6
    V = Vr{n+1};
    W = V;
    Z = Zr{n+1};
    tensSSr = tensSS.PetrovGalerkinLPV(W,V,Z,constantTerm);

    [tensSSMM, TxMM] = tensSS.lpvTensMMStates(n, 'R');

    [~,~,yrMM,prMM] = tensSSMM.simulateSS(utrain(t),t,TxMM'*x02,1); 
    [~,~,yr,pr] = tensSSr.simulateSS(utrain(t),t,V'*x02,constantTerm); 
    
    TxrMM{n} = TxMM;
    tensSSrMM{n} = tensSSMM;

    tensSSrN{n} = tensSSr;
    yrN{n} = yr;
    erN{n} = norm(y-yr)./norm(y);

    yrMMN{n} = yrMM;
    erMMN{n} = norm(y-yrMM)./norm(y);

end
%%
set(0, 'DefaultLineLineWidth', 1.2);

clrs = parula(6);
FigDiffN_MM
FigDiffN_MM = figure(3),
clf(FigDiffN_MM);
plot(t,y,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yrMMN{1},'Color',clrs(2,:))
plot(t,yrMMN{2},'Color',clrs(3,:))
plot(t,yrMMN{3},'Color',clrs(4,:))
plot(t,yrMMN{4},'Color',clrs(5,:))
% plot(t,yrN{5},'Color',clrs(6,:))
legend("FOM",sprintf("N=1 (Rx = %d, Rp = %d)",tensSSrMM{1}.Nx,tensSSrMM{1}.Np), ...
    sprintf("N=2 (Rx = %d, Rp = %d)",tensSSrMM{2}.Nx,tensSSrMM{2}.Np), ...
    sprintf("N=3 (Rx = %d, Rp = %d)",tensSSrMM{3}.Nx,tensSSrMM{3}.Np), ...
    sprintf("N=4 (Rx = %d, Rp = %d)",tensSSrMM{4}.Nx,tensSSrMM{4}.Np))
sgtitle("Different $N$ - Moment Matching state reduction")
xlabel("Time (s)"); ylabel("y (m)");


set(0, 'DefaultLineLineWidth', 1.2);

clrs = parula(6);
FigDiffN = figure(4),
clf(FigDiffN);
plot(t,y,'LineWidth',1.5,'LineStyle','--','Color',clrs(1,:)); hold on; grid on;
plot(t,yrN{1},'Color',clrs(2,:))
plot(t,yrN{2},'Color',clrs(3,:))
plot(t,yrN{3},'Color',clrs(4,:))
plot(t,yrN{4},'Color',clrs(5,:))
% plot(t,yrN{5},'Color',clrs(6,:))
legend("FOM",sprintf("N=1 (Rx = %d, Rp = %d)",tensSSrN{1}.Nx,tensSSrN{1}.Np), ...
    sprintf("N=2 (Rx = %d, Rp = %d)",tensSSrN{2}.Nx,tensSSrN{2}.Np), ...
    sprintf("N=3 (Rx = %d, Rp = %d)",tensSSrN{3}.Nx,tensSSrN{3}.Np), ...
    sprintf("N=4 (Rx = %d, Rp = %d)",tensSSrN{4}.Nx,tensSSrN{4}.Np))
sgtitle("Different $N$ - Moment Matching joint reduction")
xlabel("Time (s)"); ylabel("y (m)");
%%
FigEdiffN = figure(5);
plot(1:6, [erN{:}],'ro','LineWidth',2); grid on;  hold on;
plot(1:6,[erMMN{:}],'bo','LineWidth',2)
xlabel("Horizon N (-)"); ylabel("Relative L2 error e (-)")
legend("Joint MM","State MM");
%%
exportgraphics(FigDiffN,'FigDiffNJointRed.pdf');
exportgraphics(FigEdiffN,'FigEdiffN.pdf');
exportgraphics(FigDiffN_MM,'FigDiffNStateRed.pdf');

%% 
TxJointRed = cell(N,1);
TxStateRed = cell(N,1);

for n =1:N
TxJointRed{n} = Vr{n+1}
TxStateRed{n} = TxrMM{n}
end



%%
Wnred = ReachabilityTensors(tensSSrN{2}.A,tensSSrN{2}.B,N);
Qnred = ObservabilityTensors(tensSSrN{2}.A,tensSSrN{2}.C,N);
%%
for n =1:N
q = SameSubspace(TxJointRed{n},TxStateRed{n})

end

%% Lets answer the question: if it is a lossless compression, transform back from rom to fom

% FOM
tensSS

% ROM - sysr


Rn = Vr{5};
Nn = Vq{5};
[N1,R] = qr(Nn'*Rn);
V = Rn/R;
W = Nn'\N1;
Z = Z2;

constantTerm = 0;
sysr = tensSS.PetrovGalerkinLPV(W,V,Z,constantTerm);
[~,xr2opt,yr2opt,pr2opt] = sysr.simulateSS(utrain(t),t,V'*x02,constantTerm);

%  
% Vr = V';
% Wr = W';
% Zr = Z';
% sysReconstructed = sysr.PetrovGalerkinLPV(Wr,Vr,Zr,constantTerm)
% [~,~,yreconstr2opt,prreconstr2opt] = sysReconstructed.simulateSS(utrain(t),t,Vr'*V'*x02,1);


%%
[sysr_statesMM, T_statesMM] = tensSS.lpvTensMMStates(3, 'T');

subspace(T_statesMM,TxJointRed{4})

%%
Wred = ReachabilityTensors(sysr.A,sysr.B,5)
T = hosvd(tensor(Wred{4 }),tol)
%%
for n=1:N
angl(n) = subspace(TxJointRed{n},TxStateRed{n})
end
function q = SameSubspace(A,B)
    A = orth(A);
    B = orth(B);
    if rank(A)~= rank(B)
        q = 0;
        return
    end

    R = [A, B];
    if rank(R) ~= rank(A)
        q = 0;
        return;
    end
    q = 1;
end
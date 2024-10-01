%% Analyzing tensor-based moment matching method (using benchmark MSD1)

% Questions:
% 1) What happens if we keep the number of states and schedluing the same
% as of the original model? 
% 2) How to choose Rx and Rp?
% To answer question 2) we will provide singular values, ||T-Tr||_F/||T||_F
% relative frobenious norm of the error tensor, and ||y-yr||/||y|| relative
% L2 norm of the response error singal

lpvss = lpvss_msd1;
tensSS = msd1;
eta_map = eta_msd1;

x0 = zeros(tensSS.Nx,1);

[~,~,yTrain,pTrain] = tensSS.simulateSS(umsd_train,t,x0,1);
[~,~,yTest,pTest] = tensSS.simulateSS(umsd_test,t,x0,1);
%%

N = 3;  % this cant be too high (for msd3) because of memory
        % issues (max N = 3); for other models it is ok to be higher

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
t_hosvd_r = cell(N,1);
t_hosvd_q = cell(N,1);
t_tsvd_r = cell(N,1);
t_tsvd_q = cell(N,1);
%%

Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);

%%
decomp = 'hosvd';

tol = 1e-6;
for n = 1:N
    disp(n);
    if strcmp(decomp,"hosvd") == 1
        tic
        T = hosvd(tensor(Wn{n}),tol);
        t_hosvd_r{n} = toc;
        Wn_hosvd{n} = T;
        Tstates{n} = [Tstates{n}, T{1}];
        for j = 2:n+1
            Tsched{n} = [Tsched{n}, T{j}];
        end
        tic
        Tq = hosvd(tensor(Qn{n}),2*sqrt(eps));
        t_hosvd_q{n} = toc;
        Qn_hosvd{n} = Tq;
        Tstatesq{n} = [Tstatesq{n}, Tq{n+2}];
        for j = 2:n+1
            Tschedq{n} = [Tschedq{n}, Tq{j}];
        end
    elseif strcmp(decomp,"tsvd") == 1
         
        t_tsvd_r{n} = toc;
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
        tic
        [x,gains,iterations,Xa] = TSVDND_SW(Qn{n},K);
        t_tsvd_q{n} = toc;
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
n=1
norm(double(Wn{n})-double(Wn_tsvd{n}.Xa),'fro')
norm(double(Qn{n})-(Qn_tsvd{n}.Xa),'fro');

%%
Vr = cell(N,1);
Zr = cell(N,1);
Vq = cell(N,1);
Zq = cell(N,1);
Wr = cell(N,1);
Wq = cell(N,1);

Vr{1} = orth(Tstates{1});
Wr{1} = Vr{1};
Zr{1} = orth(Tsched{1});
Vq{1} = orth(Tstatesq{1});
Zq{1} = orth(Tschedq{1});
Wq{1} = Vq{1};
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
mode = 'R';
if mode == 'R'
    V = Vr;N
    W = Wr;
    Z = Zr;
end
if mode == 'Q'
    V = Vq;
    W = Wq;
    Z = Zq;
end
if mode == 'H'
    V = Vh;
    W = Wh;
    Z = Zh;
end

   
e = cell(N,1);
for i = 2:N
    if mode == 'H' && opt_case(i) == 0
        continue
    end
    sysr = tensSS.PetrovGalerkinLPV(W{i},V{i},Z{i},0);
    [~,~,yr_train,~] = sysr.simulateSS(umsd_train,t,V{i}'*x0,constantTerm);
    e{i} = nrmse(yTrain,yr_train);
end

% e_r_hosvd = e;
% e_q_hosvd = e;
e_h_hosvd = e;
% e_r_tsvd = e;
% e_q_tsvd = e;
% e_h_tsvd = e;


%%
yr = cell(N,1);
yrTest = cell(N,1);
sysr = cell(N,1);
for i = 1:4
    sysri = tensSS.PetrovGalerkinLPV(Wr{i},Vr{i},Zr{i},0);
    sysr{i} = sysri;
    [~,~,yri,~] = sysri.simulateSS(umsd_train,t,Vr{i}'*x0,constantTerm);
    yr{i} = yri;

    [~,~,yrTesti,~] = sysr{i}.simulateSS(umsd_test,t,Vr{i}'*x0,constantTerm);
    yrTest{i}= yrTesti;

end

%%

FigWnTSVDMSD1 = figure(141);
clrs = lines(5);
tiledlayout(2,1,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,yTrain,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr{1},'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr{2},'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yr{3},'Color',clrs(4,:),'LineWidth',1.2);
plot(t,yr{4},'Color',clrs(5,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend(sprintf("FOM, $n_x = %d, n_p = %d$",tensSS.Nx,tensSS.Np), ...
       sprintf("$N=0, r_x = %d, r_p = %d$", sysr{1}.Nx,sysr{1}.Np), ...
       sprintf("$N=1, r_x = %d, r_p = %d$", sysr{2}.Nx,sysr{2}.Np), ...
       sprintf("$N=2, r_x = %d, r_p = %d$", sysr{3}.Nx,sysr{3}.Np), ...
       sprintf("$N=3, r_x = %d, r_p = %d$", sysr{4}.Nx,sysr{4}.Np))

nexttile;
plot(t,yTest,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yrTest{1},'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yrTest{2},'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yrTest{3},'Color',clrs(4,:),'LineWidth',1.2);
plot(t,yrTest{4},'Color',clrs(5,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend("FOM","$N=0$","$N=1$","$N=2$","$N=3$")

%%

erhosvd = [e_r_hosvd{:}];
eqhosvd = [e_q_hosvd{:}];
ertsvd = [e_r_tsvd{:}];
eqtsvd = [e_q_tsvd{:}];
%%
FigErrorsMSD1 = figure(145);
clf(FigErrorsMSD1);
tiledlayout(1,1,"TileSpacing",'tight','padding','compact');
plot(0:6,erhosvd,'rx','LineWidth',2); hold on; grid on;
plot(0:6,eqhosvd,'bx','LineWidth',2);
plot(0:6,ertsvd,'ro','LineWidth',2);
plot(0:6,eqtsvd,'bo','LineWidth',2);

%%
FigErrorsMSD1 = figure(145);
clf(FigErrorsMSD1);
tiledlayout(1,1,"TileSpacing",'tight','padding','compact');

% Main plot
plot(0:6,erhosvd,'rx','MarkerSize',12); hold on; grid on;
plot(0:6,eqhosvd,'bx','MarkerSize',12);
plot(0:6,ertsvd,'ro','MarkerSize',12);
plot(0:6,eqtsvd,'bo','MarkerSize',12);
xticks(0:6);
legend('$HOSVD, W_n$','$HOSVD, Q_n$','$TSVD, W_n$','$TSVD, Q_n$')
xlabel("Horizon N"); ylabel("NRMSE")
% Create zoomed-in inset
axes('Position',[0.6 0.6 0.25 0.25]); % Adjust position as necessary
box on; % Add a border around the inset plot
hold on;
plot(0:6,erhosvd,'rx','MarkerSize',12); % Re-plot data for zoomed-in section
plot(0:6,eqhosvd,'bx','MarkerSize',12);
plot(0:6,ertsvd,'ro','MarkerSize',12);
plot(0:6,eqtsvd,'bo','MarkerSize',12);

% Set axis limits for zoomed-in plot
xlim([2.5 6.5]); % Zoom into x-axis from 2 to 5
ylim([-0.1 1.5]); % Adjust y-axis from 0 to 5
grid on;
title('Zoomed');

%% Lets create non minimal representation of msd1

V1 = [eye(tensSS.Nx), eye(tensSS.Nx)];
Z1 = eye(tensSS.Np+1);
W1 = V1;
% 
Ar = double(ttm(tensSS.A,{pinv(W1'*V1)*W1', V1', Z1'},[1,2,3]));
Br = double(ttm(tensSS.B,{pinv(W1'*V1)*W1',Z1'},[1,3]));
Cr = double(ttm(tensSS.C,{V1',Z1'},[2,3]));
Dr = double(ttm(tensSS.D,{Z1'},3));
eta_red = @(xr,u) pinv(Z1)*[1; tensSS.eta_map(V1*xr,u)];
tensMSDNonMinimal = tensorSS(Ar,Br,Cr,Dr,eta_red,tensSS.Ts)

%%
[tout,xNonMin,yNonMin] = tensMSDNonMinimal.simulateSS(umsd_train,t,V1'*x0,0)

%%
n=4;
sysrEq = tensMSDNonMinimal.lpvTensMM(6,'tsvd','R')

%%
% sysrEq = tensSS.PetrovGalerkinLPV(Wr{n},Vr{n},Zr{n},0);
[~,~,yeq_train,pr1] = sysrEq.simulateSS(umsd_train,t,zeros(sysrEq.Nx,1),constantTerm);
[~,~,yeq_test,pr2] = sysrEq.simulateSS(umsd_test,t,zeros(sysrEq.Nx,1),constantTerm);
%% tensMM vs MM+PCA
tic
tensSSr = msd1.PetrovGalerkinLPV(Wr{3},Vr{3},Zr{3},0);
t1 = toc
% [lpvmmSS, Tx1] = tensSS.lpvTensMMStates(2,'R')

tic
[lpvssR,Tx1] = lpvmmred(msd1_dt, 2, 'R');
eta1_mm = @(x,u) eta_msd1.map(Tx1*x,u);

pTraj = p1'; Rp = tensSSr.Np;
[lpvssPCA1, mapPCA1] = lpvpcared(lpvssR, Rp, pTraj,'trajectory');
etar1 = @(x,u) mapPCA1(eta1_mm(x,u).');
t2 = toc
[~, ~, yrmmpca_train] = affineLpvSim(lpvssPCA1, etar1, umsd_train, t, Tx1'*x0);
[~, ~, yrmmpca_test] = affineLpvSim(lpvssPCA1, etar1, umsd_test, t, Tx1'*x0);
% [lpvssR,Tx1] = lpvmmred(lpvss_dt{1}, 2, 'R');

%%
[~,~,yrtmm_train,~] = tensSSr.simulateSS(umsd_train,t,Vr{3}'*x0,constantTerm);
[~,~,yrtmm_test,~] = tensSSr.simulateSS(umsd_test,t,Vr{3}'*x0,constantTerm);


%%

etmm_train = nrmse(yTrain,yrtmm_train)
etmm_test = nrmse(yTest,yrtmm_test)

emmpca_train = nrmse(yTrain,yrmmpca_train)
emmpca_test = nrmse(yTest,yrmmpca_test)

%%
Fig_CompMSD1 = figure(312);
clrs = lines(3);
tiledlayout(2,1,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,yTrain,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yrtmm_train,'Color',clrs(2,:),'LineWidth',2,'LineStyle','-.');
plot(t,yrmmpca_train,'Color',clrs(3,:),'LineWidth',1.2);
legend('FOM',"TMM","MM+PCA",'FontSize',10);
xlabel("Time (s)"); ylabel("$y (m)$");
title("Training simulation")


nexttile;
plot(t,yTest,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yrtmm_test,'Color',clrs(2,:),'LineWidth',2,'LineStyle','-.');
plot(t,yrmmpca_test,'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend('FOM',"TMM","MM+PCA",'FontSize',10);
title("Testing simulation")



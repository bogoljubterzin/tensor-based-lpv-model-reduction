%% Very careful result obtainer
% today's episode: data-based tensor model reduction method which I will
% call POD method
% We should say that there are two types of POD: Galerkin and
% Petrov-Galerkin

% todays star: msd1 model

msd1; % tensorSS representation
msd1_dt; % lpvss representation which we will use to simulate cause its faster

% input signals
figure(12414), plot(t,umsd_train,t,umsd_test);
x0 = zeros(msd1.Nx,1);

[~, xTrain, yTrain, pTrain] = affineLpvSim(msd1_dt, msd1.eta_map, umsd_train, t, x0);
[~, xTest, yTest,pTest] = affineLpvSim(msd1_dt, msd1.eta_map, umsd_test, t, x0);

%% use less data cause you can
Ts = 20;
ts = 1:Ts:length(t);
N = length(ts);

%% Petrov POD - 3D POD - W = V
tol = 1e-6;

dataP = zeros(msd1.Nx, msd1.Np, N);
for i = 1:length(ts)
    dataP(:,:,i) = xTrain(:,ts(i))*pTrain(:,ts(i))';
end

tic
Thosvd = hosvd(tensor(dataP),tol);
CPUtime_HOSVD_POD = toc;
V1 = Thosvd{1};
W1 = V1;
Z1 = Thosvd{2};

K = min(modrank(dataP));
tic
[x,gains,iterations,Xa] = TSVDND_SW(dataP,K);
CPUtime_TSVD_POD = toc;
V2 = x{1};
W2 = x{1};
Z2 = x{2};
%%
EtsvdPOD = norm(dataP - Xa,'fro')
%% values from TMM
% V1= V2; W1 = W2; Z1 = Z2;
rx = 4;
rp = 3;
msd1_pod1 = msd1.PetrovGalerkinLPV(W1(:,1:rx),V1(:,1:rx),Z1(:,1:rp),1);
lpvss_msd1pod = msd1_pod1.tensSS2lpvss();

[~, ~, ypod_Train,ppod_Train] = affineLpvSim(lpvss_msd1pod, msd1_pod1.eta_map, umsd_train, t, V1(:,1:rx)'*x0);
[~, ~, ypod_Test,ppod_Test]   = affineLpvSim(lpvss_msd1pod, msd1_pod1.eta_map, umsd_test, t, V1(:,1:rx)'*x0);


e_pod_train = nrmse(yTrain,ypod_Train)
e_pod_test = nrmse(yTest,ypod_Test)
%% lower values
rx = 2;
rp = 1;
msd1_pod2 = msd1.PetrovGalerkinLPV(W1(:,1:rx),V1(:,1:rx),Z1(:,1:rp),1);
lpvss_msd1pod2 = msd1_pod2.tensSS2lpvss();

[~, ~, ypod_Train2,~] = affineLpvSim(lpvss_msd1pod2, msd1_pod2.eta_map, umsd_train, t, V1(:,1:rx)'*x0);
[~, ~, ypod_Test2,~]   = affineLpvSim(lpvss_msd1pod2, msd1_pod2.eta_map, umsd_test, t, V1(:,1:rx)'*x0);


e_pod_train2 = nrmse(yTrain,ypod_Train2)
e_pod_test2 = nrmse(yTest,ypod_Test2)
%% Petrov-Galerkin POD - 4D POD - W independent


xdot = xTrain(:,2:end) - xTrain(:,1:end-1);
dataPG = zeros(msd1.Nx, msd1.Nx, msd1.Np, length(ind));
for i = 1:length(ind)
    t1 = ttt(tensor(xTrain(:,ind(i))),tensor(xdot(:,ind(i))));
    t1 = squeeze(t1);
    t2 = ttt(t1,tensor(pTrain(:,ind(i))));
    t2 = squeeze(t2);
    dataPG(:,:,:,i) = t2;
end

tic
Thosvd = hosvd(tensor(dataPG),tol);
CPUtime_HOSVD_PODPG = toc;
V3 = Thosvd{1};
W3 = Thosvd{2};
Z3 = Thosvd{3};

K = min(modrank(dataPG));
tic
[x,gains,iterations,Xa] = TSVDND_SW(dataPG,K);
CPUtime_TSVD_PODPG = toc;
V4 = x{1};
W4 = x{2};
Z4 = x{3};

%%
ETSVD_PODPG = norm(dataPG-Xa,'fro')
%%
%% values from TMM
% V4= V3; W4 = W3; Z4 = Z3;
% V3 = V4; W3 = W4; Z3 = Z4;
rx = 4;
rp = 3;
msd1_pod3 = msd1.PetrovGalerkinLPV(W3(:,1:rx),V3(:,1:rx),Z3(:,1:rp),1);
lpvss_msd1pod3 = msd1_pod3.tensSS2lpvss();

[~, ~, ypod_Train3,~] = affineLpvSim(lpvss_msd1pod3, msd1_pod3.eta_map, umsd_train, t, V3(:,1:rx)'*x0);
[~, ~, ypod_Test3,~]   = affineLpvSim(lpvss_msd1pod3, msd1_pod3.eta_map, umsd_test, t, V3(:,1:rx)'*x0);


e_pod_train3 = nrmse(yTrain,ypod_Train3)
e_pod_test3 = nrmse(yTest,ypod_Test3)
%% lower values
rx = 2;
rp = 1;
msd1_pod4 = msd1.PetrovGalerkinLPV(W3(:,1:rx),V3(:,1:rx),Z3(:,1:rp),1);
lpvss_msd1pod4 = msd1_pod4.tensSS2lpvss();

[~, ~, ypod_Train4,~] = affineLpvSim(lpvss_msd1pod4, msd1_pod4.eta_map, umsd_train, t, V3(:,1:rx)'*x0);
[~, ~, ypod_Test4,~]   = affineLpvSim(lpvss_msd1pod4, msd1_pod4.eta_map, umsd_test, t, V3(:,1:rx)'*x0);


e_pod_train4 = nrmse(yTrain,ypod_Train4)
e_pod_test4 = nrmse(yTest,ypod_Test4)

%%

FigPODMSD1 = figure(141);
clrs = lines(5);
tiledlayout(2,1,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,yTrain,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,ypod_Train,'Color',clrs(2,:),'LineWidth',1.2);
plot(t,ypod_Train3,'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend("FOM",'Galerkin','Petrov-Galerkin','FontSize',10);
title(sprintf("ROM $r_x = %d, r_p = %d$",msd1_pod.Nx,msd1_pod.Np));

nexttile;
plot(t,yTest,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,ypod_Test,'Color',clrs(2,:),'LineWidth',1.2);
plot(t,ypod_Test3,'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend("FOM",'Galerkin','Petrov-Galerkin','FontSize',10);

%%
FigPODMSD11 = figure(142);
clrs = lines(5);
tiledlayout(2,1,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,yTrain,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,ypod_Train2,'Color',clrs(2,:),'LineWidth',1.2);
plot(t,ypod_Train4,'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend("FOM",'Galerkin','Petrov-Galerkin','FontSize',10);
title(sprintf("ROM $r_x = %d, r_p = %d$",msd1_pod4.Nx,msd1_pod4.Np),'FontSize',10);
nexttile;
plot(t,yTest,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,ypod_Test2,'Color',clrs(2,:),'LineWidth',1.2);
plot(t,ypod_Test4,'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend("FOM",'Galerkin','Petrov-Galerkin','FontSize',10);

%%
exportgraphics(FigPODMSD1,'FigPODMSD1.pdf')
exportgraphics(FigPODMSD11,'FigPODMSD11.pdf')



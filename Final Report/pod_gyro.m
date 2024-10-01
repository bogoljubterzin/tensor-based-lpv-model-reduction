gyro
x0Gyro
uGyro_train;
uGyro_test;

[~, xTrain, yTrain,pTrain] = affineLpvSim(gyro_dt, gyro.eta_map, uGyro_train, t, x0Gyro);
[~, xTest, yTest,pTest] = affineLpvSim(gyro_dt, gyro.eta_map, uGyro_test, t, x0Gyro);

%% (1) POD gyro;

Ts = 40
ind = 1:Ts:length(t);
data = zeros(gyro.Nx,gyro.Np,length(ind));
for i = 1:length(ind)
    data(:,:,i) = xTrain(:,ind(i))*pTrain(:,ind(i)).';
end
%%

tol = 1e-6;

tic
T = hosvd(tensor(data),tol);
tcomp = toc;

%%
rx = 4;
rp = 4;

V = T{1}(:,1:rx);
W = V;
Z = T{2}(:,1:rp);

gyro_pod = gyro.PetrovGalerkinLPV(W,V,Z,1);
lpvss_gyropod = gyro_pod.tensSS2lpvss();

[~, xpod_Train, ypod_Train,ppod_Train] = affineLpvSim(lpvss_gyropod, gyro_pod.eta_map, uGyro_train, t, V'*x0Gyro);
[~, xpod_Test, ypod_Test,ppod_Test] = affineLpvSim(lpvss_gyropod, gyro_pod.eta_map, uGyro_test, t, V'*x0Gyro);

e1 = nrmse(yTrain,ypod_Train)
e1test = nrmse(yTest,ypod_Test)
%%

rx1 = 5;
rp1 = 8;

V = T{1}(:,1:rx1);
W = V;
Z = T{2}(:,1:rp1);

gyro_pod2 = gyro.PetrovGalerkinLPV(W,V,Z,1);
lpvss_gyropod2 = gyro_pod2.tensSS2lpvss();

[~, ~, ypod_Train2,~] = affineLpvSim(lpvss_gyropod2, gyro_pod2.eta_map, uGyro_train, t, V'*x0Gyro);
[~, ~, ypod_Test2,~] = affineLpvSim(lpvss_gyropod2, gyro_pod2.eta_map, uGyro_test, t, V'*x0Gyro);
%%
e2 = nrmse(yTrain,ypod_Train2)
e2test = nrmse(yTest,ypod_Test2)
%%

%% (2) POD4D gyro;

ind = 1:20:length(t);

xdot = xTrain(:,2:end) - xTrain(:,1:end-1);
data = zeros(gyro.Nx, gyro.Nx, gyro.Np, length(ind));
for i = 1:length(ind)
    t1 = ttt(tensor(xTrain(:,ind(i))),tensor(xdot(:,ind(i))));
    t1 = squeeze(t1);
    t2 = ttt(t1,tensor(pTrain(:,ind(i))));
    t2 = squeeze(t2);
    data(:,:,:,i) = t2;
end
%%
tic
T = hosvd(tensor(data),tol);
tcomp = toc;

%%
rx = 4;
rp = 4;

V = T{1}(:,1:rx);
W = T{2}(:,1:rx);
Z = T{3}(:,1:rp);

gyro_pod3 = gyro.PetrovGalerkinLPV(W,V,Z,1);
lpvss_gyropod3 = gyro_pod3.tensSS2lpvss();

[~, ~, ypod_Train3,~] = affineLpvSim(lpvss_gyropod3, gyro_pod3.eta_map, uGyro_train, t, V'*x0Gyro);
[~, ~, ypod_Test3,~] = affineLpvSim(lpvss_gyropod3, gyro_pod3.eta_map, uGyro_test, t, V'*x0Gyro);

e3 = nrmse(yTrain,ypod_Train3)
e3test = nrmse(yTest,ypod_Test3)

% rx1, rp1;

V = T{1}(:,1:rx1);
W = T{2}(:,1:rx1);
Z = T{3}(:,1:rp1);

gyro_pod4 = gyro.PetrovGalerkinLPV(W,V,Z,1);
lpvss_gyropod4 = gyro_pod4.tensSS2lpvss();

[~, ~, ypod_Train4,~] = affineLpvSim(lpvss_gyropod4, gyro_pod4.eta_map, uGyro_train, t, V'*x0Gyro);
[~, ~, ypod_Test4,~] = affineLpvSim(lpvss_gyropod4, gyro_pod4.eta_map, uGyro_test, t, V'*x0Gyro);
e4 = nrmse(yTrain,ypod_Train4)
e4test = nrmse(yTest,ypod_Test4)

%%
% Galerkin:
% rx = 4 rp = 4
ypod_Train;
ypod_Test;
% rx1 = 5; rp1 = 8;
ypod_Train2;
ypod_Test2;

% Petrov-Galerkin:
% rx = 4 rp = 4
ypod_Train3;
ypod_Test3;
% rx1 = 5; rp1 = 8;
ypod_Train4;
ypod_Test4;

e1 = nrmse(yTrain,ypod_Train)
e14d = nrmse(yTrain,ypod_Train3)

e2 = nrmse(yTest,ypod_Test)
e24d = nrmse(yTest,ypod_Test3)

nrmse(yTrain,ypod_Train2)
nrmse(yTrain,ypod_Train4)
nrmse(yTest,ypod_Test2)
nrmse(yTest,ypod_Test4)
%%
gyropod_train = figure(3),
clrs = lines(3);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t,yTrain(1,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train(1,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train3(1,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_1(t)$');


nexttile;
plot(t,yTrain(2,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train(2,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train3(2,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_2(t)$');

nexttile;
plot(t,yTrain(3,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train(3,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train3(3,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_3(t)$');
legend(sprintf("FOM"),sprintf("Galerkin, $r_x = %d, r_p = %d$", gyro_pod.Nx,gyro_pod.Np), ...
       sprintf("Petrov-Galerkin, $r_x = %d, r_p = %d$", gyro_pod.Nx,gyro_pod.Np),'FontSize',11)
sgtitle("Galerkin vs Petrov-Galerkin on CMG ($r_x = 4, r_p = 4)$ - Training data",'Interpreter','Latex');


%%
gyropod_Test = figure(4),
clrs = lines(3);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t,yTest(1,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test(1,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test3(1,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_1(t)$');

nexttile;
plot(t,yTest(2,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test(2,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test3(2,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_2(t)$');

nexttile;
plot(t,yTest(3,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test(3,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test3(3,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_3(t)$');
legend(sprintf("FOM"),sprintf("Galerkin, $r_x = %d, r_p = %d$", gyro_pod.Nx,gyro_pod.Np), ...
       sprintf("Petrov-Galerkin, $r_x = %d, r_p = %d$", gyro_pod.Nx,gyro_pod.Np),'FontSize',11)
sgtitle("Galerkin vs Petrov-Galerkin on CMG ($r_x = 4, r_p = 4)$ - Testing data",'Interpreter','Latex');


%%

gyropod_train2 = figure(5),
clrs = lines(3);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t,yTrain(1,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train2(1,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train4(1,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_1(t)$');



nexttile;
plot(t,yTrain(2,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train2(2,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train4(2,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_2(t)$');

nexttile;
plot(t,yTrain(3,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Train2(3,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Train4(3,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_3(t)$');
legend(sprintf("FOM"),sprintf("Galerkin, $r_x = %d, r_p = %d$", gyro_pod4.Nx,gyro_pod4.Np), ...
       sprintf("Petrov-Galerkin, $r_x = %d, r_p = %d$", gyro_pod4.Nx,gyro_pod4.Np),'FontSize',11)
sgtitle("Galerkin vs Petrov-Galerkin on CMG ($r_x = 5, r_p = 8)$ - Training data",'Interpreter','Latex');


gyropod_test2 = figure(7),
clrs = lines(3);
tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t,yTest(1,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test2(1,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test4(1,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_1(t)$');



nexttile;
plot(t,yTest(2,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test2(2,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test4(2,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_2(t)$');

nexttile;
plot(t,yTest(3,:),'LineWidth',3,'LineStyle','--'); hold on; grid on;
plot(t,ypod_Test2(3,:),'LineWidth',1.2,'Color',clrs(2,:))
plot(t,ypod_Test4(3,:),'LineWidth',1.2,'Color',clrs(3,:))
xlabel('Time (s)'); ylabel('$y_3(t)$');
legend(sprintf("FOM"),sprintf("Galerkin, $r_x = %d, r_p = %d$", gyro_pod4.Nx,gyro_pod4.Np), ...
       sprintf("Petrov-Galerkin, $r_x = %d, r_p = %d$", gyro_pod4.Nx,gyro_pod4.Np),'FontSize',11)
sgtitle("Galerkin vs Petrov-Galerkin on CMG ($r_x = 5, r_p = 8)$ - Testing data",'Interpreter','Latex');

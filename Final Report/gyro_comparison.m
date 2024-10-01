set(groot, 'DefaultTextInterpreter', 'latex');        % For text objects
set(groot, 'DefaultLegendInterpreter', 'latex');      % For legends
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % For axis tick labels
set(groot, 'DefaultColorbarTickLabelInterpreter', 'latex'); % For colorbar tick labels


%% Gyroscope
% comparing TMM vs LPVMM + PCA
x0Gyro = rand(gyro.Nx,1);
uGyro_train = [sin(t); cos(t); 0.5*sin(t)-0.1];
uGyro_test = [3*cos(t); 2*sin(t); sin(t)+0.4];
[~,~,ygyro_train,pgyro_train] = gyro.simulateSS(uGyro_train,t,x0Gyro,1); 
[~,~,ygyro_test,~] = gyro.simulateSS(uGyro_test,t,x0Gyro,1); 


%%
tensSS = gyro;

mode =  'R';
decomp = 'hosvd';
N = 0;
tic
[tensSSr,Tx,Tp,info] = tensSS.lpvTensMM(N,decomp,mode);
t1 = toc;
%%
[~,~,yrtmm_train,~] = tensSSr.simulateSS(uGyro_train,t,Tx'*x0Gyro,0); 
[~,~,yrtmm_test,~] = tensSSr.simulateSS(uGyro_test,t,Tx'*x0Gyro,0); 

%%
% [~,~,yrtmm_train0] = affineLpvSim(lpv1, tensSSr.eta_map, uGyro_train, t, Tx'*x0Gyro);
%%
tic
[lpvssR,Tx1] = lpvmmred(gyro_dt,0,'R')
eta1_mm = @(x,u) eta_gyro.map(Tx1*x,u);

pTraj = pgyro_train'; Rp = tensSSr.Np;
[lpvssMMPCA, mapPCA1] = lpvpcared(lpvssR, Rp, pTraj,'trajectory');
etar1 = @(x,u) mapPCA1(eta1_mm(x,u).');
t2 = toc;
%%
[~, ~, yrmmpca_train] = affineLpvSim(lpvssMMPCA, etar1, uGyro_train, t, Tx1'*x0Gyro);
[~, ~, yrmmpca_test] = affineLpvSim(lpvssMMPCA, etar1, uGyro_test, t, Tx1'*x0Gyro);

%% Again for N=1
N = 1;
tic
[tensSSr1,Tx2,Tp2,info2] = tensSS.lpvTensMM(N,decomp,mode);
t3 = toc;

%
[~,~,yrtmm_train1,~] = tensSSr1.simulateSS(uGyro_train,t,Tx2'*x0Gyro,0); 
[~,~,yrtmm_test1,~] = tensSSr1.simulateSS(uGyro_test,t,Tx2'*x0Gyro,0); 

%
tic
[lpvssR1,Tx3] = lpvmmred(gyro_dt,N,'R')
eta3_mm = @(x,u) eta_gyro.map(Tx3*x,u);

pTraj = pgyro_train'; Rp = tensSSr1.Np;
[lpvssMMPCA1, mapPCA11] = lpvpcared(lpvssR1, Rp, pTraj,'trajectory');
etar3 = @(x,u) mapPCA11(eta3_mm(x,u).');
t4 = toc;
%
[~, ~, yrmmpca_train1] = affineLpvSim(lpvssMMPCA1, etar3, uGyro_train, t, Tx3'*x0Gyro);
[~, ~, yrmmpca_test1] = affineLpvSim(lpvssMMPCA1, etar3, uGyro_test, t, Tx3'*x0Gyro);

%%

nrmse(ygyro_train,yrtmm_train)
nrmse(ygyro_test,yrtmm_test)

nrmse(ygyro_train,yrtmm_train1)
nrmse(ygyro_test,yrtmm_test1)


nrmse(ygyro_train,yrmmpca_train)
nrmse(ygyro_test,yrmmpca_test)

nrmse(ygyro_train,yrmmpca_train1)
nrmse(ygyro_test,yrmmpca_test1)
%%

FigGyro = figure(1214),
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')
nexttile;
plot(t,ygyro_train(1,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train(1,:))
plot(t,yrmmpca_train(1,:))
legend("FOM","TMM","MM+PCA",'FontSize',10);
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',9);
xlabel("Time (s)"); ylabel('$y_1(t)$')
nexttile;
plot(t,ygyro_train(2,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train(2,:))
plot(t,yrmmpca_train(2,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_2(t)$')
nexttile;
plot(t,ygyro_train(3,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train(3,:))
plot(t,yrmmpca_train(3,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_3(t)$')

%%

FigGyroTest = figure(1215),
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')
nexttile;
plot(t,ygyro_test(1,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_test(1,:))
plot(t,yrmmpca_test(1,:))
legend("FOM","TMM","MM+PCA",'FontSize',10);
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',9);
xlabel("Time (s)"); ylabel('$y_1(t)$')
nexttile;
plot(t,ygyro_test(2,:),'LineStyle','--','LineWidth',2); hold on;grid on;
plot(t,yrtmm_test(2,:))
plot(t,yrmmpca_test(2,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_2(t)$')
nexttile;
plot(t,ygyro_test(3,:),'LineStyle','--','LineWidth',2); hold on;grid on;
plot(t,yrtmm_test(3,:))
plot(t,yrmmpca_test(3,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_3(t)$')
%%

[lpvssR_gyro,Tx1] = lpvmmred(gyro_dt, 1, 'R');
eta1_mm = @(x,u) eta_gyro.map(Tx1*x,u);
%%
pTraj = p1'; Rp = tensSSr.Np;
[lpvssPCA1, mapPCA1] = lpvpcared(lpvssR_gyro, Rp, pTraj,'trajectory');
etar1 = @(x,u) mapPCA1(eta1_mm(x,u).');

[~, ~, yrmmpca_train] = affineLpvSim(lpvssPCA1, etar1, umsd_train, t, Tx1'*x0);
[~, ~, yrmmpca_test] = affineLpvSim(lpvssPCA1, etar1, umsd_test, t, Tx1'*x0);
%% %%%%%%%%%%%%%%%%%%%%%%55

FigGyro1 = figure(1216),
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')
nexttile;
plot(t,ygyro_train(1,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train1(1,:))
plot(t,yrmmpca_train1(1,:))
legend("FOM","TMM","MM+PCA",'FontSize',10);
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',9);
xlabel("Time (s)"); ylabel('$y_1(t)$')
nexttile;
plot(t,ygyro_train(2,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train1(2,:))
plot(t,yrmmpca_train1(2,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_2(t)$')
nexttile;
plot(t,ygyro_train(3,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_train1(3,:))
plot(t,yrmmpca_train1(3,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_3(t)$')

%%

FigGyroTest1 = figure(1217),
tiledlayout(3,1,'TileSpacing','tight','Padding','compact')
nexttile;
plot(t,ygyro_test(1,:),'LineStyle','--','LineWidth',2); hold on; grid on;
plot(t,yrtmm_test1(1,:))
plot(t,yrmmpca_test1(1,:))
legend("FOM","TMM","MM+PCA",'FontSize',10);
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',9);
xlabel("Time (s)"); ylabel('$y_1(t)$')
nexttile;
plot(t,ygyro_test(2,:),'LineStyle','--','LineWidth',2); hold on;grid on;
plot(t,yrtmm_test1(2,:))
plot(t,yrmmpca_test1(2,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_2(t)$')
nexttile;
plot(t,ygyro_test(3,:),'LineStyle','--','LineWidth',2); hold on;grid on;
plot(t,yrtmm_test1(3,:))
plot(t,yrmmpca_test1(3,:))
% legend(sprintf("$FOM, n_x = %d, n_p = %d$", gyro.Np,gyro.Np), ...
%     sprintf("$TMM, r_x = %d, r_p = %d$", tensSSr.Np,tensSSr.Np), ...
%     sprintf("$MM+PCA, r_x = %d, r_p = %d$", lpvssMMPCA.Np,lpvssMMPCA.Np),'FontSize',11);
xlabel("Time (s)"); ylabel('$y_3(t)$')
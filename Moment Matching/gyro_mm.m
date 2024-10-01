[nlss_gyro, lpvsys_gyro, eta_gyro] = gen_gyro
%%
uGyro = [sin(2*t); cos(2*t); zeros(1,length(t))];
Ts = 0.005;
lpvsysdt = c2d(lpvsys_gyro,Ts);
tensGyro = tensorSS(lpvsysdt,eta_gyro);

x0 = rand(tensGyro.Nx,1);
[tout,x,y] = tensGyro.simulateSS(uGyro,t,x0,1)

%%
rGyro = cell(3,1);
%%
decomp = 'hosvd';
mode = 'R';
info_decomposition = cell(3,1);
for N = 1:3
    [redGyro,Vr,Vp,info] = tensGyro.lpvTensMM(N, decomp, mode)
    rGyro{N} = redGyro;
    info_decomposition{N} = info{1};
end

%%
% save('rGyro_TSVD.mat','rGyro_TSVD','infoTSVD')
% save('rGyro_HOSVD.mat','rGyro_HOSVD','infoHOSVD')

%%
[tout,xr2,yr] = redGyro.simulateSS(uGyro,t,Vr'*x0,0)

%%
y1 = y(1,:);
y2 = y(2,:);
y3 = y(3,:);

yr1 = yr(1,:);
yr2 = yr(2,:);
yr3 = yr(3,:);

figure(10),
plot(t,y1,t,yr1);
figure(11),
plot(t,y2,t,yr2);
figure(12),
plot(t,y3,t,yr3);


%%
yhosvd = cell(3,1);
yhosvd{2} = yr;

%%

GyroOutput1 = figure(1);
clf(GyroOutput1)
plot(t, y1,'LineWidth',2,'LineStyle','--','Color','k');
hold on; grid on;
plot(t,yhosvd{1}(1,:),'LineWidth',1);
plot(t,yhosvd{2}(1,:),'LineWidth',1);
plot(t,yhosvd{3}(1,:),'LineWidth',1);
legend(sprintf("FOM, N_x = %d, N_p = %d", tensGyro.Nx, tensGyro.Np),sprintf("N=1, R_x = %d, R_p = %d", rGyro{1}.Nx,rGyro{1}.Np), ...
    sprintf("N=2, R_x = %d, R_p = %d", rGyro{2}.Nx,rGyro{2}.Np), ...
    sprintf("N=3, R_x = %d, R_p = %d", rGyro{3}.Nx,rGyro{3}.Np));
xlabel("Time (s)"); ylabel("y1")

GyroOutput2 = figure(2);
clf(GyroOutput2)
plot(t, y2,'LineWidth',2,'LineStyle','--','Color','k');
hold on; grid on;
plot(t,yhosvd{1}(2,:),'LineWidth',1.2);
plot(t,yhosvd{2}(2,:),'LineWidth',1.2);
plot(t,yhosvd{3}(2,:),'LineWidth',1.2);
legend(sprintf("FOM, N_x = %d, N_p = %d", tensGyro.Nx, tensGyro.Np),sprintf("N=1, R_x = %d, R_p = %d", rGyro{1}.Nx,rGyro{1}.Np), ...
    sprintf("N=2, R_x = %d, R_p = %d", rGyro{2}.Nx,rGyro{2}.Np), ...
    sprintf("N=3, R_x = %d, R_p = %d", rGyro{3}.Nx,rGyro{3}.Np));
xlabel("Time (s)"); ylabel("y2")

GyroOutput3 = figure(3);
clf(GyroOutput3)
plot(t, y3,'LineWidth',2,'LineStyle','--','Color','k');
hold on; grid on;
plot(t,yhosvd{1}(3,:),'LineWidth',1.2);
plot(t,yhosvd{2}(3,:),'LineWidth',1.2);
plot(t,yhosvd{3}(3,:),'LineWidth',1.2);
legend(sprintf("FOM, N_x = %d, N_p = %d", tensGyro.Nx, tensGyro.Np),sprintf("N=1, R_x = %d, R_p = %d", rGyro{1}.Nx,rGyro{1}.Np), ...
    sprintf("N=2, R_x = %d, R_p = %d", rGyro{2}.Nx,rGyro{2}.Np), ...
    sprintf("N=3, R_x = %d, R_p = %d", rGyro{3}.Nx,rGyro{3}.Np));
xlabel("Time (s)"); ylabel("y3")
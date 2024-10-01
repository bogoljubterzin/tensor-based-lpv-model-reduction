[nlss_gyro, lpvss_gyro, eta_gyro] = gen_gyro;
[nlss_robot, lpvss_robot, eta_robot] = gen_robot;
[nlss_msd1, lpvss_msd1, eta_msd1] = gen_MSD(5,9);
[nlss_msd3, lpvss_msd3, eta_msd3] = gen_MSD(50,99);
%% discretization

Td = 1e-3;
msd1_dt = c2d(lpvss_msd1,Td);
msd3_dt = c2d(lpvss_msd3,Td);
gyro_dt = c2d(lpvss_gyro,Td);
robot_dt = c2d(lpvss_robot,Td);

%% Convert to tensorSS models

msd1 = tensorSS(msd1_dt,eta_msd1);
msd3 = tensorSS(msd3_dt,eta_msd3);
gyro = tensorSS(gyro_dt,eta_gyro);
robot = tensorSS(robot_dt,eta_robot);

%% Gather previous input signals
load trainref_100Hz_25s.mat
t = 0:0.01:20;
ind = 1:length(t);
umsd_train = ref_100Hz_25s{1,1,1}(ind).';
umsd_test = 0.2*ref_100Hz_25s{3,1,1}(ind).';

uGyro_train = ref_100Hz_25s{1,3,1}(ind).';
uGyro_test = ref_100Hz_25s{3,3,1}(ind).';

uRobot_train = ref_100Hz_25s{1,2,1}(ind).';
uRobot_test = ref_100Hz_25s{3,2,1}(ind).';

constantTerm = 1;
[~,xmsd1_train,ymsd1_train,~] = msd1.simulateSS(umsd_train,t,zeros(msd1.Nx,1),constantTerm);
[~,xmsd1_test,ymsd1_test,~] = msd1.simulateSS(umsd_test,t,zeros(msd1.Nx,1),constantTerm); 
[~,~,ymsd3_train,~] = msd3.simulateSS(umsd_train,t,zeros(msd3.Nx,1),constantTerm); 
[~,~,ymsd3_test,~] = msd3.simulateSS(umsd_test,t,zeros(msd3.Nx,1),constantTerm); 
%%
x0Gyro = rand(gyro.Nx,1);
uGyro = [sin(t); cos(t); 0.5*sin(t)-0.1];
uRobot = [sin(t); cos(t)];
[~,~,ygyro,~] = gyro.simulateSS(uGyro,t,x0Gyro,constantTerm); 
[~,~,yrobot,~] = robot.simulateSS(uRobot,t,zeros(robot.Nx,1),constantTerm); 

%%
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.2);
%%
Fig1 = figure(1);
tiledlayout(3,1,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t, umsd_train); hold on; grid on;
plot(t, umsd_test);
xlabel("Time (s)"); ylabel("Amplitude");
title('MSD');
legend("Training data","Testing data");

nexttile;
plot(t,uGyro); grid on;
xlabel("Time (s)"); ylabel("Amplitude");
title('Gyroscope');
legend("$u_1(t)$","$u_2(t)$","$u_3(t)$");

nexttile;
plot(t,uRobot); grid on;
xlabel("Time (s)"); ylabel("Amplitude");
title('Robotic Arm');
legend("$u_1(t)$","$u_2(t)$");

%%

Fig2 = figure(2);
tiledlayout(2,1,'TileSpacing','tight','Padding','compact');

nexttile;
plot(t,ymsd1_train); hold on; grid on;
title('Training data (MSD models)');
xlabel('Time (s)'); ylabel('Position (m)');
nexttile;
plot(t, ymsd1_test);grid on; hold on
title('Testing data (MSD models)');
xlabel('Time (s)'); ylabel('Position (m)');

% sgtitle('MSD');

Fig3 = figure(3);
tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
nexttile;
plot(t, ygyro); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Gyroscope'); legend('$y_1(t)$','$y_2(t)$','$y_3(t)$');
nexttile;
plot(t, yrobot);grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Robotic arm'); legend('$y_1(t)$','$y_2(t)$');


%% Space exploration plot

[U1,S1,V1] = svd(xmsd1_train);
x1 = U1(:,1:3)'*xmsd1_train;

[U2,S2,V2] = svd(xmsd1_test);
x2 = U2(:,1:3)'*xmsd1_test;
%%
x1 = xmsd1_train(4,:);
x2 = xmsd1_train(5,:);
x3 = xmsd1_train(10,:);


x4 = xmsd1_test(4,:);
x5 = xmsd1_test(5,:);
x6 = xmsd1_test(10,:);
%%
state_space_exploration = figure(4); clf(state_space_exploration);
ind1 = 1:4:length(t);
tiledlayout(1,1,"TileSpacing","tight","Padding","tight")
cmp = parula(5);
ax(1) = nexttile; 
pval(2)=plot3(x1(ind1),x2(ind1),x3(ind1)); hold on;
pval(3)=plot3(x4(ind1),x5(ind1),x6(ind1)); hold on;
title("$\mathrm{MSD}_1$");

grid(ax(1),'on');
set(pval(2),'Marker','o'); set(pval(2),'LineStyle','none'); set(pval(2),'Color',cmp(3,:)) 
set(pval(3),'Marker','o'); set(pval(3),'LineStyle','none'); set(pval(3),'Color',cmp(4,:)) 

set([pval(2)],'MarkerEdgeColor',cmp(3,:));
set([pval(3)],'MarkerEdgeColor',cmp(4,:));
legend("Training data",'Testing data')
xlabel('$x_4$'); ylabel('$x_5$'); zlabel('$x_{10}$');
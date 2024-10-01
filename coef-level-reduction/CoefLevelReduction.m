%% Coefficient Level Reduction - Tensor-based LPV model reduction method
% Testing it on benchmarks: MSD1, Gyroscope, Robotic Arm
% Used decompositions: TSVD, succR1, HOSVD, (CP)
% Type 1: tensor A; Type 2: tensors A,B,C,D

addpath(genpath("..\tools\benchmarks-inputs"));
addpath(genpath("..\tensorSS"))
addpath(genpath("..\tools\tensor_toolbox-master"));
addpath(genpath("..\siep-functions"));
set(groot, 'defaultTextInterpreter', 'latex');
% bm_systems = getBenchmarkSystems

load trainref_10Hz_25s.mat
data = ref_10Hz_25s;

msd_train = data{1,1,1}.';
robot_train = data{1,2,1}.';
% Gyroscope is misssing one input signal! I added manually u3 = u1
gyro_train = data{1,3,1}.'; gyro_train(3,:) = gyro_train(1,:);

% Another input signals for simulation
% msd_train = data{2,1,1}.';
% robot_train = data{2,2,1}.';
% % Gyroscope is misssing one input signal! I added manually u3 = u1
% gyro_train = data{2,3,1}.'; gyro_train(3,:) = gyro_train(1,:);

Ts = 1/10;
tin = 0:Ts:(length(msd_train)-1)*Ts;
[~,~,y1,~] = bm_systems{1}.simulateSS(msd_train,tin,zeros(bm_systems{1}.Nx,1),'FOM');
[~,~,y2,~] = bm_systems{2}.simulateSS(gyro_train,tin,zeros(bm_systems{2}.Nx,1),'FOM');
[~,~,y3,~] = bm_systems{3}.simulateSS(msd_train,tin,zeros(bm_systems{3}.Nx,1),'FOM');
[~,~,y4,~] = bm_systems{4}.simulateSS(robot_train,tin,zeros(bm_systems{4}.Nx,1),'FOM');

%% BM1 - msd1 - applying methods
rom = 'ROM';
bm = 1;
tensSS = bm_systems{bm};

rx = 8;
rp = 4;

% TSVD
type = 1;
decomp = "TSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t11 = toc
[~,~,yr1_tsvd_type1] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);


type = 2;
decomp = "TSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t12 = toc
[~,~,yr1_tsvd_type2] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);

% succR1
type = 1;
decomp = "succR1";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t21 = toc
[~,~,yr1_succR1_type1] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);


type = 2;
decomp = "succR1";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t22 = toc
[~,~,yr1_succR1_type2] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);


% HOSVD
type = 1;
decomp = "HOSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t31 = toc
[~,~,yr1_hosvd_type1] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);


type = 2;
decomp = "HOSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp)
t32 = toc
[~,~,yr1_hosvd_type2] = sysr.simulateSS(msd_train,tin,zeros(sysr.Nx,1),rom);


time_msd1 = [t11,t12;t21 t22;t31,t32];
 
%% msd1 - plotting
set(groot, 'defaultTextInterpreter', 'latex');

FigMSD1CoefRed_Type1 = figure(111);
tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y1,'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr1_tsvd_type1,'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr1_succR1_type1,'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr1_hosvd_type1,'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output y");

t(2) = nexttile;
semilogy(tin, abs(y1 - yr1_tsvd_type1),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y1 - yr1_succR1_type1),'Color',cmp(130,:));
semilogy(tin, abs(y1 - yr1_hosvd_type1),'Color',cmp(200,:));
xlabel("Time (s)"); ylabel("Absolute error e");

legend("TSVD","SUCC R1","HOSVD")
sgtitle(sprintf("CoefRed 1 MSD1 (States %d to %d) (Sched. %d to %d)",bm_systems{1}.Nx,rx,bm_systems{1}.Np,rp));


FigMSD1CoefRed_Type2 = figure(112);
tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y1,'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr1_tsvd_type2,'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr1_succR1_type2,'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr1_hosvd_type2,'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output y");

t(2) = nexttile;
semilogy(tin, abs(y1 - yr1_tsvd_type2),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y1 - yr1_succR1_type2),'Color',cmp(130,:));
semilogy(tin, abs(y1 - yr1_hosvd_type2),'Color',cmp(200,:));
xlabel("Time (s)"); ylabel("Absolute error e");

legend("TSVD","SUCC R1","HOSVD")
sgtitle(sprintf("CoefRed 2 MSD1 (States %d to %d) (Sched. %d to %d)",bm_systems{1}.Nx,rx,bm_systems{1}.Np,rp));


exportgraphics(FigMSD1CoefRed_Type1,"CoefRed1MSD1rx8rp4.pdf");
exportgraphics(FigMSD1CoefRed_Type2,"CoefRed2MSD1rx8rp4.pdf");


%% %% BM2 - Gyro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bm_systems = getBenchmarkSystems
rom = 'ROM';
bm = 2;
tensSS = bm_systems{bm};
rx = 3;
rp = 4;
% TSVD
type = 1;
decomp = "TSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t11 = toc

[~,~,yr2_tsvd_type1] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "TSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t12 = toc

[~,~,yr2_tsvd_type2] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);

% succR1
type = 1;
decomp = "succR1";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t21 = toc
[~,~,yr2_succR1_type1] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "succR1";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t22 = toc
[~,~,yr2_succR1_type2] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);

% HOSVD
type = 1;
decomp = "HOSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t31 = toc

[~,~,yr2_hosvd_type1] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "HOSVD";
tic
sysr = tensSS.coef_level_reduction(decomp, type, rx, rp);
t32 = toc

[~,~,yr2_hosvd_type2] = sysr.simulateSS(gyro_train,tin,zeros(sysr.Nx,1),rom);


time_gyro = [t11,t12;t21,t22;t31,t32]

%% gyroscope - plotting
set(groot, 'defaultTextInterpreter', 'latex');

FigGyroCoefRed_Type1 = figure(121);
tiledlayout(3,2);
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y2_train(1,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type1(1,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type1(1,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type1(1,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_1$");


t(2) = nexttile;
semilogy(tin, abs(y2_train(1,:) - yr2_tsvd_type1(1,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(1,:) - yr2_succR1_type1(1,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(1,:) - yr2_hosvd_type1(1,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_1$");


t(3) = nexttile;
plot(tin, y2_train(2,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type1(2,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type1(2,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type1(2,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_2$");

t(4) = nexttile;
semilogy(tin, abs(y2_train(2,:) - yr2_tsvd_type1(2,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(2,:) - yr2_succR1_type1(2,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(2,:) - yr2_hosvd_type1(2,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_2$");

t(5) = nexttile;
plot(tin, y2_train(3,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type1(3,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type1(3,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type1(3,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_3$");

t(6) = nexttile;
semilogy(tin, abs(y2_train(3,:) - yr2_tsvd_type1(3,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(3,:) - yr2_succR1_type1(3,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(3,:) - yr2_hosvd_type1(3,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_3$");

sgtitle(sprintf("CoefRed 1 Gyro (States %d to %d) (Sched. %d to %d)",tensSS.Nx,rx,tensSS.Np,rp));

%%
FigGyroCoefRed_Type2 = figure(122);
tiledlayout(3,2);
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y2_train(1,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type2(1,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type2(1,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type2(1,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_1$");


t(2) = nexttile;
semilogy(tin, abs(y2_train(1,:) - yr2_tsvd_type2(1,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(1,:) - yr2_succR1_type2(1,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(1,:) - yr2_hosvd_type2(1,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_1$");


t(3) = nexttile;
plot(tin, y2_train(2,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type2(2,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type2(2,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type2(2,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_2$");

t(4) = nexttile;
semilogy(tin, abs(y2_train(2,:) - yr2_tsvd_type2(2,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(2,:) - yr2_succR1_type2(2,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(2,:) - yr2_hosvd_type2(2,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_2$");

t(5) = nexttile;
plot(tin, y2_train(3,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr2_tsvd_type2(3,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr2_succR1_type2(3,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr2_hosvd_type2(3,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_3$");

t(6) = nexttile;
semilogy(tin, abs(y2_train(3,:) - yr2_tsvd_type2(3,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y2_train(3,:) - yr2_succR1_type2(3,:)),'Color',cmp(130,:));
semilogy(tin, abs(y2_train(3,:) - yr2_hosvd_type2(3,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_3$");

sgtitle(sprintf("CoefRed 2 Gyro (States %d to %d) (Sched. %d to %d)",tensSS.Nx,rx,tensSS.Np,rp));



%%

exportgraphics(FigGyroCoefRed_Type1,"CoefRed1Gyrorx3rp4.pdf");
exportgraphics(FigGyroCoefRed_Type2,"CoefRed2Gyrorx3rp4.pdf");  

%% Robotic ARm - applying method

rom = 'ROM';
bm = 4;

rx = 4; rp = 5;

% TSVD
type = 1;
decomp = "TSVD";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t11 = toc

[~,~,yr4_tsvd_type1] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "TSVD";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t12 = toc

[~,~,yr4_tsvd_type2] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);

% succR1
type = 1;
decomp = "succR1";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t21 = toc
[~,~,yr4_succR1_type1] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "succR1";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t22 = toc
[~,~,yr4_succR1_type2] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);

% HOSVD
type = 1;
decomp = "HOSVD";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t31 = toc

[~,~,yr4_hosvd_type1] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);

type = 2;
decomp = "HOSVD";
tic
sysr = bm_systems{bm}.coef_level_reduction(decomp, type, rx, rp)
t32 = toc

[~,~,yr4_hosvd_type2] = sysr.simulateSS(robot_train,tin,zeros(sysr.Nx,1),rom);
time_robot = [t11,t12;t21,t22;t31,t32]


%% robotic arm - plotting results
set(groot, 'defaultTextInterpreter', 'latex');


FigRobotCoefRed_Type1 = figure(141);
tiledlayout(2,2);
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y4(1,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr4_tsvd_type1(1,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr4_succR1_type1(1,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr4_hosvd_type1(1,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_1$");


t(2) = nexttile;
semilogy(tin, abs(y4(1,:) - yr4_tsvd_type1(1,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y4(1,:) - yr4_succR1_type1(1,:)),'Color',cmp(130,:));
semilogy(tin, abs(y4(1,:) - yr4_hosvd_type1(1,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_1$");


t(3) = nexttile;
plot(tin, y4(2,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr4_tsvd_type1(2,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr4_succR1_type1(2,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr4_hosvd_type1(2,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_2$");

t(4) = nexttile;
semilogy(tin, abs(y4(2,:) - yr4_tsvd_type1(2,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y4(2,:) - yr4_succR1_type1(2,:)),'Color',cmp(130,:));
semilogy(tin, abs(y4(2,:) - yr4_hosvd_type1(2,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_2$");

sgtitle(sprintf("CoefRed 1 Robot Arm (States %d to %d) (Sched. %d to %d)",bm_systems{4}.Nx,rx,bm_systems{4}.Np,rp));

%%
FigRobotCoefRed_Type2 = figure(142);
tiledlayout(2,2);
cmp = colormap('parula');

t(1) = nexttile;
plot(tin, y4(1,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr4_tsvd_type2(1,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr4_succR1_type2(1,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr4_hosvd_type2(1,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_1$");

t(2) = nexttile;
semilogy(tin, abs(y4(1,:) - yr4_tsvd_type2(1,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y4(1,:) - yr4_succR1_type2(1,:)),'Color',cmp(130,:));
semilogy(tin, abs(y4(1,:) - yr4_hosvd_type2(1,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_1$");

t(3) = nexttile;
plot(tin, y4(2,:),'LineWidth',2,'Color',cmp(1,:)); grid on; hold on;
plot(tin, yr4_tsvd_type2(2,:),'LineWidth',1,'Color',cmp(80,:));
plot(tin, yr4_succR1_type2(2,:),'LineWidth',1,'Color',cmp(130,:));
plot(tin, yr4_hosvd_type2(2,:),'LineWidth',1,'Color',cmp(200,:));
legend("Original","TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Output $y_2$");

t(4) = nexttile;
semilogy(tin, abs(y4(2,:) - yr4_tsvd_type2(2,:)),'Color',cmp(80,:)); hold on; grid on;
semilogy(tin, abs(y4(2,:) - yr4_succR1_type2(2,:)),'Color',cmp(130,:));
semilogy(tin, abs(y4(2,:) - yr4_hosvd_type2(2,:)),'Color',cmp(200,:));
legend("TSVD","SUCC R1","HOSVD")
xlabel("Time (s)"); ylabel("Error $e_2$");

sgtitle(sprintf("CoefRed 2 Robot Arm (States %d to %d) (Sched. %d to %d)",bm_systems{4}.Nx,rx,bm_systems{4}.Np,rp));




%%
exportgraphics(FigRobotCoefRed_Type1,"CoefRed1Robotrx4rp9.pdf");
exportgraphics(FigRobotCoefRed_Type2,"CoefRed2Robotrx4rp9.pdf");    


% Testing benchmark generator of NL MSD interconnected system
clear; close all; clc;
addpath(genpath("../../Tools/")) %Custom branch
%% Test 1: custom MSD benchmark
% Input parameters:
M = 5; % Number of mass-blockz
%Structure containing the indexes where we want nonlin stuff
springs = [1:1:(M-1)];
springs_wall = [1:1:M];
dampers = [1:1:(M-1)];
dampers_wall = [1:1:M];
nl_index.springs = springs;
nl_index.springs_wall = springs_wall;
nl_index.dampers = dampers;
nl_index.dampers_wall = dampers_wall;

%Define nonlinearities with the nl function
k = @(x) x^2; % Spring connecting mass blocks
% k = @(x) 0.5*x; % Spring connecting mass blocks
d = @(x,xdot)  xdot^3; % Damper connecting mass blocks

k_wall = @(x) x^2; % Spring connecting each mass with the wall
% k_wall = @(x) 0.5*x; % Spring connecting each mass with the wall
d_wall = @(x,xdot) xdot^3; % Damper connecting each mass with the wall

param = [1 0.5 1 0.5 1];

%% Generate custom MSD NL equations
tic
[fMSD,state,input] = generalized_msd_Javi(M,nl_index,k,d,k_wall,d_wall);
toc

% Nonlinear state-space model
f = @(x,u) fMSD(x,u,param);
h = @(x,u) x(M,:);
nx = length(state);
nu = length(input);
ny = 1;

sys = LPVcore.nlss(f,h,nx,nu,ny,0,true);
xx = sym('x',[nx,1]); uu = sym('u',[nu,1]);
f(xx,uu)
% SIM our custom MSD
tin = [0 1000];
x0 = zeros(nx,1);
u = @(tin) 0.5*cos(0.01*tin)+0.5;
[y,tout,x] = sim(sys,u,tin',x0,'ode45');
figure(1)
plot(tout,x(:,3:4))
%% Test lpv conversion
[lpvsys, eta] = LPVcore.nlss2lpvss(sys, 'analytical', 'element');
lpvsys.SchedulingTimeMap.RateBound(:) = {[-1, 1]};
lpvsys.SchedulingTimeMap.Range(:) = {[-1, 1]};
% gam = lpvnorm(lpvsys,'h2');
%% Test conversion to lpvgridss
cellArray = cell(nx + nu,1);
cellArray(1:nx/2) = {[0 0]}; %Position grid
cellArray(nx/2 + 1: nx) = {[-0.5 0.5]}; %Velocity grid
cellArray(nx + 1 : end) = {[-1 1]}; %Input grid
fields = string([sym('x',[nx 1]); input]);
grid = cell2struct(cellArray,fields);

lpvgridsys = LPVcore.nlss2lpvgridss(sys, grid);
% gam = lpvnorm(lpvgridsys,'h2');
%% Now compare against Mathworks example
M = 100;        % Number of masses
m = 1;          % Mass of individual blocks, kg
k1 = 0.5;       % Linear spring constant, N/m
k2 = 1;         % Cubic spring constant, N/m^3
b1 = 1;         % Linear damping constant, N/(m/s)
b2 = 0.2;       % Cubic damping constant, N/(m/s)^3

Tf = 50;
dt = 0.001;
tin = 0:dt:Tf;
Nt = numel(tin);

% Specify sinusoidal parameter (F) trajectory
amp = -1;
bias = 1;
freq = 0.5;
Ft = amp*cos(freq*tin)+bias;
dFt = [tin(:) Ft(:)];
F = 0;

% Specify input Force
dut = [tin(:) zeros(Nt,1)];
% dut(tin>=25,2) = 0;
U0 = 0;

% Simulate the full nonlinear model from Mathworks
simout_NL_MDOF = sim('MSD_NL_Javi_fix.slx');
tout_mathworks = simout_NL_MDOF.tout;
simout_NL_MDOF_data = squeeze(simout_NL_MDOF.simout.Data);

y_mathworks = simout_NL_MDOF_data(M,:)';
 
%% Generate Mathworks example with our much cooler MSD gen
% Input arguments
param = [1 0.5 1 0.5 1];

k = @(x) k2*x^3;
d = @(x,xdot) b2*xdot^3;
k_wall = @(x) k2*x^3; %They consider k_wall = k.
d_wall = @(x,xdot) b2*xdot^3; %They consider d_wall = d.

nl_index.springs = 1:M-1; nl_index.springs_wall = 1:M;
nl_index.dampers = 1:M-1; nl_index.dampers_wall = 1:M;

%MSD gen
tic
[fMSD,state,input] = generalized_msd_Javi(M,nl_index,k,d,k_wall,d_wall);
toc

%NLSS gen
f = @(x,u) fMSD(x,u,param);
h = @(x,u) x(M,:);
nx = length(state);
nu = length(input);
ny = 1;

sys = LPVcore.nlss(f,h,nx,nu,ny,0,true);

%SIM our custom MSD
x0 = zeros(nx,1);
u = @(tin) amp*cos(freq*tin)+bias;
[y,tout,x] = sim(sys,u,tin',x0,'ode45');

%% Plot Mathworks vs. Custom
outputFig = figure(1);
tiledlayout(2, 1, "TileSpacing","tight","Padding","tight");
nexttile; plot(tout,y,tout_mathworks,y_mathworks); grid on; legend("Generalized MSD","Mathworks"); xlabel("time [s]"); ylabel("Position q_M [m]");
nexttile; semilogy(tout, abs(y_mathworks-y)); title("Absolute error"); grid on;xlabel("time [s]");

error = y_mathworks-y;
thresh = 1e-8;
if norm(error) < thresh
    disp("We cool")
else
    fprintf("Minor bug alert. Error: %d",norm(error))
end




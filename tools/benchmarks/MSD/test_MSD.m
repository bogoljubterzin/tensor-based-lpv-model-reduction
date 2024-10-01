% Testing benchmark generator of NL MSD interconnected system

clear; close all; clc;

%% Add tools, such as LPVcore (custom branch), to path
addpath(genpath("../../Tools/"))
%%

% Input parameters:
M = 5; % Number of blocks: Nx = 2M
nl_springs = [4];
nl_springs_wall = [4,5];
nl_dampers = [];
nl_dampers_wall = [];
nl_index.nl_springs = nl_springs;
nl_index.nl_springs_wall = nl_springs_wall;
nl_index.nl_dampers = nl_dampers;
nl_index.nl_dampers_wall = nl_dampers_wall;

k = @(x) 0.5*x+3*x^3; % Spring nonlinearity function handle
d = @(x,xdot) xdot; % Damper nonlinearity function handle


load("MSDParameters.mat")
param = [params.m params.k0 params.d0];

%% Generate custom MSD NL equations
tic
[fMSD,state,input] = generalized_msd(M,nl_index,k,d);
toc

%% Nonlinear state-space model
f = @(x,u) fMSD(x,u,param);
h = @(x,u) x(M,:);

nx = length(state);
nu = length(input);
ny = 1;

sys = LPVcore.nlss(f,h,nx,nu,ny,0,true);

%% 

xx = sym('x',[nx,1]); uu = sym('u',[nu,1]);
f(xx,uu)

%% Convert nlss 2 lpvss
tic
[lpvsys, eta] = LPVcore.nlss2lpvss(sys, 'numerical', 'element');
toc

%% Mathworks example

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
% dut(tin>=25,2) = 0.1;
U0 = 0;

% Simulate the high-fidelity nonlinear model
simout_NL_MDOF = sim('MSD_NL');
tout = simout_NL_MDOF.tout;
simout_NL_MDOF_data = squeeze(simout_NL_MDOF.simout.Data)';

y_mathworks = simout_NL_MDOF_data(:,M);
 
%%

param = [1 0.5 1 0.5 1];

k = @(x) k1*x + k2*x^3;
d = @(x,xdot) b1*xdot + b2*xdot^3;

k_wall = @(x) k2*x^3;
d_wall = @(x,xdot) b2*xdot^3;

nl_springs = 1:M-1; nl_springs_wall = 1:M;
nl_dampers = 1:M-1; nl_dampers_wall = 1:M;
% nl_dampers = 1:M-1; nl_dampers_wall = 1:M;

nl_index.nl_springs = nl_springs;
nl_index.nl_springs_wall = nl_springs_wall;
nl_index.nl_dampers = nl_dampers;
nl_index.nl_dampers_wall = nl_dampers_wall;

tic
[fMSD,state,input] = generalized_msd(M,nl_index,k,d);
toc

f = @(x,u) fMSD(x,u,param);
h = @(x,u) x(M,:);

nx = length(state);
nu = length(input);
ny = 1;

x0 = zeros(nx,1);
%%
sys = LPVcore.nlss(f,h,nx,nu,ny,0,true);

u = @(tin) amp*cos(freq*tin) + bias;
[y,tout,x] = sim(sys,u,tin',x0,'ode45');

%%
outputFig = figure(1);
tiledlayout(2, 1, "TileSpacing","tight","Padding","tight");
nexttile; plot(tout,y,tout,y_mathworks); grid on; legend("Generalized MSD","Mathworks"); xlabel("time [s]"); ylabel("Position q_M [m]");
nexttile; semilogy(tout, abs(y_mathworks-y)); title("Absolute error"); grid on;xlabel("time [s]");

error = y_mathworks-y;
thresh = 1e-8;
if norm(error) < thresh
    disp("Error insignificant")
else
    fprintf("Minor bug alert. Error: %d",norm(error))
end




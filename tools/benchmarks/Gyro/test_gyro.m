clear; close all;

%% Add tools, such as LPVcore (custom branch), to path
addpath(genpath("../../Tools/"))

%% Load gyro physical parameters
load("GyroParameters.mat");

param = [params.Ka params.Ib params.Jb ...
         params.Kb params.Ic params.Jc ...
         params.Kc params.Id params.Jd...
         params.Km1 params.Km2 params.Km3 params.Km4 ...
         params.fv1 params.fv2 params.fv3 params.fv4]';

%% Generate Custom gyro equations
actuators = [1;1;0;0]; locked = [0;0;1;0];
tic
[fGyro, f_dec_, e_, state, input] = generalized_gyro(actuators, locked);
toc

%% Nonlinear state-space model
f = @(x,u) fGyro(x, u, param);
h = @(x,u) [x(1,:);x(2,:);x(4,:)];

nx = length(state);
nu = length(input);
ny = 3;
sys = LPVcore.nlss(f, h, nx, nu, ny, 0, true);

%% Convert nlss2lpvss
tic
[lpvsys, eta] = LPVcore.nlss2lpvss(sys, 'numerical', 'element');
toc

%It works!

%% Now just check if the generalized gyro EoM generator provides the expected outputs
% Discretize the system
Ts = 1e-3;
sysd = c2d(sys, Ts, "rk4");

%% Simulation parameters
x0 = zeros(nx,1);
x0(end-2) = 40;         % intial condition
Tend = 10;                  % simulation end time
t = 0:Ts:Tend;
N = numel(t);

u1 = .5 * sin(t) + 0.3;            % i1 input signal
u2 = 1 * sin(2 * pi * t - pi);     % i2 input signal
% u3 = .5 * sin(2 * pi * t - pi);     % i3 input signal
u = [u1;u2];

%% Simulation nl (DT)
x = zeros(sysd.Nx, numel(t));
x(:,1) = x0;

for i = 1:numel(t)-1
    x(:,i+1) = sysd.f(x(:,i), u(:,i));
end

%% Plot
stateFig = figure;
tiledlayout(nx, 1, "TileSpacing","tight","Padding","tight");
for i = 1:nx
    nexttile
    plot(t, x(i,:));
    ylabel(sprintf('$x_%i$', i), 'Interpreter', 'latex');
    grid on;
end
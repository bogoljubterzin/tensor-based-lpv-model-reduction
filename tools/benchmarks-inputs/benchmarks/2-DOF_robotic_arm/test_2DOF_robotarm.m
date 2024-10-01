clear; close all;

%% Add tools, such as LPVcore (custom branch), to path
addpath(genpath("../../Tools/"))

%% Load 2-DOF robotic arm parametes
load("2-DOF_robotic_arm_parameters.mat");
param = [params.a params.b params.c params.d params.e params.f params.n];

%% Generate Custom 2-DOF robotic arm NL equations
tic
[f_, f_dec_, e_, state, input] = generalized_2DOF_robotarm;
toc

%% NL-SS model
f = @(x,u) f_(x, u, param);
h = @(x,u) [x(1,:);x(2,:)]; %Just positions for example

nx = length(state);
nu = length(input);
ny = 2;
sys = LPVcore.nlss(f, h, nx, nu, ny, 0, true);

%% Convert nlss 2 lpvss
xx = sym('x', [nx,1]); uu = sym('u', [nu,1]);
tic
[lpvsys, eta] = LPVcore.nlss2lpvss(sys, 'numerical', 'element');
toc
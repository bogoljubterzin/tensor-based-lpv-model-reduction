%% Comparison with internship
addpath(genpath("../"))
%%
Ts = 1e-2; Tend = 10;
time = [0:Ts:Tend]';

input = 1*sin(3*time) + 0.5*cos(2.5*time);
M = 5;
system = example_lpv_model(M);
data = simulate_example(system, input, time,M);
y_internship = data.yNL;
%%

M = 5; % Number of blocks: Nx = 2M
nl_springs = 1:4;
nl_springs_wall = 1:5;
nl_dampers = [];
nl_dampers_wall = [];
nonlinearities.nl_springs = nl_springs;
nonlinearities.nl_springs_wall = nl_springs_wall;
nonlinearities.nl_dampers = nl_dampers;
nonlinearities.nl_dampers_wall = nl_dampers_wall;

k = @(x) 0.5*x+x^3; % Spring nonlinearity function handle
d = @(x,xdot) xdot; % Damper nonlinearity function handle


load("MSDParameters.mat")
param = [params.m params.k0 params.d0];

tic
[fMSD,state,input1] = generalized_msd(M,nonlinearities,k,d);
toc
f = @(x,u) fMSD(x,u,param);
h = @(x,u) x(M,:);

nx = length(state);
nu = length(input1);
ny = 1;

sys = LPVcore.nlss(f,h,nx,nu,ny,0,true);
[lpvsys, eta] = LPVcore.nlss2lpvss(sys, 'numerical', 'element');

x0 = zeros(nx,1);
[y,tout,x] = sim(sys,input,time,x0);

%%
outputFig = figure(1);
tiledlayout(2, 1, "TileSpacing","tight","Padding","tight");
nexttile; plot(tout,y,tout,y_internship); grid on; legend("Generalized MSD","Internship");
nexttile; semilogy(tout, abs(y_internship-y)); title("Absolute error"); grid on;

error = y_internship-y;
thresh = 1e-8;
if norm(error) < thresh
    disp("Error insignificant")
else
    fprintf("Minor bug alert. Error: %d",norm(error))
end

%%

xx = sym('x',[nx,1]); uu = sym('u',[nu,1]);
f(xx,uu)



% This file contains the implementation of tensor-based moment matching
% with additional analysis in terms of singular values and singular vector,
% here we do not necesserily take all singular vectors, but we can choose
% to discard some of them to obtain lower-order models

tensSS = tensRobot;
% tensSS = tensMSD1;

horizon = 4;
N = horizon+1;
decomp = 'hosvd';
mode = 'T';
tol = 1e-6;

Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);

%% Lets see what we get for N=0

W0 = Wn{1};

Kw = min(modrank(W0));
[xW,gainsW,iterW,XaW] = TSVDND_SW(W0,Kw);
Ew0_tsvd = norm(double(W0)-XaW,'fro');

V0t = xW{1};
W0t = xW{1};
Z0t = xW{2};

%%

W0_hosvd = hosvd(W0,tol);
svW0 = hosvd_sv(W0_hosvd);
Ew0_hosvd = norm(double(W0)-double(W0_hosvd),'fro');

V0h = W0_hosvd{1};
W0h = V0h;
Z0h = W0_hosvd{2};

%%
constantTerm = 0;
redRobot0t = tensRobot.PetrovGalerkinLPV(W0t,V0t,Z0t,constantTerm);
redRobot0h = tensRobot.PetrovGalerkinLPV(W0h,V0h,Z0h,constantTerm);

%%
[~,x,y] = tensRobot.simulateSS(uRobot,t,zeros(tensRobot.Nx,1),1);
[~,xr0t,yr0t] = redRobot0t.simulateSS(uRobot,t,V0t'*zeros(tensRobot.Nx,1),constantTerm);
[~,xr0h,yr0h] = redRobot0h.simulateSS(uRobot,t,V0h'*zeros(tensRobot.Nx,1),constantTerm);
%%
x0h = V0h*xr0h;
x0t = V0t*xr0t;
%%
e0h = norm(y-yr0h)
e0t = norm(y-yr0t)
e0hx = norm(x-x0h)/norm(x)
x0tx = norm(x-x0t)/norm(x)
%% output plot
FigRobotRed0 = figure(1);
tiledlayout(2,1,'TileSpacing','tight');
nexttile;
plot(t, y(1,:)); grid on; hold on;
plot(t,yr0t(1,:));
plot(t,yr0h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, y(2,:)); grid on; hold on;
plot(t,yr0t(2,:));
plot(t,yr0h(2,:));
legend("FOM", "TSVD","HOSVD");

%% reconstructed state plot
FigRobotRedx0 = figure(2);
tiledlayout(2,2,'TileSpacing','tight');
nexttile;
plot(t, x(1,:)); grid on; hold on;
plot(t,x0t(1,:));
plot(t,x0h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(2,:)); grid on; hold on;
plot(t,x0t(2,:));
plot(t,x0h(2,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(3,:)); grid on; hold on;
plot(t,x0t(3,:));
plot(t,x0h(3,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(4,:)); grid on; hold on;
plot(t,x0t(4,:));
plot(t,x0h(4,:));
legend("FOM", "TSVD","HOSVD");


%% %%%%%%%%%%%%%%%%%%%%%%%%%% N = 1
W1 = Wn{2};

Kw = min(modrank(W1));
[xW1,gainsW1,iterW1,XaW1] = TSVDND_SW(W1,Kw);
Ew1_tsvd = norm(double(W1)-XaW1,'fro');

V1t0 = xW1{1};
W1t0 = xW1{1};
Z1t0 = [xW1{2},xW1{3}];

%%

W1_hosvd = hosvd(W1,tol);
svW1 = hosvd_sv(W1_hosvd);
Ew1_hosvd = norm(double(W1)-double(W1_hosvd),'fro');

V1h0 = W1_hosvd{1};
W1h0 = V1h;
Z1h0 = [W1_hosvd{2}, W1_hosvd{3}];



%%

V1t = orth([V0t,V1t0]);
W1t = V1t;
Z1t = orth([Z0t,Z1t0]);

V1h = orth([V0h,V1h0]);
W1h = V1h;
Z1h = orth([Z0h,Z1h0]);

%%
constantTerm = 0;
redRobot1t = tensRobot.PetrovGalerkinLPV(W1t,V1t,Z1t,constantTerm);
redRobot1h = tensRobot.PetrovGalerkinLPV(W1h,V1h,Z1h,constantTerm);

%%
[~,xr1t,yr1t] = redRobot1t.simulateSS(uRobot,t,V1t'*zeros(tensRobot.Nx,1),constantTerm);
[~,xr1h,yr1h] = redRobot1h.simulateSS(uRobot,t,V1h'*zeros(tensRobot.Nx,1),constantTerm);
%%
x1h = V1h*xr1h;
x1t = V1t*xr1t;
%%
e1h = norm(y-yr1h)
e1t = norm(y-yr1t)
e1hx = norm(x-x1h)/norm(x)
x1tx = norm(x-x1t)/norm(x)
%% output plot
FigRobotRed1 = figure(3);
tiledlayout(2,1,'TileSpacing','tight');
nexttile;
plot(t, y(1,:)); grid on; hold on;
plot(t,yr1t(1,:));
plot(t,yr1h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, y(2,:)); grid on; hold on;
plot(t,yr1t(2,:));
plot(t,yr1h(2,:));
legend("FOM", "TSVD","HOSVD");

%% reconstructed state plot
FigRobotRedx1 = figure(4);
tiledlayout(2,2,'TileSpacing','tight');
nexttile;
plot(t, x(1,:)); grid on; hold on;
plot(t,x1t(1,:));
plot(t,x1h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(2,:)); grid on; hold on;
plot(t,x1t(2,:));
plot(t,x1h(2,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(3,:)); grid on; hold on;
plot(t,x1t(3,:));
plot(t,x1h(3,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(4,:)); grid on; hold on;
plot(t,x1t(4,:));
plot(t,x1h(4,:));
legend("FOM", "TSVD","HOSVD");

%% %%%%%%%%%%%%%%%%%%%%%%%%%% N = 2
W2 = Wn{3};

Kw = min(modrank(W2));
[xW2,gainsW2,iterW2,XaW2] = TSVDND_SW(W2,Kw);
Ew2_tsvd = norm(double(W2)-XaW2,'fro');

V2t0 = xW2{1};
W2t0 = xW2{1};
Z2t0 = [xW2{2},xW2{3},xW2{4}];

%%

W2_hosvd = hosvd(W2,tol);
svW2 = hosvd_sv(W2_hosvd);
Ew2_hosvd = norm(double(W2)-double(W2_hosvd),'fro');

V2h0 = W2_hosvd{1};
W2h0 = V2h0;
Z2h0 = [W2_hosvd{2}, W2_hosvd{3},W2_hosvd{4}];

%%

V2t = orth([V0t,V1t,V2t0]);
W2t = V2t;
Z2t = orth([Z0t,Z1t,Z2t0]);

V2h = orth([V0h,V1h,V2h0]);
W2h = V2h;
Z2h = orth([Z0h,Z1h,Z2h0]);

%%
constantTerm = 0;
redRobot2t = tensRobot.PetrovGalerkinLPV(W2t,V2t,Z2t,constantTerm);
redRobot2h = tensRobot.PetrovGalerkinLPV(W2h,V2h,Z2h,constantTerm);

%%
[~,xr2t,yr2t] = redRobot2t.simulateSS(uRobot,t,V2t'*zeros(tensRobot.Nx,1),constantTerm);
[~,xr2h,yr2h] = redRobot2h.simulateSS(uRobot,t,V2h'*zeros(tensRobot.Nx,1),constantTerm);
%%
x2h = V2h*xr2h;
x2t = V2t*xr2t;
%%
e2h = norm(y-yr2h)
e2t = norm(y-yr2t)
e2hx = norm(x-x2h)/norm(x)
x2tx = norm(x-x2t)/norm(x)
%% output plot
FigRobotRed2 = figure(5);
tiledlayout(2,1,'TileSpacing','tight');
nexttile;
plot(t, y(1,:)); grid on; hold on;
plot(t,yr2t(1,:));
plot(t,yr2h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, y(2,:)); grid on; hold on;
plot(t,yr2t(2,:));
plot(t,yr2h(2,:));
legend("FOM", "TSVD","HOSVD");

%% reconstructed state plot
FigRobotRedx2 = figure(6);
tiledlayout(2,2,'TileSpacing','tight');
nexttile;
plot(t, x(1,:)); grid on; hold on;
plot(t,x2t(1,:));
plot(t,x2h(1,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(2,:)); grid on; hold on;
plot(t,x2t(2,:));
plot(t,x2h(2,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(3,:)); grid on; hold on;
plot(t,x2t(3,:));
plot(t,x2h(3,:));
legend("FOM", "TSVD","HOSVD");
nexttile;
plot(t, x(4,:)); grid on; hold on;
plot(t,x2t(4,:));
plot(t,x2h(4,:));
legend("FOM", "TSVD","HOSVD");

%% %%%%%%%%%%%%%%%%%%%%%%%%%% N = 3
W3 = Wn{4};

Kw = min(modrank(W3));
[xW3,gainsW3,iterW3,XaW3] = TSVDND_SW(W3,Kw);
Ew3_tsvd = norm(double(W3)-XaW3,'fro');

V3t0 = xW3{1};
W3t0 = xW3{1};
Z3t0 = [xW3{2},xW3{3},xW3{4},xW3{5}];

%%

W3_hosvd = hosvd(W3,tol);
svW3 = hosvd_sv(W3_hosvd);
Ew3_hosvd = norm(double(W3)-double(W3_hosvd),'fro');

V3h0 = W3_hosvd{1};
W3h0 = V3h0;
Z3h0 = [W3_hosvd{2}, W3_hosvd{3},W3_hosvd{4},W3_hosvd{5}];

%%

V3t = orth([V0t,V1t,V2t,V3t0]);
W3t = V3t;
Z3t = orth([Z0t,Z1t,Z2t,Z3t0]);

V3h = orth([V0h,V1h,V2h,V3h0]);
W3h = V3h;
Z3h = orth([Z0h,Z1h,Z2h,Z3h0]);

%%
constantTerm = 0;
redRobot3t = tensRobot.PetrovGalerkinLPV(W3t,V3t,Z3t,constantTerm);
redRobot3h = tensRobot.PetrovGalerkinLPV(W3h,V3h,Z3h,constantTerm);

%%
[~,xr3t,yr3t] = redRobot3t.simulateSS(uRobot,t,V3t'*zeros(tensRobot.Nx,1),constantTerm);
[~,xr3h,yr3h] = redRobot3h.simulateSS(uRobot,t,V3h'*zeros(tensRobot.Nx,1),constantTerm);
%%
x3h = V3h*xr3h;
x3t = V3t*xr3t;
%%
e3h = norm(y-yr3h)
e3t = norm(y-yr3t)
e3hx = norm(x-x3h)/norm(x)
e3tx = norm(x-x3t)/norm(x)

%%

eh = [e0h, e1h, e2h, e3h]/norm(y)
et = [e0t, e1t, e2t, e3t]/norm(y)

figure(123);
plot(0:3,eh,'rx','LineWidth',2); hold on; grid on;
plot(0:3,et,'bo','LineWidth',2)
%% lower order approximation
constantTerm = 0;
rx = 3;
rp = 7;
redRobot1 = tensRobot.PetrovGalerkinLPV(W3t(:,1:rx),V3t(:,1:rx),Z3t(:,1:rp),constantTerm);
redRobot2 = tensRobot.PetrovGalerkinLPV(W3h(:,1:rx),V3h(:,1:rx),Z3h(:,1:rp),constantTerm);

%%
[~,xrt,yrt] = redRobot1.simulateSS(uRobot,t,V3t(:,1:rx)'*zeros(tensRobot.Nx,1),constantTerm);
[~,xrh,yrh] = redRobot2.simulateSS(uRobot,t,V3h(:,1:rx)'*zeros(tensRobot.Nx,1),constantTerm);
%%
figure,
plot(t,y(1,:)); hold on; grid on;
plot(t, yrt(1,:));
plot(t,yrh(1,:))

figure,
plot(t,y(2,:)); hold on; grid on;
plot(t, yrt(2,:));
plot(t,yrh(2,:))
%%
set(groot, 'DefaultTextInterpreter', 'latex');        % For text objects
set(groot, 'DefaultLegendInterpreter', 'latex');      % For legends
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % For axis tick labels
set(groot, 'DefaultColorbarTickLabelInterpreter', 'latex'); % For colorbar tick labels

%%
FigRoboticArmReachabilityY1 = figure(101);
clrs = lines(5);
tiledlayout(2,2,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,y(1,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr0t(1,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr0h(1,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_1$");
legend('FOM','TSVD','HOSVD');
title("$N = 0$")
nexttile;
plot(t,y(1,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr1t(1,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr1h(1,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_1$");
legend('FOM','TSVD','HOSVD');
title("$N = 1$")

nexttile;
plot(t,y(1,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr2t(1,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr2h(1,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_1$");
legend('FOM','TSVD','HOSVD');
title("$N = 2$")

nexttile;
plot(t,y(1,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr3t(1,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr3h(1,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_1$");
legend('FOM','TSVD','HOSVD');
title("$N = 3$")

%%
FigRoboticArmReachabilityY2 = figure(102);
clrs = lines(5);
tiledlayout(2,2,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr0t(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr0h(2,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_2$");
legend('FOM','TSVD','HOSVD');
title("$N = 0$")
nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr1t(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr1h(2,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_2$");
legend('FOM','TSVD','HOSVD');
title("$N = 1$")

nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr2t(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr2h(2,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_2$");
legend('FOM','TSVD','HOSVD');
title("$N = 2$")

nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr3t(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr3h(2,:),'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y_2$");
legend('FOM','TSVD','HOSVD');
title("$N = 3$")
%%
plot(t,yr1t(1,:),'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yr2t(1,:),'Color',clrs(4,:),'LineWidth',1.2);
plot(t,yr3t(1,:),'Color',clrs(5,:),'LineWidth',1.2);
nexttile;
plot(t,y(1,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr0h(1,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr1h(1,:),'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yr2h(1,:),'Color',clrs(4,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("y_1");
legend('FOM','TSVD','HOSVD');
title("N = 1")
nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr0t(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr1t(2,:),'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yr2t(2,:),'Color',clrs(4,:),'LineWidth',1.2);
plot(t,yr3t(2,:),'Color',clrs(5,:),'LineWidth',1.2);
nexttile;
plot(t,y(2,:),'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr0h(2,:),'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr1h(2,:),'Color',clrs(3,:),'LineWidth',1.2);
plot(t,yr2h(2,:),'Color',clrs(4,:),'LineWidth',1.2);
%%


Vr_tsvd = cell(N,1);
Zr_tsvd = cell(N,1);

Vq_tsvd = cell(N,1);
Zq_tsvd = cell(N,1);

gainsVr_tsvd = cell(N,1);
gainsZr_tsvd = cell(N,1);
gainsVq_tsvd = cell(N,1);
gainsZq_tsvd = cell(N,1);



sv_sched = cell(N,1);
sv_states = cell(N,1);

for n = 1:N
    K = min(modrank(Wn{n}));
    [x,gains,iterations,Xa] = TSVDND_SW(Wn{n},K);
    % separate singular vectors: 

    % for the state projection V, we gather the
    % vectors from first dimension x{1}
    Vr_tsvd{n} = x{1};

    % for the scheduling projection Z, we gather vectors from dimensions
    % x{2},..., x{n}
    for j = 2:n+1
        Zr_tsvd{n} = [Zr_tsvd{n}, x{j}];
        % Tsched{j-1} = null(x{j}.')
    end

    % save info from decomposition from every level
    tsvd_info{n}.x = x;
    tsvd_info{n}.gains = gains;
    tsvd_info{n}.iterations = iterations;


    % Wn_hosvd = hosvd(tensor(Wn{n}),tol);
    % gainsWn_hosvd = hosvd_sv(Wn_hosvd);
    % 
    % gainsVr_tsvd{n} = gainsWn_hosvd{1}(1:size(Wn_hosvd{1},2));
    % Vr_tsvd{n} = Wn_hosvd{1};
    % gainsZr_tsvd{n} = [];
    % for j = 2:n+1
    %     Zr_tsvd{n} = [Zr_tsvd{n}, Wn_hosvd{j}];
    %     gainsZr_tsvd{n} = [gainsZr_tsvd{n}; gainsWn_hosvd{j}(1:size(Wn_hosvd{j},2))];
    % end
    % 
    % Qn_hosvd = hosvd(tensor(Qn{n}),tol);
    % gainsQn_hosvd = hosvd_sv(Qn_hosvd);
    % 
    % gainsVq_tsvd{n} = gainsQn_hosvd{n+2}(1:size(Qn_hosvd{n+2},2));
    % Vq_tsvd{n} = Qn_hosvd{n+2};
    % 
    % gainsVq_tsvd{n} = [];
    % 
    % for j = 2:n+1
    %     Zq_tsvd{n} = [Zq_tsvd{n}, Qn_hosvd{j}];
    %     gainsZq_tsvd{n} = [gainsZq_tsvd{n}; gainsQn_hosvd{j}(1:size(Qn_hosvd{j},2))];
    % end
end
%%
Vr = [];
Zr = [];

svVr = [];
svZr = [];
for n = 1:N
    Vr = [Vr, Vr_tsvd{n}];
    svVr = [svVr, gainsVr_tsvd{n}.'];

    Zr = [Zr, Zr_tsvd{n}];
    svZr = [svZr, gainsZr_tsvd{n}.'];
end

Vr1 = orth(Vr);
Zr1 = orth(Zr);

% now i want out of Vr_hosvd and gainsVr_hosvd to obtain Vr
% Zr_hosvd and gainsZr_hosvd -> Zr, etc.

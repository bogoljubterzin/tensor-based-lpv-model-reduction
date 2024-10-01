[nlss_robot, lpvsys_robot, eta_robot] = gen_robot
%%
Ts = 0.005;
lpvsysdt_robot = c2d(lpvsys_robot,Ts);
tensRobot = tensorSS(lpvsysdt_robot,eta_robot);
%%
x0 = zeros(tensRobot.Nx,1);
t = 0:0.01:20;
uRobot = [sin(t); cos(t)];
[tout,x_robot,y_robot] = tensRobot.simulateSS(uRobot,t,x0,1)
%%
N = 3;
decomp = 'hosvd';
mode = 'R';
[redRobot, Tx, Tp, info] = tensRobot.lpvTensMM(N, decomp, mode)

%%
constantTerm = 0;
[tout,xr1,yr1] = redRobot.simulateSS(uRobot,t,Tx'*x0,constantTerm);
%%
figure(1);
clf(1)
plot(t,y_robot(1,:)); hold on; grid on;
plot(t, yr1(1,:));

figure(2);
clf(2)
plot(t,y_robot(2,:)); hold on; grid on;
plot(t,yr1(2,:));


%% Define non-minimal representation
V = [eye(tensRobot.Nx), eye(tensRobot.Nx)];
Z = eye(tensRobot.Np+1);
W = V;
% 
Ar = double(ttm(tensRobot.A,{pinv(W'*V)*W', V', Z'},[1,2,3]));
Br = double(ttm(tensRobot.B,{pinv(W'*V)*W',Z'},[1,3]));
Cr = double(ttm(tensRobot.C,{V',Z'},[2,3]));
Dr = double(ttm(tensRobot.D,{Z'},3));
eta_red = @(xr,u) pinv(Z)*[1; tensRobot.eta_map(V*xr,u)];
tensRobotNonMinimal = tensorSS(Ar,Br,Cr,Dr,eta_red,tensRobot.Ts)

%%
[tout,xNonMin,yNonMin] = tensRobotNonMinimal.simulateSS(uRobot,t,V'*x0,0)


%%
figure(1);
plot(t,y_robot(1,:)); hold on; grid on;
plot(t, yNonMin(1,:));

figure(2);
plot(t,y_robot(2,:)); hold on; grid on;
plot(t,yNonMin(2,:));

%%
horizon = 2;
N = horizon+1;
decomp = 'hosvd';
mode = 'T';
tol = 1e-6;

Wn = ReachabilityTensors(tensRobotNonMinimal.A,tensRobotNonMinimal.B,N);
Qn = ObservabilityTensors(tensRobotNonMinimal.A,tensRobotNonMinimal.C,N);


Vr_hosvd = cell(N,1);
Zr_hosvd = cell(N,1);
Vq_hosvd = cell(N,1);
Zq_hosvd = cell(N,1);
gainsWn_hosvd = cell(N,1);
gainsQn_hosvd = cell(N,1);

sv_sched = cell(N,1);
sv_states = cell(N,1);
for n = 1:N
    Wn_hosvd = hosvd(tensor(Wn{n}),tol);
    gainsWn_hosvd{n} = hosvd_sv(Wn_hosvd);
    sv_states{n} = gainsWn_hosvd{n}{1};

    Vr_hosvd{n} = Wn_hosvd{1};
    sv_sched{n} = [];
    for j = 2:n+1
        Zr_hosvd{n} = [Zr_hosvd{n}, Wn_hosvd{j}];
        sv_sched{n} = [sv_sched{n}; gainsWn_hosvd{n}{j}]
    end
    
    Qn_hosvd = hosvd(tensor(Qn{n}),tol);
    Vq_hosvd{n} = Qn_hosvd{n+2};
    for j = 2:n+1
        Zq_hosvd{n} = [Zq_hosvd{n}, Qn_hosvd{j}];
    end
    gainsQn_hosvd{n} = hosvd_sv(Qn_hosvd);
end

%%
tol = 1e-4;
sv = [];
Tx = [];
for n = 1:N
    ni = size(Vr_hosvd{n},2);
    sv = [sv, gainsWn_hosvd{n}{1}(1:ni)'];
    Tx = [Tx, Vr_hosvd{n}];
end
sv_new = [];
Tx_new = [];
sv_imp = [];
cnt = 1;
for i = 1:length(sv)
    if sv(i)/sum(sv) > tol
        sv_imp(cnt) = sv(i)/sum(sv)
        sv_new = [sv_new, sv(i)];
        Tx_new = [Tx_new, Tx(:,i)];
        cnt = cnt+1;
    end
end


%%

Tp =[];
sv_all = [];
for n = 1:N
    ni = size(Zr_hosvd{n},2);
    for i = 2:n+1
        sv_all = [sv_all,gainsWn_hosvd{n}{i}'];
    end
    Tp = [Tp, Zr_hosvd{n}];
end

Tp_new = [];
sv_new = [];
for i = 1:length(sv_all)
    if sv_all(i)/sum(sv_all) > tol
        Tp_new = [Tp_new, Tp(:,i)];
        sv_new = [sv_new, sv_all(i)];
    end
end
%     sv = [sv, gainsWn_hosvd{n}{1}(1:ni)'];
%     Tx = [Tx, Vr_hosvd{n}];
% end
% sv_new = [];
% Tx_new = [];
% sv_imp = [];
% cnt = 1;
% for i = 1:length(sv)
%     if sv(i)/sum(sv) > tol
%         sv_imp(cnt) = sv(i)/sum(sv)
%         sv_new = [sv_new, sv(i)];
%         Tx_new = [Tx_new, Tx(:,i)];
%         cnt = cnt+1;
%     end
% end
%%
[redRobot, Tx, Tp, info] = tensRobotNonMinimal.lpvTensMM(N, decomp, mode)
constantTerm = 0;
x0 = zeros(tensRobotNonMinimal.Nx,1);
[tout,xr1,yr1] = redRobot.simulateSS(uRobot,t,Tx'*x0,constantTerm);

%%
close all
figure(3);
plot(t,yNonMin(1,:)); hold on; grid on;
plot(t, yr1(1,:));

figure(4);
plot(t,yNonMin(2,:)); hold on; grid on;
plot(t,yr1(2,:));
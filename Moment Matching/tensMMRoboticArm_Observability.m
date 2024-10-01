%% robot arm observability
%% N=0

Q0 = Qn{1};

Kw = min(modrank(Q0));
[xWq,gainsWq,iterWq,XaWq] = TSVDND_SW(Q0,Kw);
Ew0_tsvd = norm(double(Q0)-XaWq,'fro');

V0tq = xWq{3};
W0tq = xWq{3};
Z0tq = xWq{2};


Q0_hosvd = hosvd(Q0,tol);
svQ0 = hosvd_sv(Q0_hosvd);
EQ0_hosvd = norm(double(Q0)-double(Q0_hosvd),'fro');

V0hq = Q0_hosvd{3};
W0hq = V0hq;
Z0hq = Q0_hosvd{2};

%%
constantTerm = 0;
redRobot0tq = tensRobot.PetrovGalerkinLPV(W0tq,V0tq,Z0tq,constantTerm);
redRobot0hq = tensRobot.PetrovGalerkinLPV(W0hq,V0hq,Z0hq,constantTerm);

%%
% [~,x,y] = tensRobot.simulateSS(uRobot,t,zeros(tensRobot.Nx,1),1);
[~,xr0tq,yr0tq] = redRobot0tq.simulateSS(uRobot,t,V0tq'*zeros(tensRobot.Nx,1),constantTerm);
[~,xr0hq,yr0hq] = redRobot0hq.simulateSS(uRobot,t,V0hq'*zeros(tensRobot.Nx,1),constantTerm);
%%
x0hq = V0hq*xr0hq;
x0tq = V0tq*xr0tq;
%%
e0h = norm(y-yr0hq)
e0t = norm(y-yr0tq)
e0hx = norm(x-x0hq)/norm(x)
x0tx = norm(x-x0tq)/norm(x)

%% N=0 ...
[redRobot0tq,V0tq, Z0tq, info0tq] = tensRobot.lpvTensMM(0,'tsvd','T')
[~,x0tq,y0tq] = redRobot0tq.simulateSS(uRobot,t,V0tq'*zeros(tensRobot.Nx,1),constantTerm);
e0tq = norm(y-y0tq)

[redRobot0hq,V0hq, Z0hq, info0hq] = tensRobot.lpvTensMM(0,'hosvd','T')
[~,x0hq,y0hq] = redRobot0hq.simulateSS(uRobot,t,V0hq'*zeros(tensRobot.Nx,1),constantTerm);
e0hq = norm(y-y0hq)
%% N=1
[redRobot1tq,V1tq, Z1tq, info1tq] = tensRobot.lpvTensMM(1,'tsvd','T')
[~,x1tq,y1tq] = redRobot1tq.simulateSS(uRobot,t,V1tq'*zeros(tensRobot.Nx,1),constantTerm);
e1tq = norm(y-y1tq)

[redRobot1hq,V1hq, Z1hq, info1hq] = tensRobot.lpvTensMM(1,'hosvd','T')
[~,x1hq,y1hq] = redRobot1hq.simulateSS(uRobot,t,V1hq'*zeros(tensRobot.Nx,1),constantTerm);
e1hq = norm(y-y1hq)
%% N=2
[redRobot2tq,V2tq, Z2tq, info2tq] = tensRobot.lpvTensMM(2,'tsvd','T')
[~,x2tq,y2tq] = redRobot2tq.simulateSS(uRobot,t,V2tq'*zeros(tensRobot.Nx,1),constantTerm);
e2tq = norm(y-y2tq)

[redRobot2hq,V2hq, Z2hq, info2hq] = tensRobot.lpvTensMM(2,'hosvd','T')
[~,x2hq,y2hq] = redRobot2hq.simulateSS(uRobot,t,V2hq'*zeros(tensRobot.Nx,1),constantTerm);
e2hq = norm(y-y2hq)

%% N=3
[redRobot3tq,V3tq, Z3tq, info3tq] = tensRobot.lpvTensMM(3,'tsvd','T')
[~,x3tq,y3tq] = redRobot3tq.simulateSS(uRobot,t,V3tq'*zeros(tensRobot.Nx,1),constantTerm);
e3tq = norm(y-y3tq)

%%
[redRobot3hq,V3hq, Z3hq, info3hq] = tensRobot.lpvTensMM(4,'hosvd','T')
[~,x3hq,y3hq] = redRobot3hq.simulateSS(uRobot,t,V3hq'*zeros(tensRobot.Nx,1),constantTerm);
e3hq = norm(y-y3hq)


%%

eHhosvd = []

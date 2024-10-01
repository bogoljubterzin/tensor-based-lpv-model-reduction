%% simulate LPV model and Nonlinear function. Gather data 
% Models has to contain all nonlinear springs and created with
% example_lpv_model.m
function data_structure = simulate_example(system_lpv, input, time,M)

    d = 1; % Linear damping coefficient
    k1 = 0.5; % Linear spring coefficient
    k2 = 1; % Cubic spring coefficient
    m = 1; % Mass
    parameters = [d,k1,k2,m,M];

    data_structure.u = input;

    nx = 2*M; 
    nu = 1; % Force applied to the sybsystem M
    ny = 1; % Position qM

    % Simulating the nonlinear system

    x0 = zeros(nx,1);
    % using fixed-step solver ode4
    [tout,xNL] = ode45(@(t,x)f_nl(t,x,input,time,parameters),time,x0); 
    yNL = xNL(:,M);

    q = xNL(:,1:M);
    p = sched_map(q);
    [yLPV,t,xLPV] = lsim(system_lpv,p,input,time);

    data_structure.yNL = yNL;
    data_structure.p = p;
    data_structure.yLPV = yLPV;
    data_structure.xLPV = xLPV;

end


% Nonlinear system equation
function xdot = f_nl(t,x,u,time,parameters)
    u_int = interp1(time, u, t);
    d = parameters(1); k1 = parameters(2); k2 = parameters(3); m = parameters(4); M = parameters(5);

    dx1 = x(M+1:2*M);
    dx2 = zeros(M,1);
    for i = 2:M-1
        dx2(i) = -1/m*(Fij(x,i,i-1,parameters) + Fij(x,i,-1,parameters) + Fij(x,i,i+1,parameters));
    end
    dx2(1) = -1/m*(Fij(x,1,-1,parameters) + Fij(x,1,2,parameters));
    dx2(M) = -1/m*(Fij(x,M,-1,parameters) + Fij(x,M,M-1,parameters) - u_int);
    xdot = [dx1;dx2];
end

% Calculating forces
function y = Fij(x,i,j,parameters)
    M = parameters(5);
    if j == -1
        qj = 0; qjdot = 0;
    else
        qj = x(j);
        qjdot = x(M+j);
    end
    qi = x(i); qidot = x(M+i); 
    d = parameters(1); k1 = parameters(2); k2 = parameters(3);
    y = d*(qidot-qjdot) + k1*(qi-qj) + k2*(qi-qj)^3;
end


% trajectory of NL system -> scheduling variables
function p = sched_map(q)
    N = size(q,1); M = size(q,2);
    p = zeros(N,2*M-1);
    a = []; b = [];

    for i = 1:M
        qi = q(:,i);
        a = [a qi.^2];
    end

    for i = 2:M
        qi = q(:,i);
        qiminus1 = q(:,i-1);
        b = [b (qi - qiminus1).^2];
    end    
    p(:,1:M) = a;
    p(:,M+1:2*M-1) = b;
end
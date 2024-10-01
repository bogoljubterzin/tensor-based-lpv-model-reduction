function [nlss, lpvsys, eta] = gen_robot_manual
    %Use this for the paper, as it is the one with minimal amount of sched
    %vars!
    %Using the exact same model as used in LPVDNNRED!
    a = 5.6794;     b = 1.473;      c = 1.7985;
    d = 4e-1;       e = 4e-1;       f = 2;
    n = 1;
    param = [a,b,c,d,e,f,n];
    [f_, ~, ~, state, input] = generalized_2DOF_robotarm;
    f = @(x,u) f_(x, u, param);
    h = @(x,u) [x(1,:);x(2,:)]; %Just positions for example
    nx = length(state);
    nu = length(input);
    ny = 2;
    nlss = LPVcore.nlss(f, h, nx, nu, ny, 0, true);
    
    %Create LPV model (not exactly same as the one produced in nlss2lpv)
    np = 10;
    p = cell(np,1);
    for i = 1:np
        p{i} = preal(sprintf(sprintf('p(%%0%id)', numel(num2str(np))), i));
    end
    A = [ 0,         0,        1,    0;
        0,         0,        0,    1;
        c*d*p{3}, -b*e*p{4}, p{5}, b*p{6};
        -b*d*p{7},  a*e*p{8}, p{9}, p{10}];

    B = [ 0,         0;
        0,         0;
        c*n*p{1}, -b*n*p{2};
        -b*n*p{2},  a*n*p{1}];

    C = [eye(2),zeros(2)];
    D = zeros(2);
    eta = @(x) schedMap(x,param);
    lpvsys = LPVcore.lpvss(A,B,C,D);
end

function p = schedMap(x, par)
    a = par(1);
    b = par(2);
    c = par(3);
    f = par(6);

    q1 = x(1,:);
    q2 = x(2,:);
    q1d = x(3,:);
    q2d = x(4,:);

    cosd = cos(q1-q2);
    sind = sin(q1-q2);

    h = 1 ./ (a * c - b^2 * cosd.^2);
    
    p = [h;
         cosd .* h;
         sincFun(q1) .* h;
         cosd .* sincFun(q2) .* h;
         (-b^2 * cosd .* sind .* q1d - (c + b * cosd) * f) .* h;
         (-c * sind .* q2d + cosd * f) .* h;
         cosd .* sincFun(q1) .* h;
         sincFun(q2) .* h;
         (a * b * sind .* q1d + f * (a + b * cosd)) .* h;
         (b^2 * sind .* cosd .* q2d - a * f) .* h];
end
% sinc function
function  y = sincFun(x)
    y = sin(x)./x;
    y(x == 0) = 1;
end
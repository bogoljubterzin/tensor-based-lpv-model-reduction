function [nlss, lpvsys, eta] = gen_robot
    load("2-DOF_robotic_arm_parameters.mat","params");
    param = [params.a params.b params.c params.d params.e params.f params.n];
    [f_, ~, ~, state, input] = generalized_2DOF_robotarm;
    f = @(x,u) f_(x, u, param);
    h = @(x,u) [x(1,:);x(2,:)]; %Just positions for example
    nx = length(state);
    nu = length(input);
    ny = 2;
    nlss = LPVcore.nlss(f, h, nx, nu, ny, 0, true);
    [lpvsys, eta] = LPVcore.nlss2lpvss(nlss, 'numerical', 'factor');
end
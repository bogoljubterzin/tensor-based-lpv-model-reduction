function [nlss, lpvsys, eta] = gen_gyro
    %Function to generate the gyro used for the paper
    actuators = [1;1;0;1]; locked = [0;0;1;0];
    load("GyroParameters.mat","params");
    param = [params.Ka params.Ib params.Jb ...
        params.Kb params.Ic params.Jc ...
        params.Kc params.Id params.Jd...
        params.Km1 params.Km2 params.Km3 params.Km4 ...
        params.fv1 params.fv2 params.fv3 params.fv4]';
    [fGyro, ~, ~, state, input] = generalized_gyro(actuators, locked);
    f = @(x,u) fGyro(x, u, param);
    h = @(x,u) [x(4,:);x(5,:);x(6,:)];
    nx = length(state);
    nu = length(input);
    ny = 3;
    nlss = LPVcore.nlss(f, h, nx, nu, ny, 0, true);
    [lpvsys, eta] = LPVcore.nlss2lpvss(nlss, 'numerical', 'element');
end
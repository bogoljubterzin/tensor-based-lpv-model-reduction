%% Symbolically construct the nonlinear state-space equations of the gyro
% Provides either the classical nonlinear form
%           \dot{x} = f(x,u)
% Or a the descriptor nonlinear form
%    e(x,u) \dot{x} = f_decoupled(x,u)
%
% Example: a giro with actuators in the 1st and 2nd component, and the 3rd
% gimbal locked
%    actuators = [1; 1; 0; 0];
%    locked = [0; 0; 1; 0];
%    [f, f_decoupled, e] = generalized_gyro(actuators, locked);


function [f, f_decoupled, e, x, u] = generalized_gyro(actuators, locked)
    arguments
        actuators (4,1) double {mustBeNumeric,mustBeFinite}
        locked (4,1) double {mustBeNumeric,mustBeFinite}
    end
        if any(actuators & locked)
            eid = 'Settings:notCompatible';
            msg = 'Cannot actuate a locked component. Please, change the position of the actuators or locks.';
            error(eid,msg)
        end

    %% Variables
    syms q1 q2 q3 q4 qd1 qd2 qd3 qd4 real
    syms Ka Ib Jb Kb Ic Jc Kc Id Jd real
    syms i1 i2 i3 i4 real
    syms Km1 Km2 Km3 Km4 real
    syms fv1 fv2 fv3 fv4 real
    
    actuators = boolean(actuators);
    locked = boolean(locked);

    q = [q1;q2;q3;q4];
    qd = [qd1;qd2;qd3;qd4];
    i = [i1;i2;i3;i4];

    s1 = sin(q1); c1 = cos(q1);
    s2 = sin(q2); c2 = cos(q2);
    s3 = sin(q3); c3 = cos(q3);
    s4 = sin(q4); c4 = cos(q4);

    a1 = Jc - Kc;
    a2 = Jd - Id;
    a3 = Id - Jc - Jd + Kc;
    a4 = Ic + Id;
    a5 = Ib + Ic - Kb - Kc;

    param = [Ka Ib Jb Kb Ic Jc Kc Id Jd Km1 Km2 Km3 Km4 fv1 fv2 fv3 fv4]';

    %% M - Mass matrix
    Ma = blkdiag(0,0,0,Ka);
    Mb = blkdiag(0,0,Jb,Ib*s3^2+Kb*c3^2);
    Mc = [0,0,0,0;
        0,Ic,0,-Ic*s3;
        0,0,Jc*c2^2+Kc*s2^2,a1*s2*c2*c3;
        0,-Ic*s3,a1*s2*c2*c3,Ic*s3^2+(Jc*s2^2+Kc*c2^2)*c3^2];
    Md = [Jd,0,Jd*c2,Jd*s2*c3;
        0,Id,0,-Id*s3;
        Jd*c2,0,Id*s2^2+Jd*c2^2,a2*s2*c2*c3;
        Jd*s2*c3,-Id*s3,a2*s2*c2*c3,Id*s3^2+(Id*c2^2+Jd*s2^2)*c3^2];

    M = Ma+Mb+Mc+Md;
    M = simplify(M,100);

    %% C - Coriolis matrix
    Qd = blkdiag(qd',qd',qd',qd');

    Gamma1 = .5*[0, 0,        0,         0;
        0, 0,        -Jd*s2,     Jd*c2*c3;
        0, -Jd*s2,   0,         -Jd*s2*s3;
        0, Jd*c2*c3, -Jd*s2*s3,     0];
    % simplify(Gamma1-Gamma1'==0)
    Gamma2 = .5*[0,         0, Jd*s2,                      -Jd*c2*c3;
        0,         0, 0,                          0;
        Jd*s2,     0, -2*a3*s2*c2,                a3*(c2^2*c3-s2^2*c3)-a4*c3;
        -Jd*c2*c3, 0, a3*(c2^2*c3-s2^2*c3)-a4*c3, 2*a3*c2*c3^2*s2];
    % simplify(Gamma2-Gamma2'==0)
    Gamma3 = .5*[0,        -Jd*s2,                     0,          Jd*s2*s3;
        -Jd*s2,   0,                          2*a3*s2*c2, a4*c3+a3*(c3*s2^2-c2^2*c3);
        0,        2*a3*s2*c2,                 0,          0;
        Jd*s2*s3, a4*c3+a3*(c3*s2^2-c2^2*c3), 0,          -2*(a5+a3*s2^2)*c3*s3];
    % simplify(Gamma3-Gamma3'==0)
    Gamma4 = .5*[0,         Jd*c2*c3,                   -Jd*s2*s3,                  0;
        Jd*c2*c3,  0,                          a3*(c3*s2^2-c2^2*c3)-a4*c3, -2*a3*c2*c3^2*s2;
        -Jd*s2*s3, a3*(c3*s2^2-c2^2*c3)-a4*c3, 2*a3*c2*s2*s3,              2*(a5+a3*s2^2)*c3*s3;
        0,         -2*a3*c2*c3^2*s2,           2*(a5+a3*s2^2)*c3*s3,       0];
    % simplify(Gamma4-Gamma4'==0)

    C = Qd*[Gamma1;Gamma2;Gamma3;Gamma4];

    C = simplify(C,100);

    %% Damping
    F = blkdiag(fv1,fv2,fv3,fv4);

    %% Motor constant
    Km = blkdiag(Km1,Km2,Km3,Km4);

    %% Remove rows/columns corresponding to q3 and qd3
    % and set q3, qd3 to zero

    qd = [qd1;qd2;qd3;qd4];

    qind = [1 2 3 4]; qind(locked) = [];
    iind = [1 2 3 4]; iind(~actuators) = [];

    M = M(qind, qind);
    M = simplify(subs(M, q(locked), zeros(length(q(locked)),1)),100);

    C = C(qind,qind);
    C = subs(C, [q(locked); qd(locked)], [zeros(length(q(locked)),1);zeros(length(qd(locked)),1)]);

    qd = qd(qind);
    i = i(iind);

    F = F(qind, qind);

    Km = Km(qind, iind);

    %% Total

    Mi = inv(M);
    Mi = simplify(Mi,100);

    qdd = Mi*(Km*i-C*qd-F*qd);
    if all(locked == 0)
        warning('No locked components, so long EoM will take long time (~100 sec) to generate :/')
        qdd  = simplify(qdd, 1);
    elseif all(locked == [1; 0; 0; 0])
        qdd  = simplify(qdd, 20);
    else
        qdd  = simplify(qdd, 100);
    end
    
    qdd_decoupled = (Km*i-C*qd-F*qd);
    qdd_decoupled = simplify(qdd_decoupled,100);

    f_ = [qd;qdd];
    f_decoupled = [qd;qdd_decoupled];
    e_ = blkdiag(eye(length(qd)), M); %the first states are integrating velocity to position. Identity in this case

    x = [q(qind);qd];
    u = i;

    f = matlabFunction(f_, 'vars', {x, u, param});
    f_decoupled = matlabFunction(f_decoupled, 'vars', {x, u, param});
    e = matlabFunction(e_, 'vars', {x, u, param});
end
function [f1, f_decoupled, e, x, u] = generalized_2DOF_robotarm
    %% Variables
    syms q1 q2 qd1 qd2 real
    syms a b c1 d1 e f1 n real
    syms i1 i2 real
    
    q = [q1; q2];
    qd = [qd1; qd2];
    i = [i1; i2];
    
    param = [a b c1 d1 e f1 n];

    %% M - Mass matrix
    M = [a,                b * cos(q1 - q2);
        b * cos(q1 - q2), c1];
    %% C - Coriolis matrix
    C = [b * sin(q1 - q2) * qd2^2 + f1 * qd1;
        -b * sin(q1 - q2) * qd1^2 + f1 * (qd2 - qd1)];
    %%
    g = [-d1 * sin(q1);
        -e * sin(q2)];

    %% Total
    Mi = inv(M);
    Mi = simplify(Mi);

    qdd = Mi * (n * i - C - g);
    qdd_decoupled = (n * i - C - g);

    f_ = [qd; qdd];
    f_decoupled_ = [qd; qdd_decoupled];
    e_ = blkdiag(eye(length(qd)), M);
    
    x = [q; qd];
    u = i;

    f1 = matlabFunction(f_, 'vars',{x, u, param});
    f_decoupled = matlabFunction(f_decoupled_, 'vars',{x, u, param});
    e = matlabFunction(e_, 'vars',{x, u, param});
end
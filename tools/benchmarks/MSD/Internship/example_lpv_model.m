%% example_lpv_model.m
% Function that creates LPVSS model of MSD example 
% provided the number of masses M:
% Equations of motion, F being forces
% Nonlinear springs: k(q) = k1*q + k2*q^3
% mddot(q1) = -F1 - F1,2
% mddot(qi) = -Fi - Fi,i-1 - Fi,i+1 ; i = 2,..M-1
% mddot(qM) = -FM - FM,M-1 + u

function system_lpv = example_lpv_model(M)

    % M = 50; % Number of masses, subsystems
    d = 1; % Linear damping coefficient
    k1 = 0.5; % Linear spring coefficient
    k2 = 1; % Cubic spring coefficient
    m = 1; % Mass

    sys = {};
    A_sc = {}; % contains scheduling variables ai = xi^2
    B_sc = {}; % contains scheduling variables bi = (xi-xi-1)^2
    for i = 1:M
        A_sc{i} = preal(strcat('a',num2str(i)));
        B_sc{i} = preal(strcat('b',num2str(i)));
    end
    for i = 2:M-1
    
        % For every subsystem inputs are ui = [qi-1 qdoti-1 qi+1 qdoti+1]
    
        % defining 3 scheduling variables for each sub-system
        p1 = A_sc{i};
        p2 = B_sc{i};
        p3 = B_sc{i+1};
        % Parameter-varying SS matrices

        A = [0 1; -1/m*(3*k1 +k2*(p1+p2+p3)) -1/m*(3*d)];
        B = [0 0 0 0;
            1/m*(k1+k2*p2), 1/m*d, 1/m*(k1+k2*p3), 1/m*d];
        C = eye(2);
        D = zeros(2,4);
    
        sys{i} = LPVcore.lpvss(A,B,C,D);
    
        % Properly defining input and output names so the interconnected system
        % can be created
        sys{i}.InputName = [{strcat('q',num2str(i-1))}, {strcat('qdot',num2str(i-1))},...
            {strcat('q',num2str(i+1))}, {strcat('qdot',num2str(i+1))}];
        sys{i}.OutputName = [{strcat('q',num2str(i))},{strcat('qdot',num2str(i))}];
            
    end
    
    % Edge cases: First subsystem
    % input u1 = [q2, qdot2]
    p1 = A_sc{1};
    p2 = B_sc{2};
    
    A = [0 1; -1/m*(2*k1+k2*(p1+p2)) -1/m*(2*d)];
    B = [0 0; 1/m*(k1+k2*p2) d/m];
    C = eye(2);
    D = zeros(2,2);
    sys{1} = LPVcore.lpvss(A,B,C,D);
    sys{1}.InputName = {'q2','qdot2'};
    sys{1}.OutputName = {'q1','qdot1'};
    
    % Edge cases: Last subsystem
    % input uM = [qM-1, qdotM-1, u]
    p1 = A_sc{M};
    p2 = B_sc{M};

    A = [0 1; -1/m*(2*k1 + k2*(p1+p2)) -1/m*(2*d)];
    B = [0 0 0; 1/m*(k1+k2*p2) d/m 1/m];
    C = eye(2); D = zeros(2,3);
    sys{M} = LPVcore.lpvss(A,B,C,D);
    sys{M}.InputName = {strcat('q',num2str(M-1)),strcat('qdot',num2str(M-1)),'u'};
    sys{M}.OutputName = {strcat('q',num2str(M)),strcat('qdot',num2str(M))};

    sysM =  connect(sys{:},'u',sprintf('q%i',M));

    system_lpv = sysM;

end



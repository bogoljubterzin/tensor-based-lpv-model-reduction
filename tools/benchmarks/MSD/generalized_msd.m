%% Symbolically construct the nonlinear state-space equations of generalized MSD
% Provides the classical nonlinear form
%           \dot{x} = f(x,u)
%                 y = h(x,u)
%
% 
% Input arguments: 
% M - Number of blocks
% nl_springs - array of indices of nonlinear springs (between masses) in
% {1,...,M-1}
% nl_springs_wall - array of indices of nonlinear springs (between mass and
% wall) in {1,...,M}
% nl_dampers - array of indices of nonlinear dampers (between masses) in
% {1,...,M-1}
% nl_dampers_wall - array of indices of nonlinear springs (between mass and
% wall) in {1,...,M}
% k - Spring nonlinearity function handle (e.g. k = @(x) 0.5*x + x^3)
% d - Damper nonlinearity function handle (e.g. d = @(x,xdot) sin(xdot) + x^2)

% parameters: m, k0, d0
% m - Mass of each block
% k0 - Linear spring coefficient
% d0 - Linear damper coefficient

% example: nl_springs = [3] denotes nonlinear spring between m3 and m4
%          nl_damper_wall = [5] denotes nonlinear damper between m5 and
%          wall
function [f,x,u] = generalized_msd(M,nonlinearities,k,d)
    % Validate input arguments
    arguments
        M (1,1) {mustBeNumeric,mustBeFinite,mustBePositive}
        nonlinearities
        k function_handle
        d function_handle
    end
    nl_springs = nonlinearities.nl_springs;
    nl_springs_wall = nonlinearities.nl_springs_wall;
    nl_dampers = nonlinearities.nl_dampers;
    nl_dampers_wall = nonlinearities.nl_dampers_wall;
    if nargin < 4
        error("Not enough input arguments.");
    end
    if ~isempty(nl_springs) && (max(nl_springs) >= M || min(nl_springs) <= 0)
        error("Invalid choice of nonlinear springs indices (between blocks). Values must be between 1 and M-1.");

    end
    if ~isempty(nl_springs_wall) && (max(nl_springs_wall) > M || min(nl_springs_wall) <= 0)
        error("Invalid choice of nonlinear spring indices (connected to wall). Values must be between 1 and M.")
    end

    if ~isempty(nl_dampers) && (max(nl_dampers) >= M || min(nl_dampers) <= 0)
        error("Invalid choice of nonlinear damper indices (between blocks). Values must be between 1 and M-1.")
    end
    if ~isempty(nl_dampers_wall) && (max(nl_dampers_wall) > M || min(nl_dampers_wall) <= 0)
        error("Invalid choice of nonlinear damper indices (connected to wall). Values must be between 1 and M.")
    end    
    %% Variables
    syms k0 d0 m real    
    syms u real

    q = sym('q',[M,1],'real');
    qd = sym('qd',[M,1],'real');
    qdd = sym(NaN(M,1));
    param = [m k0 d0];

    dlin = @(x,xd) d0*xd;
    klin = @(x) k0*x;

    for i = 1:M
        
        
        kwall = klin;
        kright = klin;
        kleft = klin;
        
        dwall = dlin;
        dright = dlin;
        dleft = dlin;

        if ismember(i,nl_springs_wall)
            kwall = k;
        end
        if ismember(i,nl_dampers_wall)
            dwall = d;
        end
        if ismember(i,nl_springs)
            kright = k;
        end
        if ismember(i,nl_dampers)
            dright = d;
        end
        if ismember(i-1,nl_springs)
            kleft = k;
        end
        if ismember(i-1,nl_dampers)
            dleft = d;
        end

        qdd(i,1) = -1/m*(kwall(q(i)) + dwall(q(i),qd(i)));
        if i ~= 1
            qdd(i,1) = qdd(i,1) - 1/m*(kleft(q(i)-q(i-1)) + dleft(q(i)-q(i-1),qd(i)-qd(i-1)));
        end
        if i ~= M
            qdd(i,1) = qdd(i,1) - 1/m*(kright(q(i)-q(i+1)) + dright(q(i)-q(i+1),qd(i)-qd(i+1)));
        end


        % Fwall = - kwall(q(i)) - dwall(q(i),qd(i));
        % if i == 1 && M == 1
            % Ftot = u;
        % elseif i == 1
       
            % Ftot = kright(q(i+1)-q(i)) + dright(q(i+1)-q(i),qd(i+1)-qd(i));
        % elseif i == M
            % Ftot = u - kleft(q(i)-q(i-1)) - dleft(q(i)-q(i-1),qd(i)-qd(i-1));
        % if i ~= M
            % qdd(i,1) = qdd(i,1) - 1/m*(kright(q(i)-q(i+1)) + dright(q(i)-q(i+1),qd(i)-qd(i+1)));
        % else
        %     Ftot = - kleft(q(i)-q(i-1)) - dleft(q(i)-q(i-1),qd(i)-qd(i-1)) +...
        %         kright(q(i+1)-q(i)) + dright(q(i+1)-q(i),qd(i+1)-qd(i));
        % end
            % Ftot = Ftot + Fwall;
            % qdd(i) = Ftot/m;
    end
    qdd(M) = qdd(M) + 1/m*u;
    x = [q; qd];
    f_ = [qd; qdd];
    
    f = matlabFunction(f_,'Vars',{x,u,param});

end



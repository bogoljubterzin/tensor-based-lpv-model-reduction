function [tout, x, y] = nonAffineLpvSim(obj, map, ref, t, x0)
    Np = obj.Np; Nx = obj.Nx; Ny = obj.Ny; Nu = obj.Nu;
    % A_(abs(A_)<=1e-10) = 0;  B_(abs(B_)<=1e-10) = 0; %remove numeric noise
    % C_(abs(C_)<=1e-10) = 0;  D_(abs(D_)<=1e-10) = 0; %remove numeric noise
    if obj.isct == 1
        % Simulate
        ufun = @(tau) interp1(t, ref, tau, 'previous')';
        xdot = @(t, x) xfun(x,ufun(t),obj,map,Np);
        [tout,x] = ode45(xdot,t,x0); %x \in Nt x Nx!
        x = x'; %x \in Nx x Nt!
        %now get output
        y = hfun(x, ufun(tout), obj,map,Np,Ny);
    else
        N = length(t);
        x = [x0, zeros(Nx,N)];
        y = zeros(Ny,N);
        for i = 1:N
            if mod(i,20) == 0
                disp(i);
            end
            p = map(x(:,i),ref(:,i)).';
            Ai = feval(obj.A,p);
            Bi = feval(obj.B,p);
            Ci = feval(obj.C,p);
            Di = feval(obj.D,p);
            x(:,i+1) = [Ai, Bi]*[x(:,i);ref(:,i)];
            y(:,i) = [Ci, Di]*[x(:,i);ref(:,i)];
        end
        tout = t;
    end
end

%%%%%%%%%%%%%% LOCAL FXs
function xdot = xfun(x,u,A,A_,B,B_,map,Np)
    p = map(x,u).';
    Ai = feval(obj.A,p);
    Bi = feval(obj.B,p);
    xdot = [Ai, Bi]*[x;u];
end

function y = hfun(x, u, C,C_,D,D_,map,Np,Ny)
    y = zeros(Ny, size(x,2));
    for i = 1:size(x,2) % Need for-loop here such that 'sim' works properly
        %NOTE: x \in Nx x Nt, u \in Nu x Nt
        p = map(x(:,i),u(:,i)); 
        Ci = feval(obj.C,p);
        Di = feval(obj.D,p);
        y(:,i) = [Ci, Di] * [x(:,i);u(:,i)];
    end
end
function [tout, x, y, pout] = affineLpvSim(obj, map, ref, t, x0)
    Np = obj.Np; Nx = obj.Nx; Ny = obj.Ny; Nu = obj.Nu;
    pout = [];
    %Gather affine matrices
    A_ = obj.A.matrices; A = zeros(Nx,Nx,Np);
    B_ = obj.B.matrices; B = zeros(Nx,Nu,Np);
    C_ = obj.C.matrices; C = zeros(Ny,Nx,Np);
    D_ = obj.D.matrices; D = zeros(Ny,Nu,Np);
    % A_(abs(A_)<=1e-10) = 0;  B_(abs(B_)<=1e-10) = 0; %remove numeric noise
    % C_(abs(C_)<=1e-10) = 0;  D_(abs(D_)<=1e-10) = 0; %remove numeric noise
    A(:,:,1) = A_(:,:,1); B(:,:,1) = B_(:,:,1); %non-dependent term
    C(:,:,1) = C_(:,:,1); D(:,:,1) = D_(:,:,1); %non-dependent term
    if obj.isct == 1
        % Simulate
        ufun = @(tau) interp1(t, ref, tau, 'previous')';
        xdot = @(t, x) xfun(x,ufun(t),A,A_,B,B_,map,Np);
        [tout,x] = ode45(xdot,t,x0); %x \in Nt x Nx!
        x = x'; %x \in Nx x Nt!
        %now get output
        y = hfun(x, ufun(tout), C,C_,D,D_,map,Np,Ny);
    else
        N = length(t);
        x = [x0, zeros(Nx,N)];
        y = zeros(Ny,N);
        pout = zeros(Np,N);
        for i = 1:N
            p = map(x(:,i),ref(:,i));            
            for j = 1:Np
                A(:,:,j+1) = A_(:,:,j+1)*p(j);
                B(:,:,j+1) = B_(:,:,j+1)*p(j);
                C(:,:,j+1) = C_(:,:,j+1)*p(j);
                D(:,:,j+1) = D_(:,:,j+1)*p(j);

            end
            x(:,i+1) = [sum(A,3), sum(B,3)]*[x(:,i);ref(:,i)];
            y(:,i) = [sum(C,3), sum(D,3)]*[x(:,i);ref(:,i)];
            pout(:,i) = p;
        end
        tout = t;
    end
end

%%%%%%%%%%%%%% LOCAL FXs
function xdot = xfun(x,u,A,A_,B,B_,map,Np)
    p = map(x,u);
    for i=1:Np
        A(:,:,i+1) = A_(:,:,i+1)*p(i);
        B(:,:,i+1) = B_(:,:,i+1)*p(i);
    end
    xdot = [sum(A,3), sum(B,3)]*[x;u];
end

function y = hfun(x, u, C,C_,D,D_,map,Np,Ny)
    y = zeros(Ny, size(x,2));
    for i = 1:size(x,2) % Need for-loop here such that 'sim' works properly
        %NOTE: x \in Nx x Nt, u \in Nu x Nt
        p = map(x(:,i),u(:,i)); 
        for j=1:Np
            C(:,:,j+1) = C_(:,:,j+1)*p(j);
            D(:,:,j+1) = D_(:,:,j+1)*p(j);
        end
        y(:,i) = [sum(C,3), sum(D,3)] * [x(:,i);u(:,i)];
    end
end
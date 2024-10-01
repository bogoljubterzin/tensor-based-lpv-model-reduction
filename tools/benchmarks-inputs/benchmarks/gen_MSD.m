function [nlss, lpvsys, eta] = gen_MSD(M,np)
    %Reproducing the MSD system from Bogoljub's Internship
    nl_index.springs = [M-fix(np/2):M-1];
    nl_index.springs_wall = [M-fix(np/2):M];
    nl_index.dampers = [];
    nl_index.dampers_wall = [];
    %Define nonlinearities with the nl function
    k = @(x) x^3; % Spring connecting mass blocks
    d = @(x,xdot) 0; % Damper connecting mass blocks
    k_wall = @(x) x^3; % Spring connecting each mass with the wall
    d_wall = @(x,xdot) 0; % Damper connecting each mass with the wall
    param = [1 0.5 1 0.5 1]; % [M k d kwall dwall]
    [fMSD,state,input] = generalized_msd_Javi(M,nl_index,k,d,k_wall,d_wall);
    f = @(x,u) fMSD(x,u,param); h = @(x,u) x(M,:);
    nx = length(state); nu = length(input); ny = 1;
    nlss = LPVcore.nlss(f,h,nx,nu,ny,0,true);
    [lpvsys, eta] = LPVcore.nlss2lpvss(nlss,"analytical","factor");
end
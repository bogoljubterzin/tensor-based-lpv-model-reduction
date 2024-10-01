% Developing tensor-based POD methods
l2err = @(y,ybar) norm(y-ybar)/norm(y-mean(y-ybar))*100;

load('benchmark_models.mat')
lpvss_dt
tensSS_dt

msd = 1; % MSD1: Nx = 10; Np = 9;

lpvss = lpvss_dt{msd};
tensSS = tensSS_dt{msd};
eta_map = eta{msd}.map;

[~,x,~,p] = tensSS.simulateSS(utrain(t),t,x0{msd});
xdot = x(:,2:end)-x(:,1:end-1);
%%

% decomp = 'tsvd';
decomp = 'hosvd';

data = zeros(tensSS.Nx, tensSS.Nx, tensSS.Np, length(t)-1);
for i = 1:length(t)-1
    t1 = ttt(tensor(x(:,i)),tensor(xdot(:,i)));
    t1 = squeeze(t1);
    t2 = ttt(t1,tensor(p(:,i)));
    t2 = squeeze(t2);
    data(:,:,:,i) = t2;
end

time_horizon = 1:5:200;
data = data(:,:,:,time_horizon);
R = min(modrank(data));
% rx = 6; rp = 6;
%%
clear decomp_info
switch decomp
    case "tsvd"
        tic
        [x,gains,iterations,Xa]=TSVDND_SW(data,R);
        tcomp = toc;
        V = x{1}(:,1:rx);
        W = x{2}(:,1:rx);
        Z = x{3}(:,1:rp);

        decomp_info.gains = gains;
        decomp_info.iterations = iterations;
        decomp_info.low_rank_app = Xa;
    case "hosvd"
        tic
        T = hosvd(tensor(data),2*sqrt(eps));
        tcomp = toc;
        % V = T{1}(:,1:rx);
        % W = T{2}(:,1:rx);
        % Z = T{3}(:,1:rp);
        V = T{1};
        W = T{2};
        Z = T{3};
        % info

        decomp_info.S = T.core;
        decomp_info.T = T;
        
end
%%

normxsqr = collapse(data.^2);
diffnormsqr = collapse((data-full(T)).^2);
relnorm = sqrt(diffnormsqr/normxsqr);


% tensSSrPOD4d_tsvd = tensSS.PetrovGalerkinLPV(W,V,Z,1);
tensSSrPOD4d_hosvd = tensSS.PetrovGalerkinLPV(W,V,Z,1);


%%
function [tensSSr, V, info] = POD4D(obj, train_data, decomp, rx, rp)
            decomp = lower(decomp);
            u_train = train_data.input;
            x0 = train_data.x0;
            tin = 0:obj.Ts:(length(u_train)-1)*obj.Ts;
            [~,x,~,p] = obj.simulateSS(u_train,tin,x0);
            xdot = x(:,2:end) - x(:,1:end-1);

            data = zeros(obj.Nx, obj.Nx, obj.Np, length(tin)-1);
            for i = 1:length(tin)-1
                t1 = ttt(tensor(x(:,i)),tensor(xdot(:,i)));
                t1 = squeeze(t1);
                t2 = ttt(t1,tensor(p(:,i)));
                t2 = squeeze(t2);
                data(:,:,:,i) = t2;
            end

            R = max(rx,rp);
            switch decomp
                case "tsvd"
                    tic
                    [x,gains,iterations,Xa]=TSVDND_SW(data,R);
                    tcomp = toc;
                    V = x{1}(:,1:rx);
                    W = x{2}(:,1:rx);
                    Z = x{3}(:,1:rp);

                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "succr1"
                    tic
                    [x1,x2,x3,gains,iterations,Xa]=succR1_SW(data,R);
                    tcomp = toc;
                    V = x1(:,1:rx);
                    W = V;
                    Z = x2(:,1:rp);

                    decomp_info.gains = gains;
                    decomp_info.iterations = iterations;
                    decomp_info.low_rank_app = Xa;
                case "hosvd"
                    tic
                    T = hosvd(tensor(data),2*sqrt(eps));
                    tcomp = toc;
                    V = T{1}(:,1:rx);
                    W = T{2}(:,1:rx);
                    Z = T{3}(:,1:rp);
                    
                    % info
        
                    decomp_info.S = T.core;
                    decomp_info.T = T;
                    
            end

            info.tcomp = tcomp;
            info.decomp = decomp_info;

            tensSSr = obj.PetrovGalerkinLPV(W,V,Z,1);
        end
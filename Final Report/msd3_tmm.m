lpvss = lpvss_msd3;
tensSS = msd3;
eta_map = eta_msd3;

x0 = zeros(tensSS.Nx,1);

[~,~,yTrain,p1] = tensSS.simulateSS(umsd_train,t,x0,1);
[~,~,yTest,p2] = tensSS.simulateSS(umsd_test,t,x0,1);
%%

N = 3;  % this cant be too high (for msd3) because of memory
        % issues (max N = 3); for other models it is ok to be higher

Wn_hosvd = cell(N,1);
Qn_hosvd = cell(N,1);

Wn_tsvd = cell(N,1);
Qn_tsvd = cell(N,1);
Tstates = cell(N,1);
Tsched = cell(N,1);
Tstatesq = cell(N,1);
Tschedq = cell(N,1);

for i = 1:N
    Tstates{i} = [];
    Tstates{i} = [];
    Tsched{i} = [];
    Tschedq{i} = [];
end
t_hosvd_r = cell(N,1);
t_hosvd_q = cell(N,1);
t_tsvd_r = cell(N,1);
t_tsvd_q = cell(N,1);

Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);

%%
decomp = 'hosvd';

tol = 1e-6;
for n = 1:N
    disp(n);
    if strcmp(decomp,"hosvd") == 1
        tic
        T = hosvd(tensor(Wn{n}),tol);
        t_hosvd_r{n} = toc;
        Wn_hosvd{n} = T;
        Tstates{n} = [Tstates{n}, T{1}];
        for j = 2:n+1
            Tsched{n} = [Tsched{n}, T{j}];
        end
        tic
        Tq = hosvd(tensor(Qn{n}),2*sqrt(eps));
        t_hosvd_q{n} = toc;
        Qn_hosvd{n} = Tq;
        Tstatesq{n} = [Tstatesq{n}, Tq{n+2}];
        for j = 2:n+1
            Tschedq{n} = [Tschedq{n}, Tq{j}];
        end
    elseif strcmp(decomp,"tsvd") == 1
        K = min(modrank(Wn{n}));
        tic
        [x,gains,iterations,Xa] = TSVDND_SW(Wn{n},K);
        t_tsvd_r{n} = toc;
        Wn_tsvd{n}.x = x;
        Wn_tsvd{n}.gains = gains;
        Wn_tsvd{n}.iterations = iterations;
        Wn_tsvd{n}.Xa = Xa;

        % for the state projection V, we gather the
        % vectors from first dimension x{1}
        Tstates{n} = [Tstates{n}, x{1}];
    
        % for the scheduling projection Z, we gather vectors from dimensions
        % x{2},..., x{n}
        for j = 2:n+1
            Tsched{n} = [Tsched{n}, x{j}];
            % Tsched{j-1} = null(x{j}.')
        end

        K = min(modrank(Qn{n}));
        tic
        [x,gains,iterations,Xa] = TSVDND_SW(Qn{n},K);
        t_tsvd_q{n} = toc;
        Qn_tsvd{n}.x = x;
        Qn_tsvd{n}.gains = gains;
        Qn_tsvd{n}.iterations = iterations;
        Qn_tsvd{n}.Xa = Xa;
        
        Tstatesq{n} = [Tstatesq{n}, x{n+2}];
        for j = 2:n+1
            Tschedq{n} = [Tschedq{n}, x{j}];
        end
    end
end
%%
n=1
norm(double(Wn{n})-double(Wn_tsvd{n}.Xa),'fro')
norm(double(Qn{n})-(Qn_tsvd{n}.Xa),'fro');

%%
Vr = cell(N,1);
Zr = cell(N,1);
Vq = cell(N,1);
Zq = cell(N,1);
Wr = cell(N,1);
Wq = cell(N,1);

Vr{1} = orth(Tstates{1});
Wr{1} = Vr{1};
Zr{1} = orth(Tsched{1});
Vq{1} = orth(Tstatesq{1});
Zq{1} = orth(Tschedq{1});
Wq{1} = Vq{1};
for i = 2:N
    Vr{i} = orth([Vr{i-1}, Tstates{i}]);
    Wr{i} = Vr{i};
    Zr{i} = orth([Zr{i-1}, Tsched{i}]);
    Vq{i} = orth([Vq{i-1}, Tstatesq{i}]);
    Wq{i} = Vq{i};
    Zq{i} = orth([Zq{i-1}, Tschedq{i}]);
end

% Now we have the matrices Tx and Tp for every horizon (for both
% reachability and observability case)
% Let's see if for some N we can use the optimal case
opt_case = zeros(N,1);
for i = 1:N
    Rn = Vr{i};
    Nn = Vq{i};
    rx = size(Rn,2);
    if rank(Nn'*Rn) == rx
        opt_case(i) = true;
        [N1,R] = qr(Nn'*Rn);
        Vh{i} = Rn/R;
        Wh{i} = Nn'\N1;
        Zh{i} = orth([Zr{i},Zq{i}]);
    end
end


%%
er_hosvd = zeros(N,1);
sysr_hosvd = cell(N,1);
yr_hosvd = cell(N,1);
%%
er_tsvd = zeros(N,1);
sysr_tsvd = cell(N,1);
yr_tsvd = cell(N,1);
%%
yr_hosvd_test = cell(N,1);
if strcmp(decomp,'hosvd')
    for i = 1:N
        sysr = tensSS.PetrovGalerkinLPV(Wr{i},Vr{i},Zr{i},0);
        [~,~,yri_hosvd,~] = sysr.simulateSS(umsd_train,t,Vr{i}'*x0,constantTerm);
        [~,~,yri_hosvd_test,~] = sysr.simulateSS(umsd_test,t,Vr{i}'*x0,constantTerm);

        sysr_hosvd{i} = sysr;
        yr_hosvd{i} = yri_hosvd;
        er_hosvd(i) = nrmse(yTrain,yri_hosvd);
        yr_hosvd_test{i} = yri_hosvd_test;
    end

else
    for i = 1:N
        sysr = tensSS.PetrovGalerkinLPV(Wr{i},Vr{i},Zr{i},0);
        [~,~,yri_tsvd,~] = sysr.simulateSS(umsd_train,t,Vr{i}'*x0,constantTerm);
        sysr_tsvd{i} = sysr;
        yr_tsvd{i} = yri_tsvd;
        er_tsvd(i) = nrmse(yTrain,yri_tsvd);
    end

end


%%

FigMSD3_output = figure(1245);
clf(FigMSD3_output)
clrs = lines(3);
tiledlayout(2,1,"TileSpacing",'tight','padding','compact');
nexttile;
plot(t,yTrain,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr_hosvd{2},'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr_hosvd{3},'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");
legend(sprintf("FOM, $n_x = %d, n_p = %d$",tensSS.Nx,tensSS.Np), ...
       sprintf("$N=1, r_x = %d, r_p = %d$", sysr_hosvd{2}.Nx,sysr_hosvd{2}.Np), ...
       sprintf("$N=2, r_x = %d, r_p = %d$", sysr_hosvd{3}.Nx,sysr_hosvd{3}.Np))

nexttile;
plot(t,yTest,'Color',clrs(1,:),'LineStyle','--','LineWidth',3);
hold on; grid on;
plot(t,yr_hosvd_test{2},'Color',clrs(2,:),'LineWidth',1.2);
plot(t,yr_hosvd_test{3},'Color',clrs(3,:),'LineWidth',1.2);
xlabel("Time (s)"); ylabel("$y (m)$");

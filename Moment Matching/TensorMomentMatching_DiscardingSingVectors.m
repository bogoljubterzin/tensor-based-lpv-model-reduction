% This file contains the implementation of tensor-based moment matching
% with additional analysis in terms of singular values and singular vector,
% here we do not necesserily take all singular vectors, but we can choose
% to discard some of them to obtain lower-order models
tensSS = tensRobot;

horizon = 0;
N = horizon+1;
decomp = 'hosvd';
mode = 'T';
tol = 1e-6;

Wn = ReachabilityTensors(tensSS.A,tensSS.B,N);
Qn = ObservabilityTensors(tensSS.A,tensSS.C,N);


Vr_hosvd = cell(N,1);
Zr_hosvd = cell(N,1);

Vq_hosvd = cell(N,1);
Zq_hosvd = cell(N,1);

gainsVr_hosvd = cell(N,1);
gainsZr_hosvd = cell(N,1);
gainsVq_hosvd = cell(N,1);
gainsZq_hosvd = cell(N,1);



sv_sched = cell(N,1);
sv_states = cell(N,1);

for n = 1:N
    Wn_hosvd = hosvd(tensor(Wn{n}),tol);
    gainsWn_hosvd = hosvd_sv(Wn_hosvd);

    gainsVr_hosvd{n} = gainsWn_hosvd{1}(1:size(Wn_hosvd{1},2));
    Vr_hosvd{n} = Wn_hosvd{1};
    gainsZr_hosvd{n} = [];
    for j = 2:n+1
        Zr_hosvd{n} = [Zr_hosvd{n}, Wn_hosvd{j}];
        gainsZr_hosvd{n} = [gainsZr_hosvd{n}; gainsWn_hosvd{j}(1:size(Wn_hosvd{j},2))];
    end

    Qn_hosvd = hosvd(tensor(Qn{n}),tol);
    gainsQn_hosvd = hosvd_sv(Qn_hosvd);

    gainsVq_hosvd{n} = gainsQn_hosvd{n+2}(1:size(Qn_hosvd{n+2},2));
    Vq_hosvd{n} = Qn_hosvd{n+2};

    gainsVq_hosvd{n} = [];

    for j = 2:n+1
        Zq_hosvd{n} = [Zq_hosvd{n}, Qn_hosvd{j}];
        gainsZq_hosvd{n} = [gainsZq_hosvd{n}; gainsQn_hosvd{j}(1:size(Qn_hosvd{j},2))];
    end
end
%%
Vr = [];
Zr = [];

svVr = [];
svZr = [];
for n = 1:N
    Vr = [Vr, Vr_hosvd{n}];
    svVr = [svVr, gainsVr_hosvd{n}.'];

    Zr = [Zr, Zr_hosvd{n}];
    svZr = [svZr, gainsZr_hosvd{n}.'];
end

Vr1 = orth(Vr);
Zr1 = orth(Zr);

% now i want out of Vr_hosvd and gainsVr_hosvd to obtain Vr
% Zr_hosvd and gainsZr_hosvd -> Zr, etc.





% 
%     Qn_hosvd = hosvd(tensor(Qn{n}),tol);
%     Vq_hosvd{n} = Qn_hosvd{n+2};
%     for j = 2:n+1
%         Zq_hosvd{n} = [Zq_hosvd{n}, Qn_hosvd{j}];
%     end
%     gainsQn_hosvd{n} = hosvd_sv(Qn_hosvd);
% end

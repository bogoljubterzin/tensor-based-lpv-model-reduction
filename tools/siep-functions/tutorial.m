clear
close all

addpath('tensor_toolbox-master/')

%% Create Tensor
N1 = 10;
N2 = 10;
N3 = 1e3;

% first rank 1 tensor components
w1 = rand(1,N1);
v1 = rand(1,N2);
x = randn(N2,N3); % N2 inputs to multivariate polynomial function, N3 operating points
z1 = v1*x; % argument of h -> h(z)
h1 = z1.^3; % smooth h when sorted on z
% h1 = rand(1,N3); % random h

% second rank 1 tensor components
w2 = rand(1,N1);
v2 = rand(1,N2);
z2 = v2*x; % argument of h -> h(z), using same inputs 'x' as in the first rank 1 tensor
h2 = z2.^2; % smooth h when sorted on z
% h1 = rand(1,N3); % random h

figure
subplot(121)
plot(z1,h1,'.')
subplot(122)
plot(z2,h2,'.')

% outer product -> rank 1 tensor
T1 = zeros(N1,N2,N3);
T2 = zeros(N1,N2,N3);

for i=1:N1
    for j=1:N2
        for k=1:N3
            T1(i,j,k) = w1(i)*v1(j)*h1(k);
            T2(i,j,k) = w2(i)*v2(j)*h2(k);
        end
    end
end

% rank 2 tensor
T = T1+T2;

% reshuffle T and x %
s = RandStream('mlfg6331_64'); 
Index = datasample(s,1:N3,N3,'Replace',false); % draw random sequence

T_shuffle = zeros(size(T));
for shuffle = 1: N3
    T_shuffle(:,:,shuffle) = T(:,:,Index(shuffle)); % reshuffle the tensor in the third dimension given the random sequence
end
x_shuffle = x(:,Index); % reshuffle the inputs accoringly

%% Compute decomposition

% define filter operators
F1 = eye(N1);
F2 = eye(N2);

% second order central difference
F3 = toeplitz([-2 1 zeros(1,N3-2)],[-2 1 zeros(1,N3-2)]);
% F3 = F3*(-1); % to ensure positive definite
%F3 = eye(N3);


maxiter = 100;
r = 2;

[W,V,H,gains,iterations,Approx,W1,W2,W3,store_k_loops] = succR1_SW_dynamic_W(T,r,F1,F2,F3,x,maxiter);
% [W,V,H,gains,iterations,Approx] = succR1_SW_JD(T,r);
H = H.*repmat(gains',N3,1); % scale one of the singular vectors with the singular value

Z = V.'*x;

figure
subplot(121)
plot(Z(1,:),H(:,1),'r.')
subplot(122)
plot(Z(2,:),H(:,2),'r.')

T_approx = zeros(N1,N2,N3);

for i=1:N1
    for j=1:N2
        for k=1:N3
            T_approx_1(i,j,k) = W(i,1)*V(j,1)*H(k,1);
            T_approx_2(i,j,k) = W(i,2)*V(j,2)*H(k,2);
        end
    end
end
T_approx = T_approx_1+T_approx_2;

rel_err = frob(T_approx-T)/frob(T)
figure, plot(db(store_k_loops))
% shuffled Tensor %
[W_sh,V_sh,H_sh,gains,iterations,Approx] = succR1_SW_dynamic_W(T_shuffle,r,F1,F2,F3,x_shuffle,maxiter);
H_sh = H_sh.*repmat(gains.',N3,1);

Z_sh = V_sh.'*x_shuffle;

figure
subplot(121)
plot(Z_sh(1,:),H_sh(:,1),'r.')
subplot(122)
plot(Z_sh(2,:),H_sh(:,2),'r.')







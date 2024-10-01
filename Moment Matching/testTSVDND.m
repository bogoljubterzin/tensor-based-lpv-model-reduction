%% Compare it with 3D tensor decomposition algorithm

T = (rand(4,6,5));
[x,gains,iterations,Xa] = TSVDND_SW(T,num_levels)
[x1,x2,x3,gains1,iterations1,Xa1] = TSVD3D_SW(T,num_levels)

%%
fprintf("Error norm ||Xa - Xa1||: %d \n", norm(Xa-Xa1,'fro'))
disp("Gains 3D decomposition: "); disp(gains);
disp("Gains ND decomposition: "); disp(gains1);

%% Check the decomposition of 4D tensor
T = rand(5,5,4,4);
num_levels = 3;
[x,gains,iterations,Xa] = TSVDND_SW(T,num_levels)

%% Compute singular value core and see if the gains match

% Calculate core to check if the diagonal (singular values) are equal 
% to the gains (also singular values) from algorithm

R1 = size(x{1},2);
R2 = size(x{2},2);
R3 = size(x{3},2);
R4 = size(x{4},2);
S = zeros(R1,R2,R3,R4);
for i = 1:size(x{1},2)
    for j = 1:size(x{2},2)
        for k = 1:size(x{3},2)
            for l = 1:size(x{4},2)
                t1 = tensor(x{1}(:,i));
                t2 = tensor(x{2}(:,j));
                t3 = tensor(x{3}(:,k));
                t4 = tensor(x{4}(:,l));
    
                Trank1 = squeeze(ttt(t1,t2));
                Trank1 = squeeze(ttt(Trank1,t3));
                Trank1 = squeeze(ttt(Trank1,t4));
                S(i,j,k,l) = innerprod(tensor(T),Trank1);
            end
        end
    end
end


fprintf("sigma_1 = %f \n", S(1,1,1,1))
fprintf("sigma_2 = %f \n", S(2,2,2,2))
fprintf("sigma_3 = %f \n", S(3,3,3,3))


fprintf("Gains: %f \n", gains)

% This confirms that the decompositions are the same for 3D and ND (when
% N=3) algorithms. The error I was getting is tied to the fact that one
% dimension of tensor is nu = 1

%%
W = zeros(2,2,2);
W(1,1,1) = 4;
W(1,2,2) = 2;
W(2,2,2) = 1;
W = tensor(W);
[x,gains,iterations,Xa] = TSVDND_SW(W,2)


%%

e0 = zeros(5,1);
e1 = e0; e1(1) = 1;
e2 = e0; e2(2) = 1;
e3 = e0; e3(3) = 1;
e4 = e0; e4(4) = 1;
e5 = e0; e5(5) = 1;

T = zeros(I1,I2,I3);
I1 = 4; I2 = 5; I3 = 2;
v = rand(I1,3);
w = rand(I2,3);

M = v*w';
T(:,:,1) = M;
T(:,:,2) = -M;

R = modrank(T);

[x,gains,iterations,Xa] = TSVDND_SW(T,min(R));

M1 = Xa(:,:,1);
M2 = Xa(:,:,2);
rank(M - M1)
rank(M - M2)

E = T(:,:,1);
[x,gains,iterations,Xa] = TSVDND_SW(E,5);

hosvd(tensor(T),eps)
%%
T = zeros(I,I,I);
M = zeros(I,I);
M = 5*e1*e1' + 3*e2*e4' + 2*e4*e2' + e5*e4';
T(:,:,1) = M;
T(:,:,2) = M;
T(:,:,3) = M;
T(:,:,4) = M;
T(:,:,5) = M;

T = permute(T,[1,3,2])
R = modrank(T);
[rankval, minrankindex] = min(R);
[x,gains,iterations,Xa] = TSVDND_SW(T,3);
E1 = T-Xa;

E1c = zeros(I,I);
for i = 1:I
    ei = reshape(E1(:,i,:),[I,I]);
    E1c = E1c+ ei;
end
% Quick test on 3x3x3 tensor
X=zeros(3,3,3);
X(1,1,1)=2; 
X(2,2,2)=-6; 
X(3,3,3)=-4; 
X(1,2,3)=2; 
X(3,1,2)=-1; 
X(2,1,1)=-1; 
X(2,1,3)=4;


W1=eye(3); W2=eye(3); W3=toeplitz([2 -1 0]);
max_level = 2;

[x1,x2,x3,gains,iterations,Xa1] = succR1_SW(X,max_level);
disp('... Successive rank 1 completed ');

[x1,x2,x3,gains,iterations,Xa2]=TSVD3D_SW(X,max_level);
disp('... TSVD completed ');

[x1,x2,x3,gains,iterations,Xa3] = succR1_SW_W(X,max_level,W1,W2,W3);
disp('... Weighted successive rank 1 completed ');

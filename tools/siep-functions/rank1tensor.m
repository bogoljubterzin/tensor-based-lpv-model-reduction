function result = rank1tensor(sigmav,k1,k2,k3,x1,x2,x3)

result = (ttt(tensor(x1(:,k1)),tensor(x2(:,k2))));
result = squeeze(result);
result = ttt(result,tensor(x3(:,k3)));
result = squeeze(result);
result = sigmav.*double(result);

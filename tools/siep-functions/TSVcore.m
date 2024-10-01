function S = TSVcore(T,max_level)
    if class(T) == 'tensor'
        T = double(T);
    end
    [w1,w2,w3,sv,iter,Tapp] = TSVD3D_SW(T,max_level);
    S = zeros(size(T));
    
    for i = 1:size(w1,2)
        for j = 1:size(w2,2)
            for k = 1:size(w3,2)
                t1 = tensor(w1(:,i));
                t2 = tensor(w2(:,j));
                t3 = tensor(w3(:,k));
                Trank1 = squeeze(ttt(t1,t2));
                Trank1 = squeeze(ttt(Trank1,t3));
                S(i,j,k) = innerprod(tensor(T),Trank1);
            end
        end
    end
end
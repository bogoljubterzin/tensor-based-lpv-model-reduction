function Wn = ReachabilityTensors(A,B,n)

    Wn = cell(n,1);
    
    Wn{1} = permute(B,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Wn{i} = ttt(Ap,Wn{i-1},3,1);
    end

end


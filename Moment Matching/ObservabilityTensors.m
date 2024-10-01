function Qn = ObservabilityTensors(A,C,n)
    Qn = cell(n,1);
    Qn{1} = permute(C,[1,3,2]);
    Ap = permute(A,[1,3,2]);
    for i = 2:n
        Qn{i} = ttt(Qn{i-1},Ap,i+1,1);
    end
end
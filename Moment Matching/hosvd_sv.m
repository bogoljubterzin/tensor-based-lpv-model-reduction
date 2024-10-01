function sv = hosvd_sv(T)
    N = numel(size(T));
    sv = cell(N,1);    
    for i = 1:N
        Ti = tens2mat(T,i);
        % ri = rank(Ti,tol);
        svi = svd(Ti);
        sv{i} = svi;
    end
end
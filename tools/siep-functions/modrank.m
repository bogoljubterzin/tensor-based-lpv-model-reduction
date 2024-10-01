function R = modrank(T)
    % created by Bogoljub Terzin; 
    % still an issue if last dimension is 1 it gets deleted
    T = double(T);
    dims = size(T);
    N = length(dims);
    R = zeros(N,1);
    for i = 1:N
        W = tens2mat(T,i);
        R(i) = rank(W);
    end
end

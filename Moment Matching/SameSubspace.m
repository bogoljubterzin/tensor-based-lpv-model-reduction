function q = SameSubspace(A,B)
    A = orth(A);
    B = orth(B);
    if rank(A)~= rank(B)
        q = 0;
        return
    end

    R = [A, B];
    if rank(R) ~= rank(A)
        q = 0
    end
    q = 1;
end
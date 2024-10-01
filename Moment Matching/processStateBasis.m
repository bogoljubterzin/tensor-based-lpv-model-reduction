function [Tx_new,sv_new,sv_imp] = processStateBasis(Vr_hosvd,gainsWn_hosvd,tol)
    N = numel(Vr_hosvd);
    sv = [];
    Tx = [];
    for n = 1:N
        ni = size(Vr_hosvd{n},2);
        sv = [sv, gainsWn_hosvd{n}{1}(1:ni)'];
        Tx = [Tx, Vr_hosvd{n}];
    end
    sv_new = [];
    Tx_new = [];
    sv_imp = [];
    cnt = 1;
    for i = 1:length(sv)
        if sv(i)/sum(sv) > tol
            sv_imp(cnt) = sv(i)/sum(sv)
            sv_new = [sv_new, sv(i)];
            Tx_new = [Tx_new, Tx(:,i)];
            cnt = cnt+1;
        end
    end

end
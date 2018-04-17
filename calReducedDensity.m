function reducedDensity = calReducedDensity(mps,Vmat,para,k)
%   DEPRECATED: Use getObservable()
%
% calculates the reduced density matrix of any single site, mps{k} for size(mps{k})=a_{k-1} x a_k x n_k
%
%   created 24/05/2014 by Florian Schroeder @ University of Cambridge
%

% copied from prepare.m:
% does l -> r sweep to create state in local picture of k

for i = 1:k-1
    if para.useVmat==1
        [Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
        mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
    end
    [mps{i}, U] = prepare_onesite(mps{i}, 'lr',para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
    mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
    para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
    para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
end

% now in form: Al{1}...Al{k-1} M{k} Ar{k+1}...Ar{L}
%   with Al = left-normalized, Ar: right-normalized.

reducedDensity = contracttensors(mps{k},3,[1,2],conj(mps{k}),3,[1,2]);  % contract rD_nm = Mk_abn Mk*_abm
reducedDensity = contracttensors(reducedDensity,2,2,conj(Vmat{k}),2,2);       % contract rD_nj = rD_nm Vmat*_jm
reducedDensity = contracttensors(Vmat{k},2,2,reducedDensity,2,1);       % contract rD_ij = Vmat_in rd_nj

% unnecessary, as does not save mps back:
% % bring mps onto right-normalized form again
% for i = k:-1:2
% %     if para.useVmat==1
% %         [Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);
% %         mps{i} = contracttensors(mps{i},3,3,V,2,2);
% %     end
%     [mps{i}, U] = prepare_onesite(mps{i}, 'rl',para,i);
%     mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
%     mps{i-1} = permute(mps{i-1}, [1, 3, 2]);
%         para=gennonzeroindex(mps,Vmat,para,i);
%         para=gennonzeroindex(mps,Vmat,para,i-1);
% end

end
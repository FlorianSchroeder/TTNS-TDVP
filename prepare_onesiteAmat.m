function [Amat, V, vNE] = prepare_onesiteAmat(Amat,para,sitej)
% Prepare the MPS matrices A for the subsequent updating of Ub
% left normalize A, SVD in n to set focus on Vmat
% TODO: replace SVD by QR, if nargout == 2.

[D1, D2, d_opt] = size(Amat);


if para.parity~='n'
    s=Amat(para.Anzi{sitej});
    sdim=length(s);
    Amat_odd=reshape(s(1:sdim/2),D1*D2/2,d_opt/2);
    Amat_even=reshape(s(sdim/2+1:end),D1*D2/2,d_opt/2);
    [Amat_odd, S_odd, V_odd] = svd2(Amat_odd);
    [Amat_even, S_even, V_even] = svd2(Amat_even);

    vNE = vonNeumannEntropy(S_odd)+1i*vonNeumannEntropy(S_even);
    DB = size(S_odd, 1)+size(S_even,1);


    if D1~=1 && D2~=1;
        Amat=zeros(D1,D2,DB);
        for m=1:DB/2
            submatsize=D1*D2/4;
            Amat(D1/2+1:end,1:D2/2,m)=reshape(Amat_odd(1:submatsize,m),D1/2,D2/2);
            Amat(1:D1/2,D2/2+1:end,m)=reshape(Amat_odd(submatsize+1:end,m),D1/2,D2/2);

            Amat(1:D1/2,1:D2/2,m+DB/2)=reshape(Amat_even(1:submatsize,m),D1/2,D2/2);
            Amat(D1/2+1:end,D2/2+1:end,m+DB/2)=reshape(Amat_even(submatsize+1:end,m),D1/2,D2/2);
        end
    else
        s=cat(1,Amat_odd(:),Amat_even(:));
        Amat=zeros(D1,D2,DB);
        Amat(para.Anzi{sitej})=s;
    end
    V = blkdiag(S_odd* V_odd,S_even*V_even);

else
% for para.parity = 'n'
    Amat = reshape(Amat, [D1 * D2, d_opt]);		% A = A_{(ab),n}
    [Amat, S, V] = svd2(Amat);
    vNE = vonNeumannEntropy(S);
    d_opt = size(S, 1);
    Amat = reshape(Amat, [D1, D2, d_opt]);
    V = S * V;
end
end
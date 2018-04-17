function [B, U, para,results,sv,vNE] = prepare_onesite(A,para,sitej,results)
% Does SVD on the A-MPS matrix
% 1. prepare_onesite(A,para,sitej)
%       is enough for right-norm A
%       sitej only for logging purposes
% 2. prepare_onesite(A,para,sitej,results)
%       writes SV and vNE into results.**{sitej}
% Usage for 'r' sweep:
%   [mps{sitej}, Cn, para,results] = prepare_onesite(mps{sitej}, para,sitej,results);
% -> Saves new Bond dimension in para.D(sitej)!
% Usage for 'l' sweep:
%   [mps{sitej}, Cn, para,results] = prepare_onesite(mps{sitej},para,sitej,results);
% -> Saves new Bond dimension in para.D(sitej-1)!
%
%   Changed:
%       - FS 20/10/2014: replaced direction by para.sweepto
%
% TODO: Use QR instead of SVD to speed up if SV and nVE are not needed
% 		Use [Q, R] = qr(A,0) for vertical matrix to get only Q1 and R from (Q1 Q2) (R 0).

[D1, D2, d] = size(A);					% d is site dimension
% method = 'qr';
saveNewD = 0;
if nargin > 2 && ~isempty(sitej)
	saveNewD = 1;
end
saveResults = 0;
if nargin > 3
	saveResults = 1;
else
	results = [];
end

switch para.sweepto
    case 'r'
        if para.parity == 'n'
            %% normal parity
            A = permute(A, [3, 1, 2]);
            A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
            [B, S, U] = svd2(A);			% Could also use QR decomposition if nargin !=5
%           if nargin ~=5
%               [B, U] = qr(A);                 % Doesn't need wrapping
%               DB = size(U,1);
%
%           else
%               end
            DB = size(S, 1);				% new a2 dimension of B
			if saveNewD
	            para.D(sitej)=DB;
			end
            B = reshape(B, [d, D1, DB]);
            B = permute(B, [2, 3, 1]);
            U = S * U;                      % only for SVD

        else
            %% special parity
            A=permute(A,[3,1,2]);
            s=A(para.Anzilr{sitej});
            %s=nonzeros(A);
            sdim=length(s);
            assert(sdim==D1*D2*d/2);
            A_odd=reshape(s(1:sdim/2),[D1*d/2,D2/2]);
            A_even=reshape(s(sdim/2+1:end),[D1*d/2,D2/2]);
            [B_odd,S_odd,U_odd]=svd2(A_odd);
            [B_even,S_even,U_even]=svd2(A_even);

            vNE_odd=vonNeumannEntropy(S_odd);
            vNE_even=vonNeumannEntropy(S_even);
            vNE=vNE_odd+1i*vNE_even;
            DB_odd=size(S_odd,1);
            DB_even=size(S_even,1);
            DB=DB_odd+DB_even;
            sv_odd=diag(S_odd);
            sv_even=diag(S_even);
            sv=sv_odd+1i*sv_even;

            B=zeros(d,D1,DB);
            if D1~=1;
                for m=1:DB/2
                    submatsize=D1*d/4;
                    B(d/2+1:end,1:D1/2,m)=reshape(B_odd(1:submatsize,m),d/2,D1/2);
                    B(1:d/2,D1/2+1:end,m)=reshape(B_odd(submatsize+1:end,m),d/2,D1/2);

                    B(1:d/2,1:D1/2,m+DB/2)=reshape(B_even(1:submatsize,m),d/2,D1/2);
                    B(d/2+1:end,D1/2+1:end,m+DB/2)=reshape(B_even(submatsize+1:end,m),d/2,D1/2);
                end
            else %Suppose the left most parity is even
                assert(DB==2);
                for m=1:DB/2
                    if para.parity=='e'
                        B(1:d/2,1,m)=B_odd(:,m);
                        B(d/2+1:end,1,m+DB/2)=B_even(:,m);
                    else
                        B(2,1,1)=B_odd;
                        B(1,1,2)=B_even;
                    end
                end
            end
            B=permute(B,[2,3,1]);
            U=blkdiag(S_odd*U_odd,S_even*U_even);

        end
    case 'l'
        if para.parity == 'n'
            %% normal parity
			A = permute(A, [1, 3, 2]);
			A = reshape(A, [D1, d * D2]);
% 			A = reshape(A, [D1, D2 * d]);								% Do not use!

            % decompose: A_(l,r*n) = U_(l,a')*S_(a',a)*B_(a,r*n)
            [U, S, B] = svd2(A);            % since m<n: B = new A
            DB = size(S, 1);
			if saveNewD && sitej ~= 1
	            para.D(sitej-1) = DB;
			end
% 			B = reshape(B, [DB, D2, d]);								% Do not use!
			B = reshape(B, [DB, d, D2]); B = permute(B, [1, 3, 2]);
            % create focused center: C(n)_(l,a) = U_(l,a')*S_(a',a)
            U = U * S;
			if saveResults
	            sitej = sitej - 1;				% needed for saving into the right position of results.Amat_sv
			end
        else
            %% special parity
            A=permute(A,[2,3,1]);
            s=A(para.Anzirl{sitej});
            %s=nonzeros(A);
            sdim=length(s);
            assert(sdim==D1*D2*d/2);
            A_odd=reshape(s(1:sdim/2),[D2*d/2,D1/2]);
            A_even=reshape(s(sdim/2+1:end),[D2*d/2,D1/2]);
            [B_odd,S_odd,U_odd]=svd2(A_odd);
            [B_even,S_even,U_even]=svd2(A_even);

            vNE_odd=vonNeumannEntropy(S_odd);
            vNE_even=vonNeumannEntropy(S_even);
            vNE=vNE_odd+1i*vNE_even;
            DB_odd=size(S_odd,1);
            DB_even=size(S_even,1);
            DB=DB_odd+DB_even;
            sv_odd=diag(S_odd);
            sv_even=diag(S_even);
            sv=sv_odd+1i*sv_even;

            B=zeros(D2,d,DB);
            if D2~=1;
                for m=1:DB/2
                    submatsize=D2*d/4;
                    B(D2/2+1:end,1:d/2,m)=reshape(B_odd(1:submatsize,m),D2/2,d/2);
                    B(1:D2/2,d/2+1:end,m)=reshape(B_odd(submatsize+1:end,m),D2/2,d/2);

                    B(1:D2/2,1:d/2,m+DB/2)=reshape(B_even(1:submatsize,m),D2/2,d/2);
                    B(D2/2+1:end,d/2+1:end,m+DB/2)=reshape(B_even(submatsize+1:end,m),D2/2,d/2);
                end
            else %Suppose the right most parity is even
                for m=1:DB/2
                    B(1,1:d/2,m)=B_odd(:,m);
                    B(1,d/2+1:end,m+DB/2)=B_even(:,m);
                end
            end
            B=permute(B,[3,1,2]);
            U=blkdiag(S_odd*U_odd,S_even*U_even);
            U=transpose(U);

        end
end

sv  = diag(S);
vNE = vonNeumannEntropy(S);

if saveResults && sitej >= 1
    results.Amat_vNE(sitej) =vNE;
    if length(sv)==para.D(sitej) || (length(sv)==para.D(sitej)/2 && (para.parity~='n'))
        results.Amat_sv{sitej}=sv;
    end
end

end

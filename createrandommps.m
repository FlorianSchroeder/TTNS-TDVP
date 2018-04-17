function [mps] = createrandommps(para)
% Create random MPS as the starting point for variation
% Using optimal boson basis --> A-matrix dimension: D x D x d_opt  ; where d_opt is dimension of physical index.
L	  = para.L;
D	  = para.D(end);
d_opt = para.d_opt;

mps = cell(1, L);

if para.parity=='n'
    mps{1} = randn(1,D,d_opt(end,1))./sqrt(D*d_opt(end,1));
    mps{L} = randn(D,1,d_opt(end,L))./sqrt(D*d_opt(end,L));
	for i = 2:(L - 1)
        mps{i} = randn(D,D,d_opt(end,i));
		mps{i} = mps{i}./sqrt(numel(mps{i}));
	end
else
    mps{1} = zeros(1, D, d_opt(1));
    if para.parity=='e'
        mps{1}(1,1:D/2,1:d_opt(1)/2) = randn(1,D/2,d_opt(1)/2)./(D/2);
        mps{1}(1,D/2+1:D,d_opt(1)/2+1:end) = randn(1,D/2,d_opt(1)/2)./(D/2);
    else
        mps{1}(1,D/2+1:end,1:d_opt(1)/2) = randn(1,D/2,d_opt(1)/2)./(D/2);
        mps{1}(1,1:D/2,d_opt(1)/2+1:end) = randn(1,D/2,d_opt(1)/2)./(D/2);
    end

    mps{L} = zeros(D, 1, d_opt(L));
    mps{L}(1:D/2,1,1:d_opt(L)/2) = randn(D/2,1,d_opt(L)/2)./(D/2);
    mps{L}(D/2+1:D,1,d_opt(L)/2+1:end) = randn(D/2,1,d_opt(L)/2)./(D/2);

    for i = 2:(L - 1)
        mps{i} = zeros(D, D, d_opt(i));
        for s=1:d_opt(i)/2
            mps{i}(1:D/2,D/2+1:D,s)=randn(D/2,D/2)./(D/2);
            mps{i}(D/2+1:D,1:D/2,s)=randn(D/2,D/2)./(D/2);
        end
        for s=d_opt(i)/2+1:d_opt(i)
            mps{i}(1:D/2,1:D/2,s)=randn(D/2,D/2)./(D/2);
            mps{i}(D/2+1:D,D/2+1:D,s)=randn(D/2,D/2)./(D/2);
        end
    end
end


end

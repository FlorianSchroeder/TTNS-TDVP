function c=calboson1siteCreationCorrelator_SBM1(mps0,Vmat0,para)
% Calculate bosonic 1-site correlators: c_i = <a_i>
% The operator on the spin site is set to zero
% modify this to get left and right chain occupation! perhaps in calbosonocc_SBM2.m
% copied from calbosonocc_SBM1
% Modified:
%	FS 16/03/2014:	created


c_op = cell(1,para.L);			% here: M = 1

for j=1:para.L
    if j~=para.spinposition
        if para.foldedChain == 0
            [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);
            c_op{1,j} = bp;				% for c_j
        %Modification for folded 2chain model!! Not perfect or right yet!
        elseif para.foldedChain == 1
            if para.parity == 'n'
                [bp,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);     % gives [bp,bm,n]
                idm = eye(size(n));
                bpr = kron(bp,idm);
            else
                [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);   % Why without sqrt??
                [bpr,~,~,~,~,~]=paritykron(bp,para.bosonparity);
            end
            c_op{1,j} = bpr;
        else

        end
    else
        c_op{1,para.spinposition}=zeros(para.dk(j));
    end
end

c = correlator_allsites(c_op,mps0,Vmat0);

end
function nx=calbosonocc_SBM1(mps0,Vmat0,para,results)
%   DEPRECATED: Use getObservable()
%
%Calculate the boson occupation on x chain
%The operator on the spin site is set to zero
%modify this to get left and right chain occupation! perhaps in calbosonocc_SBM2.m
% Modified:
%	FS 23/01/2014:	- Introduced '~' to ignore unused returned values
%					- support for folded Chain models
%   FS 10/03/2014:  - Using correlator_allsites which is more general
%

n_op = cell(1,para.L);

for j=1:para.L
    if j~=para.spinposition
        if para.foldedChain == 0
            [~,~,n] = bosonop(para.dk(j),para.shift(j),para.parity);
            n_op{j} = n;
        %Modification for 2chain model!! Not perfect or right yet!
        elseif para.foldedChain == 1
            if para.parity == 'n'
                [~,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
                idm = eye(size(n));
                nr = kron(n,idm);
            else
                [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);   % Why without sqrt??
                [~,~,nr,~,~,~]=paritykron(bp,para.bosonparity);
            end
            n_op{j} = nr;
        else

        end
    else
        n_op{para.spinposition}=zeros(para.dk(j));
    end
end

%
nx = correlator_allsites(n_op,mps0,Vmat0);

end
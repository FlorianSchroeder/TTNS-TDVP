function bosonshift=calbosonshiftperSpin_SBM1(mps0,Vmat0,para,results)
% Calculate boson shift x, x^2, var(x); for each spin direction
%   (bp + bm) (1+- sz)/2
%The operator on the spin site is set to zero
% Modified:
%	FS 22/01/2014:	- changed to using para.foldedChain.
%
% expectation_allsites is old and should be replaced by correlator_allsites

x_opx=cell(1,para.L);
x2_opx=cell(1,para.L);
[~,~,sigmaZ]=spinop(para.spinbase);
for j=1:para.L
    if j~=para.spinposition
        if para.foldedChain == 1            % Constructs Supersite Operators
            [bp,bm,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
            idm=eye(size(n));
            bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
            x_opx{j}=sqrt(2)/2.*(bpx+bmx);
            x2_opx{j}=x_opx{j}*x_opx{j};
        elseif para.foldedChain == 0
            [bp,bm,n] = bosonop(para.dk(j),para.shift(j),para.parity);
            x_opx{j} = sqrt(2)/2.*(bp+bm);
            x2_opx{j} = x_opx{j}*x_opx{j};
        end
    else
        x_opx{j}=zeros(para.dk(j));
        x2_opx{j}=zeros(para.dk(j));
    end
end
% from here: copied and modified correlator_allsites

bosonshift.xUp =correlator_allsites(x_opx,mps0,Vmat0,(eye(2)+sigmaZ)/2);
bosonshift.xDown =correlator_allsites(x_opx,mps0,Vmat0,(eye(2)-sigmaZ)/2);
bosonshift.xsquareUp = correlator_allsites(x2_opx,mps0,Vmat0,(eye(2)+sigmaZ)/2);
bosonshift.xsquareDown = correlator_allsites(x2_opx,mps0,Vmat0,(eye(2)-sigmaZ)/2);
%bosonshift.xvariant = sqrt(bosonshift.xsqure-bosonshift.x.^2);
%bosonshift.xerror=mean(abs(para.shift-bosonshift.x));
end

function [c]=correlator_allsites(c_op,mps,Vmat,SpinOp)
% Calculate the expectation value of the multi-site correlator "c_op" for all bath sites.
% Only for successive operators. Computes any number M! M counts < a_1...a_m > number of sites being correlated at once
% Can be used instead of expectation_allsites if M = 1;
% Format of c_op: M x N cell
%   for i x j: in j-th column: all operators to compute correlator(k).
%       First row i=1 : operator for site j
%       2nd   row i=2 : site j+1
%       ....
% c(j)=<\psi|c_op{j}|\psi>
% e.g.: c_op{j} = {bp{j} , bp{j+1}} == bp{j}*bp{j+1}*....*bp{j+m-1}
%	Modified:
%		FS 10/03/2014
%   Commented:
%       FS 07/05/2014

% N=length(c_op);		%old definition, DELETE
[M, N] = size(c_op);											% N != Length, M = number of
assert(N==length(mps) && N==length(Vmat));

c=zeros(1,N);
cdset=cell(1,N);					% Used for complete representation including unity operators for one c_op{j};

for j=1:N
        cdset{1,j}=eye(size(Vmat{j},1));
end

%% IMPORTANT, Only part modified
% sets spinOperator on site 1 to get spin-sensitive shifts
cdset{1,1} = SpinOp;


%%
for j=1:N
% This should be a bit faster than for-loop
    nterms = min(N-j+1,M);                                          % caps last terms for the end of the chain
    temp = cdset(1,j:min(j+M-1,N));                                 % save unity ops
    cdset(1,j:min(j+M-1,N)) = c_op(1:nterms,j)';                    % use c_ops, produce: 1 .... 1 n 1 .... 1
    c(j) = expectationvalue(cdset,mps,Vmat,mps,Vmat);               % calculate correlator
    cdset(1,j:min(j+M-1,N)) = temp;                                 % load back unity ops
end
end
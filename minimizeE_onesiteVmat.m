% ******************** one-site optimization ***********************
function [B, E] = minimizeE_onesiteVmat(op, A, Blaststep,para)
% A is the local MPS A-matrix
% B is the optimized local optimal phonon basis transform matrix V
% Blaststep is the to be optimized Vmat matrix.

% Commented by Florian Schroeder 31/01/2014

tol = para.eigs_tol;
M   = para.M;
[Dblast,dlast]=size(Blaststep);						% dk x d_opt

[op] = H_Eff(A, []  , 'V' , op, para);				% multiply A into Hleft, Hright, Opleft, Opright

d = size(op.HleftA, 1);								% = d_opt
if iscell(op.h1j)
% 	Db = prod(cell2mat(cellfun(@(x) size(x,1),op.h1j, 'UniformOutput',false)));		% = dk for multi-chain Hamiltonians
	Db = prod(para.dk(:,para.sitej));				% = dk for Multi-Chain
else
	Db = size(op.h2j{1,1}, 1);                      % = dk
end

% projection on orthogonal subspace
%if ~isempty(P), Heff = P' * Heff * P; end

% optimization
opts.disp = 0;
opts.tol = tol;
%opts.p=12;
opts.issym =1;
sigma='sa';
if para.complex==1
    opts.isreal = false;
    sigma='sr';
end
assert(Db*d == Dblast*dlast);

if para.parity~='n'
    opts.v0=vertcat(reshape(Blaststep(1:Db/2,1:d/2),Db*d/4,1),reshape(Blaststep(Db/2+1:end,d/2+1:end),Db*d/4,1));
    [Bvec, E]=eigs(@(x) HmultVmat(x, op, Db, d, M, para.parity), Db*d/2,1,sigma,opts);
    B = zeros(Db, d);
    B(1:Db/2,1:d/2)=reshape(Bvec(1:Db*d/4),Db/2,d/2);
    B(Db/2+1:end,d/2+1:end)=reshape(Bvec(Db*d/4+1:end),Db/2,d/2);
else
    opts.v0   = reshape(Blaststep,numel(Blaststep),1);			% = Db * d
    [Bvec, E] = eigs(@(x) HmultVmat(x, op, Db,d, M,para.parity), Db*d,1,sigma,opts);
    B         = reshape(Bvec, [Db, d]);
end
% fprintf('\n Vmat E: %g',E);				% debug


end
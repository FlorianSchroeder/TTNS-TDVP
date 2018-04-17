function [A, E, op] = minimizeE_onesiteA(op, B, Alaststep,para)
% B can be Vmat{sitej}; for spinsite: = eye()

% Commented by Florian Schroeder 31/01/2014
% Modified:
%	- FS 05/08/2015:	- Moved OBB trafo of h1j, h2j into H_Eff; return op

% ******************** one-site optimization ***********************
tol = para.eigs_tol;
DAl = size(op.Hleft, 1);					% Bond dimension Dleft of A to left
DAr = size(op.Hright, 1);					% Bond dimension Dright of A to right
d = size(B, 2);								% = d_opt of Vmat
[ll,rr,dd] = size(Alaststep);				% = Dleft, Dright, d_opt	from A
if para.parity~='n'
    nzAlast=Alaststep(para.Anzi{para.sitej});
end

% apply OBB to h1 & h2 terms -> h1jOBB, h2jOBB
op = H_Eff([]  , B, 'A' , op, para);		% will be returned as well!

% projection on orthogonal subspace
%if ~isempty(P), Heff = P' * Heff * P; end

%% optimization

% options for eigs()
opts.disp = 0;
opts.tol = tol;
%opts.p=min(24,DAl*DAr*d);
opts.issym =1;
sigma='sa';							% smallest algebraic eigenvalue
if para.complex==1
    opts.isreal = false;
    sigma='sr';						% smallest real part eigenvalue
end
assert(DAl*DAr*d == ll*rr*dd)

if para.parity=='n'
    opts.v0   = reshape(Alaststep,numel(Alaststep),1);						% starting guess, old MPS A-matrix
    [Avec, E] = eigs(@(x) HmultA(x, op, DAl, DAr, d, para.M,para.parity),DAl*DAr*d,1,sigma,opts);
    A         = reshape(Avec, [DAl, DAr, d]);
else
    opts.v0=nzAlast;
    [Avec, E]=eigs(@(x) HmultA(x, op, DAl, DAr, d, para.M,para.parity,para.Anzi{para.sitej}),DAl*DAr*d/2,1,sigma,opts);
    A = zeros(DAl, DAr, d);
    A(para.Anzi{para.sitej})=Avec;
end
% fprintf('\n Amat E: %g',E);			% Debug

end

function [n] = expectation_allsites(n_op,mps,Vmat,Cleft)
% Calculate the expectation value of the 1-site local operator "n_op" for all sites.
% performs 2*3*N contractions
% assumes precise right-normalization
% output: n(i,j)=<\psi|n_op{i,j}|\psi>
%
% Usage:
%   [n] = expectation_allsites(n_op, mps, Vmat)
%           standard
%   [n] = expectation_allsites(n_op, mps, cell(1,L)  )
%           n_op is already brought onto OBB
%           can be used for Multi-Chain code
%
%   n_op: i x L cell, where i can stand for subchains


[M, L] = size(n_op);
% assert(L == length(mps));			% would also work for N < dim(mps)

n = zeros(M,L);

Cl = [];						% contains left part in effective j-basis
if nargin == 4
	Cl = Cleft;
end

for j = 1:L
	% j == 1 -> Cl = [];
	for m = 1:M
		% calculate site operator
		n(m,j) = trace(updateCleft(Cl,mps{j},Vmat{j},n_op{m,j},mps{j},Vmat{j}));	% this is the result!
	end
	% take Cleft to next site
	Cl = updateCleft(Cl,mps{j},Vmat{j},[],mps{j},Vmat{j});
end
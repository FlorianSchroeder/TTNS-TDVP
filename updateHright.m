function [Hright] = updateHright(Hright, h1jOBB, Opstorage2, B, BUb, h2jOBB, A, AUb, M, varargin)
% called as op.Hlrstorage{j} = updateHright(op.Hlrstorage{j + 1},  op.h1term{j}, op.Opstorage(:,2,j+1), mps{j}, Vmat{j},   op.h2term(:,1,j),  mps{j},  Vmat{j},M);	for middle sites	in initstorage.m

% As called from initstorage.m:
% Hright = op.Hlrstorage{j + 1}; == 0 for j=L
% h1j = op.h1term{j};	current on-site energy
% Opstorage2{}  = op.Opstorage(:,2,j+1); per site and  contains M matrices for each term in sum
% B = mps{j};
% BUb = Vmat{j};
% h2j = op.h2term(:,1,j);
% A = mps{j};
% AUb = Vmat{j}

	Hright = updateCright(Hright, B, [], [], A, []);						% transform old Hright from effective basis j+1 into j

	if ~isempty(h1jOBB)
		Hright = Hright + updateCright([],B,[],h1jOBB,A,[]);				% transform H1{j} into effective right basis of j-1
	end

	for m=1:M
		if ~isempty(h2jOBB{m}) && ~isempty(Opstorage2{m})
			Hright = Hright + updateCright(Opstorage2{m},B,[],h2jOBB{m},A,[]);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
		end
	end
end
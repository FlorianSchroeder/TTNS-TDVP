function [Cright] = updateCright(Cright, B, BUb, X, A, AUb)
%% function [Cright] = updateCright(Cright, B, BUb, X, A, AUb)
%
% Computes contraction to put operator X in effective basis of B and A.
% BUb and AUb are the V matrices to express A in optimal boson basis.
% B and A are left and right MPS A-matrices
%
% 1. transforms Cright from effective right basis of site k into effective right basis of k-1 if:
%		X = []	; B = A = mps{k}	; BUb, AUb = Vmat{k}
% 2. transforms operator X of site k into effective right basis of site k-1 if:
%		Cleft = []	; B = A = mps{k}	; BUb = AUb = Vmat{k}
% 3. Does 1 and 2 altogether if everything is given.
%
% Called as updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});				--> express in effective basis
% Called as updateCright(Hright, mps{j}, Vmat{j}, [], mps{j}, Vmat{j}); from updateHright.m		-->
%		Hright = op.Hlrstorage{j + 1} which == 0 for j=L
%
% if Cright input == 0 --> output == 0
% Commented by Florian Schroeder 29/01/2014
%
% Modified
%	- 21/12/14 FS: replaced OBB contraction by faster matrix product

skipX = 0;
if isempty(X) && isempty(BUb)
	skipX = 1;
elseif ~isempty(X) && isempty(BUb)							% if no Vmat given or X is in OBB already: XOBB
	newX = X;
elseif isempty(X)
	newX = BUb' * AUb;
else
	% transform X into OBB using Vmat, do (Vmat^†) . X . Vmat.
	newX = (BUb' * X) * AUb;
end

% if Cright = eye: contract X (in OBB)  with A matrices to transform into effective basis representation.
% do contraction:  C_fa = B*_fde  A_abc  Xnew_ec  C_db 	where 3rd indices are running over n_k
if ~skipX
	A      = contracttensors(A, 3, 3, newX.', 2, 1);			% Anew_abe = A_abc  Xnew_ec
end
if ~isempty(Cright)
	A = contracttensors(A, 3, 2, Cright.', 2, 1);				% Anew_aed = A_abe C_db
	Cright = contracttensors(conj(B), 3, [2, 3], A, 3, [3, 2]);	% Contracts 2-3 and 3-2
else
	%have Anew_ade since C_db = 1_db -> neglected
	Cright = contracttensors(conj(B), 3, [2, 3], A, 3, [2, 3]);	% Contracts 2-2 and 3-3
end
end
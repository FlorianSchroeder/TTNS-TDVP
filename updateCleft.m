function [Cleft] = updateCleft(Cleft, B, BUb, X, A, AUb)
% Lookup: Schollwoeck 2011 4.2.1: Contracts Left contraction residue of chain with current site A-matrices and Operators.
% 1. transforms Cleft from effective left basis of site k into effective left basis of k+1 if:
%		X = []	; B = A = mps{k}	; BUb = AUb = Vmat{k} or []
% 2. transforms operator X of site k into effective left basis of site k+1 if:
%		Cleft = []	; B = A = mps{k}	; BUb = AUb = Vmat{k}
% 3. Does 1 and 2 altogether if everything is given.

% Commented by Florian Schroeder 29/01/2014
% Modified:
%	FS 05/08/2015: better order to minimize permutations!
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

% do contraction:  C_fb = B*_dfe  (C_da  A_abc  newX_ec)_dbe	where 3rd indices are running over n_k
if ~skipX
	A = contracttensors(A,3,3, newX.',2,1);					% Anew_abe = A_abc newX_ec
end
if ~isempty(Cleft)
	A = contracttensors(Cleft,2,2, A,3,1);					% Anew_dbe = C_da A_abe
end

Cleft = contracttensors(conj(B),3,[1 3], A,3,[1 3]);		% Cleft_fb = B*_dfe Anew_dbe
end
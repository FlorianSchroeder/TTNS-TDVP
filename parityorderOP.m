function newmat=parityorderOP(mat)
% parityorderOP - Reorders Boson Operators
%	B = parityorderOP(A)
%
%	A: Boson op in local basis [m:0]
%	B: Boson op in reordered parity basis [o,e]
%
%	Reordering the local boson operators into 2*2 block matrices according to
%	parity odd and even (number). The shape of the block matrix is :
%
%	 oo oe
%	 eo ee
%
%	 So operator
%	 0 sqrt(3)   0       0
%	 0   0     sqrt(2)   0
%	 0   0       0     sqrt(1)
%	 0   0       0       0
%
%	 will become:
%	 0   0     sqrt(3)   0
%	 0   0       0     sqrt(1)
%	 0 sqrt(2)   0       0
%	 0   0       0       0
%
% Modified:
%		FS 07/07/15: 43% speedup using rotation operator

% assert input
[m,n]=size(mat);
assert(mod(m,2)==0 && m==n, 'parityorderOP needs square matrix inputs of even dimensions');

% create reorder operator, perhaps export into toolbox?
U = sparse(1:m, [1:2:m,2:2:m],1,m,n);			% nicer
% U has 1s in (row,column) = (1:m, [odd,even])

% Apply reordering
newmat = (U*mat*U');

end
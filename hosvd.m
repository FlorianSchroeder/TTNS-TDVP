function [S, U, sv, para] = hosvd(T, para, dim)
%% [S U sv, para] = hosvd(T, para, dim)
% Sequentially truncated High Order SVD of a Tensor
%	[S U sv] = HOSVD(T)
%	[S U sv, para] = HOSVD(T, para)
%	[S U sv, para] = HOSVD(T, para, dim)
%
%	T    : Tensor
%	para : VMPS parameters for truncation and dimension tracking: enables truncation
%	dim  : if given, final "focus" in S will be on this Dim
%			U{dim} will be sv*U'
%
%	S    : sequentially decomposed core tensor
% 	U    : matrices for each dimension, [] for excluded dims
%	sv   : singular values
%
%	eg [S, U, sv, para] = hosvd(randn(5,5,5), para, [0 1 1])
% is only some experimental code!


dIn = size(T); dOut = dIn;
o = length(dIn);							% order of tensor
svmintol = 10^-6.5; svmaxtol = 10^-6;
dmin = 2;
adaptST = 1;								% adaptive SVD-based Truncation.
if nargin == 3 && dim ~= o					% for dim == o do adaptive truncation, then take focus to A
	adaptST = 0;							% do only dim-fixed truncation to dOut
end
if ~isempty(para)
	svmintol = para.svmintol;
	svmaxtol = para.svmaxtol;
	dmin     = para.d_opt_min;
	dOut     = para.d_opt(:,para.sitej);	% only used if adaptST == 0
end

U = cell(1,o);   sv = cell(1,o);
S = T;            d = dIn;
if ~isempty(dim)
	unFoldDim = mod(dim,o) + 1;						% first dim to operate on
end
for i = 1:o			% always do cyclic HOSVD!			% operate on dimension mod(dim+i-1,o)+1
	[S, d] = tensShape(S,'unfold',unFoldDim,d);
	% SVD in dimension (i)
	dim_i = mod(dim+i-1,o)+1;							% current working dim
	[U{dim_i}, sv{dim_i}, S] = svdTruncate(S);
	d(1) = size(S,1);
	S = reshape(S, d);
	unFoldDim = 2;
end
if ~isempty(para)
	para.d_opt(:,para.sitej) = cellfun('length' , sv)';	% perhaps do outside of hosvd!
end

	function [U,sv,V] = svdTruncate(A)
		%% [U,sv,V] = svdTruncate(A)
		%	computes SVD of A and truncates to
		%	 - fixed dimension
		%	 - adaptive dimension wrt singular values
		%	 - computes V = U'*A. U*V = A
		[U,sv] = svd2(A);
		sv = diag(sv);
		if ~adaptST
			keepdims = 1:dOut(dim_i);						% dim-fixed truncation
		elseif length(sv) <= dmin
			V = U' * A;									% only apply trafo without
			return;
		else
			%% adaptively Truncate A dims
			keepdims = find(sv > svmintol);
			% If smallest SV too large, keep 1 more
			if (sv(keepdims(end)) > svmaxtol)
				if keepdims(end) < length(sv)
					keepdims = [keepdims;keepdims(end)+1];
				else
					return;		% no truncation possible
				end
			end
			if length(keepdims) < dmin					% keep at least Dmin bonds
				keepdims = 1:dmin;
			end
		end
		if length(keepdims) < length(sv)
            U  = U(:,keepdims);							% keep columns
            sv = sv(keepdims);
		end
		V  = U' * A;									% this sequentially truncates the tensor
	end
end
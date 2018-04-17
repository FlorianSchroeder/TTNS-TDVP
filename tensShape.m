function [T, dOut] = tensShape(T,select,i,d)
% Use as:
%	tensShape(T, 'unfold', i_unfold, [d1,d2,d2])
%	tensShape(T, 'fold', i_unfold, [d1,d2,d2])
%	tensShape(T, 'unfoldiso', i_unfold, [d1,d2,d2])
%	tensShape(T, 'foldiso', i_unfold, [d1,d2,d2])
%	tensShape(T, 'foldrotunfoldiso', i_unfold, [d1,d2,d2])
switch lower(select)
	case 'unfold'
		T = tensUnfold(T,i,d);

	case 'fold'
		T = tensFold(T,i,d);

	case 'unfoldiso'
		T = tensUnfoldIso(T,i,d);

	case 'foldiso'
		T = tensFoldIso(T,i,d);

	case 'foldrotunfoldiso'
		T = tensFoldRotUnfoldIso(T,i,d);
end

	function T = tensUnfold(T,i,d)
	%% T = tensUnfold(T,i.d)
	%	unfolds tensor along i by circular permutation
	% 	d gives initial size(T)
		if i ~= 1
			r = length(d);
			T = permute(T,[i:r,1:i-1]);					% permutes only i-1 times
		end
		dOut = size(T);
		T = reshape(T,[d(i),prod(d)/d(i)]);
	end

	function T = tensFold(T,i,d)
	%% T = tensFold(T,i,d)
	%	folds a tensor which was unfolded along dimension i
	%	d: vector indicates the original dimensions [n1,n2,n3,..];
		r = length(d);
		T = reshape(T,[d(i:r),d(1:i-1)]);
		if i ~= 1
			i = r-(i-1)+1;								% need +1 since [i:r,1:i-1] only shifts by i-1
			T = permute(T,[i:r,1:i-1]);
		end
		dOut = d;
	end

	function T = tensUnfoldIso(T,i,d)
		%% T = tensUnfoldIso(T,i,d)
		%	unfolds tensor by transposition such that T: n(not in i) x n(in i)
		%	Input:
		%		i: array indicating the dimensions along the columns
		%		T: tensor n_1 x ... n_k
		%	Output:
		%		T: prod(n_j not in i) x prod(n_i(1)*n_i(2)*...)
		%
		j = 1:length(d);
		j(i) = [];
		T = permute(T,[j,i]);						% i to the end
		dOut = [d(j),d(i)];
		T = reshape(T,[prod(d(j)), prod(d(i))]);

	end

	function T = tensFoldRotUnfoldIso(T,i,d)
		%% T = tensFoldRotUnfoldIso(T,i,d)
		%	Folds tensor with d (if necessary), circularly transposes dimensions 1:end-1 by i steps while keeping last dimension fixed,
		%    then Unfolds again
		%	Input:
		%		i: number of circular transpositions
		%		d: T dimensions of current transposition
		%		T: prod(n_j not in i) x prod(n_i(1)*n_i(2)*...)
		%	Output:
		%		T: prod() x prod(n(i)*n_k)
		%		dOut: new dimensions
		%
		if ndims(T) ~= length(d)
			T = reshape(T,d);
		end

		if i ~= 0
			r = length(d);
			T = permute(T,[(i+1):r-1,1:i,r]);					% permutes only i-1 times
		end
		dOut = size(T);
		T = reshape(T,[prod(dOut(1:r-2)),prod(dOut(end-1:end))]);
	end


end
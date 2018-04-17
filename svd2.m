function [U, S, V] = svd2(T)
% SVD and returns U, S, V' for any matrix.
% Transposes T in case m < n:  to make code faster!
% Comparing performance for m < n:
%	svd(T', 0) < svd(T', 'econ') < svd(T, 'econ') <<<< svd(T)

[m, n] = size(T);
if nargout == 3
	if m >= n
		[U, S, V] = svd(full(T), 0);
	else
		[V, S, U] = svd(full(T'), 0);
	end

	V = V';
else
	% useful for hosvd() -> speedup because no storage of V
	if m >= n
		[U, S, ~] = svd(T, 0);
	else
		[~, S, U] = svd(T', 0);
	end
end
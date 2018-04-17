function [Q,R,err] = rrQR(A,r,epsilon)
% function [Q,R,err] = rrQR(A,r,epsilon)
%	Computes the randomized reduced rank QR of m x n matrix A, m>= n
%	approximates to rank r (r <= m) with error limit epsilon
%		err gives final error
%	Q: m x r
%	R: r x n

% outline of algorithm can be found in Phys. Rev. E 91, 063306 (2015), Tamascelli et al.

% 1. check input for suitability
d = size(A);
assert( d(1) >= d(2) && r <= d(2), 'Please use matrix with m >= n >= r.');
p = floor(min(0.15*r, d(2)-r));			% oversampling needed to capture rank r; governs precision
q = 10;							% maximum power iteration steps

% Initialise R:
R = zeros(r+p,r+p);				% As A = Q_q * (R_i R~_i)^q * R_0

% 2. construct randomized Omega
Omega = randn(d(2),r+p);		% set of (r+p) lin indep. vectors
Y0 = A*Omega;					% Y0 : m x r+p
[Q0, ~] = qr(Y0,0);				% expensive
for ii = 1:q
% 	Y1 = A* (A'*Q0);			% Yj  : m x r+p
	Y1 = A'*Q0;					% Yj~ : n x r+p
	[Q1,~] = qr(Y1,0);			% cheap
	Y1 = A*Q1;					% Yj  : m x r+p
	[Q1,~] = qr(Y1,0);			% expensive
	err = r-sum(sum(conj(Q0).*Q1));		% = r - trace(Q0'*Q1);
	if err < r/10
		fprintf('q = %d, err = %g, r = %d\n',ii, err, r);
		break;
	end
	Q0 = Q1;
end

% R = Q0\A;
R = Q0'*A;
Q = Q0;


end
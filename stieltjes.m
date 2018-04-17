% STIELTJES Discretized Stieltjes procedure.
%
%    Given the discrete inner product whose nodes are contained
%    in the first column, and whose weights are contained in the
%    second column, of the nx2 array xw, the call ab=STIELTJES(n,xw)
%    generates the first n recurrence coefficients ab of the
%    corresponding discrete orthogonal polynomials. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab.
%
% Modified by Florian Schroeder @ University of Cambridge
%	- output is ab = [alpha, sqrt(beta_i+1)]
%				   = [omega_n/g, t_n/g]
%			as described in Chin2010
%			usually: g=1;

function ab = stieltjes(N,xw)
tiny=10*realmin;
huge=.1*realmax;
%
% Remove data with zero weights
%
xw   = xw(xw(:,2)~=0,:);
xw   = sortrows(xw,1);
Ncap = size(xw,1);				% number of datapoints
%
if(N<=0|N>Ncap), error('N in sti out of range'), end
ab = zeros(N,2);

s0 = ones(1,Ncap)*xw(:,2);						% <pi_0,pi_0>_mu
ab(1,1)=xw(:,2)'*xw(:,1)/s0;					% alpha_0 = <x pi_0,pi_0>_mu / <pi_0,pi_0>_mu
ab(1,2)=s0;										% store ||pi_o||^2; unused!

if N==1, return, end
p1=zeros(Ncap,1); p2=ones(Ncap,1);				% containing polynomials: p1 = pi_-1; p2 = pi_0
%
% The following scaling has to be adopted in case the
% orthogonal polynomials underflow or overflow. In the
% former case, c has to be chosen sufficiently large
% positive, in the latter case sufficiently small
% positive. (The example below relates to Table2_9.m.)
%
%if N==320
%  c=1e50;
%  p2=c*p2;
%  s0=c^2*s0;
%end
%
for k=1:N
	p0 = p1; p1 = p2;								% pi_i -> pi_i-1; pi_i+1 -> pi_i

	p2 = (xw(:,1)-ab(k,1)).*p1-ab(k,2)*p0;			% pi_i+1
	s1 = xw(:,2)'*(p2.^2);							% <pi_i+1,pi_i+1>_mu
	s2 = xw(:,2)'*(xw(:,1).*(p2.^2));				% <x pi_i+1,pi_i+1>_mu

	if(max(abs(p2))>huge)||(abs(s2)>huge)
		error(sprintf('impending overflow in stieltjes for k=%3.0f',k))
	end
	if abs(s1)<tiny
		error(sprintf('impending underflow in stieltjes for k=%3.0f',k))
	end

	ab(k+1,1) = s2/s1;								% alpha_i+1
	ab(k+1,2) = s1/s0;								% beta_i+1 = ||pi_i+1||^2/||pi_i||^2
	s0 = s1;										% ||pi_i+1||^2 -> ||pi_i||^2
end

ab(:,2) = [sqrt(ab(2:end,2));0];
ab = ab(1:end-1,:);
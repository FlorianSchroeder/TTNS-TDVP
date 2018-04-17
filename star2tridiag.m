function [alpha,beta,v]=star2tridiag(indiag,inrow)
% indiag(1) = inrow(1) = 0
% Modified:
%	FS 29/05/2014: replaced first loop to fill A
%	FS 24/09/2015: sparse A to speedup A*v

reps=0.01*min([min(indiag),min(inrow)]);
dim=length(indiag);

%A=zeros(dim);

A = diag(indiag);
A(1,:) = inrow;
A(:,1) = inrow;
A = sparse(A);						% speedup!

alpha=zeros(dim,1);
beta=zeros(dim,1);

v=zeros(dim);
r=zeros(dim,1);
r(1)=1;
beta(1)=norm(r);
for j=1:dim
  % Basic recursion, Equivalent to Orthogonal Polynomials
  v(:,j)=r/beta(j);						% Normalize vec
  r=A*v(:,j);							% Get next vec, r_j+1
  if j>1, r=r-v(:,j-1)*beta(j); end;
  alpha(j)=r'*v(:,j);
  r=r-v(:,j)*alpha(j);
  beta(j+1)=norm(r);
  % Reorthogonalize
  h1=v'*r;								% standard Gram-Schmidt
  r=r-v*h1;
  nh1=norm(h1);
  if nh1>reps,
     r=r-v*(v'*r);
  end
end

% cut first element
alpha=alpha(2:end);
beta=beta(2:end-1);
function [bp,bm,n]=bosonop(dim,shift,parity)
% Generate local boson operators
% shift sets only offset to diagonal entries
% Modified:
%       FS 24/05/2014: replaced slow for-loop by gallery() to create bp
%       FS 25/10/2014: replaced slow gallery by faster direct sparse

% bp = gallery('tridiag',zeros(1,dim-1),zeros(1,dim),sqrt(dim-1:-1:1));
% bp = sparse(1:dim-1, 2:dim, sqrt(dim-1:-1:1), dim, dim);		% n-1:-1:1 inverted
bp = sparse(2:dim, 1:dim-1, sqrt(1:dim-1), dim, dim);			% 1:n-1
if nargout == 1 && shift == 0, return, end

if parity~='n'
    bp=parityorderOP(bp);
end

bm = bp';
n  = bp*bm;

if shift~=0
    bp   = full(bp); bm = full(bm); n = full(n);
    x    = sqrt(2)/2.*(bp+bm);
    iden = eye(dim);
    bp   = bp+sqrt(2)/2.*shift*iden;
    bm   = bm+sqrt(2)/2.*shift*iden;
    n    = n+shift.*x+shift^2/2.*iden;
end

end

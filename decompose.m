function [Q,V,vNE,sv] = decompose(A, i,para)
% decomposes any tensor A_abcd along index i by rearranging index i to the end.
% Uses faster QR if called without vNE and sv in output arguments. Otherwise SVD
%
% USE:  [A,V,vNE] = decompose(A,2)  % with A level 3 tensor
%                   --> permute(A,[3,1,2])
%                   --> reshape
%                   --> [Q,S,V] = svd2(A)
% TODO: - implement odd parity svd!
%       - sorting in QR needed?
%
% added by Florian Schroeder 16/01/2014
%%
%A = randn(3,3);
Asize = size(A);                % get dimensions of tensor A

indA = 1:length(Asize);         % create indices for dimensions.
indA(i) = [];                   % remove index i for decomposition
%%
if isempty(indA) || para.parity ~='n'           % only implemented for parity = 'n'
    return                                      % A was vector
end

A = permute(A,[indA,i]);                        % put dim i at the end
A = reshape(A,[prod(Asize(indA)),Asize(i)]);    % make A 2-dimensional
%%
if nargout > 2 || strcmp(para.SVDmethod,'svd')
    %disp('SVD')
    [Q, S, V] = svd2(A);                        % sorts diag(S) decreasing
    vNE = vonNeumannEntropy(S);
    sv = diag(S);
%    d_opt = size(S,1);                         % might be needed??
    V = S*V;
elseif strcmp(para.SVDmethod,'qr')
    %disp('QR')
    [Q, V] = qr(A,0);                    % no sorting by diagonal values of V
end

Q = reshape(Q, [Asize(indA),Asize(i)]);         % reshape to tensor, works only if dim(x)>dim(y) in Q(x,y)
Q = permute(Q, [1:(i-1),length(Asize), i:(length(Asize)-1)]);

end

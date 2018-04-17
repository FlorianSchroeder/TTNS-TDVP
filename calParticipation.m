function participation = calParticipation(rdm)
%   DEPRECATED: Use getObservable()
%
% rdm: a reduced density matrix of a single site
%   e.g. rdm = calReducedDensity(mps,Vmat,para,k)
%
% takes a rdm and calculates the participation ratio

participation = 1/sum(diag(rdm).^2);

end
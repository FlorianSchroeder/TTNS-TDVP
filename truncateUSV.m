function [U, S, V, newDim, err] = truncateUSV(u,sv,v, para, minDim)
%% [U, S, V, newDim, epsilon] = truncateUSV(U,S,V, para, minDim)
% Truncates Bond Dimensions according to Singular Values
% keeps always one SV in range [para.svmaxtol, para.svmintol]
% keeps always at least minDim dimensions
%
%	u : truncated in columns
%	s : vector of SV
%	v : is V', truncated in rows
%
% by Florian Schroeder 08/08/2015
% Modified
%	FS 28/08/2015: - truncate after sv normalisation!
%	FS 11/01/2016: - return truncation error

if ~isempty(para)
	svmintol = para.svmintol;
	svmaxtol = para.svmaxtol;
else
	svmintol = 10^-4.5;
	svmaxtol = 10^1;
end
if para.tdvp.imagT
	sv = sv./norm(sv);
end
err = 0;
%% Truncate A dims
% fprintf('\n SV norm: %g\n',sum(sv.^2));
keepdims = find(sv > svmintol);
if length(keepdims) < minDim                     % keep at least Dmin bonds
	keepdims = (1:minDim)';
end
% If smallest SV too large (also after normalisation), keep 1 more
while (sv(keepdims(end))/norm(sv(keepdims)) > svmaxtol)
	if keepdims(end)+1 < length(sv)
		keepdims = [keepdims;keepdims(end)+1];
	else
		% expand Bond dims! Happens later?
		U = u; S = sv; V = v; newDim = length(sv);
		return;
	end
end
if length(keepdims) < length(sv)
	U = u(:,keepdims);                          % keep columns
	S = sv(keepdims)./norm(sv(keepdims));		% normalise to keep ||psi|| = 1
	V = v(keepdims,:);                          % keep rows
    newDim = length(keepdims);

	% truncation error estimation
	truncdims = 1:length(sv);
	truncdims(keepdims) = [];
	err = sum(sv(truncdims).^2);				% all disregarded SV^2
else
	U = u; S = sv; V = v; newDim = length(sv);
end

function conjGradient()
	%% how to optimise the truncated result?
end


end

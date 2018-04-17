function [Hleft] = updateHleft(Hleft, h1jOBB, Opstorage1, B, BUb, h2jOBB, A, AUb, M, varargin)
%% function [Hleft] = updateHleft(Hleft, h1jOBB, Opstorage1, B, BUb, h2jOBB, A, AUb, M, varargin)
%
%  Calculates the left effective Hamiltonian before moving to the next site
%  j+1.
%
% Modified
%	FS 17/07/2015: - added support for multi-chain on OBB level

% 	multiChain = 0;
% 	if nargin > 9
% 		para = varargin{1};
% 		if para.nChains > 1
% 			multiChain = 1;						% multiChain at OBB level
% 		end
% 	end

% 	if isempty(Opstorage1{1})
% 		for m=1:M
% 			Opstorage1 = zeros(size(B,1));
% 		end
% 	end
%
% 	if multiChain == 0	&& any(cellfun('isempty',h2jOBB))			% initialize matrices as 0 if empty (only use for single chains), should be fine!
% 		for m = find(cellfun('isempty',h2jOBB))
% 			h2jOBB{m} = zeros(size(BUb,2));							% zeros(para.d_opt(j)) since empty h2j in updateCright -> eye -> wrong! need explicit 0!
% 		end															% multiChain == 1 does not use updateCright -> leave empty
% 	end

	Hleft = updateCleft(Hleft, B, [], [], A, []);                   % removed Vmat since should cancel anyways!

	if ~isempty(h1jOBB)
		Hleft = Hleft + updateCleft([],B,[],h1jOBB,A,[]);
	end

	for m=1:M
		if ~isempty(h2jOBB{m}) && ~isempty(Opstorage1{m})
			Hleft = Hleft + updateCleft(Opstorage1{m},B,[],h2jOBB{m},A,[]);
		end
	end


end
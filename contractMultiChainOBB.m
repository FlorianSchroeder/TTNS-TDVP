function X = contractMultiChainOBB(Vmat, McOp, para)
%% function X = contractMultiChainOBB(Vmat, McOp, para)
% Contracts a single chain operator in mcOp with Vmat
%
%  Vmat : either Vmat with MC at OBB level, or V-tensor-network
%  McOp : 1 x NC cell with operator for subchain nc at McOp{nc}
%         empty cells indicate eye()

NC = para.nChains;
Empty  = cellfun('isempty',McOp);                    % logical vector with pos of empty cells

if NC == 1 && length(McOp) == 1 % standard OBB trafo
	X = Vmat' * McOp{1} * Vmat;
	return;
end

% Do not accept all-empty mcOp
if Empty                                             % if all are empty
% 	fprintf('ContractMultiChainOBB: Empty mcOp in loop %d at site %d', para.loop, para.sitej);
	X = [];
	return
end

if ~iscell(Vmat) && ismatrix(Vmat)
	% reshape into tensor
	Vmat = reshape(Vmat,[para.dk(:,para.sitej)',size(Vmat,2)]);
end

if iscell(Vmat)
	% Vtens code, put mcOp into OBB described by Vmat{1:end-1}
	if sum(~Empty) > 1, error('Only input 1 operator per contraction'); end;
	Ind       = find(~Empty);                        % get position, has to be single number!
	McOp{Ind} = Vmat{Ind}' * McOp{Ind} * Vmat{Ind};
	Vmat = Vmat{end};                                % take the VS out for compatibility with next step
end

% Old routine, inefficient
% % contract all empty parts
% Ind	 = find(Empty);                                                 % finds all empty indices
% Ind  = reshape(Ind, [], numel(Ind));                                % need row vector!
% X	 = contracttensors(conj(Vmat), NC+1, Ind, Vmat, NC+1, Ind);		% dk' x d_opt' x dk x d_opt
% X	 = contracttensors(X,4,[1,3],McOp{~Empty},2,[1,2]);

Ind  = find(~Empty);		% get index of chain
ord = 1:NC; ord(Ind) = [];
X = contracttensors(full(McOp{~Empty}),2,2, Vmat, NC+1, Ind);
X = contracttensors(conj(Vmat),NC+1, [Ind, ord] ,X, NC+1, 1:NC);

end
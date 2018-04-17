function op = h1h2toOBB(V, para, op)
%% DEPRECATED, use [op] = H_Eff([], Vmat, 'A', op, para) instead
%%	function op = h1h2toOBB(V, para, op)
% transform all bare H terms of current sitej into OBB
% works with multi-chain at OBB level
%
% h1term to OBB, h1j can be rescaled
% h1j_(n~',n~) = V*_(n',n~') [h1j_(n',n) V_(n,n~)]_(n',n~)

% h2term to OBB, h2j can be rescaled
% h2j_(n~',n~) = V*_(n',n~') [h2j_(n',n) V_(n,n~)]_(n',n~)

if para.nChains == 1
	op.h1j = V' * (op.h1j * V);									% faster and more accurate
	for i=1:para.M
		op.h2j{i,1} = V' * (op.h2j{i,1} * V);									% faster and more accurate
		op.h2j{i,2} = V' * (op.h2j{i,2} * V);									% faster and more accurate
	end
else  % nChains > 1
	h1jnew = 0;
	for i = find(~cellfun('isempty',op.h1j'))
		H1 = cell(para.nChains,1);
		H1(i) = op.h1j(i);
		h1jnew = h1jnew + contractMultiChainOBB(V, H1, para);
	end
	op.h1j = h1jnew;

	h2jnew = cell(para.M,2);
	for i = 1:para.M
		h2jnew{i,1} = contractMultiChainOBB(V, op.h2j(i,1,:), para);
		h2jnew{i,2} = contractMultiChainOBB(V, op.h2j(i,2,:), para);
	end
	op.h2j = h2jnew;
end

end
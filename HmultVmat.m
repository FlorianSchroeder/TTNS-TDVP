function outvec = HmultVmat(invec,op,dkj,d_opt,M,parity)
% invec = V = optimising Vmat OBB matrix.
%
% Modified:
%	- 20/12/14 FS: faster matrix product for 2-2 contractions
%	- 04/08/15 FS: replaced HrOPB writing by more intuitive HrightA
%					they are all of form: op.HrightA(n~',n~)
%					for consistency
if parity~='n'
    V=zeros(dkj,d_opt);
    V(1:dkj/2,1:d_opt/2)=reshape(invec(1:dkj*d_opt/4),dkj/2,d_opt/2);
    V(dkj/2+1:end,d_opt/2+1:end)=reshape(invec(dkj*d_opt/4+1:end),dkj/2,d_opt/2);
else
    V=reshape(invec,[dkj,d_opt]);
end

out = V * (op.HleftA.' + op.HrightA.');

if ~iscell(op.h1j)
	% original single-chain code
	if ~isempty(op.h1j)
		out = out + (op.h1j * V);
	end
	for m = 1:M
		if isempty(op.h2j{m,2}) && isempty(op.h2j{m,1}), continue; end;
		outSl = op.h2j{m,2} * (V * op.OpleftA{m}.');
		outSr = op.h2j{m,1} * (V * op.OprightA{m}.');

		out=out+outSr+outSl;
	end
else	% multi-chain!
	dkSC = cell2mat(cellfun(@(x) size(x,1),op.h1j', 'UniformOutput',false));	% dk of each subchain
	nc = length(dkSC);															% # subchains

	for i = 1:nc
		% reshape inB to contract with h1j{i} of subChain i
		outTemp1 = reshape(V,[dkSC,d_opt]);
		outTemp1 = contracttensors(op.h1j{i},2,2,outTemp1,nc+1,i);				% equiv to: out = out + (op.h1j * inB);
		if i ~= 1
			outTemp1 = permute(outTemp1,[2:i,1,(i+1):(nc+1)]);
		end
		out = out + reshape(outTemp1,[prod(dkSC),d_opt]);

		for m = 1:M
			if isempty(op.h2j{m,2,i})			% m according to other subchain
				continue;						% go to next m
			end
			outTemp1 = V * op.OpleftA{m}.';
			outTemp1 = reshape(outTemp1,[dkSC,d_opt]);
			outTemp1 = contracttensors(op.h2j{m,2,i},2,2,outTemp1,nc+1,i);
			if i ~= 1
				outTemp1 = permute(outTemp1,[2:i,1,(i+1):(nc+1)]);
			end
			out = out + reshape(outTemp1,[prod(dkSC),d_opt]);					% equiv to: outSl = op.h2j{m,2} * (inB * op.OpleftOPB{m}.');

 			outTemp1 = V * op.OprightA{m}.';
			outTemp1 = reshape(outTemp1,[dkSC,d_opt]);
			outTemp1 = contracttensors(op.h2j{m,1,i},2,2,outTemp1,nc+1,i);
			if i ~= 1
				outTemp1 = permute(outTemp1,[2:i,1,(i+1):(nc+1)]);
			end
			out = out + reshape(outTemp1,[prod(dkSC),d_opt]);					% equiv to: outSr = op.h2j{m,1} * (inB * op.OprightOPB{m}.');
		end
	end
end



if parity~='n'
    outvec=vertcat(reshape(out(1:dkj/2,1:d_opt/2),dkj*d_opt/4,1),reshape(out(dkj/2+1:end,d_opt/2+1:end),dkj*d_opt/4,1));
else
    outvec=reshape(out,[numel(out),1]);
end

end
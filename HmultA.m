function outvec=HmultA(invec, op, Dl, Dr, d, M,parity,nonzeroindex)
% calculates H|MPS> for energy minimization.
% invec & outvec are vectorized forms of the optimized MPS-matrix A
%
% Modified
%	- FS 23/11/14: changed last reshape to numel() since faster
%	- FS 05/08/15: changed from h1j, h2j -> h1jOBB, h2jOBB
if parity~='n'
    inA=zeros(Dl,Dr,d);
    inA(nonzeroindex)=invec;
else
    inA=reshape(invec,[Dl,Dr,d]);
end


out  = 0;		% output of form Dl x Dr x d
out2 = 0;		% output of form (Dl * d) x Dr

% collect non-interacting parts of Hamiltonian:
out   = reshape(op.Hleft*reshape(inA,Dl,[]),[Dl,Dr,d]);
Atemp = reshape(permute(inA,[1,3,2]),[],Dr);					% (Dl * d) x Dr
out2  = Atemp * op.Hright.';

	% from now: contract local operators first -> bring inA into suitable form
inA = reshape(inA,[],d);										% (Dl * Dr) x d

if ~isempty(op.h1jOBB)
    out = out + reshape(inA * op.h1jOBB.', [Dl,Dr,d]);
end

% collect interacting parts:
for m=1:M
	if ~isempty(op.h2jOBB{m,2}) && ~isempty(op.Opleft{m})
		Atemp = inA * op.h2jOBB{m,2}.';
		Atemp = op.Opleft{m} * reshape(Atemp,Dl,[]);
		out   = out + reshape(Atemp, [Dl,Dr,d]);
	end
	if ~isempty(op.h2jOBB{m,1}) && ~isempty(op.Opright{m})
		Atemp = inA * op.h2jOBB{m,1}.';
		Atemp = reshape(permute(reshape(Atemp, [Dl,Dr,d]),[1,3,2]),[],Dr);
		out2  = out2 + Atemp * op.Opright{m}.';
	end
end

out = out + permute(reshape(out2,[Dl,d,Dr]),[1,3,2]);


if parity~='n'
    outvec=out(nonzeroindex);
else
    outvec = reshape(out,[numel(out),1]);
end

end
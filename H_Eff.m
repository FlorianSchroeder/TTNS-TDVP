function [op] = H_Eff(mps, V, target, op, para)
%% [op] = H_Eff(mps, Vtarget, op, para)
%	creates the effective operators for 'target' to allow fast contractions
%	in matrix exponentials and eigs()
%
%	mps		: single-site A matrix, mps{sitej}
%	Vmat	: single-site V matrix, Vmat{sitej}
% e.g.
%	[op] = H_Eff(Amat, []  , 'V' , op, para);
%	[op] = H_Eff([]  , Vmat, 'A' , op, para);	% transforms op.h1j, h2j to OBB
%	[op] = H_Eff(Amat, []  , 'CA', op, para);	% relies on previous OBB trafo
%
%   Only MC-TDVP code (HOSVD): 'MC-OBB', 'MC-CV', 'MC-V', 'MC-VS', 'MC-A'
%
%   [op] = H_Eff([]  ,Vtens, 'MC-OBB', op, para);  % transforms all chain h1j/h2j into MCOBB, calls MC-CV
%   [op] = H_Eff([]  ,Vtens, 'MC-CV' , op, para);  % transforms single chain into OBB; defined by para.currentChain
%   [op] = H_Eff([]  ,Vtens, 'MC-V'  , op, para);  % create HnonInt and Hleft/rightAS for single-chain V evolution
%   [op] = H_Eff([]  ,[]   , 'MC-VS' , op, para);  % create h12jAV terms for VS evolution
%   [op] = H_Eff([]  ,Vtens, 'MC-A'  , op, para);  % create h1/2jOBB terms, depends on MC-OBB
%
%	Only in Star-MPS code:
%
%	[op] = H_Eff(mps{1}, [], 'ST-CA', op, para);   % transforms System + all other chains into Hleft, Opleft defined by para.currentChain
%
%	Only in Tree-MPS code:
%
%	[op] = H_Eff(treeMPS, [], 'TR-A' , [], para);		% transforms node h1h/h2j into OBB; for leaf/chain use standard 'A'
%	[op] = H_Eff(treeMPS, [], 'TR-CA', [], para);   % transforms Node + all other chains into H/Opleft defined by para.currentChain
%	[op] = H_Eff(treeMPS, [], 'TR-CP', [], para);   % transforms Node + all chains into H/Opright for parent node
%
%	Created 04/08/2015 by FS
M = para.M;
switch target
	case 'V'
		%% multiply MPS into op.Hright, Hleft, Opright, Opleft
		% HleftA_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
		op.HleftA = contracttensors(op.Hleft,2,2,mps,3,1);
		op.HleftA = contracttensors(conj(mps),3,[1 2], op.HleftA,3,[1 2]);
		% Hright_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
		op.HrightA = contracttensors(op.Hright,2,2,mps,3,2);
		op.HrightA = contracttensors(conj(mps),3,[1 2], op.HrightA,3,[2 1]);

		op.OpleftA  = cell(M,1);
		op.OprightA = cell(M,1);

		for k=1:M
			% same contractions as above for Hleft/Hright
			if isempty(op.Opleft{k})
				op.OpleftA{k} = [];
			else
				op.OpleftA{k}= contracttensors(op.Opleft{k}, 2,2, mps,3,1);
				op.OpleftA{k}= contracttensors(conj(mps),3,[1 2], op.OpleftA{k},3,[1 2]);
			end

			if isempty(op.Opleft{k})
				op.OpleftA{k} = [];
			else
				op.OprightA{k} = contracttensors(op.Opright{k},2,2, mps,3,2);
				op.OprightA{k} = contracttensors(conj(mps),3,[1 2], op.OprightA{k},3,[2 1]);
			end
		end
	case 'A'
		%% multiply Vmat into op.h1j, h2j
		% transform all bare H terms of current sitej into OBB
		% works with multi-chain at OBB level
		if ~isfield(op,'h1j') || isempty(op.h1j)
			op = gen_sitej_h1h2(op,para,para.sitej);
		end
		if para.nChains == 1 || (~iscell(V) && ~iscell(op.h1j))
			op.h1jOBB = V' * (op.h1j * V);									% faster and more accurate
			op.h2jOBB = cell(M,2);
			for i=1:M
				op.h2jOBB{i,1} = V' * (op.h2j{i,1} * V);									% faster and more accurate
				op.h2jOBB{i,2} = V' * (op.h2j{i,2} * V);									% faster and more accurate
			end
		else  % nChains > 1 && iscell(V)
			h1jnew = 0;
			for i = find(~cellfun('isempty',op.h1j'))
				H1 = cell(para.nChains,1);
				H1(i) = op.h1j(i);
				h1jnew = h1jnew + contractMultiChainOBB(V, H1, para);
			end
			op.h1jOBB = h1jnew;

			op.h2jOBB = cell(M,2);
			for i = 1:M
				op.h2jOBB{i,1} = contractMultiChainOBB(V, op.h2j(i,1,:), para);
				op.h2jOBB{i,2} = contractMultiChainOBB(V, op.h2j(i,2,:), para);
			end
		end
	case 'CA'
		%% contract op.Hleft,Hright,Opleft,Opright with op.h1jOBB,h2jOBB and MPS matrix
		% save into op.HleftAV,HrightAV,OprightAV,h2jAV
		switch para.sweepto
			case 'r'	% MPS and V contracted into left operators only
				% HleftAV_(r~',r~) = A*_(l',r~',n~) [Hl_(l',l) * A_(l,r~,n~)]_(l',r~,n~)
				op.HleftAV = updateCleft(op.Hleft, mps, [], [], mps, []);

				% h1jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~), h1j_(n~',n~)]_(l,r~,n~')
				op.HleftAV = op.HleftAV + updateCleft([],mps,[],op.h1jOBB,mps,[]);

				% op.Opleft will be summed over and added to HleftAV, since
				% it does not interact across CA
				op.h2jAV	= cell(M,1);
				for k = 1:M
					% OpleftAV_(r~',r~) = A*_(l',r~',n~')[[OPl_(l',l) A_(l,r~,n~)]_(l',r~,n~) h2j_(n~',n~)]_(l',r~,n~')
					op.HleftAV  = op.HleftAV + updateCleft(op.Opleft{k}, mps, [], op.h2jOBB{k,2}, mps, []);

					% h2jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~) h2j_(n~',n~)]_(l,r~,n~')
					op.h2jAV{k} = updateCleft([], mps, [], op.h2jOBB{k,1}, mps, []);
				end
			case 'l'	% MPS and V contracted into right operators only
				% HHrightAV_(l~',l~) = A*_(l~',r',n~) [Hr_(r',r) * A_(l~,r,n~)]_(r',l~,n~)
				op.HrightAV = updateCright(op.Hright, mps, [], [], mps, []);

				% h1jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~), h1j_(n~',n~)]_(l~,r,n~')
				op.HrightAV = op.HrightAV + updateCright([],mps,[],op.h1jOBB,mps,[]);

% 				op.OprightAV = cell(M,1);		% add to op.HrightAV directly
				op.h2jAV	 = cell(M,1);
				for k = 1:M
					% OprightAV_(l~',l~) = A*_(l~',r',n~')[[A_(l~,r,n~) OPr_(r',r)]_(l~,n~,r') h2j_(n~',n~)]_(l~,r',n~')
					op.HrightAV = op.HrightAV + updateCright(op.Opright{k},mps,[],op.h2jOBB{k,1},mps,[]);

					% h2jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~) h2j_(n~',n~)]_(l~,r,n~')
					op.h2jAV{k} = updateCright([],mps,[],op.h2jOBB{k,2},mps,[]);
				end
		end
	case 'MC-V'
		%% contract Tensors of Multi-Chain HOSVD-TDVP for evolution of V{k}
		% builds upon 'V' preparation. -> [op] = H_Eff(Amat, []  , 'V' , op, para);
		% V has to be cell with V{end} = VS the tensor-center, unpermuted: n1 x n2 x ... x nNC x n~
		% Saves in [H/op][left/right]AS
		% leave out ~ for clarity
		NC = para.nChains;
		nc = para.currentChain;                                                                        % # chain I am working on
		nTerms = para.M/NC;                                                                            % # interacting terms per chain
		% 1. contract VS into 4th-order tensor
		IndContract = 1:NC; IndContract(nc) = [];
		VSfullContract = contracttensors(conj(V{end}),NC+1, IndContract, V{end}, NC+1, IndContract);   % VS_(nk', n', nk, n)

		op.HleftAS  = contracttensors(op.HleftA, 2,[1,2], VSfullContract,4,[2,4]);                     % H_(nk',nk)
		op.HrightAS = contracttensors(op.HrightA,2,[1,2], VSfullContract,4,[2,4]);

		% this could collect all terms non-interacting with current chain!
		op.HnonInt  = op.HleftAS + op.HrightAS;
		% h1j terms of other chains
		for mc = 1:NC
			if mc == nc || isempty(op.h1jMCOBB{mc}), continue; end;
			IndContract = 1:NC+1; IndContract([mc,nc]) = [];
			VS = contracttensors(conj(V{end}),NC+1, IndContract, V{end}, NC+1, IndContract);
			if mc < nc          % -> VS_(nmc',nk',nmc,nk)
				order = [1,3];
			else                % -> VS_(nk',nmc',nk,nmc)
				order = [2,4];
			end
			op.HnonInt = op.HnonInt + contracttensors(op.h1jMCOBB{mc},2,[1,2], VS,4,order);
		end

		for m = 1:para.M
			mc = ceil(m/nTerms);                        % the chain number for current m
			if nc == mc
				op.OpleftAS{m}  = contracttensors(op.OpleftA{m},  2,[1,2], VSfullContract,4,[2,4]);
				op.OprightAS{m} = contracttensors(op.OprightA{m}, 2,[1,2], VSfullContract,4,[2,4]);
			else % sum into op.HnonInt
				if isempty(op.h2jMCOBB{m,2,mc}) && isempty(op.h2jMCOBB{m,1,mc}), continue; end;
				if NC > 2
					% better contraction scheme?:
					[idxA, idxB] = getIdxTensChain(NC+1,mc,nc);
% 					idxA = 1:NC+1; idxB = idxA;
% 					idxA([mc,nc]) = [];	idxA = [idxA,mc];	% Do not contract nc
% 					if mc < nc
% 						idxB(nc-1) = [];					% since nc shifted
% 					else
% 						idxB(nc)   = [];
% 					end

					OpTemp = contracttensors(V{end}, NC+1, NC+1, op.OpleftA{m}.',2,1);						%_(ni..,nk)
					OpTemp = contracttensors(OpTemp, NC+1, mc, op.h2jMCOBB{m,2,mc}.',2,1);					%_(ni,...,nk,mc)
					op.HnonInt = op.HnonInt + contracttensors(conj(V{end}),NC+1, idxA, OpTemp, NC+1, idxB);

					OpTemp = contracttensors(V{end}, NC+1, NC+1, op.OprightA{m}.',2,1);						% _(ni..,nk)
					OpTemp = contracttensors(OpTemp, NC+1, mc, op.h2jMCOBB{m,1,mc}.',2,1);					% _(nc,nk,mc)
					op.HnonInt = op.HnonInt + contracttensors(conj(V{end}),NC+1, idxA, OpTemp, NC+1, idxB);
				else % need more efficient code! NC+1 = 3
					assert(NC+1 == 3, 'this code only works under this assumption!');
					OpTemp = contracttensors(V{end}, NC+1, NC+1, op.OpleftA{m}.',2,1);						% _(n1,n2,nk)
					OpTemp = contracttensors(OpTemp, NC+1, mc, op.h2jMCOBB{m,2,mc}.',2,1);					% _(nc,nk,mc)
					op.HnonInt = op.HnonInt + contracttensors(conj(V{end}),NC+1, [mc,3], OpTemp, NC+1, [3,2]);

					OpTemp = contracttensors(V{end}, NC+1, NC+1, op.OprightA{m}.',2,1);						% _(n1,n2,nk)
					OpTemp = contracttensors(OpTemp, NC+1, mc, op.h2jMCOBB{m,1,mc}.',2,1);					% _(nc,nk,mc)
					op.HnonInt = op.HnonInt + contracttensors(conj(V{end}),NC+1, [mc,3], OpTemp, NC+1, [3,2]);
				end
			end
		end
		% finish with terms: op.HnonInt, op.OpleftAS{m}, op.OprightAS{m}
	case 'MC-CV'
		%% contract to get each single chain h1j/h2j into its OBB.
		% equivalent to creating the H_Eff for VS
		if isempty(op.h1j)
			op = gen_sitej_h1h2(op,para,para.sitej);
		end

		nTerms = para.M/para.nChains;
		nc     = para.currentChain;
		M      = 1:para.M;
		M(ceil(M/nTerms) ~=nc) = [];                         % find which m to operate on
		if ~isempty(op.h1j{nc})
			op.h1jMCOBB{nc} = V{nc}' * (op.h1j{nc} * V{nc});
		else
			op.h1jMCOBB{nc} = [];
		end
		for m = M
			if ~isempty(op.h2j{m,1,nc})
				op.h2jMCOBB{m,1,nc} = V{nc}' * (op.h2j{m,1,nc} * V{nc});
				op.h2jMCOBB{m,2,nc} = V{nc}' * (op.h2j{m,2,nc} * V{nc});
			else
				op.h2jMCOBB{m,1,nc} = [];
				op.h2jMCOBB{m,2,nc} = [];
			end
		end
	case 'MC-OBB'
		%% Create all OBB terms for each chain
		for i = 1:para.nChains
			para.currentChain = i;
			[op] = H_Eff([]  , V, 'MC-CV', op, para);
		end
	case 'MC-VS'
		%% Create the kronecker product terms for faster VS exponential
		d = para.d_opt(:,para.sitej);
		nTerms = para.M/para.nChains;
		for k = 1:para.nChains
			M      = 1:para.M;
			M(ceil(M/nTerms) ~=k) = [];                         % find which m to operate on
			op.h12jAV{k} = kron(eye(d(end)), op.h1jMCOBB{k});
			for m = M
				op.h12jAV{k} = op.h12jAV{k} + kron(op.OpleftA{m} , op.h2jMCOBB{m,2,k});
				op.h12jAV{k} = op.h12jAV{k} + kron(op.OprightA{m}, op.h2jMCOBB{m,1,k});
			end
		end
	case 'MC-A'
		%% Create the OBB terms for A, gathering all chains
		% needs all h1jMCOBB to be up-to-date (MC-OBB)
		% use contractMultiChainOBB with objects from here! Copied from 'A'
		h1jnew = 0;
		for i = find(~cellfun('isempty',op.h1jMCOBB))
			H1 = cell(para.nChains,1);
			H1(i) = op.h1jMCOBB(i);
			h1jnew = h1jnew + contractMultiChainOBB(V{end}, H1, para);
		end
		op.h1jOBB = h1jnew;

		op.h2jOBB = cell(M,2);
		for i = 1:M
			op.h2jOBB{i,1} = contractMultiChainOBB(V{end}, op.h2jMCOBB(i,1,:), para);
			op.h2jOBB{i,2} = contractMultiChainOBB(V{end}, op.h2jMCOBB(i,2,:), para);
		end

	case 'ST-CA'
		%% Create the effective Chain-System terms for STAR-MPS
		% Store in op.chain(para.currentChain).H/Opleft
		% Only used once for entering a chain
		% similar to MC-V in structure
		nc = para.currentChain;
		NC = para.nChains;
		M = para.M / NC;			% what is M for each single chain?
		d  = size(mps);
		% MPS is only mps{1}: 1 x D1 x D2 x D3 x ... X D(NC) x dk
		%   chain index in mps is shifted by 1 due to first singleton
		%   dimension!! -> Compatibility with 1-chain models and future Boundary Conditions?

		op.chain(nc).Hleft = 0;						% contains effective H of all other chains
		op.chain(nc).Opleft = cell(para.M/NC,1);	% contains only terms interacting with Chain #nc

		for mc = 1:NC
			if mc == nc, continue; end;
			% always make sure that Hlrstorage{1} and Opstorage{:,2,1} are up-to-date!
			if para.tdvp.HEffSplitIsometry == 0
				% 1. Contract to non-interacting parts -> Hleft
				[idxA, idxB] = getIdxTensChain(NC+2,mc+1,nc+1);
				OpTemp = contracttensors(mps, NC+2, mc+1, op.chain(mc).Hlrstorage{1}.',2,1);						%_(ni..,nk)
				op.chain(nc).Hleft = op.chain(nc).Hleft + contracttensors(conj(mps),NC+2, idxA, OpTemp, NC+2, idxB);

				% 2. Contract the other chain-system interacting parts -> Hleft
				for m = 1:M
					% this chain's m
					systemM = M*(mc-1) + m;
					OpTemp = contracttensors(mps   , NC+2, NC+2, op.h2jOBB{systemM,1}.'   ,2,1);
					OpTemp = contracttensors(OpTemp, NC+2, mc+1, op.chain(mc).Opstorage{m,2,1}.',2,1);					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
					op.chain(nc).Hleft = op.chain(nc).Hleft + contracttensors(conj(mps),NC+2, idxA, OpTemp, NC+2, idxB);
				end
			else
				% 1. Split-off isometry
				%    but need D(mc), D(nc) and dk in Atens!
% 				para.currentChain = mc;
				[Atens,dOut] = tensShape(mps, 'unfoldiso', [mc+1,nc+1,NC+2], d);	% +1 for leading singleton

				[Iso, Atens] = qr(Atens,0);			% Iso is isometry with all unused chains

				% TODO: comment if not testing!!
% 				[Iso2,A2,err]= rrQR(Atens, floor(0.1*prod(dOut(end-2:end))),0);		% low rank QR approximation
% 				fprintf('H_eff: rank(A) = %d, dim(A,2) = %d, rrQR error = %g; %g\n', rank(Atens), prod(dOut(end-2:end)),norm(Atens - Iso2*A2), err);


				Atens = reshape(Atens,[size(Atens,1),dOut(end-2:end)]);		% D(mc)*D(nc)*dk x D(mc) x D(nc) x dk
				idxA = [1, 2, 4]; idxB = [1, 4, 3];
				% 2. Contract to non-interacting parts -> Hleft
				OpTemp = contracttensors(Atens, 4, 2, op.chain(mc).Hlrstorage{1}.',2,1);						%_(ni..,nk)
				op.chain(nc).Hleft = op.chain(nc).Hleft + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);

				% 3. Contract the other chain-system interacting parts -> Hleft
				for m = 1:M
					% this chain's m
					systemM = M*(mc-1) + m;
					OpTemp = contracttensors(Atens , 4, 4, op.h2jOBB{systemM,1}.'         ,2,1);
					OpTemp = contracttensors(OpTemp, 4, 2, op.chain(mc).Opstorage{m,2,1}.',2,1);					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
					op.chain(nc).Hleft = op.chain(nc).Hleft + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);
				end
			end
		end

		% 3. Contract h1jOBB -> Hleft
		idx = 1:NC+2; idx(nc+1) = [];
		OpTemp = contracttensors(mps, NC+2, NC+2, op.h1jOBB.', 2, 1);
		op.chain(nc).Hleft = op.chain(nc).Hleft + contracttensors(conj(mps), NC+2, idx, OpTemp,NC+2, idx);

		% 4. Contract h2jOBB which will interact with chain nc
		for m = 1:M
			systemM = M*(nc-1) + m;
			OpTemp = contracttensors(mps, NC+2, NC+2, op.h2jOBB{systemM,1}.', 2, 1);
			op.chain(nc).Opleft{m} = contracttensors(conj(mps), NC+2, idx, OpTemp,NC+2, idx);
		end
		% overwrite main vars for evolve Kn
% 		op.Hleft   = op.chain(nc).Hleft;
% 		op.Opleft  = op.chain(nc).Opleft;
		op.Hright  = op.chain(nc).Hlrstorage{1};
		op.Opright = op.chain(nc).Opstorage(:,2,1);
	
	case 'TR-A'
		%% multiply Vmat into h1j/h2j for node
		% transform all bare H terms of current node into OBB
		% only for treeMPS
		
		op = H_Eff_TR_A(mps);
	case 'TR-CA'
		%% Create the effective Hleft terms for the children of TREE-MPS
		% Store in treeMPS.child(para.currentChain).op.H/Opleft
		%	since Hlrstorage{1} still contains Hright for CA
		% Only used once for entering a child
		% similar to MC-V and ST-CA in structure
		%
		%	mps:  treeMPS - contains node-specific para
		%	V:    []
		%	op:   []
		%	para: para - contains calculation-specific para
		
		% do call to subfunction to translate variable names -> more consistent
		
		op = H_Eff_TR_CA(mps,para);
		
	case 'TR-CP'
		%% Create the effective Hright terms for the parent of TREE-MPS
		%	TRee - CenterParent
		% Store in treeMPS.op.H/Opright
		%	since Hlrstorage{1} still has Hleft for CA
		% similar to MC-V and ST-CA in structure
		% 
		%	mps:  treeMPS - contains node-specific para
		%	V:    []
		%	op:   []
		%	para: para - contains calculation-specific para
		if mps.hasSite
			op = H_Eff_TR_CP_Site(mps,para);
		else
			op = H_Eff_TR_CP_NoSite(mps,para);			% special version to treat chain-combination nodes
		end
	
end
end

function op = H_Eff_TR_A(treeMPS)
%% function op = H_Eff_TR_A(treeMPS)
%
%	fetches h1j/h2j if empty and transforms into OBB
%	Only for internal nodes of tree

op = treeMPS.op;
if ~isfield(op,'h1j') || isempty(op.h1j)
	op = gen_sitej_h1h2(op,treeMPS);
end

V = treeMPS.Vmat{1};

op.h1jOBB = V' * (op.h1j * V);

op.h2jOBB = cell(size(op.h2j));
for i = 1:size(op.h2j,1)
	for j = 1:size(op.h2j,3)
		op.h2jOBB{i,1,j} = V' * (op.h2j{i,1,j} * V);				% faster than cellfun and more accurate
		op.h2jOBB{i,2,j} = V' * (op.h2j{i,2,j} * V);
	end
end

end

function op = H_Eff_TR_CA(treeMPS,para)
%% function op = H_Eff_TR_CA(treeMPS,para)
%
%	Creates Hleft for entering a child/chain = walking down the tree
%	Child has no site and is an entanglement renormalisation tensor!
%
%	created by FS 23/02/2016

nc = treeMPS.currentChild;				% = child to prepare operators for
NC = treeMPS.degree;					% = nChild for node
M  = treeMPS.M;
d  = size(treeMPS.mps{1});				% numel(d) = NC+1 or NC+2
Nd = numel(d);
mps = treeMPS.mps{1};					% get handle, should not copy object!
% MPS is only mps{1}: Dl x D1 x D2 x D3 x ... X D(NC) x dk
%   chain index in mps is shifted by 1 due to first left Bond
for i = 1:length(treeMPS.child)
	NCoup(i) = length(treeMPS.child(i).chainIdx);																% much faster than  arrayfun()
end
op = treeMPS.child(nc).op;				% only return op -> needs overwriting in function call
op.Hleft  = 0;							% contains effective H of all other chains
op.Opleft = cell(M,NCoup(nc));			% contains only terms interacting with Chain #nc

% 1. construct non-interacting Hleft for child
for mc = 0:NC
	if treeMPS.hasSite
		dProd = prod(d([mc+1,nc+1,NC+2]));		% mc interacts with local site
	else
		dProd = prod(d([mc+1,nc+1,1]));			% mc interacts with parent
	end
	% mc: other chain to contract with MPS to get into left basis of nc
	if mc == nc, continue, end
	% always make sure that Hlrstorage{1} and Opstorage{:,2,1} are up-to-date!
	if para.tdvp.HEffSplitIsometry == 0 || (prod(d)/dProd^2) < 1.3
		if mc == 0 && treeMPS.isRoot
			continue
		elseif mc == 0 && ~treeMPS.isRoot				% Node-parent exists -> contract
			idx = 1:Nd; idx(nc+1) = [];
			% 1. Contract to non-interacting parts
			w   = contracttensors(treeMPS.op.Hlrstorage{1},2,2, mps, Nd, 1);									% order unchanged, Focused on Node -> Hlrstorage{1} is Hleft of Node
			
			if treeMPS.hasSite
				% 2. Contract the parent-local site interacting parts
				for m = 1:M
					OpTemp   = contracttensors(mps, Nd, Nd, treeMPS.op.h2jOBB{m,2,1}.',2,1);					% parent is 1st, node is 2nd position in H_int !!
					w   = w + contracttensors(treeMPS.op.Opstorage{m,1,1},2,2, OpTemp, Nd, 1);					% (m,1,1) should be the operator of parent in the effective left basis for node site 1
				end
			end
			op.Hleft = op.Hleft + contracttensors(conj(mps),Nd, idx, w, Nd, idx);
			continue
		end
		
		% 1. Contract to non-interacting parts -> Hleft
		w   = contracttensors(mps, Nd, mc+1, treeMPS.child(mc).op.Hlrstorage{1}.',2,1);					%_(ni..,nk)

		% 2. Contract the other chain-system interacting parts -> Hleft
		% NCoup(mc) indicates how many terms are needed: 1 for leaf/chain; more for ER-tensor
		NCoupOffset = sum(NCoup(1:mc-1));																		% how many terms belong to other children -> Offset
		if treeMPS.hasSite
			for nn = 1:NCoup(mc)
				for m = 1:M
					OpTemp   = contracttensors(mps   , Nd, Nd, treeMPS.op.h2jOBB{m,1,nn+NCoupOffset}.'    ,2,1);
					w   = w + contracttensors(OpTemp, Nd, mc+1, treeMPS.child(mc).op.Opstorage{m,2,1,nn}.',2,1);		% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				end
			end
		else % ~hasSite
			for nn = 1:NCoup(mc)
				for m = 1:M
					OpTemp   = contracttensors(treeMPS.op.Opstorage{m,1,1, nn+NCoupOffset},2,2, mps   , Nd, 1);			% Take interaction from parent
					w   = w + contracttensors(OpTemp, Nd, mc+1, treeMPS.child(mc).op.Opstorage{m,2,1,nn}.',2,1);		% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				end
			end
		end
		[idxA, idxB] = getIdxTensChain(Nd,mc+1,nc+1);		
		% contraction of mps with chain ops of mc moved its index to the end 
		% -> generate idx pair for <mps|H|MPS>, where everything except nc shall be contracted
		op.Hleft = op.Hleft + contracttensors(conj(mps),Nd, idxA, w, Nd, idxB);
		
	else
		if mc == 0 && treeMPS.isRoot
			continue
		end
		% 1. Split-off isometry
		%    but need D(mc), D(nc) and dk in Atens!
		[Atens,dOut] = tensShape(mps, 'unfoldiso', [mc+1,nc+1,NC+2], d);	% +1 for leading singleton
		[~, Atens] = qr(Atens,0);			% Iso is isometry with all unused chains
		% TODO: comment if not testing!!
% 		rEst = rank(Atens,10^-4.5);
% 		fprintf('H_eff: rank(A) = %d, dim(A,2) = %d\n',rEst,prod(dOut(end-2:end)));
% 		[Q, ~]     = qr(Atens,0);			% Iso is isometry with all unused chains
% 		A2 = Q'*Atens;
% 		[Iso2,A3,err]= rrQR(Atens, rEst,0);		% low rank QR approximation
% 				fprintf('H_eff: rank(A) = %d, dim(A,2) = %d, rrQR error = %g; %g\n', rank(Atens), prod(dOut(end-2:end)),norm(Atens - Iso2*A2), err);
% 		Atens = A1;

		if mc == 0		% TODO: if these assignments take too long, then remove this shortcut!
			Hstor  = treeMPS.op.Hlrstorage{1};					% Hleft of node
			h2jOBB = treeMPS.op.h2jOBB(:,2,1);
			Opstor = treeMPS.op.Opstorage(:,1,1);				% Opleft of node
		else
			Hstor  = treeMPS.child(mc).op.Hlrstorage{1};		% Hright of child
			h2jOBB = treeMPS.op.h2jOBB(:,1,mc);
			Opstor = treeMPS.child(mc).op.Opstorage(:,2,1);		% Opleft of child
		end

		
		Atens = reshape(Atens,[size(Atens,1),dOut(end-2:end)]);		% D(mc)*D(nc)*dk x D(mc) x D(nc) x dk
		idxA = [1, 2, 4]; idxB = [1, 4, 3];
		% 2. Contract to non-interacting parts -> Hleft
		OpTemp   = contracttensors(Atens, 4, 2, Hstor.',2,1);							%_(ni..,nk)
		op.Hleft = op.Hleft + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);
		
		% 3. Contract the other chain-system interacting parts -> Hleft
		for m = 1:M
			OpTemp   = contracttensors(Atens , 4, 4, h2jOBB{m}.', 2, 1);
			OpTemp   = contracttensors(OpTemp, 4, 2, Opstor{m}.', 2, 1);					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
			op.Hleft = op.Hleft + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);
		end
	end
end

idx = 1:Nd; idx(nc+1) = [];
NCoupOffset = sum(NCoup(1:nc-1));																		% how many terms belong to other children -> Offset
if treeMPS.hasSite
	% 3. Contract h1jOBB -> Hleft
	OpTemp = contracttensors(mps, Nd, Nd, treeMPS.op.h1jOBB.', 2, 1);
	op.Hleft = op.Hleft + contracttensors(conj(mps), Nd, idx, OpTemp,Nd, idx);

	% 4. Contract h2jOBB which will interact with chain nc
	for nn = 1:NCoup(nc)
		for m = 1:M
			OpTemp = contracttensors(mps, Nd, Nd, treeMPS.op.h2jOBB{m,1,nn+NCoupOffset}.', 2, 1);
			op.Opleft{m,nn} = contracttensors(conj(mps), Nd, idx, OpTemp,Nd, idx);
		end
	end
else % ~hasSite: carry over parent interactions
	for nn = 1:NCoup(nc)
		for m = 1:M
			OpTemp = contracttensors(treeMPS.op.Opstorage{m,1,1,nn+NCoupOffset}, 2, 2, mps, Nd, 1);
			op.Opleft{m,nn} = contracttensors(conj(mps), Nd, idx, OpTemp,Nd, idx);
		end
	end
end

op.Hright  = treeMPS.child(nc).op.Hlrstorage{1};		% copy into temporary H/Opright to allow overwriting in updateop.
op.Opright = squeeze(treeMPS.child(nc).op.Opstorage(:,2,1,:));		% (M,NCoup)

end

function op = H_Eff_TR_CP_Site(treeMPS,para)
%% function op = H_Eff_TR_CP_Site(treeMPS,para)
%
%	Creates Hright for entering the parent = walking up the tree
%	Applies if current node carries a physical site which couples to the parent and which all children are coupled to
%
%	created by FS 23/02/2016

NChild = treeMPS.degree;		% = nChildren for node
NChain = treeMPS.nChains;		% = nChains of node
M  = treeMPS.M;
d  = size(treeMPS.mps{1});		% numel(d) = NC+2 if hasSite; elseif ~hasSite: NC+1
mps = treeMPS.mps{1};			% get handle, should not copy object!
% MPS is only mps{1}: Dl x D1 x D2 x D3 x ... X D(NC) [x dk ]
%   chain index in mps is shifted by 1 due to first left Bond
op = treeMPS.op;				% only return op -> needs overwriting in function call
op.Hright  = 0;					% contains effective H of node and all chains
op.Opright = cell(M,1);			% contains only terms interacting with parent

if treeMPS.isRoot
	error('VMPS:H_Eff:NotSupported','This function shall only be used for non-root nodes');
end

% 1. construct non-interacting Hright for parent
for mc = 1:NChild
	% mc: chain to contract with MPS to get into right basis of parent
	
	% always make sure that Hlrstorage{1} and Opstorage{:,2,1} are up-to-date!
	if para.tdvp.HEffSplitIsometry == 0
		[idxA, idxB] = getIdxTensChain(NChild+2,mc+1,1);		
		% contraction of mps with chain ops of mc will move its index to the end 
		% -> generate idx pair for <mps|H|MPS>, where everything except 1 shall be contracted
		
		% 1. Contract to non-interacting parts -> Hright
		OpTemp   = contracttensors(mps, NChild+2, mc+1, treeMPS.child(mc).op.Hlrstorage{1}.',2,1);
		op.Hright = op.Hright + contracttensors(conj(mps),NChild+2, idxA, OpTemp, NChild+2, idxB);

		% 2. Contract the other chain-system interacting parts -> Hright
		for m = 1:M
			OpTemp    = contracttensors(mps   , NChild+2, NChild+2, treeMPS.op.h2jOBB{m,1,mc}.'            ,2,1);
			OpTemp    = contracttensors(OpTemp, NChild+2, mc+1, treeMPS.child(mc).op.Opstorage{m,2,1}.',2,1);		% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
			op.Hright = op.Hright + contracttensors(conj(mps),NChild+2, idxA, OpTemp, NChild+2, idxB);
		end
	else
		% 1. Split-off isometry
		%    but need D(mc), D(nc) and dk in Atens!
		[Atens,dOut] = tensShape(mps, 'unfoldiso', [mc+1,1,NChild+2], d);	% +1 for leading singleton
		
		[~, Atens] = qr(Atens,0);			% Iso is isometry with all unused chains
		
		% TODO: comment if not testing!!
% 				[Iso2,A2,err]= rrQR(Atens, floor(0.1*prod(dOut(end-2:end))),0);		% low rank QR approximation
% 				fprintf('H_eff: rank(A) = %d, dim(A,2) = %d, rrQR error = %g; %g\n', rank(Atens), prod(dOut(end-2:end)),norm(Atens - Iso2*A2), err);
		
		Hstor  = treeMPS.child(mc).op.Hlrstorage{1};		% Hright of child
		h2jOBB = treeMPS.op.h2jOBB(:,1,mc);
		Opstor = treeMPS.child(mc).op.Opstorage(:,2,1);		% Opleft of child
		
		Atens = reshape(Atens,[size(Atens,1),dOut(end-2:end)]);		% D(mc)*D(nc)*dk x D(mc) x D(nc) x dk
		idxA = [1, 2, 4]; idxB = [1, 4, 3];
		% 2. Contract to non-interacting parts -> Hright
		OpTemp    = contracttensors(Atens, 4, 2, Hstor.',2,1);
		op.Hright = op.Hright + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);
		
		% 3. Contract the other chain-system interacting parts -> Hright
		for m = 1:M
			OpTemp    = contracttensors(Atens , 4, 4, h2jOBB{m}.', 2, 1);
			OpTemp    = contracttensors(OpTemp, 4, 2, Opstor{m}.', 2, 1);					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
			op.Hright = op.Hright + contracttensors(conj(Atens),4, idxA, OpTemp, 4, idxB);
		end
	end
end

% 3. Contract h1jOBB -> Hright
idx = 2:NChild+2;
OpTemp    = contracttensors(mps, NChild+2, NChild+2, treeMPS.op.h1jOBB.', 2, 1);
op.Hright = op.Hright + contracttensors(conj(mps), NChild+2, idx, OpTemp,NChild+2, idx);

% 4. Contract h2jOBB which will interact with parent
for m = 1:M
	OpTemp = contracttensors(mps, NChild+2, NChild+2, op.h2jOBB{m,2,1}.', 2, 1);
	op.Opright{m} = contracttensors(conj(mps), NChild+2, idx, OpTemp,NChild+2, idx);
end

op.Hleft  = treeMPS.op.Hlrstorage{1};		% copy into temporary H/Opleft to allow overwriting in updateop.
op.Opleft = treeMPS.op.Opstorage(:,1,1);
end

function op = H_Eff_TR_CP_NoSite(treeMPS,para)
%% function op = H_Eff_TR_CP_NoSite(treeMPS,para)
%
%	Creates Hright for entering the parent = walking up the tree
%	Applies if current node carries NO physical site, but serves to combine several nodes/leaves for entanglement renormalisation
%	All children of this node couple to the parent, thus need to channel through treeMPS.child(:).op.Opstorage{:,2,1,:}
%
%	created by FS 30/08/2016

NChild = treeMPS.degree;					% = nChildren for node
NCoup = length(treeMPS.chainIdx);			% = number of coupling channels from children to parent node
%	if child node hasSite -> chainIdx = min(child.chainIdx)
%	elseif child ~hasSite -> chainIdx = [child.chainIdx]
%	thus correctly counts the number of couplings through the ER-Node
M  = treeMPS.M;
d  = size(treeMPS.mps{1});					% numel(d) = NC+2 if hasSite; elseif ~hasSite: NC+1
mps = treeMPS.mps{1};						% get handle, should not copy object!
% MPS is only mps{1}: Dl x D1 x D2 x D3 x ... X D(NC) [x dk ]
%   chain index in mps is shifted by 1 due to first left Bond
op = treeMPS.op;							% only return op -> needs overwriting in function call
op.Hright  = 0;								% contains effective H of node and all chains
op.Opright = cell(M,NCoup);					% if no Site, then all chains from right must be channelled through to the left, thus need all coupling terms

NCoupOffset = 0;							% necessary to correctly & easily fill op.Opright

if treeMPS.isRoot
	error('VMPS:H_Eff:NotSupported','This function shall only be used for non-root nodes');
end

for mc = 1:NChild
	% mc: chain to contract with MPS to get into right basis of parent
	
	HStor   = treeMPS.child(mc).op.Hlrstorage{1};			% Hright of child
	OpStor  = treeMPS.child(mc).op.Opstorage(:,2,1,:);		% Opright of child: (M,NCoup)
	[MM,NN] = size(OpStor);
	ratio   = (d(1)*d(mc+1))^2/prod(d);						% d(needed legs)/d(rest); Benefit from Isometry only if << 1
	
	% always make sure that Hlrstorage{1} and Opstorage{:,2,1} are up-to-date!
	if para.tdvp.HEffSplitIsometry == 0 || ratio >= 1					% TODO: need to find optimal threshold
		[idxA, idxB] = getIdxTensChain(NChild+1,mc+1,1);		
		% contraction of mps with chain ops of mc will move its index to the end 
		% -> generate idx pair for <mps|H|MPS>, where everything except 1 shall be contracted
		
		% 1. Contract to non-interacting parts -> Hright
		OpTemp   = contracttensors(mps, NChild+1, mc+1, HStor.',2,1);
		op.Hright = op.Hright + contracttensors(conj(mps),NChild+1, idxA, OpTemp, NChild+1, idxB);

		% 2. Contract the other parts for the chain-parent interacting parts -> Opright
		for n = 1:NN
			for m = 1:MM
				OpTemp   = contracttensors(mps, NChild+1, mc+1, OpStor{m,n}.',2,1);
				op.Opright{m,n+NCoupOffset} = contracttensors(conj(mps),NChild+1, idxA, OpTemp, NChild+1, idxB);
			end
		end
		NCoupOffset = NCoupOffset + NN;
	else
		% 1. Split-off isometry
		%    but need D(mc), D(nc=1) in Atens! (no dk, since no Site)
		[Atens,dOut] = tensShape(mps, 'unfoldiso', [1, mc+1], d);	% +1 for leading singleton
		
		[~, Atens] = qr(Atens,0);			% [Iso,R] = qr(A); Iso is isometry with all unused chains -> contracts to Iso'*Iso = eye
		
		% TODO: comment if not testing!!
 		%		[Iso2,A2,err]= rrQR(Atens, floor(0.1*prod(dOut(end-2:end))),0);		% low rank QR approximation
 		%		fprintf('H_eff: rank(A) = %d, dim(A,2) = %d, rrQR error = %g; %g\n', rank(Atens), prod(dOut(end-2:end)),norm(Atens - Iso2*A2), err);
		
		Atens = reshape(Atens,[size(Atens,1),dOut(end-1:end)]);		% D(1)*D(mc) x D(1) x D(mc)
		idxA = [1, 3]; idxB = [1, 3];
		% 2. Contract to non-interacting parts -> Hright
		OpTemp    = contracttensors(Atens, 3, 3, HStor.',2,1);
		op.Hright = op.Hright + contracttensors(conj(Atens),3, idxA, OpTemp, 3, idxB);
		
		% 3. Contract the other parts for the chain-parent interacting parts -> Opright
		for n = 1:NN
			for m = 1:MM
				OpTemp   = contracttensors(Atens, 3, 3, OpStor{m,n}.',2,1);
				op.Opright{m,n+NCoupOffset} = contracttensors(conj(Atens),3, idxA, OpTemp, 3, idxB);
			end
		end
		NCoupOffset = NCoupOffset + NN;
	end
end
assert(size(op.Opright,2) == NCoup, 'VMPS:H_Eff:H_EFF_TR_CP_NoSite:WrongNCOUP','The number of couplings was wrong, Opright grew! Please check chainIdx for each node for correctness!');

op.Hleft  = treeMPS.op.Hlrstorage{1};		% copy into temporary H/Opleft to allow overwriting in updateop.
op.Opleft = treeMPS.op.Opstorage(:,1,1,:);
end

function [idxA, idxB] = getIdxTensChain(N,m,n)
	% generate two arrays with indices corresponding to the desired contractions in the case:
	%	all bonds ~= n will be contracted over
	%	m : 1 position which got moved to the end in Tensor 1 (e.g. by contraction)
	%	n : 1 position which should be the output bonds
	%	
	idxA = 1:N; idxB = idxA;
	idxA([m,n]) = [];	idxA = [idxA,m];

	if m < n
		idxB(n-1) = [];
	elseif m > n
		idxB(n)   = [];
	else
		error('This is only meant to be used for m ~= n');
	end

end
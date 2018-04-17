function [op] = initstorage(mps, Vmat,op,para)
% op.Hlrstorage stores the sub hamiltonian of the (left?) and right block.
% op.Opstorage stores the operators in the (left?) and right bases.
% sweeps from right to left.
% Only use with rightnormalised MPS

% Commented by Florian Schroeder 29/01/2014
% Changed:
%   FS 20/10/2014: - use op = updateop() instead of explicit calculations
%                    to remove redundancy
%                  - moved separate treatment of L into for-loop

if (isfield(para,'useTreeMPS') && para.useTreeMPS)
	op = initstorage_Tree_MPS(mps,para);		% recursive sub-function call
	return;
end

L             = length(mps);
M             = size(op.h2term,1); 		% number of terms in sum per site
op.Hlrstorage = cell(1, L + 1);
op.Opstorage  = cell(M, 2, L+1);

% treat first and last terms seperately
op.Hlrstorage{1, 1} 	= 0;
op.Hlrstorage{1, L + 1} = 0;

% moved to updateop():
% op.Hlrstorage{L} = updateHright(op.Hlrstorage{L + 1}, op.h1term{L}, op.Opstorage(:,2,L+1), mps{L}, Vmat{L}, op.h2term(:,1,L), mps{L}, Vmat{L},M);
	% gives h1term{L} in effective basis of r_L-1 as:
	% op.Hlrstorage{l+1} = 0; op.h1term{L} = exists ; op.Opstorage(:,2,L+1) = [], op.h2term(:,1,L) = 0;
for m=1:M
	% moved to updateop():
%     op.Opstorage{m,2,L}=updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});
	% transforms interaction terms into r_L-1 ef. basis
	op.Opstorage{m,1,1}	  =0;
	op.Opstorage{m,2,L+1} =0;
end

para.sweepto = 'l';
op.h1j    = []; op.h2j    = [];
op.h1jOBB = []; op.h2jOBB = [];


% middle terms build: l<-r
for j = L:-1:2
	para.sitej = j;
    op = updateop(op,mps,Vmat,j,para);
	op.h1jOBB = []; op.h2jOBB = [];
	op.h1j = []; op.h2j = [];
end

end

function treeMPS = initstorage_Tree_MPS(treeMPS,para)
%% function treeMPS = initstorage_Tree_MPS(treeMPS,para)
%
%	recursively computes the effective Hamiltonian terms

p = para;					% create working copy of para, to temporarily modify the .use* parameters
p.sweepto = 'l';
pIdx = num2cell(treeMPS.treeIdx+1);

if treeMPS.height == 0
	% this is leaf -> call prepare(mps,Vmat,para)
	%	or normalise here in place!
	initstorageChain();			% mutating subfunction
else
	for k = 1:treeMPS.degree
		% recursive calls
		treeMPS.child(k) = initstorage_Tree_MPS(treeMPS.child(k),para);
	end
	initstorageNode();
end
	

	function initstorageChain()
		%% function initstorageChain()
		%	subfunction to create effective Hamiltonians of single chainMPS
		%
		%	Modify treeMPS only
		
		% necessary fields in para:
		% p.L = treeMPS.L;
		% p.nChains = 1;
		% p.useStarMPS = 0;
		% p.useTreeMPS = 0;
		treeMPS.sweepto = 'l';
		
		% use treeMPS instead of para to make different chain topologies possible. e.g. mixed Vtens - Vmat
		treeMPS.op = initstorage(treeMPS.mps,treeMPS.Vmat,treeMPS.op,treeMPS);		% site 2:L
		treeMPS.sitej = 1;
		% call syntax for chain
		treeMPS.op = updateop(treeMPS.op,treeMPS.mps,treeMPS.Vmat,1,treeMPS);		% treat Hright from site 1 to parent node separately
	end

	function initstorageNode()
		%% function initstorageNode()
		%	subfunction to calculate the effective right Hamiltonian for the parent node
		%
		%	Modify treeMPS only
		
		treeMPS.sweepto = 'l';
		treeMPS.sitej = 1;
		
		L = treeMPS.L(1);								% should be 1
		M = treeMPS.M;									% number of terms in sum per site
		N = length(treeMPS.chainIdx);					% number of chains to couple to parent
		
		% following will only store Hright for the parent node or Hleft of/from the parent node
		% Hright of the children for this node and Hleft of this node for the children will be stored in
		%	treeMPS.child(ii).op.Hlrstorage
		% No need for op.chain{}.*storage variables (replaced op.chain{} by treeMPS.child(ii).op)
		treeMPS.op.Hlrstorage = cell(1, L);			% No +1 since right part is inside children.
		treeMPS.op.Opstorage  = cell(M, 2, L, N);
		treeMPS.op.h1j    = []; treeMPS.op.h2j    = [];
		treeMPS.op.h1jOBB = []; treeMPS.op.h2jOBB = [];

		% call syntax for node
		treeMPS.op = updateop([],treeMPS,[],[],para);
		
		if treeMPS.isRoot
			treeMPS.op.Hlrstorage{1} = 0;
			for m = 1:M
				treeMPS.op.Opstorage{m,1,1} = 0;
				% treeMPS.op.Opstorage{m,2,1} = 0; stays empty since multiple children to the right!
			end
		end
	end
end

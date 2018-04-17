function tdvp_1site_tree(treeMPS,para,results,tresults)
%% VMPS Time-dependent variational principle
%   implementation following Haegeman et al. 2014, arXiv:1408.5056
%   using Expokit for expv()
%
%   Expands OBB on sitej upon entering mps{sitej}, before SVD to Vmat
%   Truncates OBB on sitej after focusing the next mps site
%
%
%   To continue a completed calculation use:
%       inarg: tmps instead of mps, tVmat instead of Vmat. reorganization
%              will be done in code under 0.
%
%	To resume an incomplete calculation use:
%		inarg: tmps(end)-> mps, tVmat(end) -> Vmat
%
%   Created by Florian Schroeder 07/10/2014
%
%	Modified 22/02/2016
%		- Restructured Variables to support TreeMPS
%
para.tdvp.starttime = tic;

if nargin == 0
	createTestRun();
end

%% 0. Parameter settings:
if ~isfield(para,'tdvp')
    error('VMPS:tdvp_1site_tree:MissingParameters','Please define all necessary tdvp variables in para!');
end

% open file for writing results
% if storeMPS then para.tdvp.filename will contain tMPS
assert(isfield(para.tdvp,'filename'),'Need definition of para.tdvp.filename to save results!');
assert(isfield(para.tdvp,'filenameSmall'),'Need definition of para.tdvp.filenameSmall to save results!');

if isfield(para.tdvp,'storeMPS') && para.tdvp.storeMPS == 1
	outFile = matfile(para.tdvp.filenameMPS,'Writable',true);
	% this does not yet create file.
	% File will be created upon first writing!
end

resumeTDVP = 0;					% 0 = start new, 1 = resume aborted, 2 = continue completed

if ~isempty(tresults) && para.tdvp.resume == 1 && exist(para.tdvp.filename,'file')
	%% must be resumed since no writing has happend since start!
	% detect startTime:
	%  - outFile exists, and is already synced with tresults from resume
	%  - only have tresults
% 	try
% 		load(para.tdvp.filename);	% all the rest
% 	catch
% 		disp('MAT-file corrupt, load BAK');
% 		load([para.tdvp.filename(1:end-4),'.bak'],'-mat');		% use backup if file corrupted
% 	end
	if isfield(tresults,'lastIdx')
		% there must be preexisting calculation -> resume!
		startSize = para.timeslice+1;
		resumeTDVP = 1;
	else
		startSize = 1;				% fallback value. Perhaps change?
		resumeTDVP = 0;
		clear('tresults');			% throw away old results -> Overwrite!
	end
elseif ~isempty(tresults)
	% NEEDS WORK!! TODO
	% detect continuation
	if isempty(strfind(para.tdvp.fromFilename,'results.mat'))
		resumeTDVP = 2;
		%load(para.tdvp.fromFilename, 'tresults');							% if error, load .bak file!
		% tresults will be expanded in calTimeObservables!

% 		if ~exist('tresults','var')
% 			error('tresults were not in fromFile!');
% 		end
		startSize = para.timeslice+1;
		fprintf('\nSet resumeTDVP = 2, startSize = %d', startSize);
	end
end

%% 1. sweep: condition the ground state MPS and prepare hamiltonian operator terms H(n) and K(n)
    % is already done in op.opstorage and op.hlrstorage
    %
    % also prepare other stuff
if resumeTDVP == 0
	if exist('outFile','var')
		% initialize tMPS and calc observables at the end
		fprintf('Initialize tMPS\n');
		error('VMPS:tdvp_1site_tree:NotImplemented','This needs implementation.')
% 		outFile.tmps(1,para.L) = {[]};
% 		outFile.tVmat(1,para.L) = {[]};
% 		outFile.tmps(1,:) = mps;
% 		outFile.tVmat(1,:) = Vmat;
	end
	% calc tresults now to initialize tresults
	tresults = calTimeObservables(treeMPS,[],para);				% TODO: this needs overloading..
	if para.logging
% 		results = initresultsTDVP(para, results);				% TODO: this needs Work
	end
	para.trustsite(end) = para.L;       % necessary to enable several truncation / expansion mechanisms
	% For fresh starts, define number of timeslices
	para.tdvp.slices = 1:(length(para.tdvp.t)-1);
	para.tdvp.calcTime = zeros(length(para.tdvp.t)-1,1);		% logs the elapsed time after each slice in h
elseif resumeTDVP == 1
	fprintf('Resuming aborted TDVP\n')
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);
	if ~isfield(results,'tdvp')
		% results was not stored properly before -> reinitialize
		if para.logging
% 			results = initresultsTDVP(para, results);			% TODO!!!
			if para.tdvp.expandOBB
% 				results.tdvp.d_opt(startSize,:)   = para.d_opt;		% TODO!!
			end
			if para.tdvp.truncateExpandBonds
% 				results.tdvp.D(startSize,:)       = para.D;		% TODO!!
			end
		end
	end
elseif resumeTDVP == 2
	%% If continue previously completed TDVP
	fprintf('Continuing previous TDVP\n');
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);
	if exist('outFile','var');								% TODO!!
% 		outFile.tmps = mps;
% 		outFile.tVmat = Vmat;
% 		outFile.tmps(length(para.tdvp.t),para.L) = {[]};
% 		outFile.tVmat(length(para.tdvp.t),para.L) = {[]};
	end
	% increase size of results objects:
	if para.logging
		n = length(para.tdvp.t);
		if para.tdvp.expandOBB
% 			results.tdvp.d_opt(n,end)   = 0;						% TODO!!
		end
		if para.tdvp.truncateExpandBonds
% 			results.tdvp.D(n,end)       = 0;						% TODO!!
		end
% 		results.tdvp.expError(n-1,1)    = 0;						% TODO!!
% 		results.tdvp.EvaluesLog(n-1,1)  = 0;						% TODO!!
		if para.tdvp.logSV
% 			results.tdvp.Vmat_sv{n,end} = cell(n,para.L);						% TODO!!
% 			results.tdvp.Amat_sv{n,end} = cell(n,para.L-1);						% TODO!!
		else
% 			results.tdvp.Vmat_vNE(n,end)= 0;						% TODO!!
% 			results.tdvp.Amat_vNE(n,end)= 0;						% TODO!!
		end
	end
end

fprintf('Saving.');
save(para.tdvp.filenameSmall,'para','tresults','-v7.3');
if para.tdvp.serialize
	hlp_save(para.tdvp.filename,para,[],results,tresults,[],[],treeMPS);
else
	save(para.tdvp.filename,'para','results','tresults','treeMPS','-v7.3');
	fprintf('.\n');
end
%% Asserts to do:

%% 2. time sweep

for timeslice = para.tdvp.slices
    para.sweepto = 'r'; treeMPS.sweepto = 'r';
	para.timeslice = timeslice;
% 	para.tdvp.deltaT = para.tdvp.t(timeslice+1)-para.tdvp.t(timeslice);
    fprintf('t = %g\n', para.tdvp.t(timeslice+1));
    % sweep l->r and time evolve each site
	
	%% check whether treeMPS is Tree, Star or Chain
	%	Pro-Recursive: No need for call stack, automatic assignment after calculation
	%	Con-Recursive: Function call overhead. copy on write
	%	Pro-Iterative: No Function call overhead
	%	Con-Iterative: Call Stack needed. Saving to parent problematic, as not doable by reference. 
	
	%% Start recursive TDVP
	if treeMPS.height > 0
		% treeMPS is Tree (or Star)
		treeMPS = tdvp_1site_evolveTree(treeMPS,para);
	else
		% treeMPS is a Leaf/Chain
		treeMPS = tdvp_1site_evolveChain(treeMPS,para);
	end
	
	%% Sweep finished. Logging output:
	
	% Output matrix dimensions if changed:
	if para.tdvp.truncateExpandBonds
		fprintf('\npara.D:\n');
		out = D2Str(treeMPS);
		fprintf('%s\n',out{:});
		fprintf('-log10(minSV):\n');
		out = minSV2Str(treeMPS);
		fprintf('%s\n',out{:});
	end
	fprintf('\n')
	
	% report time estimates
%	results.tdvp.EvaluesLog(timeslice) = getObservable({'energy',op},mps,Vmat,para);
	completePercent = round(timeslice./(length(para.tdvp.t)-1)*1000)./10;
	para.tdvp.calcTime(timeslice) = toc(para.tdvp.starttime)./3600;
	n = 1:timeslice;
	n = n(para.tdvp.calcTime(n,1)>0);		% compensates for zeros due to resume of aborted TDVP
	hoursTotal = fnval(fnxtr(csaps([0,n],[0;para.tdvp.calcTime(n,1)])), length(para.tdvp.t)-1);
	if log10(abs(hoursTotal-para.tdvp.calcTime(timeslice))) <= 0
		fprintf('Completed: %5.3g%%, t = %4.3g, Time elapsed: %4.2gh, Time left: %4.2gm, Time total: %4.2gh\n',...
			completePercent, para.tdvp.t(timeslice+1), para.tdvp.calcTime(timeslice),...
			(hoursTotal-para.tdvp.calcTime(timeslice))*60, hoursTotal);
	else
		fprintf('Completed: %5.3g%%, t = %4.3g, Time elapsed: %4.2gh, Time left: %4.2gh, Time total: %4.2gh\n',...
			completePercent, para.tdvp.t(timeslice+1), para.tdvp.calcTime(timeslice),...
			hoursTotal-para.tdvp.calcTime(timeslice), hoursTotal);
	end
    %% save tmps and tVmat and log parameters
	if exist('outFile','var')
		fprintf('Saving tMPS\n');
		outFile.tmps(timeslice+1,:)  = mps;					% writes to File
		outFile.tVmat(timeslice+1,:) = Vmat;				% writes to File
		outFile.currentSize = timeslice+1;
	end
	if para.logging && mod(timeslice, round(para.tdvp.extractObsInterval/para.tdvp.deltaT))==0			% at extractObsInterval
		n = 1+tresults.lastIdx;
		if para.tdvp.logSV
			if size(results.Vmat_sv,1) == para.L
				results.tdvp.Vmat_sv(n,:,:)	 = results.Vmat_sv;
			else
				results.tdvp.Vmat_sv(n,:,:)	 = results.Vmat_sv.';
			end
			results.tdvp.Amat_sv(n,:) = results.Amat_sv;
		else
			results.tdvp.Vmat_vNE(n,:,:)= results.Vmat_vNE.';
			results.tdvp.Amat_vNE(n,:,:)= results.Amat_vNE.';
		end
		if para.tdvp.expandOBB
	        results.tdvp.d_opt(n,:)   = reshape(para.d_opt.',1,[]) - sum(results.tdvp.d_opt);	% sum in 1st Dim
		end
		if para.tdvp.truncateExpandBonds
	        results.tdvp.D(n,:)       = reshape(para.D.',1,[]) - sum(results.tdvp.D);
		end
		if para.tdvp.useDkExpand
			results.tdvp.dk(n,:)      = reshape(para.dk.',1,[]) - sum(results.tdvp.dk);
		end
	end
	if ~exist('tresults','var')
		% calculate everything up to now!
		if exist('outFile','var')
			tresults = calTimeObservables(outFile.tmps(1:outFile.currentSize,:),outFile.tVmat(1:outFile.currentSize,:),para);
		else
			% Should not happen!
			error('VMPS:tdvp_1site_tree:Error','How did you get here??');
			%tresults = calTimeObservables(mps,Vmat,para);
		end
	elseif isfield(tresults,'lastIdx') % && tresults.lastIdx == timeslice not applicable anymore for extractObsInterval
		% only calculate the current slice
		tresults = calTimeObservables(treeMPS,[],para,tresults);
	end
	
	
	%% Additional saving
	if mod(timeslice,para.tdvp.saveInterval) == 0 || timeslice ==  para.tdvp.slices(end)
		%% backup results, op and para less often than tmps and tVmat
		results.tdvp.time = toc(para.tdvp.starttime);
		fprintf('Saving results.. ')
		copyfile(para.tdvp.filename,[para.tdvp.filename(1:end-4),'.bak']);			% in case saving gets interrupted
		save(para.tdvp.filenameSmall,'para','tresults','-v7.3');
		if para.tdvp.serialize
			hlp_save(para.tdvp.filename,para,[],results,tresults,[],[],treeMPS);
		else
			save(para.tdvp.filename,'para','results','tresults','treeMPS','-v7.3');
			fprintf('.\n');
		end
	end
end

delete([para.tdvp.filename(1:end-4),'.bak']);			% get rid of bak file

	function createTestRun()
		%% function createTestRun()
		%	
		%	creates a testcase for debugging and development. Will be called if no arguments are given.
		%	
		Lchains = 5;
		D = 10;			% one BondDim for all edges
		d_opt = 10;
		dk = 20;
		
		% Define tree structure in para:
		para.logging = 0;
		para.L = 7;								% total max height of tree, including length of leaves
		para.useTreeMPS = 1;
		para.useStarMPS = 0;
		para.useVmat = 1;
		para.useVtens = 0;
		para.trustsite = para.L;
		para.model = 'testTree';
		para.nChains = 4;						% = external nodes = leaves
% 		para.nEnvironments = para.nChains;		% perhaps needed for calTimeObservables
		para.alpha = 0.1;
		para.hx = 0.1;
		para.hz = 0;
		para.M = 2;
		para.spinbase = 'Z';
		para.parity = 'n';
		para.svmaxtol = 10^-4;					% keep 1 below this!
		para.svmintol = 10^-4.5;				% throw away all below
		para.timeslice = 0;

% 		para.rescaling = 0;
		para.treeMPS.nNodes = 2;								% = internal nodes = nodes
		para.treeMPS.height = 2;								% maximum height of tree
		para.treeMPS.maxDegree = [3,2]; % size(treeMPS.L);		% for each level
		para.treeMPS.leafIdx = num2cell(ones(para.nChains,para.treeMPS.height));				% indices of leaves in para; take -1 to get treeIdx
		para.treeMPS.leafIdx{1,1} = 1+1;														% map from chain number to leaf index in para
		para.treeMPS.leafIdx{2,1} = 2+1;
		para.treeMPS.leafIdx(3,1:2) = num2cell([3,1]+1);
		para.treeMPS.leafIdx(4,1:2) = num2cell([3,2]+1);
		para.treeMPS.nodeIdx = num2cell(ones(para.treeMPS.nNodes,para.treeMPS.height));			% maps node number -> nodeIdx
		para.treeMPS.nodeIdx{1,1} = 0+1;
		para.treeMPS.nodeIdx{2,1} = 3+1;
		para.treeMPS.chainIdx = {};																% maps from leafIdx -> chain number
		for ii = 1:size(para.treeMPS.leafIdx,1)
			para.treeMPS.chainIdx{para.treeMPS.leafIdx{ii,:}} = ii;
		end
		para.chain{1}.alpha   = 0.1;
		para.chain{1}.L       = Lchains;
		para.chain{1}.s       = 1;
		para.chain{1}.epsilon = ones(Lchains,1)*0.5;
		para.chain{1}.t       = ones(Lchains,1)*0.25;
		for ii = 1:para.nChains
			para.chain{ii} = para.chain{1};
		end
		
		% Define each struct of chainMPS and initialise random MPS
		structChainMPS(para.nChains) = struct();		% smallest subunit
		for ii = 1:para.nChains
			Lc = para.chain{ii}.L;
			chainMPS = cell(1,Lc);
			chainVmat = cell(1,Lc);
			for jj = 1:Lc
				if jj == Lc
					chainMPS{jj} = randn(D,1,d_opt);			% end-of-chain -> Dr = 1
				else
					chainMPS{jj} = randn(D,D,d_opt);
				end
				chainVmat{jj} = sparse(1:d_opt,1:d_opt,1,dk,d_opt);
			end
			% Build Struct Chain MPS
			structChainMPS(ii).mps = chainMPS;
			structChainMPS(ii).Vmat = chainVmat;
			structChainMPS(ii).treeIdx = [para.treeMPS.leafIdx{ii,:}]-1;
			structChainMPS(ii).chainIdx = ii;
			structChainMPS(ii).height = 0;									% this struct is a leaf
			structChainMPS(ii).degree = 0;									% leaf
			structChainMPS(ii).level = nnz(structChainMPS(ii).treeIdx)+1;	% = site number of first chain site if counted from the root node
			structChainMPS(ii).child = [];									% contains all children nodes
			structChainMPS(ii).currentChild = [];							% leaf: [], node: currentChain;
			structChainMPS(ii).currentChain = [];							% leaf: [], node: currentChain;
			structChainMPS(ii).L = getTreeLength(structChainMPS(ii));		% StarMPS can have array L
			structChainMPS(ii).D = [ones(1,Lc)*D,1];						% first is Dl of chain, last is Dr of chain, here horizontal!
			structChainMPS(ii).d_opt = ones(1,Lc)*d_opt;
			structChainMPS(ii).dk = ones(1,Lc)*dk;
			structChainMPS(ii).M = para.M;
			structChainMPS(ii).shift = zeros(1,Lc);
			structChainMPS(ii).useVmat = 1;
			structChainMPS(ii).useVtens = 0;
			structChainMPS(ii).useStarMPS = 0;
			structChainMPS(ii).useTreeMPS = 0;
			structChainMPS(ii).nChains = 1;
			structChainMPS(ii).sweepto = 'l';					% 'l' for initial preparation
			structChainMPS(ii).spinposition = [];				% pure Boson Chain
			structChainMPS(ii).sitej = 1;						% currently focused site
			structChainMPS(ii).isRoot = 0;
			structChainMPS(ii).BondCenter = [];
			structChainMPS(ii).tdvp = [];
			pIdx = para.treeMPS.leafIdx(ii,:);
			para.dk{pIdx{:}} = structChainMPS(ii).dk;
			para.D{pIdx{:}} = structChainMPS(ii).D;
			para.d_opt{pIdx{:}} = structChainMPS(ii).d_opt;
			para.shift{pIdx{:}} = structChainMPS(ii).shift;
			para.treeMPS.L(pIdx{:}) = structChainMPS(ii).L;
			
			% generate Hamiltonian terms:
			structChainMPS(ii).op.h1term = cell(1,Lc);
			structChainMPS(ii).op.h2term = cell(para.M,2,Lc);
			for jj = 1:Lc
				temp = genh1h2term_onesite(para,structChainMPS(ii).treeIdx,jj);
				structChainMPS(ii).op.h1term(jj)     = temp.h1term;
				structChainMPS(ii).op.h2term(:,:,jj) = temp.h2term;
			end
		end
		
		% Construct tree node 2
		ii = 2;
		pIdx = para.treeMPS.nodeIdx(ii,:);
		dk = 2; d_opt = 2;
		structStarMPS = structChainMPS(1);							% make a copy to inherit fields. Only modify difference
		structStarMPS.treeIdx = [para.treeMPS.nodeIdx{ii,:}]-1;
		structStarMPS.level = nnz(structStarMPS.treeIdx)+1;			% = site number of first chain site if counted from the root node
		structStarMPS.mps = {randn(D,D,D,d_opt)};					% Dl, Dc1, Dc2, dOBB
		structStarMPS.Vmat = {sparse(1:d_opt,1:d_opt,1,dk,d_opt)};
		structStarMPS.child = [structChainMPS(3);structChainMPS(4)];
		structStarMPS.currentChild = 1;
		structStarMPS.currentChain = 1;
		structStarMPS.chainIdx = min([structStarMPS.child.chainIdx]);
		structStarMPS.height = max([structStarMPS.child.height])+1;
		structStarMPS.degree = 2;
		structStarMPS.nChains = 2;									% equals degree for nodes
		structStarMPS.L = [];
		structStarMPS.L = getTreeLength(structStarMPS);				% StarMPS can have array L
		structStarMPS.d_opt = d_opt;
		structStarMPS.dk = dk;
		structStarMPS.D = ones(3,1)*D;								% [Dl;Dc1;Dc2]  here vertical!
		structStarMPS.shift = 0;
		structStarMPS.useStarMPS = 1;								% star is subset of tree!
		structStarMPS.useTreeMPS = 1;
		structStarMPS.spinposition = [1];							% Spin Node
		structStarMPS.tdvp = [];
		para.dk{pIdx{:}} = structStarMPS.dk;
		para.D{pIdx{:}} = structStarMPS.D;
		para.d_opt{pIdx{:}} = structStarMPS.d_opt;
		para.shift{pIdx{:}} = structStarMPS.shift;
		para.treeMPS.L(pIdx{:}) = 1;
		% generate Hamiltonian terms:
		structStarMPS.op = genh1h2term_onesite(para,structStarMPS.treeIdx,1);

		ii = 1;
		pIdx = para.treeMPS.nodeIdx(ii,:);
		dk = 2; d_opt = 2;
		treeMPS = structStarMPS;							% make a copy to inherit fields. Only modify difference
		treeMPS.mps = {randn(1,D,D,D,d_opt)};
		treeMPS.Vmat = {sparse(1:d_opt,1:d_opt,1,dk,d_opt)};
		treeMPS.child = [structChainMPS(1); structChainMPS(2); structStarMPS];
		treeMPS.treeIdx = [para.treeMPS.nodeIdx{ii,:}]-1;
		treeMPS.chainIdx = min([treeMPS.child.chainIdx]);
		treeMPS.level = 1;									% root node always level 1
		treeMPS.currentChild = 1;
		treeMPS.height = max([treeMPS.child.height])+1;
		treeMPS.degree = 3;
		treeMPS.nChains = 3;							% equals degree for nodes
		treeMPS.L = [];
		treeMPS.L = getTreeLength(treeMPS);				% L as array
		treeMPS.d_opt = d_opt;
		treeMPS.dk = dk;
		treeMPS.D = [1;ones(3,1)*D];					% [Dl;Dc1;Dc2;Dc3]  here vertical! Leading singleton -> root of tree!
		treeMPS.useStarMPS = 0;					% star is subset of tree!
		treeMPS.useTreeMPS = 1;
		treeMPS.spinposition = [1];				% Spin Node
		treeMPS.isRoot = 1;
		treeMPS.op = genh1h2term_onesite(para,treeMPS.treeIdx,1);
		para.treeMPS.L(pIdx{:}) = 1;
		
		para.tdvp.filename = 'test.mat';
		para.tdvp.filenameSmall = 'test-small.mat';
		para.tdvp.StoreMPS = 0;
		para.tdvp.resume = 0;
		para.tdvp.tmax = 100;
		para.tdvp.deltaT = 1;
		para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
		para.tdvp.extractObsInterval = para.tdvp.deltaT;
		para.tdvp.Observables = '.n.';
		para.tdvp.serialize = 0;
		para.tdvp.HEffSplitIsometry = 0;
		para.tdvp.imagT = 0;
		para.tdvp.useDkExpand = 0;		% does not yet work!
		para.tdvp.expandOBB = 1;
		para.tdvp.maxOBBDim = 10;
		para.tdvp.maxBondDim = [20];		%
		para.tdvp.truncateExpandBonds = 1;
		para.tdvp.evolveSysTrotter = 1;
		para.tdvp.expvTol = 1e-15;
		para.tdvp.expvM   = 20;
		para.tdvp.maxExpMDim = 300;				% For Lappy: 100, OE-PC: 80, pc52: 260; E5: 300 System dependent, use benchmark!
		para.tdvp.maxExpVDim = 700;				% higher dim -> use expvCustom() if expvCustom == 1. Number from benchmarking. Lappy: 400, Haswell: 800; E5: 700 maxExpMDim < maxExpVDim
		para.tdvp.expvCustom = 1;				% 1 for Custom programmed, 0 for standard expv()
		para.tdvp.expvCustomTestAccuracy = 0;	% do expvCustom alongside expv for testing.
		para.tdvp.expvCustomTestAccuracyRMS = 0;	% display RMS of expvCustom from expv(); set only if para.tdvp.expvCustomTestAccuracy = 1;
		para.tdvp.logSV = 0;
		para.tdvp.starttime = tic;
		para.tdvp.saveInterval = 10;
		
		tresults = [];
		results = initresults(para);
% 		results = initresultsTDVP(para, results);
		
		%% Do one prepare sweep to bring MPS onto right-canonical form
		[treeMPS,~,para] = prepare(treeMPS,[],para);
		
		%% Calculate effective Hamiltonians
		treeMPS = initstorage(treeMPS,[],[],para);		% new call pattern for treeMPS
		
	end
end

function out = D2Str(treeMPS)
%% function out = D2Str(treeMPS)
%	generates string output of bond dimensions returned in a flattened cell array for each line to display
if treeMPS.height == 0
	out = {num2str(treeMPS.D,'%-5u')};		% converts a chain to string; pack into cell
	return;
end
% else: hasChild!
childStr = {};
for i = 1:treeMPS.degree
	childStr{i} = D2Str(treeMPS.child(i));
end

for i = 1:treeMPS.degree
	% str1 for the first line of each child
	if i == 1
		str1 = sprintf('%-5u%s',treeMPS.D(1));
	elseif i == treeMPS.degree
		str1 = sprintf('  \\- ');
	else
		str1 = sprintf('  |- ');
	end
	if iscell(childStr{i})
		childStr{i}{1} = sprintf('%s%s',str1,childStr{i}{1});
		for j = 2:length(childStr{i})
			if i == treeMPS.degree
				childStr{i}{j} = sprintf('     %s',childStr{i}{j});
			else
				childStr{i}{j} = sprintf('  |  %s',childStr{i}{j});
			end
		end
	end
end
childStr = [childStr{:}];		% flatten cell array
out = childStr;

end

function out = minSV2Str(treeMPS)
%% function minSV2Str(treeMPS)
%	generates string output of smallest Amat singular value returned in a flattened cell array for each line to display
if treeMPS.height == 0
	out = {num2str(cellfun(@(x) -log10(x(find(x,1,'last'))),treeMPS.Amat_sv),'%-5.1f')};		% converts a chain to string; pack into cell
	return;
end
% else: hasChild!
childStr = {};
for i = 1:treeMPS.degree
	childStr{i} = minSV2Str(treeMPS.child(i));
end

for i = 1:treeMPS.degree
	% str1 for the first line of each child
	if i == 1
		str1 = sprintf('%-5.1f%s',-log10(treeMPS.Amat_sv{1}(end)));
	elseif i == treeMPS.degree
		str1 = sprintf('  \\- ');
	else
		str1 = sprintf('  |- ');
	end
	if iscell(childStr{i})
		childStr{i}{1} = sprintf('%s%s',str1,childStr{i}{1});
		for j = 2:length(childStr{i})
			if i == treeMPS.degree
				childStr{i}{j} = sprintf('     %s',childStr{i}{j});
			else
				childStr{i}{j} = sprintf('  |  %s',childStr{i}{j});
			end
		end
	end
end
childStr = [childStr{:}];		% flatten cell array
out = childStr;

% out = num2str(cellfun(@(x) -log10(x(find(x,1,'last'))),results.Amat_sv),'%-4.2g');
end

function results = initresultsTDVP(para, results)
	% save Dimension Log only as difference! reconstruct with cumsum(A)
% 	fprintf('Initialize results.tdvp\n');
	n = round(para.tdvp.tmax/para.tdvp.extractObsInterval) +1;
	NC = para.nChains; L = para.L;
	if para.tdvp.logSV
		if para.useVtens
			results.tdvp.Vmat_sv{n,L,NC+1}	= [];
		elseif para.useStarMPS
			results.tdvp.Vmat_sv{n,L,NC}	= [];
		end
		if size(results.Vmat_sv,1) == L
			results.tdvp.Vmat_sv(1,:,:)	 = results.Vmat_sv;
		else
			results.tdvp.Vmat_sv(1,:,:)	 = results.Vmat_sv.';
		end
		results.tdvp.Amat_sv{n,L-1} = [];
		results.tdvp.Amat_sv(1,:)		 = results.Amat_sv;
	else
		if para.useVtens
			results.tdvp.Vmat_vNE(n,L, NC+1) = 0;
			results.tdvp.Amat_vNE(n,L-1)	 = 0;
		elseif para.useStarMPS
			results.tdvp.Vmat_vNE(n,L, NC)   = 0;
			if L ~= 1
				results.tdvp.Amat_vNE(n,L-1, NC) = 0;
			end
		end
		results.tdvp.Vmat_vNE(1,:,:)	= results.Vmat_vNE.';
		results.tdvp.Amat_vNE(1,:,:)	= results.Amat_vNE.';
	end
	if para.tdvp.expandOBB
		if para.useVtens
			results.tdvp.d_opt              = sparse(n,L*(NC+1));
		else
			results.tdvp.d_opt              = sparse(n,L*NC);
		end
		results.tdvp.d_opt(1,:)			= reshape(para.d_opt.',1,[]);
	end
	if para.tdvp.truncateExpandBonds
		if para.useStarMPS
			results.tdvp.D				= sparse(n,(L-1)*NC);
		else
			results.tdvp.D				= sparse(n,L-1);
		end
		results.tdvp.D(1,:)				= reshape(para.D.',1,[]);
	end
	if para.tdvp.useDkExpand
		if para.useStarMPS
			results.tdvp.dk				= sparse(n,L*NC);
		else
			results.tdvp.dk				= sparse(n,L);
		end
		results.tdvp.dk(1,:)			= reshape(para.dk.',1,[]);
	end
%	results.tdvp.dk 		  = sparse(n,para.L);
%	results.tdvp.dk(1,:)      = para.dk;
	results.tdvp.expvTime     = [];
% 	results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, (para.L-1)*4*2+3);
	results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, 1);		% otherwise file gets too large
	results.tdvp.EvaluesLog	  = zeros(length(para.tdvp.t)-1, 1);
end

function [mps,Vmat,para,results] = copyChainToStar(mps,mpsC,Vmat,VmatC,para,paraC,results,resultsC)
% only use this after sweep along a single chain
nc = para.currentChain;
NC = para.nChains;
L  = para.chain{nc}.L;

% Copy the MPS from Chain to Star!
for ii = 1:length(mpsC)
	mps{ii+1}{nc} = mpsC{ii};
	Vmat{ii+1}{nc} = VmatC{ii};
end

para.dk(nc,2:L)    = paraC.dk;
para.d_opt(nc,2:L) = paraC.d_opt;
para.D(nc,2:L-1)   = paraC.D;

results.Vmat_sv(nc,2:L)    = resultsC.Vmat_sv;
results.Amat_sv(nc,2:L-1)  = resultsC.Amat_sv;
results.Vmat_vNE(nc,2:L)   = resultsC.Vmat_vNE;
results.Amat_vNE(nc,2:L-1) = resultsC.Amat_vNE;

end

function [Atens, CA] = prepare_Tens(Atens, para)
	%% Input: focused AS = mps{1} for star-MPS
	%  sets focus onto Vtens{mc}
	%  no truncation for now!
	%
	mc = para.currentChain;								% the chain to focus on
	d  = size(Atens);

	[AS]     = tensShape(Atens, 'unfold', mc+1, d);		% D(mc) x rest

	if diff(size(AS)) >= 0
		[AS, CA] = prepare_onesiteVmat(AS.', para);			% need for transpose since m x n, m > n input needed
	else
		% AS: D x prod(rest)
		[CA, AS] = prepare_onesiteVmat(AS, para);			% means that prod(rest) < D, thus AS will be smaller than before!
		CA = CA.';
		AS = AS.';											% prod(rest) x prod(rest)
		d(mc+1) = size(AS,1);
	end

	Atens	 = tensShape(AS.', 'fold', mc+1, d);
end

function treeMPS = tdvp_1site_evolveTree(treeMPS,para)
%% function treeMPS = tdvp_1site_evolveTree(treeMPS,para)
%
%	evolves node and calls evolution of children recursively
%	treeMPS is not shared with workspace of tdvp_1site_tree. Only shared with tdvp_1site_evolveNode()
%	Created by FS 23/02/2016

fprintf('\nNode %s:\n',mat2str(treeMPS.treeIdx));
%% Perform TDVP on Node
% only use expvCustom due to complicated contractions!
fprintf('1');
para.sitej = 1; treeMPS.sitej = 1; treeMPS.sweepto = 'r';
tdvp_1site_evolveNode();				% nested function for now!
fprintf('-');

for ii = 1:treeMPS.degree
	% Iterate through children, go down each branch
	treeMPS.currentChild = ii; treeMPS.currentChain = ii; para.currentChain = ii;
	para.sweepto = 'r'; treeMPS.child(ii).sweepto = 'r';
	
	% SVD to get A_C(mc)
	[treeMPS.mps{1}, Cn] = prepare_Tens(treeMPS.mps{1}, para);

	% calculate op.Hleft, op.Opleft for the next chain! Also copies H/Opright.
	treeMPS.child(ii).op = updateop([],treeMPS,[],[],para);
	
	% evolve AC back in time
	% declare site as 0 to achieve contraction into mps{1}
	treeMPS.child(ii).mps = tdvp_1site_evolveKn(treeMPS.child(ii).mps,[],para,[],treeMPS.child(ii).op,0,Cn,[]);
% 	[mpsChain, ~, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,[],paraChain,resultsChain,op,0,Cn,[]);
	% now focused on mps{1} of child!
	
	% Start recursion
	if treeMPS.child(ii).height == 0
		treeMPS.child(ii) = tdvp_1site_evolveChain(treeMPS.child(ii),para);
	else
		treeMPS.child(ii) = tdvp_1site_evolveTree(treeMPS.child(ii),para);
	end
	% Now need to contract BondCenter of child into treeMPS.mps{1}
	nd = ndims(treeMPS.mps{1});
	treeMPS.mps{1}  = contracttensors(treeMPS.child(ii).BondCenter.', 2, 2, treeMPS.mps{1}, nd, ii+1);
	treeMPS.mps{1}  = permute(treeMPS.mps{1},[2:ii+1, 1, ii+2:nd]);				% order back
	treeMPS.D(ii+1) = size(treeMPS.child(ii).BondCenter,2);								% this is new Bond dimension in the child's edge
	
end

treeMPS.sweepto = 'l'; para.sweepto = 'l';
tdvp_1site_evolveNode();				% nested function for now!

% Shift focus to the parent!
d = size(treeMPS.mps{1});
if treeMPS.isRoot
	[treeMPS.mps{1}, Cn,~,~,sv,vNE] = prepare_onesite(reshape(treeMPS.mps{1},d(1),[],d(end)),para,1,[]);		% TODO: replace by faster 2D SVD?
else
	[treeMPS.mps{1}, Cn,~,~,sv,vNE] = prepare_onesite_truncate(reshape(treeMPS.mps{1},d(1),[],d(end)),para,1,[]);		% Do truncation with maxBondDim(1) TODO: replace by faster 2D SVD?
end
d(1) = size(Cn,2); treeMPS.D(1) = d(1);																			% can be smaller for ER-tensors -> update

treeMPS.mps{1} = reshape(treeMPS.mps{1},d);
treeMPS.Amat_sv = {sv};

treeMPS.op = updateop([],treeMPS,[],[],para);

if ~treeMPS.isRoot
	% sitej = 0:
	treeMPS.BondCenter = tdvp_1site_evolveKn(treeMPS.mps,[],para,[],treeMPS.op,0,Cn,[]);
	% [    Cn  , VmatChain, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej,Cn,Hn);
else
	treeMPS.BondCenter = Cn;		% should == 1
	if abs(1-abs(Cn)) > 10^-5 
		fprintf('\n State norm: %g', Cn);
	end
end
	
	function tdvp_1site_evolveNode()
	%% function tdvp_1site_evolveNode()
	%	evolves the Node = site 1 of a Tree-MPS
	%
	% Created 23/02/2016 by FS
	%
	t = para.tdvp.deltaT./2;
	
	if para.tdvp.imagT
		t = -1i*t;
	end
	
	d = size(treeMPS.mps{1});
	
	% skip time-evolution of Vmat for now, since only small dk!
	% add TDVP for Vmat if placing boson on node
	
	% update site 1 OBB terms h1j/h2jOBB
	if isempty(treeMPS.op.h1jOBB) || isempty(treeMPS.op.h1j)
		treeMPS.op = H_Eff(treeMPS, [], 'TR-A' , [], para);		% bring into OBB
	end
	
	%% Take matrix exponential
	% A (t+dt) = exp(-i ?? dt) * A(t)
	treeMPS.tdvp.expvTol = para.tdvp.expvTol;		% replace para by treeMPS in expvCustom call for access to children operators
	treeMPS.tdvp.expvM   = para.tdvp.expvM;
	if para.tdvp.evolveSysTrotter == 0 
		% Evolves system in one single step interacting with all chains at once
		% better for nodes with < 3-4 children
		[treeMPS.mps{1}, ~] = expvCustom(- 1i*t,'TREE-Hn1',treeMPS.mps{1}, treeMPS, treeMPS.op);

		treeMPS.mps{1} = reshape(treeMPS.mps{1},d);
	else
		% newer scheme, trotterise in system-chain interactions to decrease complexity
		% only for nodes which also have a local site
		Atens = treeMPS.mps{1};
		dIn = d;
		for mc = 0:treeMPS.nChains
			%% create isometry into chain direction
			% rotate chain bonds
			para.currentChain = mc; treeMPS.currentChain = mc;
			if mc == 0 && treeMPS.isRoot
				continue;		% nothing to evolve!
			elseif mc == 1 && treeMPS.isRoot
				[Atens,dOut] = tensShape(Atens, 'foldrotunfoldiso', 2, dIn);	% rot by 2 since leading singleton!
			else
				[Atens,dOut] = tensShape(Atens, 'foldrotunfoldiso', 1, dIn);	% rotates such that A: prod(D(1:NC without mc)) x (D(mc,1) * dk)
			end

			[Iso, A] = qr(Atens,0);			% Iso is isometry with all unused chains
			
			% TODO: comment if not testing!!
	% 		[Iso2,A2]= rrQR(Atens, floor(0.6*prod(dOut(end-1:end))),0);		% low rank QR approximation
	% 		fprintf('evolve: rank(A) = %d, dim(A,2) = %d, rrQR error = %g\n', rank(Atens), prod(dOut(end-1:end)),norm(Atens - Iso2*A2));
			
			if size(A,2) ~= size(A,1)
				warning('VMPS:tdvp_1site_tree:tdvp_1site_evolveNode:BadDimensions','Matrix has wrong shape for trotter splitting. If error occurs please use para.tdvp.evolveSysTrotter = 0.');
			end
			A = reshape(A,[],dOut(end-1),dOut(end));
			% evolve simplified mps matrix
			[A, err] = expvCustom(- 1i*t,'TREE-Hn1Trotter',A, treeMPS, treeMPS.op);
			A = reshape(A, prod(dOut(end-1:end)),[]);			% D*dk x D*dk
			% contract back together
			Atens = Iso * A;

			dIn = dOut;			% reset the new dimensions after rotation

		end

		treeMPS.mps{1} = reshape(Atens, d);
	end
	
	end

%% function [A, Iso] = prepare_IsometryDkChain(Atens, para)
% 	%% Input: focused AS = mps{1} of star-MPS
% 	%  Splits Atens into Center A, containing bonds n_k, and D(para.currentChain)
% 	%	and the Isometry containing all other chain bonds
% 	%	see PEPS techniques for Trotter gates
% 	%  no truncation!
% 	% Output:
% 	%	A: D(NC,1)*dk x D(NC,1) x dk
% 	%	I:
% 	mc = para.currentChain;								% the chain to focus on
% 	d  = size(Atens);
%
% 	[AS]     = tensShape(Atens, 'unfold', mc+1, d);		% D(mc) x rest
%
% 	[AS, CA] = prepare_onesiteVmat(AS.', para);			% need for transpose since m x n, m > n input needed
%
% 	Atens	 = tensShape(AS.', 'fold', mc+1, d);
%
% end

end

function treeMPS = tdvp_1site_evolveChain(treeMPS,para)
%% function treeMPS = tdvp_1site_evolveChain(treeMPS,para)
%
%	evolves a single chain
%	since para is not returned, copy relevant fields from treeMPS

copyPara();
results = initresults(para);		% TODO: change later!
results = initresultsTDVP(para, results);
% Quick dirty fix: want chain to node dimension == maxBondDim(2)
maxBondDim = para.tdvp.maxBondDim;
para.tdvp.maxBondDim = para.tdvp.maxBondDim(end);

idx = num2cell(treeMPS.treeIdx + 1);		% get idx in para to print chain number
fprintf('\nChain %d:\n',para.treeMPS.chainIdx{idx{:}});
%% Sweep l -> r:
for sitej = 1:para.L
	fprintf('%g', sitej); para.sitej = sitej;
	%% Update on-site Operators
	treeMPS.op = gen_sitej_op(treeMPS.op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

	%% Do the time-evolution of A and V with H(n)
	% this is symmetric for l->r and l<-r
	if ~para.useVtens
		[treeMPS.mps{sitej}, treeMPS.Vmat{sitej}, para, results, treeMPS.op, Hn] = tdvp_1site_evolveHn(treeMPS.mps,treeMPS.Vmat,para,results,treeMPS.op,sitej);
	else
		para.tdvp.expvCustomNow = 1;                                     % necessary setting!
		[treeMPS.mps{sitej}, treeMPS.Vmat{sitej}, para, results, treeMPS.op, Hn] = tdvp_1site_evolveHnMC(treeMPS.mps,treeMPS.Vmat,para,results,treeMPS.op,sitej);
	end
	% now: A and V are time-evolved, A is focused
	% OBBDim has increased by 50%. This must be truncated in the next
	% step again!
	% if sitej = L, then start lr sweep with decomposition of mps

	%% Take focus to next site, evolve with K(n)
	if sitej ~= para.L
		fprintf('-');
		%% Left-normalize A and get Center matrix C(n,t+dt)_(rl,r)
		% expand/truncate BondDimensions here?
		[treeMPS.mps{sitej}, Cn, para,results] = prepare_onesite_truncate(treeMPS.mps{sitej}, para,sitej,results);

		%% Do the time-evolution of C
		% evolve non-site center between sitej and sitej+1
		% and focus A(n+1)
		[treeMPS.mps, ~, para, results] = tdvp_1site_evolveKn(treeMPS.mps,[],para,results,treeMPS.op,sitej,Cn,Hn);
		clear('Hn','Cn');

		%% update Left Hamiltonian operators
		treeMPS.op = updateop(treeMPS.op,treeMPS.mps,treeMPS.Vmat,sitej,para);
		
	else % sitej = L
		%% Normalize with last SVD
		[treeMPS.mps{sitej}, Cn] = prepare_onesite_truncate(treeMPS.mps{sitej}, para,sitej);
		% Cn = 1 approximately. can be thrown away -> mps normalized
		% finish with focus on A after time evolution

	end
end

% 		fprintf('\n');							%Debug!
% 		fprintf('%2g-',para.d_optnew);
% 		fprintf('\n');

%% Log vNE etc ?
% from prepare_onesite() and prepare_onesiteVmat():
% results.Amat_vNE  array
% results.Amat_sv   cell
% results.Vmat_vNE
% results.Vmat_sv
%% Output matrix dimensions if changed:
fprintf('\n');
if para.tdvp.truncateExpandBonds
	fprintf('para.D:\t');
	out = strrep(mat2str(para.D),';','\n');
	fprintf([out(2:end-1),'\n']);
	out = mat2str(cellfun(@(x) x(find(x,1,'last')),results.Amat_sv),2);
	fprintf([out(2:end-1),'\n']);
end
if para.tdvp.expandOBB
	fprintf('para.d_opt:\t');
	out = strrep(mat2str(para.d_opt),';','\n');
	fprintf([out(2:end-1),'\n']);
	out = mat2str(cellfun(@(x) x(end),results.Vmat_sv(end,~cellfun('isempty',results.Vmat_sv(end,:)))),2);
	fprintf([out(2:end-1),'\n']);
end

%% SWEEP l <- r:
para.sweepto = 'l';
% now Hn = H(L)_(l'*r'*n',l*r*n)
for sitej = para.L-1:-1:0
	fprintf('-'); para.sitej = sitej;

	%% Right-normalize A(n+1) and get Center matrix C(n,t+dt)_(l,lr)
	% normalisation needed for updateop()!
	% Applies to bond n -> use sitej for result storage
	if para.tdvp.truncateExpandBonds
		[treeMPS.mps{sitej+1}, Cn, para, results, sv, vNE] = prepare_onesite_truncate(treeMPS.mps{sitej+1},para,sitej+1, results);
	else
		[treeMPS.mps{sitej+1}, Cn, para, results, sv, vNE] = prepare_onesite(treeMPS.mps{sitej+1},para,sitej+1, results);
	end

	%% Do the time-evolution of C
	% evolve non-site center between sitej and sitej+1
	% and focus A(n)
	if sitej == 0
		[    Cn     ,~, para, results] = tdvp_1site_evolveKn(treeMPS.mps,[],para,results,treeMPS.op,sitej,Cn,Hn);
		% contract Cn with mps{1} after for loop
	else
		[treeMPS.mps,~, para, results] = tdvp_1site_evolveKn(treeMPS.mps,[],para,results,treeMPS.op,sitej,Cn,Hn);
	end
	clear('Hn');

	fprintf('%g', sitej);

	%% Expand / Truncate dk if needed!
	if para.tdvp.useDkExpand && sitej+1 >= 1
		% para is needed to call genh1h2term_onesite correctly
		[treeMPS.Vmat{sitej+1}, para, treeMPS.op] = truncateExpandDk(treeMPS.Vmat{sitej+1},para,results,treeMPS.op,sitej+1);
	end

	%% update right Hamiltonian operators
	% prepare operators in (sitej+1) for time-evolution on sitej
	para.sitej = para.sitej+1;					% needed for multi-chain reshape
	treeMPS.op         = updateop(treeMPS.op,treeMPS.mps,treeMPS.Vmat,sitej+1,para);
	para.sitej = para.sitej-1;
	
	if sitej == 0
		fprintf('\n');
		continue;
	end
	%% Get on-site Operators and dimensions
	treeMPS.op = gen_sitej_op(treeMPS.op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

% 		shift = getObservable({'1siteshift'}, mps{sitej}, Vmat{sitej}, para); % calculate applicable shift
% 		fprintf('\n Shift: %s', mat2str(shift));
%
	%% Do the time-evolution of A and V
	% this is symmetric for l->r and l<-r
	if ~para.useVtens
		[treeMPS.mps{sitej}, treeMPS.Vmat{sitej}, para, results, treeMPS.op, Hn] = tdvp_1site_evolveHn(treeMPS.mps,treeMPS.Vmat,para,results,treeMPS.op,sitej);
	else
		[treeMPS.mps{sitej}, treeMPS.Vmat{sitej}, para, results, treeMPS.op, Hn] = tdvp_1site_evolveHnMC(treeMPS.mps,treeMPS.Vmat,para,results,treeMPS.op,sitej);
	end
end

% sitej = 0: Done above!
% treeMPS.BondCenter = tdvp_1site_evolveKn(treeMPS.mps,[],para,[],treeMPS.op,0,Cn,[]);
% [    Cn  , treeMPS.Vmat, para, results] = tdvp_1site_evolveKn(treeMPS.mps,treeMPS.Vmat,para,results,opChain,sitej,Cn,Hn);
treeMPS.BondCenter = Cn;

% Copy Dimension Information
treeMPS.D(2:end-1) = para.D;
treeMPS.D(1)       = size(treeMPS.mps{1},1);
treeMPS.dk         = para.dk;
treeMPS.d_opt      = para.d_opt;
treeMPS.Amat_sv    = {sv,results.Amat_sv{:}};
% can be omitted:
% results.Amat_sv{1,1}  = sv;
% results.Amat_vNE(1,1) = vNE;		% from last prepare_onesite_truncate
para.tdvp.maxBondDim = maxBondDim;		% quick fix: rest bonddim

	function copyPara()
		% Copy relevant fields from treeMPS to para
		para.sweepto      = 'r';
		para.sitej        = 1;
		para.currentChain = 1;
		para.L            = treeMPS.L;
		para.M            = treeMPS.M;
		para.D            = treeMPS.D(2:end-1);
		para.dk           = treeMPS.dk;
		para.d_opt        = treeMPS.d_opt;
		para.shift        = treeMPS.shift;
		para.useVmat      = treeMPS.useVmat;
		para.useVtens     = treeMPS.useVtens;
		para.useStarMPS   = treeMPS.useStarMPS;
		para.useTreeMPS   = treeMPS.useTreeMPS;
		para.nChains      = treeMPS.nChains;
		para.treeIdx      = treeMPS.treeIdx;
		para.spinposition = treeMPS.spinposition;
	end

	function [Vmat,para,op] = truncateExpandDk(Vmat,para,results,op,s)
		%% function [Vmat, para, op] = truncateExpandDk(Vmat,para,results,op,s)
		%	truncates or expands the local dimension dk based on the OBB usage
		%	the vectors in Vmat and their SV are the indicators
		%	
		%	For now: not for Vtens
		%	Input: Vmat = Vmat{sitej}, single-site Vmat
		if para.useVtens || ~para.useDkExpand || ~para.useVmat
			return;			% Do nothing
		end
		nc = para.currentChain;

		%% Determine whether to truncate or to expand: Copied from adjustdopt.m
		adddim = estimateDkExpand(Vmat,results.Vmat_sv{1,s}, para);
		if adddim > 0
		%% Apply the expansion
			[dk,dOBB]      = size(Vmat);
			para.dk(nc,s)  = dk+adddim;						% operators will be expanded in genh1h2term?
			addmat  = zeros(adddim,dOBB);
	% 			Vmat = cat(1,addmat,Vmat);			% N:1 ordering
			Vmat = cat(1,Vmat,addmat);			% 1:N ordering
			para.hasexpanded = 1;
			para.increasedk  = 0;							% reset this value
	% 			dispif('Increased dk',para.logging)
		end
		% expand the operators:
		temp = genh1h2term_onesite(para,para.treeIdx,s);
		op.h1term(s)     = temp.h1term;
		op.h2term(:,:,s) = temp.h2term;
		% since only zeros were added in Vmat, no further update of h1j/OBB etc is needed.
	end
end


function L = getTreeLength(treeMPS)
%% function L = getTreeLength(treeMPS)
%	extracts the length of the treeMPS, based on the length of its children
%
if isfield(treeMPS,'L') && ~isempty(treeMPS.L)
	L = treeMPS.L;
	return
elseif treeMPS.height == 0
	if isfield(treeMPS,'mps')
		L = length(treeMPS.mps);
	else
		L = 0;					% is leaf without MPS
	end
	return
end

% first determine size(L)
h = treeMPS.height;			% ndims of L > 0
d = zeros(1,h);				% this will be size(L)
d(1) = treeMPS.degree+1;

if h == 1					% this is StarMPS with leaves only
	Lchild = arrayfun(@getTreeLength, treeMPS.child);
	L = [1;reshape(Lchild,d(1)-1,1)];
	return
end
% else height >= 2

if h == 2
	d(3) = 1;		% need singleton to allow comparisons later
end
Lchild = arrayfun(@getTreeLength, treeMPS.child,'UniformOutput',false);
% find out maximum dimensions per level to zero-pad Lchild
idx = find([treeMPS.child.height] == h-1);		% only pick the ones with maximum height
for ii = idx
	d(2:end) = max(d(2:end),size(Lchild{ii}));
end
% now d contains the maximum dimensions of L -> initialise
L = zeros(d);
L(1) = 1;					% length of nodes = 1 always for now!
for ii = 1:treeMPS.degree
	idx = arrayfun(@(x) 1:x,size(Lchild{ii}),'UniformOutput',false);		% create cell array containing index ranges of Lchild{ii}
	L(ii+1,idx{:}) = Lchild{ii};
end

end


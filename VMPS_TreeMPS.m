function filename = VMPS_TreeMPS(varargin)
%% function filename = VMPS_TreeMPS(varargin)
%
%	Initialises the model and starts the TDVP computation with treeMPS.
%	Supports flexible calling syntax via input parser. Input can be Struct with appropriate field names
%
%	Created by FS 26/02/2016
%

%%	Define Input Parser
p = inputParser;					% parses for para
p.CaseSensitive = true;
p.KeepUnmatched = true;
p.FunctionName = 'VMPS_TreeMPS';

pt = inputParser;					% parses for para.tdvp
pt.CaseSensitive = false;
pt.KeepUnmatched = true;
pt.FunctionName = 'VMPS_TreeMPS';

pSBM = inputParser;					% parses for para.SpinBoson
pSBM.CaseSensitive = false;
pSBM.KeepUnmatched = true;
pSBM.FunctionName = 'VMPS_TreeMPS';

pDPMES = inputParser;					% parses for para.SpinBoson
pDPMES.CaseSensitive = false;
pDPMES.KeepUnmatched = true;
pDPMES.FunctionName = 'VMPS_TreeMPS';

%% Define local parameters
SBMInitialStates   = {'sz', '-sz', 'sx', '-sx', 'sy', '-sy', 'none'};
DPMESInitialStates = {'TT','LE+', 'LE-', 'CT+', 'CT-'};					% Position of Name = number of InitialState

%% Define validation functions for string expressions
% para input parser
p.addParameter('model'			,'SpinBoson',@(x) any(validatestring(x,...
								{'SpinBoson','SpinBoson2C', 'DPMES5-7C', 'DPMESclust7-1', 'DPMES-Tree1','DPMES-Tree2',...
								 'DPMES-Tree3','DPMES-Tree4','DPMES-Tree6','testTree'})));		% possible model inputs
p.addParameter('s'				,1	,@isnumeric);
p.addParameter('L'				,50	,@isnumeric);		% longest L between root and all leaves
p.addParameter('alpha'			,0.1,@isnumeric);
p.addParameter('delta'			,0	,@isnumeric);		% in sx
p.addParameter('epsilon'		,0.1,@isnumeric);		% in sz
p.addParameter('wc'				,1	,@isnumeric);
p.addParameter('dk_start'		,5	,@isnumeric);
p.addParameter('D_start'		,5	,@isnumeric);
p.addParameter('d_opt_start'	,5	,@isnumeric);
p.addParameter('M'				,2	,@isnumeric);
p.addParameter('useDkExpand'	,0	,@isnumeric);
p.addParameter('dkmax'			,1e3,@isnumeric);
p.addParameter('Dmin'			,2	,@isnumeric);
p.addParameter('svmaxtol'	,1e-6	,@isnumeric);
p.addParameter('svmintol'	,10^-8,@isnumeric);
p.addParameter('logging'		,0	,@isnumeric);
p.addParameter('rescaling'		,0	,@isnumeric);
p.addParameter('complex'		,0	,@isnumeric);
p.addParameter('tdvp'			,1	,@isnumeric);		% do time-evolution?
p.addParameter('spinposition'	,1	,@isnumeric);
p.addParameter('resume'			,0	,@isnumeric);
p.addParameter('precision'	,5e-15	,@isnumeric);
p.addParameter('parity'			,'n',@isstr);
p.addParameter('eigs_tol'	,1e-8	,@isnumeric);
p.addParameter('foldedChain'	,0	,@isnumeric);
p.addParameter('d_opt_min'		,2	,@isnumeric);
p.addParameter('spinbase'		,'Z',@isstr);
p.addParameter('useshift'		,0	,@isnumeric);
p.addParameter('useTreeMPS'		,1	,@isnumeric);
p.addParameter('useStarMPS'		,0	,@isnumeric);
p.addParameter('useVmat'		,1	,@isnumeric);
p.addParameter('useVtens'		,0	,@isnumeric);
p.addParameter('initChainState'	,'vac',@(x) any(validatestring(x,...
								{'rand','vac'})));		% possible initial chain states
p.addParameter('comment'		,'',@isstr);

% para.tdvp input parser
% pt.addParameter('model'			,'SpinBoson',@(x) any(validatestring(x,...
% 								{'SpinBoson','SpinBoson2C', 'DPMES5-7C','testTree'})));		% possible model inputs

pt.addParameter('imagT'				,0	,@isnumeric);
pt.addParameter('tmax'				,20	,@isnumeric);
pt.addParameter('deltaT'			,1	,@isnumeric);
pt.addParameter('resume'			,0	,@isnumeric);
pt.addParameter('saveInterval'		,10	,@isnumeric);
pt.addParameter('serialize'			,1	,@isnumeric);
pt.addParameter('logSV'				,0	,@isnumeric);
pt.addParameter('Observables'		,'.dm.n.' ,@isstr);
pt.addParameter('extractObsInterval',1	,@isnumeric);
pt.addParameter('extractStarInterval',1	,@isnumeric);
pt.addParameter('storeMPS'			,0	,@isnumeric);
pt.addParameter('evolveSysTrotter'	,1	,@isnumeric);
pt.addParameter('HEffSplitIsometry'	,1	,@isnumeric);
pt.addParameter('maxExpMDim'		,120,@isnumeric);	% E5: 160, P40: 120
pt.addParameter('maxExpVDim'		,200,@isnumeric);	% E5: 270, P40: 256
pt.addParameter('expvCustom'		,1	,@isnumeric);
pt.addParameter('useDkExpand'		,0	,@isnumeric);
pt.addParameter('expandOBB'			,1	,@isnumeric);
pt.addParameter('truncateExpandBonds',1	,@isnumeric);
pt.addParameter('maxOBBDim'			,20	,@isnumeric);
pt.addParameter('maxBondDim'		,5	,@isnumeric);
pt.addParameter('zAveraging'		,0	,@isnumeric);
pt.addParameter('logError'			,0	,@isnumeric);
pt.addParameter('logTruncError'		,0	,@isnumeric);
pt.addParameter('rescaling'			,0	,@isnumeric);
% pt.addParameter('storeMPS'			,0	,@isnumeric);

pSBM.addParameter('GroundStateMode'	,'artificial',@(x) any(validatestring(x,...
									{'decoupled','coupled', 'artificial','artTTM'})));			% GroundState preparations
pSBM.addParameter('InitialState'	,'sx',@(x) any(validatestring(x,SBMInitialStates)));		% Initial States

pDPMES.addParameter('InitialState'	,'LE+',@(x) any(validatestring(x,DPMESInitialStates)));		% Initial States
pDPMES.addParameter('CTShift'		,0	,@isnumeric);
pDPMES.addParameter('BroadenChain'  ,zeros(1,7),@isnumeric);									% logic, indicate which chains shall be broadened
pDPMES.addParameter('PeakWidth'     ,cmToeV(20),@isnumeric);									% width of 20 cm^-1 as standard?

%%
p.parse(varargin{:});
para = p.Results;
if para.tdvp == 1
	pt.parse(varargin{:});
	para.tdvp = pt.Results;
	para.complex = 1;
end

%% Defaults to the models

%% Settings for SBM
if strfind(para.model,'SpinBoson')
	pSBM.parse(varargin{:});
	para.SpinBoson = pSBM.Results;

%% Setting Hamiltonian for single-spin models
	% H0 = -para.hx./2.*sigmaX -para.hz./2.*sigmaZ
    para.hx = -para.delta;                       % Splitting with sigma_X
    para.hz = -para.epsilon;                     % Splitting with sigma_Z
%% Setting Chain for SBM
	if strcmp(para.model,'SpinBoson')
		para.nChains = 1;
	else
		out = regexp(para.model,'SpinBoson(?<nChains>\d)C','names');
		para.nChains = str2double(out.nChains);
	end
% 	para.wc = [1,5];		% only for testing purposes
	for mc = 1:para.nChains
		para.chain{mc}.mapping			= 'OrthogonalPolynomials';
		para.chain{mc}.spectralDensity	= 'Leggett_Hard';
		para.chain{mc}.discretization	= 'None';
		% para.chain{1}.discrMethod		= 'Numeric';

		para.chain{mc}.s				= para.s(min(mc,length(para.s)));			% SBM spectral function power law behaviour
		para.chain{mc}.alpha			= para.alpha(min(mc,length(para.alpha)));		% SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
		para.chain{mc}.L				= para.L(min(mc,length(para.L)));
		para.chain{mc}.w_cutoff         = para.wc(min(mc,length(para.wc)));
		para.chain{mc}.initState		= para.initChainState;						% choose: 'rand', 'vac', 
	end
	
%% Setting TreeMPS for SBM
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.nNodes    = 1;
	para.treeMPS.maxDegree = para.nChains;
	para.treeMPS.leafIdx   = num2cell((1:para.nChains)'+1);			% maps from chain number to leaf index in para
	para.treeMPS.nodeIdx   = {0+1};									% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     para.SpinBoson.GroundStateMode = 'artificial';				% input arg
        % choose: 'decoupled', 'coupled', 'artificial', 'artTTM'
		% -artificial & artTTM does no optimization! this only sets up an artificial
		%		ground state with <n> = 0 on chain and InitialState 'sz'
		% - 'artificial' for SBM2CT sets up maximally entangled state between odd and even chains in the Bath
		%	odd: thermal bath, even: ancilla
%     para.SpinBoson.InitialState = 'sz';							% input arg
        % choose: 'sz', '-sz', 'sx', '-sx', 'sy', '-sy', 'none'
		% works with all options

	if strcmp(para.SpinBoson.GroundStateMode, 'decoupled')
        para.SpinBoson.t1 = para.t(1);
        para.t(1) = 0;                              % switches off interaction with bath
	end
	if strcmp(para.SpinBoson.InitialState, 'sx') && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.SpinBoson.hx = para.hx;
        para.SpinBoson.hz = para.hz;
        para.hz = 0;
        para.hx = 10;                               % some value to do a sx splitting
    elseif strcmp(para.SpinBoson.InitialState, 'sz') && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.SpinBoson.hx = para.hx;
        para.SpinBoson.hz = para.hz;
        para.hz = 10;                               % some value to get sz+ eigenstate
        para.hx = 0;
	end
	if isempty(strfind(para.model,'TTM'))
		assert(length(para.spinposition) == 1, 'Only one spin is allowed in SBM');
	end
end

%% Settings for DPMES
if strfind(para.model,'DPMES5-7C')
 	pDPMES.parse(varargin{:});
	
%% Setting Chain for SBM
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20160129/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.chain{1}.mapping			= 'LanczosTriDiag';
	para.chain{1}.spectralDensity	= 'CoupDiscr';
	para.chain{1}.Lambda            = 1;
	para.chain{1}.initState		    = para.initChainState;
	
	para.chain{1}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-7-01.dat'));
	para.chain{1}.L					= min(length(para.chain{1}.dataPoints),para.L-1);
	
	%% chain 2 & more:
	para.chain{2}					= para.chain{1};		% simple copy
	para.chain{2}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-10-x1.dat'));
	para.chain{2}.L					= min(length(para.chain{2}.dataPoints),para.L-1);
	para.chain{3}					= para.chain{2};
	para.chain{3}.dataPoints		= cmToeV(load('DPMESdata_20160129/W14-A2-10-highv2.dat'));
	para.chain{3}.L					= min(length(para.chain{3}.dataPoints),para.L-1);
	para.chain{4}					= para.chain{2};
	para.chain{4}.dataPoints		= cmToeV(load('DPMESdata_20160129/W24-B1-17-highv2.dat'));
	para.chain{4}.L					= min(length(para.chain{4}.dataPoints),para.L-1);
	para.chain{5}					= para.chain{2};
	para.chain{5}.dataPoints		= cmToeV(load('DPMESdata_20160129/W23-B2-8-10.dat'));
	para.chain{5}.L					= min(length(para.chain{5}.dataPoints),para.L-1);
	para.chain{6}					= para.chain{2};
	para.chain{6}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1x.dat'));
	para.chain{6}.L					= min(length(para.chain{6}.dataPoints),para.L-1);
	para.chain{7}					= para.chain{2};
	para.chain{7}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1-x.dat'));
	para.chain{7}.L					= min(length(para.chain{7}.dataPoints),para.L-1);
	
%% Setting TreeMPS for SBM
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.nNodes    = 1;
	para.treeMPS.maxDegree = para.nChains;
	para.treeMPS.leafIdx   = num2cell((1:para.nChains)'+1);			% maps from chain number to leaf index in para
	para.treeMPS.nodeIdx   = {0+1};									% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

if strfind(para.model,'DPMESclust7-1')
	% This is star-model obtained using fkmeans from 08/04/2017
	% clusters under weighting by coupling strength only!
	% Works only for TreeMPS now!
 	pDPMES.parse(varargin{:});
	
%% Setting Chain for SBM
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20170408/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.BroadenChain				= pDPMES.Results.BroadenChain;								% save for reference
	
	% Chains:		1		2		3		4		5		6		7
	ChainFiles = {	'B11',	'B12',	'A11',	'A12',	'B22',  'B23',	'A2'};	% order similar to Tree2 since this is more likely to be optimal
	
	for ii = 1:length(ChainFiles)
		para.chain{ii}.mapping			= 'LanczosTriDiag';
		para.chain{ii}.spectralDensity	= 'CoupDiscr';
		para.chain{ii}.Lambda           = 1;
		para.chain{ii}.initState		= para.initChainState;
		para.chain{ii}.dataPoints		= cmToeV(load(sprintf('DPMESdata_20170408/Modes-%s.dat',ChainFiles{ii})));
		para.chain{ii}.H1				= load(sprintf('DPMESdata_20170408/H1-%s.dat',ChainFiles{ii}));
		para.chain{ii}.L				= min(length(para.chain{ii}.dataPoints),para.L-1);
		if pDPMES.Results.BroadenChain(ii)
			% allows individual broadening of each chain
			para.chain{ii}.spectralDensity	= 'CoupBroad';
			para.chain{ii}.discrMethod		= 'Direct';
			para.chain{ii}.discretization	= 'Linear';
			para.chain{ii}.mapping			= 'Stieltjes';
			para.chain{ii}.L				= max(length(para.chain{ii}.dataPoints),para.L-1);
			para.chain{ii}.peakWidth		= pDPMES.Results.PeakWidth;
			para.chain{ii}.w_cutoff			= max(para.chain{ii}.dataPoints(:,1))+pDPMES.Results.PeakWidth*10;			% 10*width above highest mode to be safe?
		end
	end
	
%% Setting TreeMPS for SBM
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.nNodes    = 1;
	para.treeMPS.maxDegree = para.nChains;
	para.treeMPS.leafIdx   = num2cell((1:para.nChains)'+1);			% maps from chain number to leaf index in para
	para.treeMPS.nodeIdx   = {0+1};									% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

%% Settings for DPMES-TreeN
if strfind(para.model,'DPMES-Tree1')
 	pDPMES.parse(varargin{:});
	
%% Setting Chain 1
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20160129/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.chain{1}.mapping			= 'LanczosTriDiag';
	para.chain{1}.spectralDensity	= 'CoupDiscr';
	para.chain{1}.Lambda            = 1;
	para.chain{1}.initState		    = para.initChainState;
	
	para.chain{1}.dataPoints		= cmToeV(load('DPMESdata_20160129/W24-B1-17-highv2.dat'));
	para.chain{1}.L					= min(length(para.chain{1}.dataPoints),para.L-1);
	
	%% chain 2 & more:
	para.chain{2}					= para.chain{1};		% simple copy
	para.chain{2}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-7-01.dat'));
	para.chain{2}.L					= min(length(para.chain{2}.dataPoints),para.L-1);
	para.chain{3}					= para.chain{2};
	para.chain{3}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-10-x1.dat'));
	para.chain{3}.L					= min(length(para.chain{3}.dataPoints),para.L-1);
	para.chain{4}					= para.chain{2};
	para.chain{4}.dataPoints		= cmToeV(load('DPMESdata_20160129/W23-B2-8-10.dat'));
	para.chain{4}.L					= min(length(para.chain{4}.dataPoints),para.L-1);
	para.chain{5}					= para.chain{2};
	para.chain{5}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1x.dat'));
	para.chain{5}.L					= min(length(para.chain{5}.dataPoints),para.L-1);
	para.chain{6}					= para.chain{2};
	para.chain{6}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1-x.dat'));
	para.chain{6}.L					= min(length(para.chain{6}.dataPoints),para.L-1);
	para.chain{7}					= para.chain{2};
	para.chain{7}.dataPoints		= cmToeV(load('DPMESdata_20160129/W14-A2-10-highv2.dat'));
	para.chain{7}.L					= min(length(para.chain{7}.dataPoints),para.L-1);
	
%% Setting TreeMPS structure
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	leafTreeIdx = [1,1,1,0;											% the treeIdx for each leaf; Order 1-7 as defined in genh1h2term_onesite
		           1,1,2,0;
				   1,2,0,0;
				   2,1,1,1;
				   2,1,1,2;
				   2,1,2,0;
				   2,2,0,0];
	para.treeMPS.leafIdx   = num2cell(leafTreeIdx+1);    			% maps from chain number to leaf index in para
	nodeTreeIdx = [0,0,0,0;                                         % the treeIdx of each node; order as defined by occurrence sweep
		           1,0,0,0;
				   1,1,0,0;
				   2,0,0,0;
				   2,1,0,0;
				   2,1,1,0];
	para.treeMPS.nNodes    = size(nodeTreeIdx,1);
	para.treeMPS.nodeIdx   = num2cell(nodeTreeIdx+1);				% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

if strfind(para.model,'DPMES-Tree2')
 	pDPMES.parse(varargin{:});
	
%% Setting Chain 1
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20160129/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states
	para.BroadenChain				= pDPMES.Results.BroadenChain;								% save for reference
	
	% Chains:			1				2				3				4			5				6				7
	ChainFiles = {'W24-B1-17-highv2','W44-A1-7-01','W44-A1-10-x1','W23-B2-8-10','W45-B2-9-1x','W45-B2-9-1-x','W14-A2-10-highv2'}
	%					B1				A11				A12				B21			B22				B23				A2
	for ii = 1:length(ChainFiles)
		para.chain{ii}.mapping			= 'LanczosTriDiag';
		para.chain{ii}.spectralDensity	= 'CoupDiscr';
		para.chain{ii}.Lambda           = 1;
		para.chain{ii}.initState		= para.initChainState;
		para.chain{ii}.dataPoints		= cmToeV(load(sprintf('DPMESdata_20160129/%s.dat',ChainFiles{ii})));
		para.chain{ii}.L				= min(length(para.chain{ii}.dataPoints),para.L-1);
		if pDPMES.Results.BroadenChain(ii)
			% allows individual broadening of each chain
			para.chain{ii}.spectralDensity	= 'CoupBroad';
			para.chain{ii}.discrMethod		= 'Direct';
			para.chain{ii}.discretization	= 'Linear';
			para.chain{ii}.mapping			= 'Stieltjes';
			para.chain{ii}.L				= max(length(para.chain{ii}.dataPoints),para.L-1);
			para.chain{ii}.peakWidth		= pDPMES.Results.PeakWidth;
			para.chain{ii}.w_cutoff			= max(para.chain{ii}.dataPoints(:,1))+pDPMES.Results.PeakWidth*10;			% 10*width above highest mode to be safe?
		end
	end
	
%% Setting TreeMPS structure
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	leafTreeIdx = [1,0,0,0;											% the treeIdx for each leaf; Order 1-7 as defined in genh1h2term_onesite
		           2,1,1,0;
				   2,1,2,0;
				   2,2,1,1;
				   2,2,1,2;
				   2,2,2,1;
				   2,2,2,2];
	para.treeMPS.leafIdx   = num2cell(leafTreeIdx+1);    			% maps from chain number to leaf index in para
	nodeTreeIdx = [0,0,0,0;                                         % the treeIdx of each node; order as defined by occurrence sweep
		           2,0,0,0;
				   2,1,0,0;
				   2,2,0,0;
				   2,2,1,0;
				   2,2,2,0];
	para.treeMPS.nNodes    = size(nodeTreeIdx,1);
	para.treeMPS.nodeIdx   = num2cell(nodeTreeIdx+1);				% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

if strfind(para.model,'DPMES-Tree3')
	% This is star-model obtained using fkmeans from 08/04/2017
	% clusters under weighting by coupling strength only!
 	pDPMES.parse(varargin{:});
	
%% Setting Chain for SBM
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20170408/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.BroadenChain				= pDPMES.Results.BroadenChain;								% save for reference
	
	% Chains:		1		2		3		4		5		6		7
	ChainFiles = {	'A2',	'A11',	'A12',	'B11',	'B12',	'B22',  'B23'	};	% order similar to Tree2 since this is more likely to be optimal
	
	for ii = 1:length(ChainFiles)
		para.chain{ii}.mapping			= 'LanczosTriDiag';
		para.chain{ii}.spectralDensity	= 'CoupDiscr';
		para.chain{ii}.Lambda           = 1;
		para.chain{ii}.initState		= para.initChainState;
		para.chain{ii}.dataPoints		= cmToeV(load(sprintf('DPMESdata_20170408/Modes-%s.dat',ChainFiles{ii})));
		para.chain{ii}.H1				= load(sprintf('DPMESdata_20170408/H1-%s.dat',ChainFiles{ii}));
		para.chain{ii}.L				= min(length(para.chain{ii}.dataPoints),para.L-1);
		if pDPMES.Results.BroadenChain(ii)
			% allows individual broadening of each chain
			para.chain{ii}.spectralDensity	= 'CoupBroad';
			para.chain{ii}.discrMethod		= 'Direct';
			para.chain{ii}.discretization	= 'Linear';
			para.chain{ii}.mapping			= 'Stieltjes';
			para.chain{ii}.L				= max(length(para.chain{ii}.dataPoints),para.L-1);
			para.chain{ii}.peakWidth		= pDPMES.Results.PeakWidth;
			para.chain{ii}.w_cutoff			= max(para.chain{ii}.dataPoints(:,1))+pDPMES.Results.PeakWidth*10;			% 10*width above highest mode to be safe?
		end
	end
	
%% Setting TreeMPS structure
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	leafTreeIdx = [1,0,0,0;											% the treeIdx for each leaf; Order 1-7 as defined in genh1h2term_onesite
		           2,1,1,0;
				   2,1,2,0;
				   2,2,1,1;
				   2,2,1,2;
				   2,2,2,1;
				   2,2,2,2];
	para.treeMPS.leafIdx   = num2cell(leafTreeIdx+1);    			% maps from chain number to leaf index in para
	nodeTreeIdx = [0,0,0,0;                                         % the treeIdx of each node; order as defined by occurrence sweep
		           2,0,0,0;
				   2,1,0,0;
				   2,2,0,0;
				   2,2,1,0;
				   2,2,2,0];
	para.treeMPS.nNodes    = size(nodeTreeIdx,1);
	para.treeMPS.nodeIdx   = num2cell(nodeTreeIdx+1);				% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

if strfind(para.model,'DPMES-Tree4')
	% This is star-model obtained using fkmeans from 08/04/2017
	% clusters under weighting by coupling strength only!
 	pDPMES.parse(varargin{:});
	
%% Setting Chain for SBM
	% Set Folders: 20170408: 252 modes; 20170626: 75 Modes
	para.nChains = 7;
	para.systemStates				= load('DPMESdata_20170408/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.BroadenChain				= pDPMES.Results.BroadenChain;								% save for reference
	
	% Chains:		1		2		3		4		5		6		7
	ChainFiles = {	'B11',	'B12',	'A11',	'A12',	'A2',	'B22',  'B23'	};
	
	for ii = 1:length(ChainFiles)
		para.chain{ii}.mapping			= 'LanczosTriDiag';
		para.chain{ii}.spectralDensity	= 'CoupDiscr';
		para.chain{ii}.Lambda           = 1;
		para.chain{ii}.initState		= para.initChainState;
		para.chain{ii}.dataPoints		= cmToeV(load(sprintf('DPMESdata_20170408/Modes-%s.dat',ChainFiles{ii})));
		para.chain{ii}.H1				= load(sprintf('DPMESdata_20170408/H1-%s.dat',ChainFiles{ii}));
		para.chain{ii}.L				= min(length(para.chain{ii}.dataPoints),para.L-1);
		if pDPMES.Results.BroadenChain(ii)
			% allows individual broadening of each chain
			para.chain{ii}.spectralDensity	= 'CoupBroad';
			para.chain{ii}.discrMethod		= 'Direct';
			para.chain{ii}.discretization	= 'Linear';
			para.chain{ii}.mapping			= 'Stieltjes';
			para.chain{ii}.L				= max(length(para.chain{ii}.dataPoints),para.L-1);
			para.chain{ii}.peakWidth		= pDPMES.Results.PeakWidth;
			para.chain{ii}.w_cutoff			= max(para.chain{ii}.dataPoints(:,1))+pDPMES.Results.PeakWidth*10;			% 10*width above highest mode to be safe?
		end
	end
	
%% Setting TreeMPS structure
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	leafTreeIdx = [1,0,0,0,0;										% the treeIdx for each leaf; Order 1-7 as defined in genh1h2term_onesite
		           2,1,0,0,0;
				   2,2,1,1,1;
				   2,2,1,1,2;
				   2,2,1,2,0;
				   2,2,2,1,0;
				   2,2,2,2,0];
	para.treeMPS.leafIdx   = num2cell(leafTreeIdx+1);    			% maps from chain number to leaf index in para
	nodeTreeIdx = [0,0,0,0,0;                                        % the treeIdx of each node; order as defined by occurrence sweep
		           2,0,0,0,0;
				   2,2,0,0,0;
				   2,2,1,0,0;
				   2,2,1,1,0;
				   2,2,2,0,0];
	para.treeMPS.nNodes    = size(nodeTreeIdx,1);
	para.treeMPS.nodeIdx   = num2cell(nodeTreeIdx+1);				% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end

if strfind(para.model,'DPMES-Tree6')
	% This is star-model obtained using fkmeans from 08/04/2017
	% clusters under weighting by coupling strength only!
 	pDPMES.parse(varargin{:});
	
%% Setting Chain for SBM
	% Set Folders: DPMESdata_20170720
	para.nChains = 9;
	para.systemStates				= load('DPMESdata_20170720/states.dat');					% [#state, E(eV)]
	para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+pDPMES.Results.CTShift);	% use delta as percentual shift of CT states!
	para.BroadenChain				= pDPMES.Results.BroadenChain;								% save for reference
	
	% Chains:		1		2		3		4		5		6		7		8		9
	ChainFiles = {	'B21',	'B22',	'B11',	'B12',	'B13',	'B14',  'A2',	'A11',	'A12'	};
	
	for ii = 1:length(ChainFiles)
		para.chain{ii}.mapping			= 'LanczosTriDiag';
		para.chain{ii}.spectralDensity	= 'CoupDiscr';
		para.chain{ii}.Lambda           = 1;
		para.chain{ii}.initState		= para.initChainState;
		para.chain{ii}.dataPoints		= cmToeV(load(sprintf('DPMESdata_20170720/Modes-%s.dat',ChainFiles{ii})));
		para.chain{ii}.H1				= load(sprintf('DPMESdata_20170720/H1-%s.dat',ChainFiles{ii}));
		para.chain{ii}.L				= min(length(para.chain{ii}.dataPoints),para.L-1);
		if pDPMES.Results.BroadenChain(ii)
			% allows individual broadening of each chain
			para.chain{ii}.spectralDensity	= 'CoupBroad';
			para.chain{ii}.discrMethod		= 'Direct';
			para.chain{ii}.discretization	= 'Linear';
			para.chain{ii}.mapping			= 'Stieltjes';
			para.chain{ii}.L				= max(length(para.chain{ii}.dataPoints),para.L-1);
			para.chain{ii}.peakWidth		= pDPMES.Results.PeakWidth;
			para.chain{ii}.w_cutoff			= max(para.chain{ii}.dataPoints(:,1))+pDPMES.Results.PeakWidth*10;			% 10*width above highest mode to be safe?
		end
	end
	
%% Setting TreeMPS structure
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	leafTreeIdx = [1,0,0,0,0;										% the treeIdx for each leaf; Order 1-7 as defined in genh1h2term_onesite
		           2,1,0,0,0;
				   2,2,1,1,1;
				   2,2,1,1,2;
				   2,2,1,2,1;
				   2,2,1,2,2;
				   2,2,2,1,0;
				   2,2,2,2,1;
				   2,2,2,2,2];
	para.treeMPS.leafIdx   = num2cell(leafTreeIdx+1);    			% maps from chain number to leaf index in para
	nodeTreeIdx = [0,0,0,0,0;                                        % the treeIdx of each node; order as defined by occurrence sweep
		           2,0,0,0,0;
				   2,2,0,0,0;
				   2,2,1,0,0;
				   2,2,1,1,0;
				   2,2,1,2,0;
				   2,2,2,0,0;
				   2,2,2,2,0];
	para.treeMPS.nNodes    = size(nodeTreeIdx,1);
	para.treeMPS.nodeIdx   = num2cell(nodeTreeIdx+1);				% maps from node number to nodeIdx
	
%% Set-up parameters for specific ground state preparation!
%     Initial States:			TT		LE+		LE-		CT+		CT-
% 			para.IntialState:	 1		  2		  3		  4		  5
	para.InitialState               = find(~cellfun('isempty', strfind(DPMESInitialStates,pDPMES.Results.InitialState)));
end
%% Defaults to para
para.loopmax		= 50;						% (minimizeE)
para.increasedk		= 0;						% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.
para.adjust			= 0;						% (minimizeE) Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
para.dimlock		= 0;						% (minimizeE) set to 0 will change D and d_dop adaptively
para.minDimChange	= 0.01;						% (minimizeE) sets dimlock = 1 if relative dimension change < minDimChange. (larger makes less loops)
para.timeslice      = 0;
para.sitej          = 0;

%%%% dk Expansion related %%%%
% only works together with OBB! (para.useVmat = 1)
% para.useDkExpand     = 1;						% input arg. Enable dk expansion, own algorithm. General switch
para.expandprecision = 1e-5;					% unused in dk expand?
para.hasexpanded	 = 0;						% for dk Expand
% Method 1:
para.useDkExpand1	 = 0;						% Expand dk if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.expandBelowSV   = 0.995;
% Method 2:
para.useDkExpand2	 = 1;						% Expand if wavefunction on site occupies the high energy dimensions
para.dkEx2_tail		 = 0.4;						% tail length [0 1] of occupation in Vmat to analyse
para.dkEx2_maxDev	 = 1.5;						% if std(log10(tail)) < maxDev --> no increase; Measures orders of magnitude in fluctuations of tail.
para.dkEx2_minExp	 = 13;						% if tail below this order than do not expand.

%%%% Boson shift related %%%%
% para.useshift=0;										% input arg
% only choose one of the following Methods
para.useFloShift = 0;                                   % shift all sites if trustsite > 0
para.useFloShift2 = 1;  para.FloShift2minSV = 0.995;    % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;                                   % shifts sites < trustsite
para.useEveryShift=0;                                   % shift in every loop

% para.shift = zeros(para.nChains,para.L);					% separate shift for each chain!
para.shift = {};								% separate shift for each chain!
para.relativeshift=zeros(para.nChains,para.L);
para.relativeshiftprecision=0.01;				% When the relative shift is below this value then stop shifting
% para = maxshift(para);						% TODO!

para.version = 'v79';

if isstruct(para.tdvp)
	%% Defaults to para for tdvp
	para.trustsite = para.L+1;	% add system manually here
	
	%% Defaults to para.tdvp
	para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
	para.tdvp.expvCustomTestAccuracy = 0;		% do expvCustom alongside expv for testing.
	para.tdvp.expvCustomTestAccuracyRMS = 0;	% display RMS of expvCustom from expv(); set only if para.tdvp.expvCustomTestAccuracy = 1;
	para.tdvp.expvTol = 1e-15;					% error tolerance of expv(); default: 1e-7
	para.tdvp.expvM   = 50;						% dim of Krylov subspace in expv(); default: 30
	para.tdvp.version = para.version;
end

%% Folder and File name definition
[~, name] = system('hostname');
para.hostname = strtrim(name);						% save hostname for later reference

Descr = para.version;
if ~strcmp(computer,'PCWIN64')
	Descr = sprintf('%sTCM%s',Descr,para.hostname(3:end));
end
if para.useTreeMPS
	MPStype = 'Tree-';
elseif para.useStarMPS
	MPStype = 'Star-';
elseif para.useVtens
	MPStype = 'VT-';
else
	MPStype = '';
end
Descr = [MPStype,Descr];
if isfield(para.chain{1},'s') && para.chain{1}.s ~= 1 && isfield(para,'SpinBoson')
	Descr = sprintf('%s-s%g',Descr,para.chain{1}.s);
end

if any(strfind(para.model,'DPMES'))
	para.folder = sprintf('%s-%s-%s%s-L%dCT%g%s',datestr(now,'yyyymmdd-HHMM-SS'), para.model, Descr, para.comment, para.L,pDPMES.Results.CTShift,pDPMES.Results.InitialState);
else
	para.folder = sprintf('%s-%s-%s-L%d',datestr(now,'yyyymmdd-HHMM-SS'), para.model, Descr, para.L);
end

para.filename      = [para.folder,'/results.mat'];

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%.2g%s'],para.tdvp.tmax,para.tdvp.deltaT,para.tdvp.version);
% para.tdvp.filename = sprintf('%s-alpha%g',para.tdvp.filename,alpha);

if para.tdvp.expandOBB
	para.tdvp.filename = sprintf('%s-OBBmax%d',para.tdvp.filename, para.tdvp.maxOBBDim);
end
if para.tdvp.truncateExpandBonds
	if length(para.tdvp.maxBondDim) > 1
		para.tdvp.filename = sprintf('%s-Dmax(%d-%d)',para.tdvp.filename,para.tdvp.maxBondDim(1),para.tdvp.maxBondDim(2));
	else
		para.tdvp.filename = sprintf('%s-Dmax%d',para.tdvp.filename,para.tdvp.maxBondDim);
	end
end
if para.tdvp.expvCustom
	para.tdvp.filename = sprintf('%s-expvCustom%d',para.tdvp.filename,para.tdvp.maxExpVDim);
end
% finish Filenames
para.tdvp.filename		= sprintf('%s-%dcore.mat',para.tdvp.filename,maxNumCompThreads);
para.tdvp.filenameSmall = sprintf('%s-small.mat',para.tdvp.filename(1:end-4));		% only for para, tresults
% Set MPS filename if needed
if para.tdvp.storeMPS == 1
	para.tdvp.filenameMPS = [para.tdvp.filename,'-MPS.mat'];		% only for tmps, tVmat
end

if ~exist(para.tdvp.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end


for k = 1:para.nChains
	if ~strcmp(para.chain{k}.mapping,'')
		[para.chain{k}]=SBM_genpara(para.chain{k});               % only need alpha, s, Lambda. Returns epsilon and t of Wilson chain. Auto choose L if L == 0
	end
end

[treeMPS, para] = initState(para);
results = initresults(para);
%% Do one prepare sweep to bring MPS onto right-canonical form
[treeMPS,~,para] = prepare(treeMPS,[],para);
		
%% Calculate effective Hamiltonians
treeMPS = initstorage(treeMPS,[],[],para);		% new call pattern for treeMPS

filename = para.tdvp.filename;
save(para.tdvp.filename,'para','results','treeMPS','-v7.3');

if strcmp(computer,'PCWIN64')
	% Only calculate if this is not the cluster. Otherwise needs moving to scratch
	tdvp_1site_tree(treeMPS,para,results,[]);
end

end

function [treeMPS, para] = initState(para)
%% function treeMPS = initState(treeMPS)
%
%	initialisation of the treeMPS.
%
%	Created by FS 29/02/2016
temp = allocateTreeStruct();
chains = repmat(temp,para.nChains,1);			% will contain all chains
nodes = repmat(temp,para.treeMPS.nNodes,1);		% will contain all nodes

dk    = para.dk_start;
d_opt = para.d_opt_start;
D     = para.D_start;
%%
for mc = 1:para.nChains
	%% initialise each chain
	L = para.chain{mc}.L;
	
	chains(mc).mps      = cell(1,L);
	chains(mc).Vmat     = cell(1,L);
	chains(mc).L        = L;
	chains(mc).treeIdx  = [para.treeMPS.leafIdx{mc,:}]-1;
	chains(mc).chainIdx = mc;
	chains(mc).level    = nnz(chains(mc).treeIdx)+1;
	chains(mc).D        = [ones(1,L)*D,1];
	chains(mc).d_opt    = ones(1,L)*d_opt;
	chains(mc).dk       = ones(1,L)*dk;
	chains(mc).shift    = zeros(1,L);
	
	%% Fill the MPS of each chain
	% default: vacuum in the chains
	for ii = 1:L
		switch para.chain{mc}.initState
			case 'rand'
				chains(mc).mps{ii}        = randn([chains(mc).D(ii:ii+1),d_opt]);			% end-of-chain -> Dr = 1
			case 'vac'
				chains(mc).mps{ii}        = zeros([chains(mc).D(ii:ii+1),d_opt]);			% end-of-chain -> Dr = 1
				chains(mc).mps{ii}(1,1,1) = 1;
		end
		chains(mc).Vmat{ii} = sparse(1:d_opt,1:d_opt,1,dk,d_opt);
	end
	
	%% Set in para
	pIdx = para.treeMPS.leafIdx(mc,:);
	para.dk{pIdx{:}}        = chains(mc).dk;
	para.D{pIdx{:}}         = chains(mc).D;
	para.d_opt{pIdx{:}}     = chains(mc).d_opt;
	para.shift{pIdx{:}}     = chains(mc).shift;
	para.treeMPS.L(pIdx{:}) = chains(mc).L;
	
	para.treeMPS.chainIdx{para.treeMPS.leafIdx{mc,:}} = mc;
	%% generate Hamiltonian terms:
	chains(mc).op.h1term = cell(1,L);
	chains(mc).op.h2term = cell(para.M,2,L);
	for ii = 1:L
		temp = genh1h2term_onesite(para,chains(mc).treeIdx,ii);
		chains(mc).op.h1term(ii)     = temp.h1term;
		chains(mc).op.h2term(:,:,ii) = temp.h2term;
	end

end
leafIndices = cell2mat(para.treeMPS.leafIdx)-1;						% nLeaves x nLevels
nodeIndices = cell2mat(para.treeMPS.nodeIdx)-1;						% nNodes x nLevels

para.L = max([chains.level]+[chains.L]-1);
%%
for mc = para.treeMPS.nNodes:-1:1
	L = 1;								% for now: only single nodes
	pIdx = para.treeMPS.nodeIdx(mc,:);
	nodes(mc).mps      = cell(1,L);
	nodes(mc).Vmat     = cell(1,L);
	nodes(mc).L        = L;
	nodes(mc).treeIdx  = [para.treeMPS.nodeIdx{mc,:}]-1;
% 	nodes(mc).chainIdx = mc;								% dummy  value?
	nodes(mc).level    = nnz(nodes(mc).treeIdx)+1;
% 	nodes(mc).D        = [ones(1,L)*D,1];					% assign after children
% 	nodes(mc).d_opt    = ones(1,L)*d_opt;					% assign after genh1h2term_onesite, derive from H
% 	nodes(mc).dk       = ones(1,L)*dk;
	nodes(mc).shift    = zeros(1,L);
	if nodes(mc).level == 1
		nodes(mc).isRoot = 1;
	end
	nodes(mc).useStarMPS = 0;
	nodes(mc).useTreeMPS = 1;
	nodes(mc).spinposition = 1;
	
	
	%% Attach children
	currentChild = 0;
	hasChild = 1;
	while hasChild
		currentChild = currentChild + 1;								% each node has at least 1 child
		% construct index of child to look for
		idx                  = nodes(mc).treeIdx;
		idx(nodes(mc).level) = currentChild;			% look for [ 4 3 3 i 0 0 0] if at current node is at level 4 of 8
		leafPos              = find(all(bsxfun(@eq,leafIndices,idx),2));
		nodePos              = find(all(bsxfun(@eq,nodeIndices,idx),2));
		if ~isempty(leafPos) && leafPos
			if isempty(nodes(mc).child)
				nodes(mc).child                 = chains(leafPos);
			else
				nodes(mc).child(currentChild,1) = chains(leafPos);
			end
		elseif ~isempty(nodePos) && nodePos
			if isempty(nodes(mc).child)
				nodes(mc).child                 = nodes(nodePos);
			else
				nodes(mc).child(currentChild,1) = nodes(nodePos);
			end
		else
			hasChild = 0;
		end
	end
	
	%% Set parameters
	nodes(mc).op       = genh1h2term_onesite(para,nodes(mc).treeIdx,1);
	nodes(mc).degree   = currentChild-1;
	nodes(mc).chainIdx = [nodes(mc).child.chainIdx];					% min(chainIdx) is chainIdx of current node
	nodes(mc).nChains  = length(nodes(mc).chainIdx);
	nodes(mc).height   = max([nodes(mc).child.height])+1;
	nodes(mc).L        = getTreeLength(nodes(mc));
	nodes(mc).D        = ones(nodes(mc).degree+1,1)*D;
	if nodes(mc).isRoot
		nodes(mc).D(1) = 1;
	end
	if isempty(nodes(mc).op.h1term{1})
		nodes(mc).hasSite = 0;
		nodes(mc).dk    = 1;			% leads to singleton
		nodes(mc).d_opt = 1;
	else
		nodes(mc).hasSite = 1;
		nodes(mc).chainIdx = min(nodes(mc).chainIdx);					% collapse chain index to indicate, that node has only one coupling to right
		nodes(mc).dk    = size(nodes(mc).op.h1term{1},1);
		nodes(mc).d_opt = nodes(mc).dk;									% assume non-boson
	end
	
	para.treeMPS.L(pIdx{:}) = 1;
	para.shift{pIdx{:}}     = 0;
	%% Fill the MPS
	dk    = nodes(mc).dk;
	d_opt = nodes(mc).d_opt;
	
	NC = nodes(mc).degree;
	nodes(mc).mps{1}        = zeros([nodes(mc).D',d_opt]);
	if ~nodes(mc).isRoot
		idx = num2cell(ones(length(nodes(mc).D),1));            % indexes all left and child bonds;
		nodes(mc).mps{1}(idx{:},1) = 1;                         % initially connect as product states
	elseif any(strfind(para.model,'DPMES'))
		idx = num2cell([ones(1,NC+1),para.InitialState]);			% start in second excited state!
		nodes(mc).mps{1}(idx{:}) = 1;
	elseif isfield(para,'SpinBoson') && strcmp(para.SpinBoson.GroundStateMode, 'artificial')
		d_opt = 2;
		if strcmp(para.SpinBoson.InitialState, 'sz')
			%% prepare +Sz eigenstate
			idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
			nodes(mc).mps{1}(idx{:},1) = 1;
		elseif strcmp(para.SpinBoson.InitialState, '-sz')
			idx = num2cell(ones(1,NC+1));
			nodes(mc).mps{1}(idx{:},2) = 1;
		elseif strcmp(para.SpinBoson.InitialState, 'sx')
			idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
			nodes(mc).mps{1}(idx{:},1) = 1/sqrt(2);
			nodes(mc).mps{1}(idx{:},2) = 1/sqrt(2);
		elseif strcmp(para.SpinBoson.InitialState, '-sx')
			idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
			nodes(mc).mps{1}(idx{:},1) = -1/sqrt(2);
			nodes(mc).mps{1}(idx{:},2) = 1/sqrt(2);
		end
	else
		error('VMPS:VMPS_TreeMPS:WrongParameters','');
	end
	nodes(mc).Vmat{1}       = sparse(1:d_opt,1:d_opt,1,dk,d_opt);
	
end

treeMPS = nodes(1);

	function treeMPS = allocateTreeStruct()
	%% function treeMPS = allocateTreeStruct()
	%
	%	initialises an empty struct template and defines all used fields
	%	Already fill in defaults for leaves
		treeMPS = struct;
		treeMPS.mps = {};
		treeMPS.Vmat = {};
		treeMPS.treeIdx = [];
		treeMPS.chainIdx = 0;								% = [# chains connected to]
		treeMPS.height = 0;									% this struct is a leaf
		treeMPS.degree = 0;									% = n(children)
		treeMPS.level = 1;									% = site number of first chain site if counted from the root node
		treeMPS.child = [];									% contains all children nodes
		treeMPS.currentChild = [];							% leaf: [], node: currentChain;
		treeMPS.currentChain = [];							% leaf: [], node: currentChain;
		treeMPS.L = [];										% StarMPS can have array L
		treeMPS.D = 0;										% first is Dl of chain, last is Dr of chain, here horizontal!
		treeMPS.d_opt = 0;
		treeMPS.dk = 0;
		treeMPS.M = para.M;
		treeMPS.shift = 0;
		treeMPS.useVmat = para.useVmat;
		treeMPS.useVtens = 0;
		treeMPS.useStarMPS = 0;
		treeMPS.useTreeMPS = 0;
		treeMPS.nChains = 1;					% n(chains connected through node) = length(chainIdx) >= degree
		treeMPS.sweepto = 'l';					% 'l' for initial preparation
		treeMPS.spinposition = [];				% pure Boson Chain
		treeMPS.sitej = 1;						% currently focused site
		treeMPS.isRoot = 0;
		treeMPS.BondCenter = [];
		treeMPS.tdvp = [];
		treeMPS.op.h1term = {};
		treeMPS.op.h2term = {};
		treeMPS.hasSite   = 1;					% whether node has real site; no Site nodes need: op.h1term{1} = []
		treeMPS.Amat_sv = {};					% Singular values for each bond within node/leaf
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
			d(3) = 1;				% need singleton to allow comparisons later
		end
		Lchild = arrayfun(@getTreeLength, treeMPS.child,'UniformOutput',false);
		% find out maximum dimensions per level to zero-pad Lchild
		Idx = find([treeMPS.child.height] == h-1);		% only pick the ones with maximum height
		for jj = Idx
			d(2:end) = max(d(2:end),size(Lchild{jj}));
		end
		% now d contains the maximum dimensions of L -> initialise
		L = zeros(d);
		L(1) = 1;					% length of nodes = 1 always for now!
		for jj = 1:treeMPS.degree
			Idx = arrayfun(@(x) 1:x,size(Lchild{jj}),'UniformOutput',false);		% create cell array containing index ranges of Lchild{ii}
			L(jj+1,Idx{:}) = Lchild{jj};
		end
	end

end

function dummy()
	%% This function just contains cells which can be executed to quickly run a certain model!

	%% DPMES 5-7C
	%	full chains, no truncation!
	VMPS_TreeMPS('model','DPMES5-7C','evolveSysTrotter',1,'L',18,'tmax',10,'deltaT',0.1,'extractObsInterval',0.1,'InitialState','LE+','useDkExpand',1,'maxBondDim',[5,10],...
		'Observables','.n.na.nd.x.xa.xd.x2.x2a.x2d.dm.ss.');
	% L: longest L of from root to all leaves

	%% DPMES Tree1
	%	full chains, no truncation!
	VMPS_TreeMPS('model','DPMES-Tree1','evolveSysTrotter',1,'L',18,'tmax',10,'deltaT',0.1,'extractObsInterval',0.1,'InitialState','LE+','useDkExpand',1,'maxBondDim',[10],...
		'Observables','.n.dm.');
	% L: longest L of from root to all leaves

	%% SpinBoson2C
	VMPS_TreeMPS('model','SpinBoson2C','evolveSysTrotter',0,'L',49,'deltaT',0.2,'wc',[1,5],'InitialState','-sx','extractObsInterval',0.2,'useDkExpand',1)
	
	%% Spin Boson, SigmaZ splitting, SigmaX coupling to chains
	% coupling operator has to be changed explicitly!
	VMPS_TreeMPS('model','SpinBoson',...
			 's',1,...
			 'alpha',4.8,...
			 'delta',0,...
			 'epsilon',8.42,...
			 'L',1,...
			 'Observables','.n.nd.na.x.xd.xa.dm.heff.pes.',...
			 'InitialState','-sz',...
			 'deltaT',0.02,...
			 'extractObsInterval',0.02,...
			 'tmax',10,...
			 'useDkExpand',1,...
			 'evolveSysTrotter',0,...
			 'maxOBBDim',40,...
			 'maxBondDim',[50,20],...
			 'D_start',2,...
			 'saveInterval',20)
	%% Spin Boson, SigmaX splitting, SigmaX tunneling and SigmaZ coupling to Chain
	% coupling operator has to be changed explicitly!
	VMPS_TreeMPS('model','SpinBoson',...
			 's',1,...
			 'alpha',0.65,...
			 'delta',0.2,...
			 'epsilon',0.5,...
			 'L',1,...
			 'Observables','.n.nd.na.x.xd.xa.dm.heff.pes.',...
			 'InitialState','sz',...
			 'deltaT',0.2,...
			 'extractObsInterval',0.2,...
			 'tmax',40,...
			 'useDkExpand',1,...
			 'evolveSysTrotter',0,...
			 'maxOBBDim',40,...
			 'maxBondDim',[50,20],...
			 'D_start',2,...
			 'serialize',0,...
			 'saveInterval',20)
end
function fileName = VMPS_FullSBM(s,alpha,delta,epsilon,L,dk,d_opt,D)
%Variational matrix product method to study spin-boson model. Rok's
%logrithimic discretization and optimal boson basis are implemented in
%this code.
% adapted to use the SBM as benchmark system

% modified to support a multi-level spin site with choosable coupling to the bath.


%Cheng Guo @Munich
% 3 Sep 2010


% mod by Florian Schroeder @Cambridge
% 23 October 2014 -
% copied from MLSBM 23/10/14


% maxNumCompThreads(1);
format short e
rng('shuffle');

starttime = tic;
if isdeployed           % take care of command line arguments
% 	if ischar(hx), hx = str2num(hx); end
% 	if ischar(hz), hz = str2num(hz); end
	if ischar(s), s = str2double(s); end
	if ischar(alpha), alpha = str2double(alpha); end
	if ischar(delta), delta = str2double(delta); end
	if ischar(epsilon), epsilon = str2double(epsilon); end
	if ischar(L), L = str2double(L); end
	if ischar(dk), dk = str2double(dk); end
	if ischar(d_opt), d_opt = str2double(d_opt); end
	if ischar(D), D = str2double(D); end
% 	if ischar(Lambda), Lambda = str2num(Lambda); end
% 	if ischar(parity), parity = str2num(parity); end
end

%% Choose model and chain mapping
para.model='MLSpinBoson';
    % choose: 'SpinBoson', 'SpinBoson2folded', 'MLSpinBoson','ImpurityQTN'
	%         '2SpinPhononModel', 'SpinBoson2C', 'SpinBosonTTM', 'SpinBoson2CT'
	%		  'SpinBoson2CXZ',
	%		  'DPMES3-4C', 'DPMES4-5C', 'DPMES5-7C'
	%		  'UniformBosonTTM'
% para.chainMapping = 'OrthogonalPolynomials';
para.nEnvironments  = 1;
	% number of different spectral functions
	% supported 1 to any
para.nChains		= 1;
	% number of chains
	% 1 for folded, can have nEnvironments = 2;
	% = nEnvironments for multi-chain models;
para.useVtens = 0;										% Enables the V-tensor-network for MultiChain models! Only Artificial State!

para.useStarMPS = 0;
para.useTreeMPS = 0;

%% System Definitions:
if ~any(strcmp(para.model,{'MLSpinBoson','2SpinPhononModel'}))
	% setting para for single-spin models
	% H0 = -para.hx./2.*sigmaX -para.hz./2.*sigmaZ
    para.hx = -delta;                       % Splitting with sigma_X
    para.hz = -epsilon;                     % Splitting with sigma_Z
end

%% Chain Definitions:
% para.chain{i}.mapping
	% choose: 'OrthogonalPolynomials','LanczosTriDiag', 'Stieltjes'
	%			- 'LanczosTriDiag' Lanczos tridiagonalization
	%			- 'Stieltjes': needs discretization, based on Orthogonal
	%						   Polynomials, less stable but faster than Lanczos
	%			- 'OrthogonalPolynomials'
% para.chain{i}.spectralDensity
	% choose: 'Leggett_Hard', 'Leggett_Soft', 'Renger'
	% Leggett for SBM;  see DOI: 10.1103/PhysRevLett.91.170601
	% Renger for MLSBM; see DOI: 10.1016/j.bbabio.2012.02.016
% para.chain{i}.discretization
	% choose: 'None', 'Linear', 'LogZ'
	%	'None': exact mapping
	%	other methods not yet supported
para.L = L;

if ~isempty(strfind(para.model,'DPMES'))
	if strcmp(para.model,'DPMES4-5C')
		%% chain 1 for DPMES:		ideally should have been struct array: para.chain(1).mapping ... this saves more memory and allows more operations
		para.systemStates				= load('DPMESdata_20151123/states.dat');		% [#state, E(eV)]
		para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+delta);			% use delta as percentual shift of CT states!
		para.InitialState				= 2;											% 1: TT, 2: LE+, 3: LE-, 4: CT+, 5: CT-
		para.chain{1}.mapping			= 'LanczosTriDiag';
		para.chain{1}.spectralDensity	= 'CoupDiscr';
		para.chain{1}.dataPoints		= cmToeV(load('DPMESdata_20151123/W44-A1-7-01.dat'));

		para.chain{1}.L					= min(length(para.chain{1}.dataPoints)+1,L);

		%% chain 2 & more:
		para.chain{2}					= para.chain{1};		% simple copy
		para.chain{2}.dataPoints		= cmToeV(load('DPMESdata_20151123/W44-A1-10-x1.dat'));
		para.chain{3}					= para.chain{2};
		para.chain{3}.dataPoints		= cmToeV(load('DPMESdata_20151123/W24-B1-highv2.dat'));
		para.chain{4}					= para.chain{2};
		para.chain{4}.dataPoints		= cmToeV(load('DPMESdata_20151123/W14-A2-highv2.dat'));
		para.chain{5}					= para.chain{2};
		para.chain{5}.dataPoints		= cmToeV(load('DPMESdata_20151123/W45-B2-all.dat'));
	elseif strcmp(para.model,'DPMES5-7C')
		% Change in DPMESdata_20160129: split B2 into more chains, 2 slow weak modes more in B1
		para.systemStates				= load('DPMESdata_20160129/states.dat');		% [#state, E(eV)]
		para.systemStates([4,5],2)      = para.systemStates([4,5],2)*(1+delta);			% use delta as percentual shift of CT states!
		para.InitialState				= floor(epsilon);								% 1: TT, 2: LE+, 3: LE-, 4: CT+, 5: CT-; only allow integer values>0!
		para.chain{1}.mapping			= 'LanczosTriDiag';
		para.chain{1}.spectralDensity	= 'CoupDiscr';
		para.chain{1}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-7-01.dat'));

		para.chain{1}.L					= min(length(para.chain{1}.dataPoints)+1,L);	% since all others are longer

		%% chain 2 & more:
		para.chain{2}					= para.chain{1};		% simple copy
		para.chain{2}.dataPoints		= cmToeV(load('DPMESdata_20160129/W44-A1-10-x1.dat'));
		para.chain{3}					= para.chain{2};
		para.chain{3}.dataPoints		= cmToeV(load('DPMESdata_20160129/W14-A2-10-highv2.dat'));
		para.chain{4}					= para.chain{2};
		para.chain{4}.dataPoints		= cmToeV(load('DPMESdata_20160129/W24-B1-17-highv2.dat'));
		para.chain{5}					= para.chain{2};
		para.chain{5}.dataPoints		= cmToeV(load('DPMESdata_20160129/W23-B2-8-10.dat'));
		para.chain{6}					= para.chain{2};
		para.chain{6}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1x.dat'));
		para.chain{7}					= para.chain{2};
		para.chain{7}.dataPoints		= cmToeV(load('DPMESdata_20160129/W45-B2-9-1-x.dat'));
	end
elseif ~isempty(strfind(para.model,'SpinBoson'))
	%% chain 1 for SBM:
	para.chain{1}.mapping			= 'OrthogonalPolynomials';
	para.chain{1}.spectralDensity	= 'Leggett_Hard';
	para.chain{1}.discretization	= 'None';
	% para.chain{1}.discrMethod		= 'Numeric';

	para.chain{1}.s					= s;			% SBM spectral function power law behaviour
	para.chain{1}.alpha				= alpha;		% SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
	para.chain{1}.L					= L;
	para.chain{1}.w_cutoff          = 1;
	if alpha == 0 && para.chain{1}.L == 0                  
		para.chain{1}.L = 10;						% otherwise encounter error
	end
% 	para.chain{2} = para.chain{1};
% 	para.chain{2}.w_cutoff = 5;
elseif ~isempty(strfind(para.model,'UniformBosonTTM'))
	% put in parameters by hand!
	para.chain{1}.epsilon = 0.5;
	para.chain{1}.t       = 0.25;
end

assert(para.nEnvironments == length(para.chain),'number of environments is wrong');		% redundant, sanity check!
%% Parameters
for k = 1:para.nEnvironments
	% add random static disorder to Bosons: (comment if not wanted!)
% 	para.chain{k}.dataPoints = arrayfun(@(x) (randn*0.1+1)*x, para.chain{k}.dataPoints);		% 10% static disorder (SDV)

	if strcmp(para.chain{k}.mapping,'OrthogonalPolynomials')
		%% now only for Leggett power law spectral densities.
		% For arbitrary J(w) use 'Stieltjes'
		% no need to define:
		% z, Lambda
		% since site energies converge to w_c/2 = 0.5, Optimum chain length can
		% not easily be determined -> give para.L
		% based on analytic results only! For numeric approach use Stieltjes

		if para.chain{k}.L == 0					% chain length override
			para.chain{k}.L = 200;				% default chain length if input L=0
		end
		para.rescaling = 0;						% only for LogZ discretization applicable, rescaling might be broken now!
		para.chain{k}.discretization = 'None';	% Only for continuous mapping
		
	elseif any(strcmp(para.chain{k}.mapping,{'LanczosTriDiag','Stieltjes'}))
		% 	if Lambda > 1 -> LogZ
		% 	if Lambda = 1 -> Linear discretization, rescaling = 0
		para.chain{k}.discrMethod = 'None';
		% choose: 'Analytic', 'Numerical', 'None', 'Direct'
		%	Sets way of evaluation of integrals
		%	'None' takes values of J(w); should be used with 'Linear'
		%	Analytic only for 'Leggett_hard' and 'LogZ'. Uses modified scheme by Žitko
		%
		para.chain{k}.discretization = 'Linear';
		% choose: 'LogZ','Linear'

		%%
		para.chain{k}.Lambda=1;							% Bath Discretization parameter
		para.chain{k}.z=1;                               % z-shift of bath; see Zitko 2009 - 10.1103/PhysRevB.79.085106
% 		para.chain{k}.L=0;                               % Length per bath; if L=0: choose optimal chain length according to para.precision;
% 		if L > 0								% chain length override
% 			para.chain{k}.L = L;
% 		end
		
	% Consistency checks
		if strcmp(para.chain{k}.discrMethod,'Analytic')
			assert(strcmp(para.chain{k}.spectralDensity,'Leggett_Hard'),'Analytic discretization only for Leggett with hard cutoff!');
		end
		if para.chain{k}.Lambda == 1
			para.chain{k}.discretization = 'Linear';
			para.rescaling = 0;                     % rescale h1term, h2term for bosonchain with \lambda^{j-2}*h1term
		elseif para.chain{k}.Lambda > 1
			para.chain{k}.discretization = 'LogZ';
			para.rescaling = 1;                     % rescale h1term, h2term for bosonchain with \lambda^{j-2}*h1term, yields faster convergence for ground state optimization. Do not use for TDVP!
		else
			error('VMPS:ModelDefinition','Lambda >= 1!')
		end
	elseif strcmp(para.chain{k}.mapping,'')
		para.rescaling = 0;
		fprintf('\nChain %d: No mapping defined -> no call of SBM_genpara',k);
	else
		error('VMPS:ModelDefinition','You need to define para.chain.mapping!');
	end
end

%% Starting MPS Dimensions
% D = 5;
% dk = 15;
% d_opt = 5;

if strcmp(para.model,'MLSpinBoson')     % definitions needed in SBM_genpara for spectral function & Wilson chain
    % Phonon Model Definition para.MLSB_mode:
    %   1:  from diagonalised, constant spacing Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
	%	3:  Define H0 and H1 as matrix from file input
	%		Needs Define: MLSB_H0, MLSB_H1
	if ~isempty(strfind(para.chain{1}.spectralDensity,'Leggett'))
		para.MLSB_mode = 3;
	elseif ~isempty(strfind(para.chain{1}.spectralDensity,'Renger'))
		para.MLSB_mode = 2;
	end
end

if isempty(strfind(para.model,'folded'))
	para.foldedChain = 0;                       % parameter to tell that Supersites for chain are used!
	para.M = 2;									% Number of 2-operator interaction terms per site in Hamiltonian.
else
	para.foldedChain = 1;
	para.M = 2*2;
	dk = dk^2;
end

para.spinposition = [1];						% This indicates all positions ~= bosonic! important for Vmat! The y chain is on the left and the z chain is on the right. (could be array !)
para.complex=0;                                 % set to 1 if any complex parameters are used.
para.resume=0;                                  % Read from saved results if available.
para.logging = 1;                               % Switch on logging and
parity = 0;										% 0: none; 1: odd, 2: even
para.precision = 5e-15;                         % was 5e-15; Determines chain length if L=0; Also E-convergence


%% %%%%%%% Calculate Wilson Chain parameters %%%%%%%%%%%%%%%%%%
% needed here: para.model, [para.MLSB_mode]
for k = 1:para.nEnvironments
	if ~strcmp(para.chain{k}.mapping,'')
		[para.chain{k}]=SBM_genpara(para.chain{k});               % only need alpha, s, Lambda. Returns epsilon and t of Wilson chain. Auto choose L if L == 0
		if (L == 0)
			L = para.chain{k}.L;                         % Take best L if not specially defined
			para.L = para.chain{k}.L;					 % not the best solution.
		end
		para.L = min(para.L, para.chain{k}.L);
	end
end
L = para.L;

%%
nc = para.nChains;											% = 1 for single chain models
para.D							= D*ones(1,L-1);			% Bond dimension; starting dimension is 2. para.D(L) is useless in this program
para.dk_start					= dk;						% local dimension per boson in bath. Will be increased effectively by oscillator shift.
para.dk							= para.dk_start*ones(nc,L);
para.dk(1,para.spinposition)	= 2;						% Impurity dimension
para.d_opt						= d_opt*ones(1,L);			% Dimension of first site is 2 (spin); Optimal Boson Basis dimension, was 16*ones
if para.useVtens
	para.d_opt = d_opt*ones(nc+1,L);
end
para.d_opt(1,para.spinposition) = 2;						% Optimal Impurity dimension
para.eigs_tol					= 1e-8;
para.loopmax					= 50;
para.increasedk					= 0;						% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.

if para.useStarMPS
	para.D     = repmat(para.D,nc,1);
	para.d_opt = repmat(para.d_opt,nc,1);					% copy for each chain
	para.d_opt(2:end,1) = 1;								% Central System only on first site!
end

if strcmp(para.model, 'SpinBosonTTM')
	para.D([1,2]) = [2,4];
end

if strcmp(para.model, 'SpinBoson2CT')
	para.d_opt(1:end-1,2:L) = para.dk(1:end,2:L);			% to allow maximal entanglement
end

if ~isempty(regexp(para.model,'SpinBoson\dC','once'))
	para.M = 2*para.nChains;				% only nearest neighbour for now!
	para.dk(2:end,para.spinposition) = 1;	% non-existent singleton!
	para.d_opt(2:end,para.spinposition) = 2;
end

if strcmp(para.model,'2SpinPhononModel')
   %para.hz = para.hx;                  % use hz as level splitting. hx = hy = 0
   %para.hx = 0;
   para.hy = 0.1;                       % coupling between the 2 spin sites = J
   para.M = 4;                          % 4 terms in sum per site
   para.dk(1) = 4;                      % As kron(site1,site2), dim=4 on first site!
   para.d_opt(1) = 4;
   para.foldedChain=1;
end

if strcmp(para.model,'DPMES3-4C')
	para.M = 4*2;
	para.dk(1,para.spinposition)	= 3;
	para.dk(2:end,para.spinposition) = 1;	% non-existent singleton!
	para.d_opt(1:end,para.spinposition) = 3;
% 	para.dk(3,2)     = 500;
elseif strcmp(para.model,'DPMES4-5C')
	para.M = 5*2;
	para.dk(1,para.spinposition)	= 4;
	para.dk(2:end,para.spinposition) = 1;	% non-existent singleton!
	para.d_opt(1:end,para.spinposition) = 4;
elseif strcmp(para.model,'DPMES5-7C')
	para.M = 7*2;
	para.dk(1,para.spinposition)	= 5;
	para.dk(2:end,para.spinposition) = 1;	% non-existent singleton!
	para.d_opt(1:end,para.spinposition) = 5;
end

%% Multi-Level Spin Boson Model for PPC

if strcmp(para.model,'MLSpinBoson')
%% Model Definition para.MLSB_mode:
    %   1:  from diagonalised, define spacings Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
    %       Can be affected with static random disorder
	%	3:  Define H0 and H1 as matrix from file input

% All modes:
% System-Bath coupling:
    
	if para.MLSB_mode ~= 3
% 		para.MLSB_t = [-1 -1 -1 1 1 1];                 % coupling of each level to bath from range [-1,1]. Used as para.t(1)*para.MLSB_t
		para.MLSB_p = period;                           % period of coupling for cosine
		para.MLSB_etaFactor = eta;                      % multiplier to increase system-bath coupling
		para.MLSB_t = @(x) para.MLSB_etaFactor.*exp(1i*2*pi/para.MLSB_p.*x);                % MLSB_t can be an anonymous function.
		para.MLSB_tOff = 1;                             % 0=diag coupling; 1 = 1. off-diagonal n,n+1 coupling
	end
	
    if para.MLSB_mode == 1
        para.MLSB_Delta = para.hx;                  % level spacings in already diagonalized model
		para.MLSB_Ls = length(para.hx)+1;           % number of levels in System, symmetric around E=0;
    end

    if para.MLSB_mode == 2
        para.MLSB_system = 'RsMolischianumB850';
            % Choose: 'RsMolischianumB850', 'RsMolischianumB800B850'
        para.MLSB_Ls = length(Hamiltonian_PPC(para));
            % Inherently defined within model
	end
	
	if para.MLSB_mode == 3
		para.MLSB_H0 = load('DPMESdata_20170408/H1-A2.dat');			% TODO: needs file input
% 		para.MLSB_H1 = load();						% TODO: needs file input, if most general coupling matrix given
		para.MLSB_H1 = diag(beta);					% e.g.: beta = [-1, 1, 3];
		para.MLSB_Ls = size(para.MLSB_H0,1);
		assert(para.MLSB_Ls == size(para.MLSB_H1,1), 'Dimension mismatch of H1 and H0');
	end

	para.MLSB_InitialState = 1;				% pick which state should be occupied initially; 0: do Ground State optimisation
	
%% random static disorder:
    para.MLSB_staticDisorder = 0;			% 0 or 1
    para.MLSB_disSDV = 200/2.3548;         % standard deviation of disorder; FWHM = 200 cm^-1; FWHM = 2.3548*SDV;
    % Vector with absolute values of disorder in cm^-1:
    para.MLSB_disDiag = para.MLSB_disSDV.*randn(para.MLSB_Ls,1);        % set to 0 if no diag disorder wanted

%% All modes:
    para.dk(1) = para.MLSB_Ls;                  % also set dk properly
    para.d_opt(1) = para.MLSB_Ls;

end

if strfind(para.model,'SpinBoson')
%% Set-up parameters for specific ground state preparation!
    para.SpinBoson.GroundStateMode = 'artificial';
        % choose: 'decoupled', 'coupled', 'artificial', 'artTTM'
		% -artificial & artTTM does no optimization! this only sets up an artificial
		%		ground state with <n> = 0 on chain and InitialState 'sz'
		% - 'artificial' for SBM2CT sets up maximally entangled state between odd and even chains in the Bath
		%	odd: thermal bath, even: ancilla
    para.SpinBoson.InitialState = 'sz';
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

%%
para.SVDmethod	  = 'qr';					% 'qr': uses QR instead of SVD wherever possible; 'svd': use SVD always (slower) (Not working now)
% smallest SV shall lie between [max, min] otherwise truncate or expand
% the smaller the higher accuracy
para.svmaxtol	  = 1e-6;
para.svmintol	  = 1e-8;					% para.svmaxtol/2; %The lower limit for the smallest Vmat singular values.
para.adjust		  = 0;						% Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
para.Dmin		  = 4;						% set a minimum Bond dimension.
para.dimlock	  = 0;						% set to 0 will change D and d_dop adaptively
para.minDimChange = 0.01;					% sets dimlock = 1 if relative dimension change < minDimChange. (larger makes less loops)

%% Parity settings
switch parity
    case 0
        para.parity='n';	% none
    case 1
        para.parity='o';	% odd
    case 2
        para.parity='e';	% even
end
para.spinbase='Z';
if para.parity~='n'
    para.Anzi=cell(para.L,1);
    para.Vnzi=cell(para.L,1);
    para.Anzilr=para.Anzi;
    para.Anzirl=para.Anzi;
end

%% %%%%%%%%%%%%%%%%%% OBB - Vmat related parameters %%%%%%%%%%%%%%%%%%%%
% Introduces Optimized Boson Basis (OBB)
para.useVmat=1;											% General switch for the use of Vmat or Vtens
if para.useVmat==0
    % then: d_opt = dk
    fprintf('Not using Vmat and OBB!\n');
    para.d_opt = para.dk;
     assert(para.dk_start==max(para.d_opt));
end
para.d_opt_min = 2;                                     % minimum d_opt dimension
if para.useVtens
	assert(para.nChains >= 2 && para.useVmat == 1, 'Use V Tensor only for more than 2 chains, and enable Vmat');
% 	para.d_opt = repmat(para.d_opt,para.nChains+1,1); done above!
end

%% %%%%%%%%%%%%%%%%%% dk Expansion - related parameters %%%%%%%
% only works together with OBB!
para.useDkExpand     = 1;		% Enable dk expansion, own algorithm. General switch

para.dkmax			 = 900;		% has to be square number for folded chains
para.expandprecision = 1e-5;	% unused?
para.hasexpanded     = 0;
% Method 1:
% needs para.useVmat = 1
para.useDkExpand1	 = 0;		% Expand dk if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.expandBelowSV   = 0.995;
% Method 2:
% needs para.useVmat = 1
para.useDkExpand2	 = 1;       % Expand if wavefunction on site occupies the high energy dimensions
para.dkEx2_tail		 = 0.4;     % tail length [0 1] of occupation in Vmat to analyse
para.dkEx2_maxDev	 = 1.5;     % if std(log10(tail)) < maxDev --> no increase; Measures orders of magnitude in fluctuations of tail.
para.dkEx2_minExp	 = 13;      % if tail below this order than do not expand.
% Method 3:
% TODO:
%  analyse based on Amat (in case para.useVmat == 0)
%  Do the SVD in the local bond to create a Vmat and use it to estimate
%  tails.
if para.useDkExpand
	assert(para.useVmat == 1, 'dk expansion only works together with OBB Vmat');
end

%% %%%%%%%%%%%%%%%%%% Boson Shift - related parameters %%%%%%%%%%%%%%%%
% Introduces shift of bosonic oscillators to increase effective dk
para.useshift=0;
% only choose one of the following Methods
para.useFloShift = 0;                                   % shift all sites if trustsite > 0
para.useFloShift2 = 1;  para.FloShift2minSV = 0.995;    % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;                                   % shifts sites < trustsite
para.useEveryShift=0;                                   % shift in every loop

para.shift=zeros(nc,para.L);							% separate shift for each chain!
para.relativeshift=zeros(nc,para.L);
para.relativeshiftprecision=0.01;        %When the relative shift is below this value then stop shifting
para=maxshift(para);

%% calculate shift analytically  (for independent sub-ohmic model)
%t = para.t; e = para.epsilon;
%A = gallery('tridiag',para.t(2:end),para.epsilon(1:end),para.t(2:end));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
%B = zeros(para.L-1,1);
%B(1) = -para.t(1)*1;
%para.shift =[0; A\B.*sqrt(2)]';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, name] = system('hostname');
para.hostname = strtrim(name);						% save hostname for later reference
para.version = 'v77';
Descr = para.version;
if ~strcmp(computer,'PCWIN64')
	Descr = sprintf('%sTCM%s',para.version,para.hostname(3:end));
end
if strcmp(para.chain{1}.mapping,'OrthogonalPolynomials')
	Descr = ['OrthPol-',Descr];
elseif strcmp(para.chain{1}.discretization,'LogZ')
	Descr = ['LogZ-',Descr];
end
if para.nChains > 1
	if para.useVtens
		Descr = ['VT-',Descr];
	elseif para.useStarMPS
		Descr = ['Star-',Descr];
	else
		Descr = ['VM-',Descr];
	end
end

if isfield(para.chain{1},'s') && para.chain{1}.s ~= 1
	Descr = sprintf('%s-s%g',Descr,para.chain{1}.s);
end

para.folder=sprintf([datestr(now,'yyyymmdd-HHMM-SS'),'-%s-%s-alpha%.10gdelta%.10gepsilon%.10gdk%.10gD%.10gdopt%gL%d'],...
    para.model,Descr,alpha,delta,epsilon,dk,D,d_opt,para.L);
if ~isempty(strfind(para.model,'DPMES'))
	para.folder=sprintf([datestr(now,'yyyymmdd-HHMM-SS'),'-%s-%s-dk%.10gD%.10gdopt%gL%dDelta%gState%d'],...
		para.model,Descr,dk,D,d_opt,para.L,delta,para.InitialState);
end
if ~isempty(strfind(para.model,'SpinBoson')) && strcmp(para.SpinBoson.GroundStateMode,'artificial')
	para.folder = sprintf('%s-art-%s',para.folder,para.SpinBoson.InitialState);
end
para.filename=strcat(para.folder,'/results.mat');
if ~exist(para.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end

%% Start Calculation
[op,para]					= genh1h2term(para);
if (strcmp(para.model,'SpinBosonTTM') && strcmp(para.SpinBoson.GroundStateMode,'artTTM')) || ~isempty(regexp(para.model,'SpinBoson\dCT','match')) ...
		|| (~isempty(regexp(para.model,'SpinBoson\dC','match')) && strcmp(para.SpinBoson.GroundStateMode,'artificial'))
	prepareArtState();
elseif ~isempty(strfind(para.model,'DPMES'))
	prepareArtState();
elseif isfield(para, 'useStarMPS') && para.useStarMPS == 1
	prepareArtState();		% for now only artificial vacuum state
elseif isfield(para, 'MLSB_InitialState') && para.MLSB_InitialState ~= 0
	prepareArtState();
else
	[mps, Vmat,para,results,op] = minimizeE(op,para);
end

if ~isempty(strfind(para.model,'SpinBoson')) && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
%% Reset original parameters after specific ground state preparation!
	if strcmp(para.SpinBoson.GroundStateMode, 'decoupled')
		% restore coupling to chain
        para.t(1) = para.SpinBoson.t1;
	end
	if any(strcmp(para.SpinBoson.InitialState,{'sx','sz'})) && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.hx = para.SpinBoson.hx;
        para.hz = para.SpinBoson.hz;
	end
	if strcmp(para.model,'SpinBosonTTM')
		%% recreate maximally entangled state on A&S
		mps{1}(:,:,1) = [1/sqrt(2) 0]; mps{1}(:,:,2) = [0 1/sqrt(2)];
		mps{2} = zeros(para.D(1),para.D(2),para.d_opt(2));
		mps{2}(1,1,1) = 1; mps{2}(2,1,2) = 1;
	end
    [op,para]   = genh1h2term(para);                % restore operators!
    [op]        = initstorage(mps, Vmat, op,para);

end


save(para.filename,'para','Vmat','mps','results','op','-v7.3');

%% Calculate some Results
results.nx         = getObservable({'occupation'},mps,Vmat,para);
% results.bosonshift = getObservable({'shift'},mps,Vmat,para);

if strcmp(para.model,'SpinBoson') || strcmp(para.model,'SpinBoson2folded')
    results.spin   = getObservable({'spin'},mps,Vmat,para);
end

if strcmp(para.model,'MLSpinBoson')
    results.participation = getObservable({'participation'},mps,Vmat,para);
    results.tunnelEnergy  = getObservable({'tunnelenergy',op},mps,Vmat,para);
end
%%
results.time = toc(starttime)
save(para.filename,'para','Vmat','mps','results','op','-v7.3');

fileName = para.filename;

	function prepareArtState()
		results = initresults(para);
		para
		NC = para.nChains;
		% system state preparation
		if isfield(para,'SpinBoson') && strcmp(para.SpinBoson.GroundStateMode,'artTTM') && strcmp(para.model, 'SpinBosonTTM')
			%% create maximally entangled state between site 1&2 TLS
			mps{1}(:,:,1) = [1/sqrt(2) 0]; mps{1}(:,:,2) = [0 1/sqrt(2)];
			mps{2}(para.D(1),para.D(2),para.d_opt(2)) = 0;
			mps{2}(1,1,1) = 1; mps{2}(2,1,2) = 1;
			Vmat{1} = eye(para.dk(1));
			Vmat{2} = eye(para.dk(2));
			nextSite = 3;
		elseif isfield(para,'SpinBoson') && strcmp(para.SpinBoson.GroundStateMode, 'artificial') && para.useStarMPS == 0
			mps = createrandommps(para);
			if strcmp(para.SpinBoson.InitialState, 'sz')
				%% prepare +Sz eigenstate
				mps{1} = reshape([1,zeros(1,numel(mps{1})-1)],[1,para.D(1),para.d_opt(1)]);
			elseif strcmp(para.SpinBoson.InitialState, '-sz')
				mps{1} = reshape([  zeros(1,numel(mps{1})/2),...
					1,zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
			elseif strcmp(para.SpinBoson.InitialState, 'sx')
				mps{1} = reshape([1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
					1/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
			elseif strcmp(para.SpinBoson.InitialState, '-sx')
				mps{1} = reshape([-1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
					1/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
			elseif strcmp(para.SpinBoson.InitialState, 'sy')
				mps{1} = reshape([ 1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
					1i/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
			elseif strcmp(para.SpinBoson.InitialState, '-sy')
				mps{1} = reshape([-1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
					1i/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
			else
				error('VMPS:minimizeE:DefineInitialState','InitialState=none is not implemented yet');
			end
			Vmat{1} = eye(para.dk(1));
			nextSite = 2;
		elseif ~isempty(strfind(para.model,'DPMES')) && ~para.useStarMPS
			mps{1} = zeros(1,para.D(1),para.dk(1,1));
			mps{1}(1,1,2) = 1;
			Vmat{1} = eye(para.dk(1,1));
			nextSite = 2;
		elseif isfield(para, 'MLSB_InitialState') && ~para.useStarMPS
			mps{1} = zeros(1,para.D(1),para.dk(1,1));
			mps{1}(1,1,para.MLSB_InitialState) = 1;
			Vmat{1} = eye(para.dk(1,1));
			nextSite = 2;
		elseif isfield(para, 'useStarMPS') && para.useStarMPS == 1
			mps{1} = zeros([1,para.D(:,1).',para.dk(1,1)]);
			if ~isempty(strfind(para.model,'DPMES'))
				idx = num2cell([ones(1,NC+1),para.InitialState]);			% start in second excited state!
				mps{1}(idx{:}) = 1;
			elseif isfield(para,'SpinBoson') && strcmp(para.SpinBoson.GroundStateMode, 'artificial')
				if strcmp(para.SpinBoson.InitialState, 'sz')
					%% prepare +Sz eigenstate
					idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
					mps{1}(idx{:},1) = 1;
				elseif strcmp(para.SpinBoson.InitialState, '-sz')
					idx = num2cell(ones(1,NC+1));
					mps{1}(idx{:},2) = 1;
				elseif strcmp(para.SpinBoson.InitialState, 'sx')
					idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
					mps{1}(idx{:},1) = 1/sqrt(2);
					mps{1}(idx{:},2) = 1/sqrt(2);
				elseif strcmp(para.SpinBoson.InitialState, '-sx')
					idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
					mps{1}(idx{:},1) = -1/sqrt(2);
					mps{1}(idx{:},2) = 1/sqrt(2);
				end
			elseif isfield(para, 'MLSB_InitialState')
				idx = num2cell(ones(1,NC+1));			% select state coupling to all first chain states
				mps{1}(idx{:},para.MLSB_InitialState) = 1;
			else
				mps{1}(1) = 1;									% just pick one random element for now!
			end
			Vmat{1} = eye(para.dk(1,1));
			nextSite = 2;
		else
			error('VMPS:prepareArtState:ParameterSettingInvalid', 'No system state prepared');
		end



		% Environment preparation
		for j = nextSite:para.L
			% Construct single-site MPS
			if para.useStarMPS
				mps{j} = cell(1,NC);
				for mc = 1:NC
					if j > para.chain{mc}.L, continue; end
					if j == para.chain{mc}.L, Dr = 1;
					else Dr = para.D(mc,j); end
					mps{j}{mc} = zeros(para.D(mc,j-1),Dr,para.d_opt(mc,j));
					mps{j}{mc}(1,1,1) = 1;
				end
			else
				if j == para.L, Dr = 1;
				else Dr = para.D(j); end
				mps{j} = zeros(para.D(j-1),Dr,para.d_opt(end,j));
				mps{j}(1,1,1) = 1;
			end

			% Construct single-site Vmat
			if para.useVtens
				Vmat{j} = cell(1,NC+1);
				for mc = 1:NC
					Vmat{j}{mc} = sparse(1:para.d_opt(mc,j),1:para.d_opt(mc,j),1,para.dk(mc,j),para.d_opt(mc,j));	% 1:n order
				end
				if ~isempty(regexp(para.model,'SpinBoson\dCT','match'))
					% prepare maximally entangled boson-ancilla state
					Vmat{j}{end} = zeros(para.d_opt(:,j)');
					if NC == 2
						Vmat{j}{end}(:,:,1) = eye(para.d_opt(1,j))/sqrt(para.d_opt(1,j));
					end
				else
					Vmat{j}{end} = sparse(1,1,1, prod(para.d_opt(1:end-1,j)),para.d_opt(end,j));
					Vmat{j}{end} = reshape(full(Vmat{j}{end}), para.d_opt(:,j)');
				end
			elseif para.useStarMPS
				Vmat{j} = cell(1,NC);
				for mc = 1:NC
					Vmat{j}{mc} = sparse(1:para.d_opt(mc,j),1:para.d_opt(mc,j),1,para.dk(mc,j),para.d_opt(mc,j));	% 1:n order
				end
			else  % no special order for MC-Vmat!
				if NC == 1
% 					Vmat{j}	= sparse(order(1:para.d_opt(j)),1:para.d_opt(j),1,para.dk(j),para.d_opt(j));	% inverted order
					Vmat{j}	= sparse(1:para.d_opt(j),1:para.d_opt(j),1,para.dk(j),para.d_opt(j));			% 1:n order
				else
					Vmat{j} = sparse(1,1,1, prod(para.dk(:,j)),para.d_opt(j));								% 1:n order
				end
			end
		end

		[mps,Vmat,para] = prepare(mps,Vmat,para);
		if ~para.useStarMPS
			[op] = initstorage(mps, Vmat, op,para);
		else
			% copied from prepare_ChainOp
			for mc = 1:NC
				cL = para.chain{mc}.L;
				cM = para.M/NC;

				mpsChain    = cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false);
				VmatChain   = cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false);

				paraChain   = para;
				paraChain.L = cL-1;
				paraChain.M = cM;
				paraChain.nChains = 1;
				paraChain.useStarMPS = 0;

				% Copy operators
				opChain.h1term = cell(1,cL-1);
				opChain.h2term = cell(cM, 2,cL-1);

				opChain.h1term = op.h1term(mc, 2:cL);
				opChain.h2term = op.h2term( cM*(mc-1)+(1:cM), :, 2:cL, mc);

				opChain = initstorage(mpsChain,VmatChain,opChain,paraChain);
				paraChain.sitej = 1;
				opChain = updateop(opChain,mpsChain,VmatChain,1,paraChain);
				op.chain(mc).Hlrstorage = opChain.Hlrstorage;
				op.chain(mc).Opstorage  = opChain.Opstorage;
			end
		end

		para.trustsite = para.L;		% needed for TDVP
		return;
	end
end
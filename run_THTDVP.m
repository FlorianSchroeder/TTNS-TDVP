function run_TDVP(alpha,maxDim)
% try
% needed if started from command line:
if isdeployed
% 	if ischar(s), s = str2double(s); end
	if ischar(alpha), alpha = str2double(alpha); end
	if ischar(maxDim), 	maxDim 	= str2double(maxDim); end
% 	if ischar(Bond), 	Bond 	= str2double(Bond); end
% 	if ischar(dk),		dk		= str2double(dk); end
% 	if ischar(L),		L		= str2double(L); end
end

fromFile = sprintf('20150907-0014-SpinBoson2CT-VT-OrthPol-v62TCMde10-alpha0.01delta0epsilon0.1dk70D5dopt5L20-art--sz/results-Till20Step0.1v62-OBBmax%d-Dmax%d-expvCustom0-1core.mat',maxDim,maxDim);

%% start ground state calculations
loadedFromFile = 0;
if isempty(fromFile)
% fileName =  VMPS_FullSBM(s,alpha,0.1,0,50,0);     % VMPS_FullSBM(s,alpha,delta,epsilon,L,rescaling)
% 	fileName =  VMPS_FullSBM(s,alpha,0,0.1,L,dk);     % iSBM(s,alpha,delta,epsilon,L,rescaling)
error('please load a file');
else
	fileName = fromFile;							% simple override!
	loadedFromFile = 1;
end
% maxNumCompThreads(1);						% safer in terms of results!				% does not work anymore, use -singleCompThread instead!

%% Define GS config

load(fileName);

% load(sprintf('20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha%gdelta0.1epsilon0dk20D5dopt5L50-artificial/results.mat',alpha));

%% Only needed if previous calc was imagT
if loadedFromFile && isfield(para.tdvp,'imagT') && para.tdvp.imagT
	% prepare -sx eigenstate!
	mps{1} = reshape([-1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
				1/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
	disp('Prepared -sx eigenstate');
end

%% Define TDVP parameters
para.tdvp.imagT = 0;					% imaginary Time = Temperature evolution?
para.tdvp.tmax = 350;
para.tdvp.deltaT = 1;					% size of timeslice in units:
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)�10^-16
para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
para.tdvp.resume = 0;					% additionally control if want to resume!
para.tdvp.saveInterval = 10;			% save '-small.mat' every n-th step
para.tdvp.serialize = 1;				% much faster I/O saving
para.tdvp.logSV = 0;					% if 1 only log SV, if 0 only log vNE (saves mem) if -1 log none!
para.tdvp.extractStarInterval = para.tdvp.deltaT;	% in [t]; for calculating star occupation! Comment if not needed!
para.tdvp.extractObsInterval  = para.tdvp.deltaT;	% in [t]; mod(extractStarInterval, extractObsInterval) = 0 !! extractObsInterval = n*deltaT
para.tdvp.Observables = '.n.s.j.sn.sx.';		% n: occupation, j: current, s: spin, sn: star n, sx: star polaron
para.tdvp.storeMPS = 0;					% save tmps or not!
para.tdvp.maxExpMDim = 00;				% For Lappy: 100, OE-PC: 80, pc52: 260; E5: 300 System dependent, use benchmark!
para.tdvp.maxExpVDim = 00;				% higher dim -> use expvCustom() if expvCustom == 1. Number from benchmarking. Lappy: 400, Haswell: 800; E5: 700 maxExpMDim < maxExpVDim
para.tdvp.expvCustom = 1;				% 1 for Custom programmed, 0 for standard expv()
para.tdvp.expvCustomTestAccuracy = 0;	% do expvCustom alongside expv for testing.
para.tdvp.expvCustomTestAccuracyRMS = 0;	% display RMS of expvCustom from expv(); set only if para.tdvp.expvCustomTestAccuracy = 1;
para.tdvp.expvTol = 1e-15;              % error tolerance of expv(); default: 1e-7
para.tdvp.expvM   = 50;                 % dim of Krylov subspace in expv(); default: 30
    % Sets threshold size for matrix exponential:
    %   if dim(A) < : use built-in expm(At)*v
    %   else        : use Expokit expv(t,A,v, expvTol, expvM)
    %   set maxExpMDim = 0 to only use expv()
% OBB settings
para.tdvp.expandOBB = 1;
% Bond-Dim settings
para.tdvp.truncateExpandBonds = 1;
% Calculate max Bond Dim: 1GB for array (l,r,n,l,r,n) with n around 20,
% 1 complex double needs 16byte. -> 20^6 * 16byte < 1GB
% para.tdvp.maxBondDim = 80;				% keep imagT maxDim
para.Dmin = 1;
% para.tdvp.maxOBBDim  = 30;
para.svmaxtol = 10^-6;					% keep 1 below this!
para.svmintol = 10^-6.5;				% throw away all below
% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization
	if ~exist('stepFrom','var')
		stepFrom = 1;					% comment if need override!
	end
end

tresults = [];						% empty variable initialization
%% comment if no new coupling needed
if loadedFromFile && para.chain{1}.alpha ~= alpha
	para.chain{1}.alpha = alpha;
	[para.chain{1}]=SBM_genpara(para.chain{1});
	[op,para]=genh1h2term(para);
	[op] = initstorage(mps, Vmat, op,para);
end

%% Format Filename
para.tdvp.version = 'v62';
if isfield(para.tdvp,'filename')
	%% Continued TDVP remember filename to load after directory change!
	% from File can be -small.mat!
	para.tdvp.fromFilename = para.tdvp.filename;		% save the reference to continued file
	% input the TDVP matrices!
% 	mps = tmps;
% 	Vmat = tVmat;
else
	para.tdvp.fromFilename = para.filename;				% or reference to VMPS Ground State File!
end

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%.2g%s'],para.tdvp.tmax,para.tdvp.deltaT,para.tdvp.version);
para.tdvp.filename = sprintf('%s-alpha%g',para.tdvp.filename,alpha);

if para.tdvp.expandOBB
	para.tdvp.filename = sprintf('%s-OBBmax%d',para.tdvp.filename, para.tdvp.maxOBBDim);
end
if para.tdvp.truncateExpandBonds
	para.tdvp.filename = sprintf('%s-Dmax%d',para.tdvp.filename,para.tdvp.maxBondDim);
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

%% Check for valid simulation arguments and declare unsettables
if para.tdvp.expvCustomTestAccuracyRMS
	assert(para.tdvp.expvCustomTestAccuracy == 1, 'Only set para.tdvp.expvCustomTestAccuracyRMS == 1 if para.tdvp.expvCustomTestAccuracy == 1');
end
if para.tdvp.zAveraging
	assert(strcmp(para.chain.discretization,'LogZ'),'Use z-Averaging only with Logarithmic Discretization!');
end
if para.tdvp.expvCustom
	assert(para.tdvp.maxExpMDim <= para.tdvp.maxExpVDim,'maxExpMDim <= maxExpVDim ! Everything else has no sense.');
end
if para.tdvp.extractObsInterval < para.tdvp.deltaT
	para.tdvp.extractObsInterval = para.tdvp.deltaT;
end

para.tdvp.expvCustomNow = 0;			% only used inside the program
para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP; needs to be off since H in exponent!
para.rescaling = para.tdvp.rescaling;
para.complex = 1;						% necessary since time-evolution is complex

%% Copy to scratch for computation
if ~strcmp(computer,'PCWIN64')
	[~, name] = system('hostname');
	para.tdvp.hostname = strtrim(name);
	para.tdvp.scratchDir = '/scratch/fayns2/TDVPtemp/'; tempFold = fileparts(para.filename);
	currentDir = pwd;
	addpath(currentDir);
	if ~exist(para.tdvp.scratchDir,'dir')
		mkdir(para.tdvp.scratchDir);
	end
	mkdir([para.tdvp.scratchDir,tempFold]);
	copyfile(para.filename,[para.tdvp.scratchDir,tempFold]);
	save(sprintf([para.tdvp.filename(1:end-4),'-incomplete-%s.mat'],para.tdvp.hostname),'para','results');
	cd(para.tdvp.scratchDir);
end

if exist(para.tdvp.filename,'file') && para.tdvp.resume && isempty(fromFile)
	% if file already exists && resume is switched on
	try
		load(para.tdvp.filename);	% all the rest
		disp('loaded para.tdvp.filename');
	catch
		disp('MAT-file corrupt, load BAK');
		load([para.tdvp.filename(1:end-4),'.bak'],'-mat');		% use backup if file corrupted
	end
	Vars = whos;				% deserialise if needed
	for ii = 1:size(Vars)
		if strcmp(Vars(ii).class,'uint8')
			eval(sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
		end
	end
	para.tdvp.resume = 1;		% need override to restore value!
elseif ~exist(para.tdvp.filename,'file') && isempty(strfind(para.tdvp.fromFilename,'results.mat')) ...
		&& exist(para.tdvp.fromFilename,'file') && ~loadedFromFile
	% Only load fromFilename if it was a TDVP file!
	% then also para and tresults are already loaded!
	% important: load mps and Vmat!
	load(para.tdvp.fromFilename, 'results','op','mps','Vmat');
	disp('loaded results,op,mps,Vmat from para.tdvp.fromFilename');
	Vars = whos;				% deserialise if needed
	for ii = 1:size(Vars)
		if strcmp(Vars(ii).class,'uint8')
			eval(sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
		end
	end
end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    para.tdvp.starttime = tic;
    tdvp_1site(mps,Vmat,para,results,op,tresults);
else
	basename = para.tdvp.filename;
    for z = stepFrom:-para.tdvp.zStep:1e-5
        %% Set new z-value, prepare chain
        para.z = z;
        [para]=SBM_genpara(para);
        [op,para]=genh1h2term(para);
        [op] = initstorage(mps, Vmat, op,para);
        para.tdvp.filename = sprintf([basename(1:end-4),'-z%.10g.mat'],para.z);

        %% Do time-evolution and save results
        para.tdvp.starttime = tic;
        tdvp_1site(mps,Vmat,para,results,op);

		if ~strcmp(computer,'PCWIN64')
			copyfile(para.tdvp.filenameSmall,[currentDir,'/',para.tdvp.filenameSmall]);
		end
        %% clear all results from sweep
        clear('tmps','tVmat','mps1','Vmat1');
        results.tdvp = struct();
    end
end
if ~strcmp(computer,'PCWIN64')
	copyfile(para.tdvp.filenameSmall,[currentDir,'/',para.tdvp.filenameSmall]);
	%delete([currentDir,'/',para.tdvp.filename(1:end-4),'-incomplete.mat']);
% 	sendmailCAM('fayns2@cam.ac.uk',...
%          'TDVP job completed',sprintf('The job \n %s\nHas successfully completed.',para.tdvp.filename));
	exit;
end

% catch err
% 	fprintf([getReport(err),'\n']);
% 	sendmailCAM('fayns2@cam.ac.uk',...
%          'TDVP job ERROR',sprintf('The job \n %s\nHas encountered an error:\n%s',para.tdvp.filename,getReport(err)));
% end
end
function tresults = calTimeObservables(tmps,tVmat,para,tresults)
%% function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%	creates tresults and computes Observables for each timeslice in tmps
%
%	tresults = calTimeObservables(tmps,tVmat,para,tresults)
%		Appends the slices in tmps to the end of tresults.
%			e.g.: tresults ranges in t = [0:3]; tmps has slices t=[4:6]
%				  then returned tresults will have t = [0:6]
%
%	tresults = calTimeObservables(treeMPS,[],para,tresults)
%		Does same as above for the treeMPS, see subfunction.
%		Only for single time step for simplicity!
%
% Modified:
%	FS 23/07/15: - changed nx to n for Multi-Chain compatibility!
%	FS 24/07/15: - added switch to allow observable selection (especially current)
%   FS 18/08/15: - changed to single precision for smaller tresults file!
%	FS 03/03/16: - removed loop over slices -> perform that in outer function
if nargin < 4 
	tresults = [];		% allow access to variable
	para.timeslice = 0;						% needed for extractObsInterval first run
end

if para.useTreeMPS
	tresults = calTimeObservables_Tree(tmps,para,tresults);
	return;
end

% switches
isNew    = 0;			% initialise variables?
skipObs  = 0;		
skipStar = 0;			% skip star observables?

% Parameters
O        = para.tdvp.Observables;		% Observables list
NC       = para.nChains;
L        = para.L;

if size(tmps,1) > 1
	% calculate Obs for many slices -> iterate over single slices
	% need correct para.timeslice to do this!
	for ii = 1:size(tmps,1)
		tresults = calTimeObservables(tmps(ii,:),tVmat(ii,:),para,tresults);
		para.timeslice = para.timeslice+1;
	end
	return;
end

% From here on: only single timeslice!
if isempty(tresults)
	isNew = 1;
	tresults = struct;
	missingN = 0;
	tresults.lastIdx = 0;
	tresults.star.lastIdx = 0;
	fprintf('Calculate Observables:\n');
end

if isfield(para.tdvp,'extractObsInterval')
	% only works with equidistant steps and single tmps slices
	if mod(para.tdvp.tmax, para.tdvp.extractObsInterval) == 0 && (para.tdvp.extractObsInterval >= para.tdvp.deltaT)
		totalN = round(para.tdvp.tmax/para.tdvp.extractObsInterval) +1;
	else
		error('VMPS:calTimeObservables:InvalidParameter','Need to define extractObsInterval so that mod(tmax,interval)=0!');
	end
	if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractObsInterval) ~= 0
		skipObs = 1;
	end
else
	para.tdvp.extractObsInterval = para.tdvp.deltaT;
	totalN = size(para.tdvp.t,2);
end

if isfield(para.tdvp,'extractStarInterval')
	% only works with equidistant steps and single tmps slices
	if mod(para.tdvp.tmax, para.tdvp.extractStarInterval) == 0 && (para.tdvp.extractStarInterval >= para.tdvp.deltaT)
		totalStarN = round(para.tdvp.tmax/para.tdvp.extractStarInterval) +1;
	else
		error('VMPS:calTimeObservables:InvalidParameter','Need to define extractStarInterval so that mod(tmax,interval)=0!');
	end
	if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractStarInterval) ~= 0
		skipStar = 1;
	end
else
	para.tdvp.extractStarInterval = para.tdvp.deltaT;
	totalStarN = size(para.tdvp.t,2);
end

% Intialise t
if isNew
	tresults.t          = single(0:para.tdvp.extractObsInterval:para.tdvp.tmax);
	tresults.star.t     = single(0:para.tdvp.extractStarInterval:para.tdvp.tmax);
	tresults.star.omega = [];
end

% Extend variables if resumed calculation
if ~isNew
	%%
	missingN = totalN - length(tresults.t);
	if missingN > 0
		% resumed calculation. All variables need to be extended
		assert(tresults.lastIdx == length(tresults.t), 'lastIdx does not agree with lenght(t)');
		% Iterate over all fields of tresults
		names = fieldnames(tresults);
		for ii = 1:length(names)
			if ~isnumeric(tresults.(names{ii})),continue,end;		% skip non-arrays
			if any(strcmp(names{ii},{'t','lastIdx'})),continue,end;	% skip t and lastIdx
			d = size(tresults.(names{ii}));
			tresults.(names{ii}) = cat(1,tresults.(names{ii}), zeros([missingN,d(2:end)],'single'));
		end
		% Iterate over all fields of tresults.star
		names = fieldnames(tresults.star);
		for ii = 1:length(names)
			if ~isnumeric(tresults.star.(names{ii})),continue,end;		% skip non-arrays
			if any(strcmp(names{ii},{'t','lastIdx'})),continue,end;	% skip t and lastIdx
			d = size(tresults.star.(names{ii}));
			tresults.star.(names{ii}) = cat(1,tresults.star.(names{ii}), zeros([missingN,d(2:end)],'single'));
		end
	end
end

i = tresults.lastIdx + 1 ;			% index in tresults
j = tresults.star.lastIdx + 1 ;		% index in tresults.star
		
fprintf('<O> for slice %d\n',i);

%% Chain Observables
% 1. Chain Occupation
% if strContains(O,'.n.') && ~skipObs
% 	if isNew
% % 		tresults.n  = zeros(totalN,L,max(NC,NE),'single');							% t x L x NC
% 		tresults.n  = zeros(totalN,L,NC,'single');									% t x L x NC
% 	end
% 	tresults.n(i,:,:) = single(getObservable({'occupation'},tmps,tVmat,para));		% (L x nChain)
% end

if strContains(O,'.n.') && ~skipObs																% Absolute non-projected occupation
	chainN = real(getObservable({'bath1correlators','n'}, tmps,tVmat,para));					% L x nChains
	if isNew
		d = size(chainN);
		tresults.n = zeros([totalN,d],'single');
	end
	tresults.n(i,:,:) = chainN;
end
if strContains(O,'.nd.') && ~skipObs															% diabatic projected occupation
	chainN = real(getObservable({'bath1correlators','n','diabatic'}, tmps,tVmat,para));			% L x nChains
	if isNew
		d = size(chainN);
		tresults.nd = zeros([totalN,d],'single');
	end
	tresults.nd(i,:,:,:) = chainN;
end
if strContains(O,'.na.') && ~skipObs															% adiabatic projected occupation
	chainN = real(getObservable({'bath1correlators','n','adiabatic'}, tmps,tVmat,para));		% L x nChains
	if isNew
		d = size(chainN);
		tresults.na = zeros([totalN,d],'single');
	end
	tresults.na(i,:,:,:) = chainN;
end

% 2. Chain displacement
if strContains(O,'.x.') && ~skipObs																% Absolute non-projected displacement
	if ~exist('chainX','var')
		chainX = real(getObservable({'bath1correlators','x'}, tmps,tVmat,para));				% L x nChains
	end
	if isNew
		d = size(chainX);
		tresults.x = zeros([totalN,d],'single');
	end
	tresults.x(i,:,:) = chainX;
end

if strContains(O,'.xd.') && ~skipObs																% diabatic projected displacement
	chainX = real(getObservable({'bath1correlators','x','diabatic'}, tmps,tVmat,para));	% L x nStates x nChains
	if isNew
		d = size(chainX);
		tresults.xd = zeros([totalN,d],'single');
	end
	tresults.xd(i,:,:,:) = chainX;
end
if strContains(O,'.xa.') && ~skipObs																% adiabatic projected displacement
	chainX = real(getObservable({'bath1correlators','x','adiabatic'}, tmps,tVmat,para));% L x nStates x nChains
	if isNew
		d = size(chainX);
		tresults.xa = zeros([totalN,d],'single');
	end
	tresults.xa(i,:,:,:) = chainX;
end

if strContains(O,'.x2.') && ~skipObs															% Absolute non-projected displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2'}, tmps,tVmat,para));					% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2 = zeros([totalN,d],'single');
	end
	tresults.x2(i,:,:) = chainX2;
end
if strContains(O,'.x2d.') && ~skipObs															% diabatic projected displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2','diabatic'}, tmps,tVmat,para));		% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2d = zeros([totalN,d],'single');
	end
	tresults.x2d(i,:,:,:) = chainX2;
end
if strContains(O,'.x2a.') && ~skipObs															% adiabatic projected displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2','adiabatic'}, tmps,tVmat,para));		% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2a = zeros([totalN,d],'single');
	end
	tresults.x2a(i,:,:,:) = chainX2;
end

if strContains(para.tdvp.Observables,'.j.','.sn.')
	% save here, to reuse later!
	AnAm = getObservable({'bath2correlators'}, tmps,tVmat,para);
end
if strfind(para.tdvp.Observables,'.j.')
	%% Calculate current along entire chain
	if ~isfield(tresults,'j') || missingN > 0
		tresults.j(totalN,para.L,nChains) = single(0);
	end
	if ~exist('AnAm','var')
		tresults.j(i,:,:) = single(getObservable({'current'},tmps,tVmat,para));
	else
		tresults.j(i,:,:) = single(getObservable({'current',AnAm},tmps,tVmat,para));
	end
end

% 3. Star Observables
if isfield(para.tdvp,'extractStarInterval') && strContains(para.tdvp.Observables,'.sn.','.sx.','.sx2.')
	Nslice = round(para.tdvp.extractStarInterval / para.tdvp.extractObsInterval);		% how often to extract Star Observables
	if mod(i-1,Nslice) == 0
		pos = ceil(i/Nslice);

		if strfind(para.tdvp.Observables,'.sn.')
			if ~exist('AnAm','var')
				occ	= getObservable({'staroccupation'}     ,tmps,tVmat,para);		% (1+1) x k x nc
			else
				occ = getObservable({'staroccupation',AnAm},tmps,tVmat,para);		% (1+1) x k x nc
			end
			starOmega = squeeze(single(occ(1,:,:)));				% get rid of leading singleton
		end

		if strContains(para.tdvp.Observables,'.sx.','.sx2.')
			if strfind(para.tdvp.Observables,'.sx.')
				polaron = getObservable({'starpolaron'},tmps,tVmat,para);				% (1+2) x k x nc, diabatic states
			elseif strfind(para.tdvp.Observables,'.sx2.')
				% not state projecting, but selecting single dominating states across the first bond!
				% kind of similar to diabatic states picture!
				polaron = getObservable({'starpolaron','adiabatic'},tmps,tVmat,para);	% (1+2) x k x nc

			end
			starOmega = squeeze(single(polaron(1,:,:)));			% get rid of leading singleton
		end

		if ~isfield(tresults, 'star')
			% initialise storage if first sweep
			nElements = para.tdvp.tmax/para.tdvp.extractStarInterval +1;
			if exist('occ','var'),     tresults.star.n = single(zeros(nElements,size(occ,2),para.nChains)); end
			if exist('polaron','var'), tresults.star.x = single(zeros(nElements,size(polaron,2),2,para.nChains)); end
			tresults.star.omega = starOmega;
			tresults.star.t     = single(zeros(1,nElements));
		end

		if strfind(para.tdvp.Observables,'.sn.')
			tresults.star.n(pos,:,:) = single(occ(2,:,:));
		end
		if strContains(para.tdvp.Observables,'.sx.','.sx2.')
			for kk = 1:size(polaron,1)-1
				tresults.star.x(pos,:,kk,:) = single(polaron(kk+1,:,:));	% SBM: 1: up-proj, 2: down-proj
			end
		end
		tresults.star.t(pos)   = single(para.tdvp.t(1,para.timeslice+1));
	end
end

if ~isempty(strfind(para.model, 'SpinBoson'))
	%% Observables for SBM
	% 1. Spin Observables
	if ~isfield(tresults,'spin')
		tresults.spin.sx = single(zeros(totalN,1));
		tresults.spin.sy = single(zeros(totalN,1));
		tresults.spin.sz = single(zeros(totalN,1));
		tresults.spin.visibility = single(zeros(totalN,1));
	elseif missingN > 0
		tresults.spin.sx = single([tresults.spin.sx; zeros(missingN,1)]);
		tresults.spin.sy = single([tresults.spin.sy; zeros(missingN,1)]);
		tresults.spin.sz = single([tresults.spin.sz; zeros(missingN,1)]);
		tresults.spin.visibility = single([tresults.spin.visibility; zeros(missingN,1)]);
	end
	temp = getObservable({'spin'},tmps,tVmat,para);
	tresults.spin.sx(i) = single(temp.sx);
	tresults.spin.sy(i) = single(temp.sy);
	tresults.spin.sz(i) = single(temp.sz);
	tresults.spin.visibility(i) = single(sqrt(temp.sx^2+temp.sy^2));
end

% 1. Density matrix
if strContains(para.tdvp.Observables,'.dm.','.dm2.')
	if ~isfield(tresults,'rho')
		tresults.rho = single(zeros(totalN,para.dk(1,1),para.dk(1,1)));
	elseif missingN > 0
		tresults.rho = single([tresults.rho; zeros(missingN,para.dk(1,1),para.dk(1,1))]);
	end
	tresults.rho(i,:,:,1) = single(getObservable({'rdm',1},tmps,tVmat,para));
	if strContains(para.tdvp.Observables,'.dm2.')
		% only for 2-lvl system for now; only calculates largest bond state.
		tresults.rho(i,:,:,2) = single(getObservable({'rdm_adiabatic',1,1},tmps,tVmat,para));  %{'rdm_adiabatic',sitej,state}
		tresults.rho(i,:,:,3) = single(getObservable({'rdm_adiabatic',1,2},tmps,tVmat,para));  %{'rdm_adiabatic',sitej,state}
	end
end

if strcmp(para.model, 'MLSBM') || ~isempty(strfind(para.model,'DPMES'))
	%% Observables for MLSBM
	% 2. PPC Wavefunction
	%    only if not extracting Density Matrix
	if ~strfind(para.tdvp.Observables,'.dm.')
		if ~isfield(tresults,'PPCWavefunction')
			tresults.PPCWavefunction = single(zeros(totalN,para.dk(1,1)));
		elseif missingN > 0
			tresults.PPCWavefunction = single([tresults.PPCWavefunction; zeros(missingN,para.dk(1,1))]);
		end
		tresults.PPCWavefunction(i,:) = single(diag(getObservable({'rdm',1},tmps,tVmat,para)));
	end

	% 3. Participation on ring
	if ~isfield(tresults,'participation')
		tresults.participation = single(zeros(totalN,1));
	elseif missingN > 0
		tresults.participation = single([tresults.participation; zeros(missingN,1)]);
	end
	tresults.participation(i) = single(getObservable({'participation'},tmps,tVmat,para));

	% 4. Hs + Hi
	if ~isfield(tresults,'hshi')
		tmp = single(getObservable({'hshi'},tmps,tVmat,para));
		tresults.hshi = single(zeros(totalN,length(tmp)));
		tresults.hshi(i,:) = tmp;
	elseif missingN > 0
		tresults.hshi = single([tresults.hshi; zeros(missingN,1)]);
	else
		tresults.hshi(i,:) = single(getObservable({'hshi'},tmps,tVmat,para));
	end

end

if strContains(para.tdvp.Observables,'.sp.')					% sp for state projection
	if ~isfield(tresults,'stateProjection')
		tresults.stateProjection = single(zeros(totalN,1));
	elseif missingN > 0
		tresults.stateProjection(totalN,1) = 0;			% does preallocation
	end
	tresults.stateProjection(i,1) = single(getObservable({'stateproject',para.InitialState,1},tmps,tVmat,para));	% project onto |IS>|0>, IS = initial state
end

if strContains(para.tdvp.Observables,'.ss.')					% ss for system state
	if ~isfield(tresults,'system') || ~isfield(tresults.system,'state')
		tresults.system.state = single(zeros(totalN,para.dk(1),para.dk(1)));		% t x dk x D (adiabatic)
	elseif missingN > 0
		tresults.system.state(totalN,para.dk(1),para.dk(1)) = 0;			% does preallocation
	end
	tresults.system.state(i,:,:) = single(getObservable({'state',1},tmps,tVmat,para));
end

if strContains(para.tdvp.Observables,'.ses.')					% ses for system-environment state
	if ~isfield(tresults,'mps')
		tresults.mps = cell(totalN,2);		% t x sites
		tresults.Vmat = cell(totalN,2);		% t x sites
	elseif missingN > 0
		tresults.mps(totalN,2) = {};			% does preallocation
		tresults.Vmat(totalN,2) = {};
	end
	out = getObservable({'sys-env-state'},tmps,tVmat,para);		% get mps([1,2]) and Vmat([1,2])
	tresults.mps(i,:) = out.mps;
	tresults.Vmat(i,:) = out.Vmat;
end

if strcmp(para.model, 'SpinBosonTTM')
	%% extract transfer tensor
	% only use for single-slice tMPS due to iterative procedure
	if length(slices) > 1, return; end;
	rdm = getObservable({'rdm',[1 2]},tmps,tVmat,para);
	% Ortho Normal Operator Basis in dxd
	d = para.dk(1,2);
	ONOB = eye(d^2); ONOB = reshape(ONOB,[d,d,d^2]);
	EAm = zeros(d,d,d^2);
	for k = 1:d^2
		EAm(:,:,k) = ncon({rdm, squeeze(ONOB(:,:,k))'},...
						  {[-1,2,-2,1], [1,2]})*d;			% apply Op, contract / trace; perhaps *d
	end
	Epsilon = reshape(EAm,[d^2,d^2]);						% Simpler: E_{(ij)k} = d* B_{(lm)k}* \rho_{i,l,j,m} can be performed as direct matrix product! TODO since much faster!
	T = Epsilon;
	for k = 3:i
		T = T - tresults.TTM.T(:,:,i+1-k)*tresults.TTM.Epsilon(:,:,k-1);          % TODO: vectorize for loop?
	end
	tresults.TTM.Epsilon(:,:,i) = Epsilon;
	if i > 1
		tresults.TTM.T(:,:,i-1) = T;
		tresults.TTM.Tnorm(i-1) = single(norm(T));
		fprintf('\n|TTM|/dt^2: %g\n',tresults.TTM.Tnorm(i-1)/para.tdvp.deltaT^2);
	end
end

if ~skipObs
	tresults.lastIdx = i;
end
if ~skipStar
	tresults.star.lastIdx = j;
end
fprintf('\n');
end

function out = strContains(str, varargin)
out = false;
for ii = 1:length(varargin)
	out = out || ~isempty(strfind(str,varargin{ii}));
end

end

function tresults = calTimeObservables_Tree(treeMPS,para,tresults)
%% function tresults = calTimeObservables_Tree(treeMPS,para,tresults)
%
%	calculates the Observables specified in para.tdvp.Observables for each time step
%	Only able to handle single time step in treeMPS for now!
%
%	para.tdvp.Observables:
%		.n.		bosonic occupation
%		.x.		bosonic displacement
%		.xd.	bosonic displacement, diabatic projection
%		.xa.	bosonic displacement, adiabatic projection
%		.sn.	bosonic star-displacement
%		.sx.	bosonic star-displacement
%		.sxd.	bosonic star-displacement, diabatic projection
%		.sxa.	bosonic star-displacement, adiabatic projection
%		.j.		bosonic current
%		.dm.	reduced density matrix of site ? (diabatic)
%		.dma.	reduced density matrix of site ? (adiabatic)
%		.ac.	autocorrelation
%		.ses.	system-environment state
%		.ss.	system state
%		.heff.  adiabatic total     energy surfaces + system state
%		.pes.   adiabatic potential energy surfaces + system state
%
%
%	Created by FS 26/02/2016

%% Initialisation
% switches
isNew    = 0;					% if constructed new tresults
skipObs  = 0;
skipStar = 0;
heffPes  = 0;					% whether heff and pes are selected

% Parameters
O        = para.tdvp.Observables;		% Observables list
NC       = para.nChains;
% NE       = para.nEnvironments;
L        = para.L;						% total max height of tree + 1 (edges+1 = #sites)

if isempty(tresults)
	% intialise tresults.
	tresults = struct;
	isNew = 1;				% switch for initialisation of each field
	missingN = 0;
	tresults.lastIdx = 0;
	tresults.star.lastIdx = 0;
	fprintf('Calculate Observables:\n');
elseif iscell(tresults)
	tresults = tresults{1};
end

if isfield(para.tdvp,'extractObsInterval')
	% only works with equidistant steps and single tmps slices
	if mod(para.tdvp.tmax, para.tdvp.extractObsInterval) == 0 && (para.tdvp.extractObsInterval >= para.tdvp.deltaT)
		totalN = round(para.tdvp.tmax/para.tdvp.extractObsInterval) +1;
	else
		error('VMPS:calTimeObservables:InvalidParameter','Need to define extractObsInterval so that mod(tmax,interval)=0!');
	end
	if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractObsInterval) ~= 0
		skipObs = 1;
	end
else
	totalN = size(para.tdvp.t,2);
end

if isfield(para.tdvp,'extractStarInterval')
	% only works with equidistant steps and single tmps slices
	if mod(para.tdvp.tmax, para.tdvp.extractStarInterval) == 0 && (para.tdvp.extractStarInterval >= para.tdvp.deltaT)
		totalStarN = round(para.tdvp.tmax/para.tdvp.extractStarInterval) +1;
	else
		error('VMPS:calTimeObservables:InvalidParameter','Need to define extractStarInterval so that mod(tmax,interval)=0!');
	end
	if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractStarInterval) ~= 0
		skipStar = 1;
	end
else
	totalStarN = size(para.tdvp.t,2);
end

if isNew
	tresults.t  = single(zeros(1,totalN));
	tresults.star.t  = single(zeros(1,totalStarN));
	tresults.star.omega = [];
end

i = tresults.lastIdx + 1;
j = tresults.star.lastIdx + 1;

%% System Observables


%% Chain Observables

% 1. Chain Occupation
if strContains(O,'.n.') && ~skipObs
	chainN = real(getObservable({'bath1correlators','n'}, treeMPS,[],para));					% L x nChains
	if isNew
		d = size(chainN);
		tresults.n  = zeros([totalN,d],'single');												% t x L x NC
	end
	tresults.n(i,:,:) = chainN;
end
if strContains(O,'.nd.') && ~skipObs															% diabatic projected occupation
	chainN = real(getObservable({'bath1correlators','n','diabatic'}, treeMPS,[],para));			% L x States x nChains
	if isNew
		d = size(chainN);
		tresults.nd = zeros([totalN,d],'single');
	end
	tresults.nd(i,:,:,:) = chainN;
end
if strContains(O,'.na.') && ~skipObs															% adiabatic projected occupation
	chainN = real(getObservable({'bath1correlators','n','adiabatic'}, treeMPS,[],para));			% L x States x nChains
	if isNew
		d = size(chainN);
		tresults.na = zeros([totalN,d],'single');
	end
	tresults.na(i,:,:,:) = chainN;
end
if strContains(O,'.nc.') && ~skipObs															% coherence projected occupation
	chainN = getObservable({'bath1correlators','n','lettcoherence'}, treeMPS,[],para);			% L x 1 x nChains
	if isNew
		d = size(chainN);
		tresults.nc = zeros([totalN,d(1),d(3)],'single');
	end
	tresults.nc(i,:,:,:) = chainN(:,1,:);		% only take one state slice, since they are all equal
end

% 2. Chain Displacement
if strContains(O,'.x.') && ~skipObs																% displacement
	chainX = real(getObservable({'bath1correlators','x'}, treeMPS,[],para));					% L x nChains
	if isNew
		d = size(chainX);
		tresults.x = zeros([totalN,d],'single');
	end
	tresults.x(i,:,:) = chainX;																	% (L x nChain)
end
% 2.1 Chain Displacement, diabatic
if strContains(O,'.xd.') && ~skipObs															% diabatic projected displacement
	chainX = real(getObservable({'bath1correlators','x','diabatic'}, treeMPS,[],para));			% L x nStates x nChains
	if isNew
		d = size(chainX);
		tresults.xd = zeros([totalN,d],'single');
	end
	tresults.xd(i,:,:,:) = chainX;																% (L x nStates x nChain)
end
% 2.2 Chain Displacement, adiabatic
if strContains(O,'.xa.') && ~skipObs															% adiabatic projected displacement
	chainX = real(getObservable({'bath1correlators','x','adiabatic'}, treeMPS,[],para));		% L x nStates x nChains
	if isNew
		d = size(chainX);
		tresults.xa = zeros([totalN,d],'single');
	end
	tresults.xa(i,:,:,:) = chainX;																% (L x nStates x nChain)
end
if strContains(O,'.xc.') && ~skipObs															% coherence projected displacement
	chainX = getObservable({'bath1correlators','x','lettcoherence'}, treeMPS,[],para);			% L x 1 x nChains
	if isNew
		d = size(chainX);
		tresults.xc = zeros([totalN,d(1),d(3)],'single');
	end
	tresults.xc(i,:,:,:) = chainX(:,1,:);		% only take first state slice, since others are 0
end

% 3. Chain spread, squared displacement
if strContains(O,'.x2.') && ~skipObs															% displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2'}, treeMPS,[],para));					% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2 = zeros([totalN,d],'single');
	end
	tresults.x2(i,:,:) = chainX2;
end
if strContains(O,'.x2d.') && ~skipObs															% diabatic projected displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2','diabatic'}, treeMPS,[],para));		% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2d = zeros([totalN,d],'single');
	end
	tresults.x2d(i,:,:,:) = chainX2;
end
if strContains(O,'.x2a.') && ~skipObs															% adiabatic projected displacement squared <x^2>
	chainX2 = real(getObservable({'bath1correlators','x^2','adiabatic'}, treeMPS,[],para));		% L x nChains
	if isNew
		d = size(chainX2);
		tresults.x2a = zeros([totalN,d],'single');
	end
	tresults.x2a(i,:,:,:) = chainX2;
end

%% Star Observables

%% Special Observables

% 1. Density matrix
if strContains(O,'.dm.','.dm2.') && ~skipObs
	if isNew
		tresults.rho = zeros(totalN,treeMPS.dk(1,1),treeMPS.dk(1,1),'single');
	end
	tresults.rho(i,:,:,1) = single(getObservable({'rdm',1},treeMPS,[],para));
	if strContains(para.tdvp.Observables,'.dm2.')
		% only for 2-lvl system for now; only calculates largest bond state.
		tresults.rho(i,:,:,2) = single(getObservable({'rdm_adiabatic',1,1},tmps,tVmat,para));  %{'rdm_adiabatic',sitej,state}
		tresults.rho(i,:,:,3) = single(getObservable({'rdm_adiabatic',1,2},tmps,tVmat,para));  %{'rdm_adiabatic',sitej,state}
	end
end

if strContains(para.tdvp.Observables,'.ss.')					% ss for system state
	if isNew
		tresults.system.state = zeros(totalN,treeMPS.dk(1),treeMPS.dk(1),'single');		% t x dk x D (adiabatic)
	end
	tresults.system.state(i,:,:) = single(getObservable({'state',1},treeMPS,[],para));
end

if strContains(para.tdvp.Observables,'.heff.')
	if strContains(para.tdvp.Observables,'.pes.')
		heffPes = 1;										% skips system.state saving in .pes.
	end
	if isNew
		D = treeMPS.D(2); dk = treeMPS.dk(1);
		tresults.Heff			= zeros([totalN, dk, dk, dk, dk],'single');
		tresults.system.state	= zeros([totalN, dk, dk],'single');
		tresults.E				= zeros([totalN,1],'single');
	end
	D = size(tresults.Heff,2);
	temp = getObservable({'sysheff',size(tresults.Heff,2)},treeMPS,[],para);
	tresults.Heff(i,:,:,:,:)    = single(temp{1}(1:D,:,1:D,:));			% t x D' x dk' x D x dk
	tresults.system.state(i,:,:)= single(temp{2}(:,1:D));					% t x dk x D
	tresults.E(i)				= single(temp{3});
end

if strContains(para.tdvp.Observables,'.pes.')
	if isNew
		D = treeMPS.D(2); dk = treeMPS.dk(1);
		tresults.pes			= zeros([totalN, dk, dk, dk, dk],'single');
		tresults.Epot			= zeros([totalN,1],'single');
		if ~heffPes
			tresults.system.state	= zeros([totalN, dk, dk],'single');
		end
	end
	temp = getObservable({'syspes',size(tresults.Heff,2)},treeMPS,[],para);
	tresults.pes(i,:,:,:,:)			= single(temp{1});			% t x D' x dk' x D x dk
	if ~heffPes
		tresults.system.state(i,:,:)= single(temp{2});			% t x dk x D
	end
	tresults.Epot(i)				= single(temp{3});
end
%% TTM Extraction



%% End
if ~skipObs
	tresults.t(i) = single(para.tdvp.t(1,para.timeslice+1));
	tresults.lastIdx = tresults.lastIdx + 1;
end

if ~skipStar
	tresults.star.t(j) = single(para.tdvp.t(1,para.timeslice+1));
	tresults.star.lastIdx = tresults.star.lastIdx + 1;
end

end
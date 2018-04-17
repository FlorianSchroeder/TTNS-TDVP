function VMPS_MLSBM(period,eta)
%Variational matrix product method to study spin-boson model. Rok's
%logrithimic discretization and optimal boson basis are implemented in
%this code.
%
% modified to support a multi-level spin site with choosable coupling to the bath.


%Cheng Guo @Munich
% 3 Sep 2010


% mod by Florian Schroeder @Cambridge
% 02 May 2014 -
% copied from FullSBM 02/05/14


maxNumCompThreads(1);
format short e

starttime = tic
if isdeployed           % take care of command line arguments
%     if ischar(hx), hx = str2num(hx); end
%     if ischar(hz), hz = str2num(hz); end
%     if ischar(s), s = str2num(s); end
%     if ischar(alpha), alpha = str2num(alpha); end
%    if ischar(delta), delta = str2num(delta); end
%    if ischar(L), L = str2num(L); end
%    if ischar(Lambda), Lambda = str2num(Lambda); end
%    if ischar(parity), parity = str2num(parity); end
end

%% Parameters
para.model='MLSpinBoson';
    % choose: 'SpinBoson', '2SpinPhononModel', 'MLSpinBoson','ImpurityQTN'

para.Lambda=2;                          % Bath log Discretization parameter
para.z=1;                               % z-shift of bath; see Zitko 2009 - 10.1103/PhysRevB.79.085106
para.L=0;                               % Length per bath; if L=0: choose optimal chain length according to para.precision;
L=para.L;

if ~strcmp(para.model,'MLSpinBoson')
    para.hx=delta;                          % Splitting with sigma_X
    para.hz=0;                              % Splitting with sigma_Z
    para.s=s;                               % SBM spectral function power law behaviour
    para.alpha=alpha;                       % SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
    if para.alpha == 0
        para.L = 50;                        % otherwise encounter error
    end
end

%% Starting MPS Dimensions
D = 5;
dk = 20;
d_opt = 5;

if strcmp(para.model,'MLSpinBoson')     % definitions needed in SBM_genpara for spectral function & Wilson chain
    % Model Definition para.MLSB_mode:
    %   1:  from diagonalised, constant spacing Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
    para.MLSB_mode = 2;
end

para.foldedChain=0;                     % parameter to tell that Supersites for chain are used!
para.spinposition=1;                    % This indicates all positions ~= bosonic! important for Vmat! The y chain is on the left and the z chain is on the right. (could be array !)
para.rescaling=1;                       % rescale h1term, h2term for bosonchain with \lambda^{j-2}*h1term
para.complex=1;							% set to 1 if any complex parameters are used.
para.resume=0;                          % Read from saved results if available.
para.logging = 1;                       % Switch on logging and
parity = 0;
para.precision = 5e-15;                 % was 5e-15; Determines chain length if L=0;

%% %%%%%%% Calculate Wilson Chain parameters %%%%%%%%%%%%%%%%%%
% needed here: para.model, [para.MLSB_mode]
[para]=SBM_genpara(para);               % only need alpha, s, Lambda. Returns epsilon and t of Wilson chain. Auto choose L if L == 0
if (L == 0)
    L = para.L;                         % Take best L if not specially defined
end

%%
para.M=2;                               % is number of terms in sum to address. Better: Number of 2-operator interaction terms per site in Hamiltonian. =2 for single chain; =4 for folded chain
para.D=D*ones(1,L-1);                   % Bond dimension; starting dimension is 2. para.D(L) is useless in this program
para.dk_start = dk;                     % local dimension per boson in bath. Will be increased effectively by oscillator shift.
para.dk=para.dk_start*ones(1,L);
para.dk(1)=2;                           % Impurity dimension
para.increasedk = 0;					% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.
para.d_opt=d_opt*ones(1,L);             % Dimension of first site is 2 (spin); Optimal Boson Basis dimension, was 16*ones
para.d_opt(1)=2;                        % Optimal Impurity dimension
para.eigs_tol=1e-8;
para.loopmax=600;

if strcmp(para.model,'2SpinPhononModelfolded')      % could be removed
   para.Delta = 0;
   para.epsilon_spin = 0;
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

%% Multi-Level Spin Boson Model for PPC

if strcmp(para.model,'MLSpinBoson')
    % Model Definition para.MLSB_mode:
    %   1:  from diagonalised, constant spacing Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
    %       Can be affected with static random disorder

% All modes:
% System-Bath coupling:
%    para.MLSB_t = [-1 -1 -1 1 1 1];                 % coupling of each level to bath from range [-1,1]. Used as para.t(1)*para.MLSB_t
    para.MLSB_p = period;                           % period of coupling for cosine
    para.MLSB_etaFactor = eta;                      % multiplier to increase system-bath coupling
    para.MLSB_t = @(x) para.MLSB_etaFactor.*exp(1i*2*pi/para.MLSB_p.*x);                % MLSB_t can be an anonymous function.
    para.MLSB_tOff = 1;                             % 0=diag coupling; 1 = 1. off-diagonal n,n+1 coupling

    if para.MLSB_mode == 1
        para.MLSB_Ls = 6;                           % number of levels in System, symmetric around E=0;
        para.MLSB_Delta = para.hx;                  % level spacing in already diagonalized model
    end

    if para.MLSB_mode == 2
        para.MLSB_system = 'RsMolischianumB850';
            % Choose: 'RsMolischianumB850', 'RsMolischianumB800B850'
        para.MLSB_Ls = length(Hamiltonian_PPC(para));
            % Inherently defined within model
    end

%% random static disorder:
    para.MLSB_staticDisorder = 1;
    para.MLSB_disSDV = 200/2.3548;         % standard deviation of disorder; FWHM = 200 cm^-1; FWHM = 2.3548*SDV;
    % Vector with absolute values of disorder in cm^-1:
    para.MLSB_disDiag = para.MLSB_disSDV.*randn(para.MLSB_Ls,1);        % set to 0 if no diag disorder wanted

%% All modes:
    para.dk(1) = para.MLSB_Ls;                  % also set dk properly
    para.d_opt(1) = para.MLSB_Ls;

end
%%

para.SVDmethod = 'qr';                      % 'qr': uses QR instead of SVD wherever possible; 'svd': use SVD always (slower)
para.svmaxtol=1e-6;
para.svmintol=1e-8;                     %para.svmaxtol/2; %The lower limit for the smallest Vmat singular values.
para.adjust=0;                          %Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
para.Dmin = 4;                          % set a minimum Bond dimension.
para.dimlock=0;                         %set to 0 will change D and d_dop adaptively
para.minDimChange = 0.01;               % sets dimlock = 1 if relative dimension change < minDimChange. (larger makes less loops)

%'e' is even; 'o' is odd; 'n' is no parity
switch parity
    case 0
        para.parity='n';
    case 1
        para.parity='o';
    case 2
        para.parity='e';
end
para.spinbase='Z';
if para.parity~='n'
    para.Anzi=cell(para.L,1);
    para.Vnzi=cell(para.L,1);
    para.Anzilr=para.Anzi;
    para.Anzirl=para.Anzi;
end
%% %%%%%%%%%%%%%%%%%%Vmat related parameters%%%%%%%%%%%%%%%%%%%%
% Introduces Optimal Bosonic Basis (OBB)
para.useVmat=1;
if para.useVmat==0
     assert(para.dk_start==max(para.d_opt));
end
para.d_opt_min = 2;                                     % minimum d_opt dimension
%% %%%%%%%%%%%%%%%%%% Expansion of Vmat related parameters %%%%%%%
para.useexpand=1;			% Enable dk expansion, own algorithm. Don't use for folded chain!!
para.dkmax=500;
para.expandprecision =1e-5; % unused?
para.hasexpanded=0;
% Method 1:
para.useDkExpand1 = 0;      % Expand dk if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.expandBelowSV = 0.995;
% Method 2:
para.useDkExpand2 = 1;        % Expand if wavefunction on site occupies the high energy dimensions
para.dkEx2_tail   = 0.4;      % tail length of occupation to analyse
para.dkEx2_maxDev = 1.5;      % if std(log10(tail)) < maxDev --> no increase; Measures orders of magnitude in fluctuations of tail.
para.dkEx2_minExp = 13;       % if tail below this order than do not expand.
%% %%%%%%%%%%%%%%%%%%Shifting related parameters%%%%%%%%%%%%%%%%
% Introduces shift of bosonic oscillators to increase effective dk
para.useshift=0;
% only choose one of the following Methods
para.useFloShift = 0;                                   % shift all sites if trustsite > 0
para.useFloShift2 = 1;  para.FloShift2minSV = 0.995;    % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;                                   % shifts sites < trustsite
para.useEveryShift=0;                                   % shift in every loop

para.shift=zeros(1,para.L);
para.relativeshift=zeros(1,para.L);
para.relativeshiftprecision=0.01;        %When the relative shift is below this value then stop shifting
para=maxshift(para);

%% calculate shift analytically  (for independent sub-ohmic model)
%t = para.t; e = para.epsilon;
%A = gallery('tridiag',para.t(2:end),para.epsilon(1:end),para.t(2:end));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
%B = zeros(para.L-1,1);
%B(1) = -para.t(1)*1;
%para.shift =[0; A\B.*sqrt(2)]';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para.folder=sprintf([datestr(now,'yyyymmdd-HHMM'),'-%s-period%.10geta%.10gdk%.10gD%.10gdopt%gL%d'],...
    para.model,period,eta,dk,D,d_opt,L);
para.filename=strcat(para.folder,'/results.mat');
if ~exist(para.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end

%% Start Calculation
[op,para]=genh1h2term(para);
[mps, Vmat,para,results,op] = minimizeE(op,para);

save(para.filename,'para','Vmat','mps','results','op');

%% Calculate some Results
results.nx         = getObservable({'occupation'},mps,Vmat,para);
results.bosonshift = getObservable({'shift'},mps,Vmat,para);

if strcmp(para.model,'MLSpinBoson')
    results.participation = getObservable({'participation'},mps,Vmat,para);
    results.tunnelEnergy  = getObservable({'tunnelenergy',op},mps,Vmat,para);
end
%%
results.time = toc(starttime)
save(para.filename,'para','Vmat','mps','results','op');

end

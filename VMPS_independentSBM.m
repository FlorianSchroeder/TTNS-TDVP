function VMPS_independentSBM(s,alpha,L)
%Variational matrix product method to study spin-boson model. Rok's
%logrithimic discretization and optimal boson basis are implemented in
%this code.
% adapted to use independent SBM as benchmark system


%Cheng Guo @Munich
% 3 Sep 2010


% mod by Florian Schroeder @Cambridge
% 13 Jan 2014 - 19 Feb 2014


maxNumCompThreads(1);
format short e

starttime = tic
if isdeployed           % take care of command line arguments
    if ischar(hx), hx = str2num(hx); end
    if ischar(hz), hz = str2num(hz); end
    if ischar(s), s = str2num(s); end
    if ischar(alpha), alpha = str2num(alpha); end
    if ischar(L), L = str2num(L); end
%    if ischar(Lambda), Lambda = str2num(Lambda); end
%    if ischar(parity), parity = str2num(parity); end
end

%% Parameters
para.hx=0;                             % Splitting with sigma_X
para.hz=0;                             % Splitting with sigma_Z
para.s=s;                               % what???
para.alpha=alpha;                       % What??
para.Lambda=2;                     % Bath log Discretization parameter
para.L=L;                               % Length per bath; if L=0: choose optimal chain length according to para.precision;
para.z=1;                               % z-shift of bath
D = 5;
dk = 20;
d_opt = 5;

para.model='SpinBoson';
para.foldedChain=0;                    % parameter to tell that Supersites for chain are used!
para.spinposition=1;                    % The y chain is on the left and the z chain is on the right.
para.rescaling=1;			% rescale h1term, h2term for bosonchain with \lambda^{j-2}*h1term
para.complex=0;
para.resume=0;                          % Read from saved results if available.
para.logging = 1;                       % Switch on logging and
parity = 0;
para.precision = 5e-15;                 % was 5e-15;

%% %%%%%%% Calculate Wilson Chain parameters %%%%%%%%%%%%%%%%%%
[para]=SBM_genpara(para);               % only need alpha, s Lambda. Returns epsilon and t of Wilson chain. Auto choose L if L == 0
if (L == 0)
    L = para.L;                         % Take best L if not specially defined
end

%%
para.M=2;                               % is number of terms in sum to address. Better: Number of 2-operator interaction terms per site in Hamiltonian. =2 for single chain; =4 for folded chain
para.D=D*ones(1,L-1);                  % Bond dimension; starting dimension is 2. para.D(L) is useless in this program
para.dk_start = dk;                     % local dimension per boson in bath. Will be increased effectively by oscillator shift.
para.dk=para.dk_start*ones(1,L);
para.dk(1)=2;                           %Impurity dimension
para.increasedk = 0;					% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.
para.d_opt=d_opt*ones(1,L);                % Dimension of first site is 2 (spin); Optimal Boson Basis dimension, was 16*ones
para.d_opt(1)=2;                        %Optimal Impurity dimension
para.eigs_tol=1e-8;
para.loopmax=600;

if strcmp(para.model,'2SpinPhononModelfolded')
   para.Delta = 0;
   para.epsilon_spin = 0;
end

if strcmp(para.model,'2SpinPhononModel')
   %para.hz = para.hx;                  % use hz as level splitting. hx = hy = 0
   %para.hx = 0;
   para.hy = 0.1;                       % coupling between the 2 spin sites = J
   para.M = 4;				% 4 terms in sum per site
   para.dk(1) = 4;                      % As kron(site1,site2), dim=4 on first site!
   para.d_opt(1) = 4;
   para.foldedChain=1;
end

para.SVDmethod = 'qr';                      % 'qr': uses QR instead of SVD wherever possible; 'svd': use SVD always (slower)
para.svmaxtol=1e-6;
para.svmintol=1e-8;                     %para.svmaxtol/2; %The lower limit for the smallest Vmat singular values.
para.adjust=0;                          %Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
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
%% %%%%%%%%%%%%%%%%%%Expansion of Vmat related parameters%%%%%%%
para.useexpand=1;			% Was unused. Now to enable dk expansion, own algorithm.
para.expandBelowSV = 0.995; % expand if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.dkmax=50000;
para.expandprecision =1e-5; % unused?
para.hasexpanded=0;
%% %%%%%%%%%%%%%%%%%%Shifting related parameters%%%%%%%%%%%%%%%%
% Introduces shift of bosonic oscillators to increase effective dk
para.useshift=1;
% only choose one of the following Methods
para.useFloShift = 0;
para.useFloShift2 = 0;  para.FloShift2minSV = 0.99;     % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;
para.useEveryShift=1;

para.shift=zeros(1,para.L);
para.relativeshift=zeros(1,para.L);
para.relativeshiftprecision=0.01; %When the relative shift is below this value then stop shifting
para=maxshift(para);

%% calculate shift analytically
t = para.t; e = para.epsilon;
A = gallery('tridiag',para.t(2:end),para.epsilon(1:end),para.t(2:end));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
B = zeros(para.L-1,1);
B(1) = -para.t(1)*1;
para.shift =[0; A\B.*sqrt(2)]';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para.folder=sprintf([datestr(now,'yyyymmdd-HHMM'),'-%s-s%.10galpha%.10gdk%.10gD%.10gdopt%gL%d'],...
    para.model,s,alpha,dk,D,d_opt,L);
para.filename=strcat(para.folder,'/results.mat');
if ~exist(para.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end

%% Start Calculation
[op,para]=genh1h2term(para);
[mps, Vmat,para,results] = minimizeE(op,para);

save(para.filename,'para','Vmat','mps','results','op');

%% Calculate Results
results.nx                = getObservable({'occupation'},mps,Vmat,para);
results.bosonshift        = getObservable({'shift'},mps,Vmat,para);
results.bosonshiftPerSpin = calbosonshiftperSpin_SBM1(mps,Vmat,para,results);
results.spin              = getObservable({'spin'},mps,Vmat,para);

results.time = toc(starttime)
save(para.filename,'para','Vmat','mps','results','op');

%if isdeployed              % don't plot if called from command line.
    return;
%end

%% Plot Results
f1=figure(1);
subplot(2,2,1);
    plot(results.nx);
    title('$$<n_x(k)>$$');
subplot(2,2,2);
    plot(para.trustsite);
    title('Trustsite')
subplot(2,2,3);
    plot(para.shift);
    title('Bosonic shift');
subplot(2,2,4);
    surf(cell2mat(results.d_opt'));
    shading interp
    set(gca,'View',[0 90]);
    title('OBB dim')
saveas(f1,[para.folder,'/nx.fig'],'fig');
%saveas(f1,['png/',para.folder,'-nx.png'],'png');
export_fig(['png/',para.folder,'-nx.png'],'-transparent',f1)    % needs package export_fig
sprintf('Sx = %.10g, Sy = %.10g, Sz = %.10g',results.spin.sx,results.spin.sy,results.spin.sz)

%% Plot d_opt, D adjustments
f2 = figure(2);
subplot(1,2,1);
surf(cell2mat(results.D'))
shading interp
title('Change in bond dimension D')
subplot(1,2,2);
surf(cell2mat(results.d_opt'))
shading interp
title('Change in $d_{opt}$')
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
%export_fig(['png/',para.folder,'-D-dopt.png'],'-transparent',f2)

end

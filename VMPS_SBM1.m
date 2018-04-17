function VMPS_SBM1(hx,hz,s,alpha,Lambda,L,parity)
%Variational matrix product method to study spin-boson model. Rok's
%logrithimic discreitzation and optimal boson basis are implemented in
%this code.

%Cheng Guo @Munich
% 3 Sep 2010

maxNumCompThreads(1);
format short e

tic
if isdeployed           % take care of command line arguments
    if ischar(hx), hx = str2num(hx); end
    if ischar(hz), hz = str2num(hz); end
    if ischar(s), s = str2num(s); end
    if ischar(alpha), alpha = str2num(alpha); end
    if ischar(L), L = str2num(L); end
    if ischar(Lambda), Lambda = str2num(Lambda); end
    if ischar(parity), parity = str2num(parity); end
end

%Parameters
para.hx=hx;
para.hz=hz;
para.s=s;
para.alpha=alpha;
para.Lambda=Lambda;
para.L=L;
para.z=1;

para.model='SpinBoson';
para.spinposition=1; %The y chain is on the left and the z chain is on the right.
para.rescaling=1;
para.complex=0;
para.resume=1;  %Read from saved results if available.


para.M=2;
para.D=10*ones(1,L-1); %starting dimension is 2. para.D(L) is useless in this program
para.dk_start=100;
para.dk=para.dk_start*ones(1,L);
para.dk(1)=2; %Impurity dimension
para.d_opt=6*ones(1,L); %Starting dimension is 2
para.d_opt(1)=2; %Optimal Imprity dimentsion
para.precision = 5e-15;
para.eigs_tol=1e-8;
para.loopmax=200;

para.svmaxtol=1e-6;
para.svmintol=1e-8;%para.svmaxtol/2; %The lower limit for the smallest Vmat singular values.
para.adjust=0;%Initialize as 0, no need the edit.
para.dimlock=0; %set to 0 will change D and d_dop adaptively

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
%%%%%%%%%%%%%%%%%%%%Vmat related parameters%%%%%%%%%%%%%%%%%%%%
para.useVmat=1;
if para.useVmat==0
     assert(para.dk_start==max(para.d_opt));
end
%%%%%%%%%%%%%%%%%%%%Expansion of Vmat related parameters%%%%%%%
para.useexpand=0;
para.dkmax=5000;
para.expandprecision =1e-5;
para.hasexpanded=0;
%%%%%%%%%%%%%%%%%%%%Shifting related parameters%%%%%%%%%%%%%%%%
para.useshift=1;
para.shift=zeros(1,para.L);
para.relativeshift=zeros(1,para.L);
para.relativeshiftprecision=0.01; %When the relative shift is below this value the stop shifting
para=maxshift(para);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para.folder=sprintf('s%.10galpha%.10ghx%.10ghz%.10gp%g',para.s,para.alpha,para.hx,para.hz,parity);
para.filename=strcat(para.folder,'/results.mat');
if ~exist(para.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end


%%%%%%%%%Calculate%%%%%%%%%%%%%%%%%%
[para]=SBM_genpara(para);
[op,para]=genh1h2term(para);
[mps, Vmat,para,results] = minimizeE(op,para);


%Calculate Results
results.nx                = getObservable({'occupation'},mps,Vmat,para);
results.bosonshift        = getObservable({'shift'},mps,Vmat,para);
results.bosonshiftPerSpin = calbosonshiftperSpin_SBM1(mps,Vmat,para,results);
results.spin              = getObservable({'spin'},mps,Vmat,para);

save(para.filename,'para','Vmat','mps','results','op');

toc
end

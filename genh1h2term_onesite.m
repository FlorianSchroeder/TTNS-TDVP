function op = genh1h2term_onesite(para,op,s)
%% function op = genh1h2term_onesite(para,op,s)
%
%Define the hamiltonian
% op.h2term{i,j,k,l}:
%	l is chain number (only used in SpinBosonMC)
%   k is site number
%   i is number of term in sum to address:  t1*a*b + t1*b*a
%                       sum:                    1       2
%   j is position nr in each hopping term:  t*a*b
%       t1*a*b:     para.t(1)*op.h2term{1,1,1}*op.h2term{1,2,2}
%       t1*b*a:     para.t(1)*op.h2term{2,1,1}*op.h2term{2,2,2}
%example: op.h2term{1,1,1} is with op.h2term{1,2,2}
% interaction terms don't seem to be normal ordered...
%
% If Model changed, modify also:
%   VMPS_SBM1.m:    para.dk(1) gives dimension of first site;
%   calspin.m:      sx,sy,sz defines spin measure operator. Change if dim(site1) changes!

if (isfield(para,'useTreeMPS') && para.useTreeMPS) || ~isa(op,'struct')
	op = genh1h2term_onesite_tree(para,op,s);			% put into separate sub-function, op = treeIdx
else
switch para.model
    case 'SpinBoson'
        %%%%%%%%%%%%%%%%%%%Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ] = spinop(para.spinbase);       % gives XYZ operators with respect to specific main base
                zm_spin			  = zeros(2);
                op.h1term{1}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n]		  = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = para.chain{1}.epsilon(para.L-1).*n;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}	  = para.chain{1}.epsilon(s-1).*n;                          % e(1) == w(0)
                op.h2term{1,1,s}  = para.chain{1}.t(s).*bp; op.h2term{1,2,s} = bm;    % t(2) == t(n=0)
                op.h2term{2,1,s}  = para.chain{1}.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBosonTTM'
        %%%%%%%%%%%%%%%%%%% Spin-boson Model for the Transfer Tensor Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
			case 1                                                  % first chain pos = ancilla system, needs to be maximally entangled with site 2
                zm_spin			  = zeros(2);
                op.h1term{s}	  = zm_spin;
                op.h2term{1,1,s}  = zm_spin; op.h2term{1,2,s} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,s}  = zm_spin; op.h2term{2,2,s} = zm_spin;
            case 2                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ] = spinop(para.spinbase);       % gives XYZ operators with respect to specific main base
                zm_spin			  = zeros(2);
                op.h1term{s}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,s}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,s} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,s}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,s} = zm_spin;
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n]		  = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = para.chain{1}.epsilon(para.L-1).*n;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}	  = para.chain{1}.epsilon(s-2).*n;                          % e(1) == w(0)
                op.h2term{1,1,s}  = para.chain{1}.t(s-1).*bp; op.h2term{1,2,s} = bm;    % t(2) == t(n=0)
                op.h2term{2,1,s}  = para.chain{1}.t(s-1).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson2folded'
        %%%%%%%%%%%%%%%%%%%2-chain Spin Boson Model - One Super-Chain%%%%%%%%%%%%%%%%%%%%%%
		% Created 07/07/15 by F.S.
        switch s
            case 1
                [sigmaX,~,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                op.h1term{1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1} = zm_spin;	% X chain
                op.h2term{2,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1} = zm_spin;
                op.h2term{3,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1} = zm_spin;	% Z chain
                op.h2term{4,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1} = zm_spin;
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);  % gives [bp, bm, n]
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
				else
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
%                     [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{para.L}     = para.chain{1}.epsilon(para.L-1).*nx + para.chain{2}.epsilon(para.L-1).*nz;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmx;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpx;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bmz;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpz;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
				else
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
%                     [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
%                 zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{s}     = para.chain{1}.epsilon(s-1).*nx + para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s} = bmx;
                op.h2term{2,1,s} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s} = bpx;
                op.h2term{3,1,s} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s} = bmz;
                op.h2term{4,1,s} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s} = bpz;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoubleBoson'
        %%%%%%%%%%%%%%%%%%%Spin double boson Model One Chain%%%%%%%%%%%%%%%%%%%%%%
		% From Cheng, DEPRECATED, para.foldedChain = 1!
		% 2 Boson chains folded into supersites
        switch s
            case 1
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                %assert(para.Delta==0);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaX./2; op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1).*sigmaX./2; op.h2term{2,2,1} = zm_spin;   %to the right
                op.h2term{3,1,1} = para.t(1).*sigmaY./2; op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1).*sigmaY./2; op.h2term{4,2,1} = zm_spin;   %to the left
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);  % gives [bp, bm, n]
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
                else
                    [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*nx+para.epsilon(para.L-1).*ny;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmx;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpx;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bmy;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpy;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
                else
                    [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
%                 zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{s}=para.epsilon(s-1).*nx+para.epsilon(s-1).*ny;
                op.h2term{1,1,s} = para.t(s).*bpx; op.h2term{1,2,s} = bmx;
                op.h2term{2,1,s} = para.t(s).*bmx; op.h2term{2,2,s} = bpx;
                op.h2term{3,1,s} = para.t(s).*bpy; op.h2term{3,2,s} = bmy;
                op.h2term{4,1,s} = para.t(s).*bmy; op.h2term{4,2,s} = bpy;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoulbeBosonFolded' %%Two chain case. Doesn't work any more.
        %%%%%%%%%%%%%%%%%%Spin Doulbe Boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % From Cheng, DEPRECATED
		switch s
            case para.spinposition      %at the middle of the chain?
                if para.parity=='n'
                    sigmaX=[0 1;1 0]; %In the non parity basis
                    sigmaY=[0 -1i;1i 0];
                    sigmaZ=[1 0;0 -1];
                else
                    sigmaX=[-1 0;0 1]; %In the parity basis
                    sigmaZ=[0 -1;-1 0];
                end
                op.h1term{s}=-para.Delta./2.*sigmaX+para.epsilon_spin./2.*sigmaZ;
                op.h2term{1,1,s} = para.tz(1).*sigmaZ./2; op.h2term{1,2,s} = para.ty(1).*sigmaY./2;
                op.h2term{2,1,s} = para.tz(1).*sigmaZ./2; op.h2term{2,2,s} = para.ty(1).*sigmaY./2;
            case 1
                [bp,bm,n] = bosonop(para.dk(1),para.shift(1),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{1}=para.epsilony(para.Ly-1).*n;
                op.h2term{1,1,1} = bm; op.h2term{1,2,1} = zm;
                op.h2term{2,1,1} = bp; op.h2term{2,2,1} = zm;
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilonz(para.Lz-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                if s<para.spinposition
                    ss=para.spinposition-s;
                    [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                    op.h1term{s}=para.epsilony(ss).*n;
                    op.h2term{1,1,s} = bm; op.h2term{1,2,s} = para.ty(ss+1).*bp;
                    op.h2term{2,1,s} = bp; op.h2term{2,2,s} = para.ty(ss+1).*bm;
                else
                    ss=s-para.spinposition;
                    [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                    op.h1term{s} = para.epsilonz(ss).*n;
                    op.h2term{1,1,s} = para.tz(ss+1).*bp; op.h2term{1,2,s} = bm;
                    op.h2term{2,1,s} = para.tz(ss+1).*bm; op.h2term{2,2,s} = bp;
                end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'DissipativeHamonicOscillator'
        %%%%%%%%%%%%%%%%%%%Dissipative Hamonic Oscillator Model%%%%%%%%%%%%%%%%%%%
        switch s
            case 1
                [bp,bm,n] = bosonop(para.dk(1),para.shift(1));
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{1}=para.Omega.*n+para.epsilon_osc./2.*(bp+bm);
                op.h2term{1,1,1} = para.t(1).*(bp+bm)./2; op.h2term{1,2,1} = zm;
                op.h2term{2,1,1} = para.t(1).*(bp+bm)./2; op.h2term{2,2,1} = zm;
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L));
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s));
                op.h1term{s}=para.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case '2SpinPhononModel'
        %%%%%%%%%%%%%%%%%%%2-State Phonon Model (folded)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copied from SpinDoubleBoson. This is indeed the folded model, as the bosonic sites have kron()
		% possibly not working!
        switch s
            case 1      % if I want to include two exciton sites with chain-ex1-ex2-chain I need to have ex1,ex2 in site 1 combined!
                        % start off with only one site and extend it for later!
                        % possibly wrong: COUPLING in Prior, Chin 10: parallel to energy splitting! so sigmaZ
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                sigmaPlus = (sigmaX+1i.*sigmaY)./2;
                sigmaMinus = (sigmaX - 1i.*sigmaY)./2;
                zm_spin=zeros(4);
                %assert(para.Delta==0);
                %op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h1term{1}= para.hx./2.*kron(sigmaZ,eye(2)) + para.hz./2.*kron(eye(2),sigmaZ) + ...
                              para.hy.*(kron(sigmaPlus,sigmaMinus) + kron(sigmaMinus,sigmaPlus));
                %here: hx=epsilon1; hz=epsilon2; hy = J for Prior 2010 with kron(site1,site2);

                op.h2term{1,1,1} = para.t(1)./4.*kron(eye(2)+sigmaZ,eye(2)); op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1)./4.*kron(eye(2)+sigmaZ,eye(2)); op.h2term{2,2,1} = zm_spin; %to the right
                op.h2term{3,1,1} = para.t(1)./4.*kron(eye(2),eye(2)+sigmaZ); op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1)./4.*kron(eye(2),eye(2)+sigmaZ); op.h2term{4,2,1} = zm_spin; %to the left
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);
                if para.parity=='n'     % always use n (might be easier?)
                    idm=eye(size(n));
                    bpr=kron(bp,idm);   bmr=bpr';   nr=kron(n,idm);     % matrices for right    (=x)
                    bpl=kron(idm,bp);   bml=bpl';   nl=kron(idm,n);     % matrices for left     (=y)
                else
                    [bpr,bmr,nr,bpl,bml,nl]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpr,1),size(bpr,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*nr+para.epsilon(para.L-1).*nl;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmr;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpr;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bml;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpl;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpr=kron(bp,idm);   bmr=bpr';   nr=kron(n,idm);
                    bpl=kron(idm,bp);   bml=bpl';   nl=kron(idm,n);
                else
                    [bpr,bmr,nr,bpl,bml,nl]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpr,1),size(bpr,1));
                op.h1term{s}=para.epsilon(s-1).*nr+para.epsilon(s-1).*nl;
                op.h2term{1,1,s} = para.t(s).*bpr; op.h2term{1,2,s} = bmr;
                op.h2term{2,1,s} = para.t(s).*bmr; op.h2term{2,2,s} = bpr;
                op.h2term{3,1,s} = para.t(s).*bpl; op.h2term{3,2,s} = bml;
                op.h2term{4,1,s} = para.t(s).*bml; op.h2term{4,2,s} = bpl;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'MLSpinBoson'
        %%%%%%%%%%%%%%%%%%%Multi-Level Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for use of this model:
        %   para.hx: de
        %   HSI is a matrix defining coupling of system to bath. Is scaled with para.t(1)
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [HS0, HSI] = MLSB_Operators(para);                  % uses local function
                zm_spin=zeros(para.MLSB_Ls);
                op.h1term{1}= HS0;
                op.h2term{1,1,1} = para.chain{1}.t(1).*HSI;  op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi) only for J_Renger!
                op.h2term{2,1,1} = para.chain{1}.t(1).*HSI'; op.h2term{2,2,1} = zm_spin;		% was sigmaX; conj() to still have hermitian Hamiltonian!
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.chain{1}.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}=para.chain{1}.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.chain{1}.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.chain{1}.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ImpurityQTN'
        %%%%%%%%%%%%%%%%%%%Impurity coupled to Quantum Telegraph Noise Model%%%%%%%%%%%%%%%%%%%
        % for use of this model:
        %   TODO: this is by far not completed!!! and ready to use..
        %
		switch s
            case 1                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ]=spinop(para.spinbase);            % gives XYZ operators with respect to specific main base
                zm_spin=zeros(2);
                op.h1term{1}=para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% was sigmaX
                op.h2term{2,1,1} = para.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;		% was sigmaX
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}=para.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.t(s).*bm; op.h2term{2,2,s} = bp;
		end

	case 'SpinBoson2C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 2-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Benchmark for multi-chain method! Should yield same as
		%		'SpinBoson2folded'
		% Spin is always in chain 1 for backward compatibility
		%
		% working (23/07/15)
		% Created 15/07/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1,1} = zm_spin;	% X chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1,1} = zm_spin;	% Z chain
                op.h2term{4,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1,1} = zm_spin;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = bpz;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s,2} = bpz;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson2CXZ'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 2-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Chain 1: X; Chain 2: Z coupled
		% To be used with StarMPS for maximum speed
		% Spin is always in chain 1 for backward compatibility
		%
		% Created 16/06/17 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1,1} = zm_spin;	% X chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1,1} = zm_spin;	% Z chain
                op.h2term{4,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1,1} = zm_spin;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = bpz;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s,2} = bpz;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson2CT'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 2-Chain Thermal Bath %%%%%%%%%%%%%%%%%%%%%%
		% SBM including ancilla chain for thermal bath states
		% Multi-chain configuration for Vtens
		%
		% Spin is always in chain 1
		%
		% working 27/08/15
		% Created 25/08/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2); zm = []; %zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1,1} = zm_spin;	% Z chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = zm; op.h2term{3,2,1,1} = zm;	% ancilla chain
                op.h2term{4,1,1,1} = zm; op.h2term{4,2,1,1} = zm;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
% 				[bpz,  ~, ~] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = []; %sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = zmz;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = zmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = zmz;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
% 				[bpz,  ~, ~] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				zmz = []; %sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = zmz;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = zmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = zmz;							% replace zmz by [] ??
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson3C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 3-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Benchmark for multi-chain method! Should be similar to
		%		'SpinBoson2C'
		% Spin is always in chain 1 for backward compatibility
		%
		% working (23/07/15)
		% Created 17/07/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1,1} = zm_spin;	% X chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1,1} = zm_spin;	% Z chain
                op.h2term{4,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1,1} = zm_spin;
				op.h2term{5,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{5,2,1,1} = zm_spin;	% Z chain
                op.h2term{6,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{6,2,1,1} = zm_spin;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				[bpy,bmy,ny] = bosonop(para.dk(3,s),para.shift(3,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = sparse(size(bpz,1),size(bpz,1));
				zmy = sparse(size(bpy,1),size(bpy,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
				op.h1term{3,s}     = para.chain{3}.epsilon(s-1).*ny;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = bpz;
				op.h2term{5,1,s,3} = zmy; op.h2term{5,2,s,3} = bmy;
                op.h2term{6,1,s,3} = zmy; op.h2term{6,2,s,3} = bpy;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				[bpy,bmy,ny] = bosonop(para.dk(3,s),para.shift(3,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
				op.h1term{3,s}     = para.chain{3}.epsilon(s-1).*ny;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s,2} = bpz;
				op.h2term{5,1,s,3} = para.chain{3}.t(s).*bpy; op.h2term{5,2,s,3} = bmy;
                op.h2term{6,1,s,3} = para.chain{3}.t(s).*bmy; op.h2term{6,2,s,3} = bpy;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson5C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 5-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Spin is always in chain 1 for backward compatibility
		%
		% working??
		% Created 17/08/15 by F.S.
        switch s
            case 1		% will be pentacene system!
				[sigmaX,~,sigmaZ]   = spinop(para.spinbase);
				zm_spin			    = zeros(2);
				op.h1term{1,1}      = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
				op.h2term{1 ,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1 ,2,1,1} = zm_spin;	% X chain
				op.h2term{2 ,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2 ,2,1,1} = zm_spin;
				op.h2term{3 ,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{4 ,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4 ,2,1,1} = zm_spin;
				op.h2term{5 ,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{5 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{6 ,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{6 ,2,1,1} = zm_spin;
				op.h2term{7 ,1,1,1} = para.chain{4}.t(1).*sigmaZ./2; op.h2term{7 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{8 ,1,1,1} = para.chain{4}.t(1).*sigmaZ./2; op.h2term{8 ,2,1,1} = zm_spin;
				op.h2term{9 ,1,1,1} = para.chain{5}.t(1).*sigmaX./2; op.h2term{9 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{10,1,1,1} = para.chain{5}.t(1).*sigmaX./2; op.h2term{10,2,1,1} = zm_spin;
            case para.L
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				for i = 1:5			% slow, but easy to modify!
					[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
					op.h2term{2*i-1,1,s,i} = zm; op.h2term{2*i-1,2,s,i} = bm;
					op.h2term{2*i  ,1,s,i} = zm; op.h2term{2*i  ,2,s,i} = bp;
				end
			otherwise
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				for i = 1:5
					[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
					op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
					op.h2term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h2term{2*i-1,2,s,i} = bm;
					op.h2term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h2term{2*i  ,2,s,i} = bp;
				end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'DPMES3-4C'
	%%%%%%%%%%%%%%%%%%% DP-MES Model - 4-Chain %%%%%%%%%%%%%%%%%%%%%%
	% Not linear, but in multi-chain configuration!
	% Spin is always in chain 1 for backward compatibility
	%
	% working
	% Created 05/10/15 by F.S.
	switch s
		case 1		% is the pentacene system!
			[H0,H1] = DPMES_Operators('3-4C',para);
			zm      = zeros(3);
			op.h1term{1,1}      = H0;
			for i = 1:4
				op.h2term{2*i-1, 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i-1, 2,1,1} = zm;
				op.h2term{2*i  , 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i  , 2,1,1} = zm;
			end
% 				op.h2term{1 ,1,1,1} = para.chain{1}.t(1).*H1{1}; op.h2term{1 ,2,1,1} = zm;	% A1 chain 1
% 				op.h2term{2 ,1,1,1} = para.chain{1}.t(1).*H1{1}; op.h2term{2 ,2,1,1} = zm;
% 				op.h2term{3 ,1,1,1} = para.chain{2}.t(1).*H1{2}; op.h2term{3 ,2,1,1} = zm;	% A1 chain 2
% 				op.h2term{4 ,1,1,1} = para.chain{2}.t(1).*H1{2}; op.h2term{4 ,2,1,1} = zm;
% 				op.h2term{5 ,1,1,1} = para.chain{3}.t(1).*H1{3}; op.h2term{5 ,2,1,1} = zm;	% B1 chain
% 				op.h2term{6 ,1,1,1} = para.chain{3}.t(1).*H1{3}; op.h2term{6 ,2,1,1} = zm;
% 				op.h2term{7 ,1,1,1} = para.chain{4}.t(1).*H1{4}; op.h2term{7 ,2,1,1} = zm;	% A2 chain
% 				op.h2term{8 ,1,1,1} = para.chain{4}.t(1).*H1{4}; op.h2term{8 ,2,1,1} = zm;
		case para.L
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:4			% slow, but easy to modify!
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				zm = sparse(size(bp,1),size(bp,1));
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = zm; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = zm; op.h2term{2*i  ,2,s,i} = bp;
			end
		otherwise
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:4
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h2term{2*i  ,2,s,i} = bp;
			end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'DPMES4-5C'
	%%%%%%%%%%%%%%%%%%% DP-MES Model - 5-Chain %%%%%%%%%%%%%%%%%%%%%%
	% Not linear, but in multi-chain configuration!
	% Spin is always in chain 1 for backward compatibility
	%
	% working
	% Created 04/10/15 by F.S.
	switch s
		case 1		% is the pentacene system!
			[H0,H1] = DPMES_Operators('4-5C',para);
			zm      = zeros(size(H0,1));
			op.h1term{1,1}      = H0;
			for i = 1:length(H1)
				op.h2term{2*i-1, 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i-1, 2,1,1} = zm;
				op.h2term{2*i  , 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i  , 2,1,1} = zm;
			end
		case para.L
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains			% slow, but easy to modify!
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				zm = sparse(size(bp,1),size(bp,1));
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = zm; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = zm; op.h2term{2*i  ,2,s,i} = bp;
			end
		otherwise
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h2term{2*i  ,2,s,i} = bp;
			end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	case 'DPMES5-7C'
	%%%%%%%%%%%%%%%%%%% DP-MES Model - 5-Chain %%%%%%%%%%%%%%%%%%%%%%
	% Not linear, but in multi-chain configuration!
	% Spin is always in chain 1 for backward compatibility
	%
	% working
	% Created 29/01/16 by F.S.
	switch s
		case 1		% is the pentacene system!
			[H0,H1] = DPMES_Operators('5-7C',para);
			zm      = zeros(size(H0,1));
			op.h1term{1,1}      = H0;
			for i = 1:length(H1)
				op.h2term{2*i-1, 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i-1, 2,1,1} = zm;
				op.h2term{2*i  , 1,1,1} = para.chain{i}.t(1).*H1{i}./sqrt(2); op.h2term{2*i  , 2,1,1} = zm;
			end
		case para.L
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains			% slow, but easy to modify!
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				zm = sparse(size(bp,1),size(bp,1));
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = zm; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = zm; op.h2term{2*i  ,2,s,i} = bp;
			end
		otherwise
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h2term{2*i  ,2,s,i} = bp;
			end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'Holstein-Cavity-16'
	%%%%%%%%%%%%%%%%%%% Holstein - Cavity Model - 16 Sites %%%%%%%%%%%%%%%%%%%%%%
	% Multi-chain configuration.
	% Cavity + Excitons are in site 1
	%
	% working
	% Created 18/01/16 by F.S.
	switch s
		case 1		% is the pentacene system!
			[H0,H1] = DPMES_Operators('4-5C',para);
			zm      = zeros(size(H0,1));
			op.h1term{1,1}      = H0;
			for i = 1:length(H1)
				op.h2term{2*i-1, 1,1,1} = para.chain{i}.t(1).*H1{i}; op.h2term{2*i-1, 2,1,1} = zm;
				op.h2term{2*i  , 1,1,1} = para.chain{i}.t(1).*H1{i}; op.h2term{2*i  , 2,1,1} = zm;
			end
		case para.L
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains			% slow, but easy to modify!
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				zm = sparse(size(bp,1),size(bp,1));
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = zm; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = zm; op.h2term{2*i  ,2,s,i} = bp;
			end
		otherwise
			if para.parity ~= 'n'
				error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			for i = 1:para.nChains
				[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
				op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
				op.h2term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h2term{2*i-1,2,s,i} = bm;
				op.h2term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h2term{2*i  ,2,s,i} = bp;
			end
	end
	
	case 'UniformBosonTTM'
	%%%%%%%%%%%%%%%%%%% Uniform Boson Chain Model for Transfer Tensor Method %%%%%%%%%%%%%%%%%%%%%%
	% site 1 is Boson ancilla
	% site 2 is actual boson
	% epsilon and t are equal for all Boson sites!
	% working ??
	% Created 20/01/16 by F.S.
        switch s
			case 1                                                  % first chain pos = ancilla system, needs to be maximally entangled with site 2
                [bp,~,~]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = zm;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,s} = zm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,s} = zm;
            case 2                                                  % second chain pos = Boson to be analysed
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));	% zero matrix
                op.h1term{s}	  = para.chain{1}.epsilon(1).*n;
                op.h2term{1,1,s}  = para.chain{1}.t(1).*bp; op.h2term{1,2,s} = zm;
                op.h2term{2,1,s}  = para.chain{1}.t(1).*bm; op.h2term{2,2,s} = zm;
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n]		  = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = para.chain{1}.epsilon(1).*n;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}	  = para.chain{1}.epsilon(1).*n;							% e(1) == w
                op.h2term{1,1,s}  = para.chain{1}.t(1).*bp; op.h2term{1,2,s} = bm;			% t(1) == t
                op.h2term{2,1,s}  = para.chain{1}.t(1).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	otherwise
		error('VMPS:genh1h2term_onesite:ModelNotFound','Could not find specified para.model')
end
end
end

function op = genh1h2term_onesite_tree(para,treeIdx,s)
%% function op = genh1h2term_onesite_tree(para,treeIdx,s)
%
%	Generates Hamiltonian terms for a treeMPS
%
%	treeIdx: 1 x N_levels array with tree-Index of node/leaf
%	s:       # site in leaf
%
% returns for treeIdx = node of tree
%	op.h1term: 1 x 1 cell
%	op.h2term: M x 2 x N_edges cell
%			op.h2term(#term,#position in term,#edge to couple to)
% returns for treeIdx = leaf of tree = chain, only single-site terms
%	op.h1term: 1 x 1 cell
%	op.h2term: M x 2 cell
%			op.h2term(#term,#position in term)
%
% #position in term:
%		1	couples to the right
%		2	couples to the left
%
% created 22/02/2016 by F.S.
%
switch para.model
	case 'SpinBoson'
		%%%%%%%%%%%%%%%%%%% Spin-boson Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	Modified 10/04/2017
		%
		op.h1term = {};
		op.h2term = cell(para.M,2);
		if treeIdx == 0
			[sigmaX,~,sigmaZ] = spinop(para.spinbase);			% gives XYZ operators with respect to specific main base
			zm_spin			  = zeros(2);
			op.h1term{1}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
			op.h2term{1,1}	  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
			op.h2term{2,1}    = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2} = zm_spin;
		else
			mc  = treeIdx;					% number of edge == chain number
			idx = treeIdx+1;				% idx = num2cell(treeIdx+1); index in para.*
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx}(s),para.shift{idx}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s
				% only need s+1 here since L got shifted by 1 in treeMPS
				case para.chain{mc}.L                                             % last chain pos: only one coupling?
					zm				  = sparse(size(bp,1),size(bp,1));
					op.h1term{1}	  = para.chain{1}.epsilon(s).*n;
					op.h2term{1,1}    = zm; op.h2term{1,2} = bm;
					op.h2term{2,1}    = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}	  = para.chain{1}.epsilon(s).*n;									% e(1) == w(0)
					op.h2term{1,1}    = para.chain{1}.t(s+1).*bp; op.h2term{1,2} = bm;					% t(2) == t(n=0)
					op.h2term{2,1}    = para.chain{1}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	case 'SpinBoson2C'
		%%%%%%%%%%%%%%%%%%% Spin-boson Model 2-chains %%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% working
		% Created 03/16 by FS
		op.h1term = {};
		op.h2term = cell(para.M,2);
		if treeIdx == 0
			[sigmaX,~,sigmaZ] = spinop(para.spinbase);			% gives XYZ operators with respect to specific main base
			zm_spin			  = zeros(2);
			op.h1term{1}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
			op.h2term{1,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;
			op.h2term{1,1,2}  = para.chain{2}.t(1).*sigmaZ./2; op.h2term{1,2,2} = zm_spin;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,2}  = para.chain{2}.t(1).*sigmaZ./2; op.h2term{2,2,2} = zm_spin;
			
		else
			mc  = treeIdx;					% number of edge == chain number
			idx = treeIdx+1;				% idx = num2cell(treeIdx+1); index in para.*
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx}(s),para.shift{idx}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s
				case para.chain{mc}.L											% last chain pos: only coupling to left
					zm				  = sparse(size(bp,1),size(bp,1));
					op.h1term{1}	  = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1}    = zm; op.h2term{1,2} = bm;
					op.h2term{2,1}    = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}	  = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1}    = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;
					op.h2term{2,1}    = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	case 'DPMES5-7C'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% working?
		% Created 22/02/16 by F.S.
		%
		op.h1term = {};
		op.h2term = cell(para.M,2);
		if treeIdx == 0
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('5-7C',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end

		else
			mc  = treeIdx;					% number of edge == chain number
			idx = treeIdx+1;				% idx = num2cell(treeIdx+1); index in para.*
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx}(s),para.shift{idx}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	case 'DPMESclust7-1'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% working?
		% Created 08/04/17 by F.S.
		%
		op.h1term = {};
		op.h2term = cell(para.M,2);
		if treeIdx == 0
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('clust7-1',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end

		else
			mc  = treeIdx;					% number of edge == chain number
			idx = treeIdx+1;				% idx = num2cell(treeIdx+1); index in para.*
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx}(s),para.shift{idx}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'DPMES-Tree1'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain Tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	Initial try to map DPMES5-7C onto an entanglement renormalisation tree
		%	Channel interaction terms through nodes without Hamiltonian terms/sites
		%   ChainIdx is numbered according to tree structure, while Number in tree diagram corresponds to number as returned by DPMES_Operators('5-7C')
		%	Structure:
		%
		%							treeIdx		Node Idx	ChainIdx		
		%	Excitons				[0,0,0,0]	1
		%	  |- Node  			    [1,0,0,0]	2
		%     | |- Node             [1,1,0,0]	3
		%     | | |- Chain 4        [1,1,1,0]				1		
		%	  | | \- Chain 1        [1,1,2,0]				2
		%	  | \- Chain 2			[1,2,0,0]				3
		%	  \- Node  				[2,0,0,0]	4
		%		|- Node  			[2,1,0,0]	5
		%       | |- Node           [2,1,1,0]	6
		%       | | |- Chain 5      [2,1,1,1]				4
		%       | | \- Chain 6      [2,1,1,2]				5
		%	    | \- Chain 7        [2,1,2,0]				6
		%		\- Chain 3			[2,2,0,0]				7
		%
		% Created 29/08/16 by F.S.
		%
		idx      = num2cell(treeIdx+1);                       % index in para.*
		mc       = para.treeMPS.chainIdx{idx{:}};             % index of chain
% 		leafPos  = find(all(bsxfun(@eq,leafIndices,idx),2));
% 		nodePos  = find(all(bsxfun(@eq,leafIndices,idx),2));
		
		op.h1term = {};
		op.h2term = cell(para.M,2);
		
		% First: list all nodes which carry a physical Hamiltonian site
		if treeIdx == [0,0,0,0]
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('Tree1',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end
% 		elseif
% 			idx = num2cell(treeIdx+1);
		% Second: list all chains
		elseif ~isempty(mc)
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		% Third: all nodes without site, only combining chains
		else
			% Do nothing! return empty, since this node only combines chain and has no Hamiltonian site!
			op.h1term = {[]};																		% put empty array to avoid errors
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	case 'DPMES-Tree2'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain Tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	Initial try to map DPMES5-7C onto an entanglement renormalisation tree
		%	Channel interaction terms through nodes without Hamiltonian terms/sites
		%   ChainIdx is numbered according to tree structure, while Number in tree diagram corresponds to number as returned by DPMES_Operators('5-7C')
		%	Structure:
		%
		%							treeIdx		Node Idx	ChainIdx		
		%	Excitons				[0,0,0,0]	1
		%	  |- Chain 1		    [1,0,0,0]				1
		%     \- Node               [2,0,0,0]	2
		%	    |- Node				[2,1,0,0]	3
		%		| |- Chain 2		[2,1,1,0]				2
		%		| \- Chain 3		[2,1,2,0]				3
		%		\- Node		        [2,2,0,0]	4
		%		  |- Node			[2,2,1,0]	5
		%		  | |- Chain 4		[2,2,1,1]				4
		%		  | \- Chain 5		[2,2,1,2]				5
		%		  \- Node			[2,2,2,0]	6
		%		    |- Chain 6		[2,2,2,1]				6
		%		    \- Chain 7		[2,2,2,2]				7
		%
		% Created 29/08/16 by F.S.
		%
		idx      = num2cell(treeIdx+1);                       % index in para.*
		mc       = para.treeMPS.chainIdx{idx{:}};             % index of chain
% 		leafPos  = find(all(bsxfun(@eq,leafIndices,idx),2));
% 		nodePos  = find(all(bsxfun(@eq,leafIndices,idx),2));
		
		op.h1term = {};
		op.h2term = cell(para.M,2);
		
		% First: list all nodes which carry a physical Hamiltonian site
		if treeIdx == [0,0,0,0]
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('Tree1',para);		% equals Tree2 operators, so no problem
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end
		% Second: list all chains
		elseif ~isempty(mc)
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		% Third: all nodes without site, only combining chains
		else
			% Do nothing! return empty, since this node only combines chain and has no Hamiltonian site!
			op.h1term = {[]};																		% put empty array to avoid errors
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	case 'DPMES-Tree3'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain Tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	First entanglement renormalisation Tree for DPMESclust7-1 parameter set.
		%	It is the best partition if the first partition is into 2 subsystems: [Sys,1,2,3] [4,5,6,7]
		%	Channel interaction terms through nodes without Hamiltonian terms/sites
		%   ChainIdx is numbered according to tree structure
		%	Structure:
		%
		%							treeIdx		Node Idx	ChainIdx		
		%	Excitons				[0,0,0,0]	1
		%	  |- Chain 1		    [1,0,0,0]				1
		%     \- Node               [2,0,0,0]	2
		%	    |- Node				[2,1,0,0]	3
		%		| |- Chain 2		[2,1,1,0]				2
		%		| \- Chain 3		[2,1,2,0]				3
		%		\- Node		        [2,2,0,0]	4
		%		  |- Node			[2,2,1,0]	5
		%		  | |- Chain 4		[2,2,1,1]				4
		%		  | \- Chain 5		[2,2,1,2]				5
		%		  \- Node			[2,2,2,0]	6
		%		    |- Chain 6		[2,2,2,1]				6
		%		    \- Chain 7		[2,2,2,2]				7
		%
		% Created 29/08/16 by F.S.
		%
		idx      = num2cell(treeIdx+1);                       % index in para.*
		mc       = para.treeMPS.chainIdx{idx{:}};             % index of chain
% 		leafPos  = find(all(bsxfun(@eq,leafIndices,idx),2));
% 		nodePos  = find(all(bsxfun(@eq,leafIndices,idx),2));
		
		op.h1term = {};
		op.h2term = cell(para.M,2);
		
		% First: list all nodes which carry a physical Hamiltonian site
		if treeIdx == [0,0,0,0]
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('Tree3',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				% Interaction terms carry 1/sqrt(2) since Davids parameters are made for x = (a+a^+)/sqrt(2)
				% Thus only relevant for DPMES!
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end
		% Second: list all chains
		elseif ~isempty(mc)
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		% Third: all nodes without site, only combining chains
		else
			% Do nothing! return empty, since this node only combines chain and has no Hamiltonian site!
			op.h1term = {[]};																		% put empty array to avoid errors
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	case 'DPMES-Tree4'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 7-Chain Tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	Second entanglement renormalisation Tree for DPMESclust7-1 parameter set.
		%	It is the best partition if the first partition is into 3 subsystems: [Sys,B11,B12] [A11,A12,A2] [B22,B23]
		%	Channel interaction terms through nodes without Hamiltonian terms/sites
		%   ChainIdx is numbered according to tree structure
		%	This code is fairly generic and can be applied to almost any tree strucures! Perhaps refactor into separate function?
		%	Structure:
		%
		%							treeIdx		Node Idx	ChainIdx		
		%	Excitons				[0,0,0,0,0]	1
		%	  |- Chain 1		    [1,0,0,0,0]				1
		%     \- Node2              [2,0,0,0,0]	2
		%	    |- Chain 2			[2,1,0,0,0]				2
		%		\- Node3	        [2,2,0,0,0]	3
		%		  |- Node4			[2,2,1,0,0]	4
		%		  | |- Node5		[2,2,1,1,0]	5
		%		  | | |- Chain 3	[2,2,1,1,1]				3
		%		  | | \- Chain 4	[2,2,1,1,2]				4
		%		  | \- Chain 5		[2,2,1,2,0]				5
		%		  \- Node6			[2,2,2,0,0]	6
		%		    |- Chain 6		[2,2,2,1,0]				6
		%		    \- Chain 7		[2,2,2,2,0]				7
		%
		% Created 09/04/17 by F.S.
		%
		idx      = num2cell(treeIdx+1);                       % index in para.*
		mc       = para.treeMPS.chainIdx{idx{:}};             % index of chain
% 		leafPos  = find(all(bsxfun(@eq,leafIndices,idx),2));
% 		nodePos  = find(all(bsxfun(@eq,leafIndices,idx),2));
		
		op.h1term = {};
		op.h2term = cell(para.M,2);
		
		% First: list all nodes which carry a physical Hamiltonian site
		if treeIdx == [0,0,0,0,0]
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('Tree4',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				% Interaction terms carry 1/sqrt(2) since Davids parameters are made for x = (a+a^+)/sqrt(2)
				% Thus only relevant for DPMES!
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end
		% Second: list all chains
		elseif ~isempty(mc)
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		% Third: all nodes without site, only combining chains
		else
			% Do nothing! return empty, since this node only combines chain and has no Hamiltonian site!
			op.h1term = {[]};																		% put empty array to avoid errors
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	case 'DPMES-Tree6'
		%%%%%%%%%%%%%%%%%%% DP-MES Model - 9-Chain Tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	Second entanglement renormalisation Tree for DPMESclust7-1 parameter set.
		%	Channel interaction terms through nodes without Hamiltonian terms/sites
		%   ChainIdx is numbered according to tree structure
		%	This code is fairly generic and can be applied to almost any tree strucures! Perhaps refactor into separate function?
		%	Structure:
		%
		%							treeIdx		Node Idx	ChainIdx		
		%	Excitons				[0,0,0,0,0]	1
		%	  |- Chain 1		    [1,0,0,0,0]				1
		%     \- Node2              [2,0,0,0,0]	2
		%	    |- Chain 2			[2,1,0,0,0]				2
		%		\- Node3	        [2,2,0,0,0]	3
		%		  |- Node4			[2,2,1,0,0]	4
		%		  | |- Node5		[2,2,1,1,0]	5
		%		  | | |- Chain 3	[2,2,1,1,1]				3
		%		  | | \- Chain 4	[2,2,1,1,2]				4
		%		  | \- Chain 5		[2,2,1,2,0]				5
		%		  \- Node6			[2,2,2,0,0]	6
		%		    |- Chain 6		[2,2,2,1,0]				6
		%		    \- Chain 7		[2,2,2,2,0]				7
		%
		% Created 09/04/17 by F.S.
		%
		idx      = num2cell(treeIdx+1);                       % index in para.*
		mc       = para.treeMPS.chainIdx{idx{:}};             % index of chain
% 		leafPos  = find(all(bsxfun(@eq,leafIndices,idx),2));
% 		nodePos  = find(all(bsxfun(@eq,leafIndices,idx),2));
		
		op.h1term = {};
		op.h2term = cell(para.M,2);
		
		% First: list all nodes which carry a physical Hamiltonian site
		if treeIdx == [0,0,0,0,0]
			% is the pentacene system!
			[H0,H1]               = DPMES_Operators('Tree6',para);
			zm                    = zeros(size(H0,1));
			op.h1term{1}          = H0;
			for mc = 1:length(H1)
				% Interaction terms carry 1/sqrt(2) since Davids parameters are made for x = (a+a^+)/sqrt(2)
				% Thus only relevant for DPMES!
				op.h2term{1,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{1,2,mc} = zm;
				op.h2term{2,1,mc} = para.chain{mc}.t(1).*H1{mc}./sqrt(2); op.h2term{2,2,mc} = zm;
			end
		% Second: list all chains
		elseif ~isempty(mc)
			if para.parity ~= 'n'
						error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
			end
			if iscell(para.dk) && iscell(para.shift)
				[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
			else
				[bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
			end
			switch s						% this is 1:L on chain
				case para.chain{mc}.L
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = zm; op.h2term{1,2} = bm;
					op.h2term{2,1} = zm; op.h2term{2,2} = bp;
				otherwise
					op.h1term{1}   = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1} = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;				% t(1) already couples to node
					op.h2term{2,1} = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		% Third: all nodes without site, only combining chains
		else
			% Do nothing! return empty, since this node only combines chain and has no Hamiltonian site!
			op.h1term = {[]};																		% put empty array to avoid errors
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	case 'testTree'
		%%%%%%%%%%%%%%%%%%% testTree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	This is a test tree with uniform boson chains.
		%
		%	Equivalent to 2 coupled spins, each coupled to 2 chains
		%
		%	Structure:
		%
		%							treeIdx
		%	spin					0
		%	  |- chain 1			1
		%	  |- chain 2			2
		%	  \- spin				3
		%		|- chain 3			[3,1]
		%		\- chain 4			[3,2]
		%
		op.h1term = {};
		op.h2term = cell(para.M,2);
		if treeIdx == 0
			% spin 1
			[sigmaX,~,sigmaZ] = spinop(para.spinbase);			% gives XYZ operators with respect to specific main base
			zm_spin			  = zeros(2);
			op.h1term{1}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
			op.h2term{1,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;
			op.h2term{1,1,2}  = para.chain{2}.t(1).*sigmaZ./2; op.h2term{1,2,2} = zm_spin;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,2}  = para.chain{2}.t(1).*sigmaZ./2; op.h2term{2,2,2} = zm_spin;
			op.h2term{1,1,3}  = para.alpha.*sigmaZ; op.h2term{1,2,3} = zm_spin;
			op.h2term{2,1,3}  = para.alpha.*sigmaX; op.h2term{2,2,3} = zm_spin;
		elseif nonzeros(treeIdx) == 3
			% spin 2
			%	op.h2term{:,2,1} couples to root node!
			[sigmaX,~,sigmaZ] = spinop(para.spinbase);			% gives XYZ operators with respect to specific main base
			zm_spin			  = zeros(2);
			op.h1term{1}	  = -para.hx./2.*sigmaZ-para.hz./2.*sigmaX;
			op.h2term{1,1,1}  = para.chain{3}.t(1).*sigmaZ./2; op.h2term{1,2,1} = sigmaZ;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,1}  = para.chain{3}.t(1).*sigmaZ./2; op.h2term{2,2,1} = sigmaX;
			op.h2term{1,1,2}  = para.chain{4}.t(1).*sigmaZ./2; op.h2term{1,2,2} = zm_spin;		% t(1) = sqrt(eta_0/pi)
			op.h2term{2,1,2}  = para.chain{4}.t(1).*sigmaZ./2; op.h2term{2,2,2} = zm_spin;
		else
			idx = num2cell(treeIdx+1);														% index in para.*
			mc  = para.treeMPS.chainIdx{idx{:}};						% chain number from linear index in para
			switch s
				case para.chain{mc}.L											% last chain pos: only coupling to left
					[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
					zm				  = sparse(size(bp,1),size(bp,1));
					op.h1term{1}	  = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1}    = zm; op.h2term{1,2} = bm;
					op.h2term{2,1}    = zm; op.h2term{2,2} = bp;
				otherwise
					[bp,bm,n]		  = bosonop(para.dk{idx{:}}(s),para.shift{idx{:}}(s),para.parity);
					op.h1term{1}	  = para.chain{mc}.epsilon(s).*n;
					op.h2term{1,1}    = para.chain{mc}.t(s+1).*bp; op.h2term{1,2} = bm;
					op.h2term{2,1}    = para.chain{mc}.t(s+1).*bm; op.h2term{2,2} = bp;
			end
		end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
end
end

function [H0, H1] = MLSB_Operators(para)
% calculates Energy levels H0 and couplings to bath HI
% define modes and parameters in VMPS_MLSBM
%
% perhaps export to file
switch para.MLSB_mode
    case 1
        %   Define spacings Delta between each level:  para.MLSB_Delta
        %   Define couplings System-Bath as vector:         para.MLSB_t
        %   Define Size of ML-System:                       para.MLSBM_Ls
        %   Energy levels will be symmetrically aligned around 0
        %   TODO: next-neighbour couping within system
        assert(length(para.MLSB_t) == para.MLSB_Ls, 'All couplings t between system and bath must be defined');
        H0 = diag(((para.MLSB_Ls:-1:1)- (para.MLSB_Ls+1)/2).*para.MLSB_Delta);
		H0 = diag([0, cumsum(para.MLSB_Delta)]);
        H1 = diag(para.MLSB_t);
    case 2
        % Hamiltonian with rotational symmetry. Get from separate function
        % untis in eV
        [H0, H1] = Hamiltonian_PPC(para);
	case 3
		% Hamiltonian from matrices
		H0 = para.MLSB_H0;
		H1 = para.MLSB_H1;
    otherwise
end

end

function [H0, H1] = DPMES_Operators(nModel,para)
%% creates the Hamiltonian terms for the DPMES molecule.
states = para.systemStates;
switch nModel
	case '3-4C'
		% TT, LE+, CT+ with 4 chains
		% chain order: 1-2: A1(1,2); 3: B1, 4: A2
		H0 = diag(states([1,2,4],2));
		n = size(H0,1);
		H1 = cell(4,1);
		H1{1} = eye(n);   H1{1}(1,1) = 0;
		H1{2} = eye(n);   H1{2}(1,1) = 2;
		H1{3} = zeros(n); H1{3}(2,3) = 1; H1{3}(3,2) = 1;
		H1{4} = zeros(n); H1{4}(1,3) = 1; H1{4}(3,1) = 1;
	case '4-5C'
		% TT, LE+, CT+, CT- with 4 chains
		% chain order: 1-2: A1(1,2); 3: B1, 4: A2, 5:B2
		H0 = diag(states([1,2,4,5],2));
		n = size(H0,1);
		H1 = cell(5,1);
		H1{1} = eye(n);   H1{1}(1,1) = 0;
		H1{2} = eye(n);   H1{2}(1,1) = 2;
		H1{3} = zeros(n); H1{3}(2,3) = 1;          H1{3}(3,2) = 1;	        % B1, W24
						  H1{3}(1,4) = -sqrt(3)/2; H1{3}(4,1) = -sqrt(3)/2;	%	, W15
		H1{4} = zeros(n); H1{4}(1,3) = 1;          H1{4}(3,1) = 1;			% A2, W14
		H1{5} = zeros(n); H1{5}(3,4) = 1;          H1{5}(4,3) = 1;			% B2, W45
	case '5-7C'
		% TT, LE+, LE-, CT+, CT- with 4 chains
		% chain order: 1-2: A1(1,2); 3: A2, 4: B1, 5-7: B2
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		H1{1} = eye(n);   H1{1}(1,1) = 0;									% A1: TT/rest = 0
		H1{2} = eye(n);   H1{2}(1,1) = 2;									% A1: TT/rest = 2
		H1{3} = zeros(n); H1{3}(1,4) = 1;          H1{3}(4,1) = 1;			% A2, W14
		H1{4} = zeros(n); H1{4}(2,4) = 1;          H1{4}(4,2) = 1;	        % B1, W24
						  H1{4}(1,5) = -sqrt(3)/2; H1{4}(5,1) = -sqrt(3)/2;	%	, W15
						  H1{4}(3,5) = -1;         H1{4}(5,3) = -1;			%   , W35
		H1{5} = zeros(n); H1{5}(2,3) = 1;          H1{5}(3,2) = 1;			% B2, W23
		H1{6} = zeros(n); H1{6}(4,5) = 1;          H1{6}(5,4) = 1;			% B2, W45
						  H1{6}(2,3) = 1.3;        H1{6}(3,2) = 1.3;		%   , W23
		H1{7} = zeros(n); H1{7}(4,5) = 1;          H1{7}(5,4) = 1;			% B2, W45
						  H1{7}(2,3) = -1.5;       H1{7}(3,2) = -1.5;		%   , W23
	case 'clust7-1'
		% TT, LE+, LE-, CT+, CT- with 7 chains
		% order: {	'B11',	'B12',	'A11',	'A12',	'B22',  'B23',	'A2'};
		% H1 can be loaded directly from para.chain{ii}.H1
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		for ii = 1:7
			H1{ii} = para.chain{ii}.H1;
		end
	case 'Tree1'
		% TT, LE+, LE-, CT+, CT- with 4 chains
		% chain order: 1-2: A1(1,2); 3: A2, 4: B1, 5-7: B2
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		H1{1} = zeros(n); H1{1}(2,4) = 1;          H1{1}(4,2) = 1;	        % B1, W24
						  H1{1}(1,5) = -sqrt(3)/2; H1{1}(5,1) = -sqrt(3)/2;	%	, W15
						  H1{1}(3,5) = -1;         H1{1}(5,3) = -1;			%   , W35
		H1{2} = eye(n);   H1{2}(1,1) = 0;									% A1: TT/rest = 0
		H1{3} = eye(n);   H1{3}(1,1) = 2;									% A1: TT/rest = 2
		H1{4} = zeros(n); H1{4}(2,3) = 1;          H1{4}(3,2) = 1;			% B2, W23
		H1{5} = zeros(n); H1{5}(4,5) = 1;          H1{5}(5,4) = 1;			% B2, W45
						  H1{5}(2,3) = 1.3;        H1{5}(3,2) = 1.3;		%   , W23
		H1{6} = zeros(n); H1{6}(4,5) = 1;          H1{6}(5,4) = 1;			% B2, W45
						  H1{6}(2,3) = -1.5;       H1{6}(3,2) = -1.5;		%   , W23
		H1{7} = zeros(n); H1{7}(1,4) = 1;          H1{7}(4,1) = 1;			% A2, W14
	case 'Tree2'
		% TT, LE+, LE-, CT+, CT- with 4 chains
		% 1: B1, 2: A11, 3: A12, 4: B21, 5: B22, 6: B23, 7: A2
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		H1{1} = zeros(n); H1{1}(2,4) = 1;          H1{1}(4,2) = 1;	        % B1, W24
						  H1{1}(1,5) = -sqrt(3)/2; H1{1}(5,1) = -sqrt(3)/2;	%	 , W15
						  H1{1}(3,5) = -1;         H1{1}(5,3) = -1;			%    , W35
		H1{2} = eye(n);   H1{2}(1,1) = 0;									% A11: TT/rest = 0
		H1{3} = eye(n);   H1{3}(1,1) = 2;									% A12: TT/rest = 2
		H1{4} = zeros(n); H1{4}(2,3) = 1;          H1{4}(3,2) = 1;			% B21, W23
		H1{5} = zeros(n); H1{5}(4,5) = 1;          H1{5}(5,4) = 1;			% B22, W45
						  H1{5}(2,3) = 1.3;        H1{5}(3,2) = 1.3;		%    , W23
		H1{6} = zeros(n); H1{6}(4,5) = 1;          H1{6}(5,4) = 1;			% B23, W45
						  H1{6}(2,3) = -1.5;       H1{6}(3,2) = -1.5;		%    , W23
		H1{7} = zeros(n); H1{7}(1,4) = 1;          H1{7}(4,1) = 1;			% A2 , W14
	case 'Tree3'
		% TT, LE+, LE-, CT+, CT- with 7 chains
		% order: {	'A2',	'A11',	'A12',	'B11',	'B12',	'B22',  'B23'};
		% H1 can be loaded directly from para.chain{ii}.H1
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		for ii = 1:7
			H1{ii} = para.chain{ii}.H1;
		end
	case 'Tree4'
		% TT, LE+, LE-, CT+, CT- with 7 chains
		% order: {	'A2',	'A11',	'A12',	'B11',	'B12',	'B22',  'B23'};
		% H1 can be loaded directly from para.chain{ii}.H1
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(7,1);								% one for each chain!
		for ii = 1:7
			H1{ii} = para.chain{ii}.H1;
		end
	case 'Tree6'
		% TT, LE+, LE-, CT+, CT- with 9 chains
		% order: {	'B22',	'B23',	'B11',	'B12',	'B13',	'B14',  'A2',	'A11',	'A12'};
		% H1 can be loaded directly from para.chain{ii}.H1
		H0 = diag(states([1,2,3,4,5],2));
		n = size(H0,1);
		H1 = cell(9,1);								% one for each chain!
		for ii = 1:9
			H1{ii} = para.chain{ii}.H1;
		end
		
end

end

function [H0, H1] = Holstein_Operators(nModel,para)
%% creates the Hamiltonian terms for the Holstein B850 model.
%	Define in k-space to reduce number of bosons by factor of 2
states = para.systemStates; %??
switch nModel
	case '4C'
		% TT, LE+, CT+ with 4 chains
		% chain order: 1-2: A1(1,2); 3: B1, 4: A2
% 		H0 = diag(states([1,2,4],2));
% 		n = size(H0,1);
% 		H1 = cell(n,1);
% 		H1{1} = eye(n);   H1{1}(1,1) = 0;
% 		H1{2} = eye(n);   H1{2}(1,1) = 2;
% 		H1{3} = zeros(n); H1{3}(2,3) = 1; H1{3}(3,2) = 1;
% 		H1{4} = zeros(n); H1{4}(1,3) = 1; H1{4}(3,1) = 1;
	case 16
		% 16 excitons with 8 boson chains
% 		H0 =
		n = size(H0,1);			% = 17
		H1 = cell(8,1);
		H1{1} = eye(n);   H1{1}(1,1) = 0;		% k = 0, diagonal coupling
		H1{2} = sparse();   H1{2}(1,1) = 2;
		H1{3} = zeros(n); H1{3}(2,3) = 1;          H1{3}(3,2) = 1;	        % B1, W24
						  H1{3}(1,4) = -sqrt(3)/2; H1{3}(4,1) = -sqrt(3)/2;	%	, W15
		H1{4} = zeros(n); H1{4}(1,3) = 1;          H1{4}(3,1) = 1;			% A2, W14
		H1{5} = zeros(n); H1{5}(3,4) = 1;          H1{5}(4,3) = 1;			% B2, W45
end

end

function op = SBM_Star_Op(para, op, s) 
	%% op = SBM_Star_Op(para, op, s) 
	%	generates the operators for a multi-environment Spin Boson Model for use with StarMPS
	
	NC = para.nChains;
	
end
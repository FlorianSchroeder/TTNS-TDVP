function [mps, Vmat, para, results, op, Hn] = tdvp_1site_evolveHnMC(mps,Vmat,para,results,op,sitej)
%% Evolves one site following Haegeman 2014
%   - Only contains A = exp(-i H(n) dt/2) A;
%   - Same procedure for l->r or l<-r
%   - Splits time-evolution into A and V part if para.useVmat == 1
%   - Specifically designed for use of HOSVD and Multiple Chains
%
% Created by Florian Schroeder @ Cambridge 17/08/2015

% for spinsites: use the older tdvp_1site_evolveHn:
if any(sitej == para.spinposition)
	[mps, Vmat, para, results, op, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej);
	return;
end
assert(para.useVmat == 1,'Code only works with Vmat enabled!');
assert(iscell(Vmat{sitej}), 'Vmat needs to be a cell');
% From now on: only for Boson sites!
M = size(op.h2j,1);
NC = para.nChains;

if sitej ~= para.L
    t = para.tdvp.deltaT./2;
else
    t = para.tdvp.deltaT;
end

if para.tdvp.imagT
	t = -1i*t;
end

% create working copy of Vmat{sitej}:
Vtens = Vmat{sitej};                    % V tensor network

%% 0. Expand OBB by p %
%% expand OBB in V by 20%, A by 40%
% since this expansion is temporarily, save change later.
% expand always BEFORE SVD
for mc = 1:NC
	dk = para.dk(mc,sitej); dOBB = para.d_opt(mc,sitej);
    if (dk ~= dOBB) && para.tdvp.expandOBB && para.tdvp.imagT == 0			% do not expand for TH since it only cools down!
		expandBy = min([ceil(dOBB*0.2),para.tdvp.maxOBBDim-dOBB,dk-dOBB]);
		if expandBy ~= 0
			Vtens{mc}(end,end+expandBy) = 0;
			d = size(Vtens{end}); d(mc) = expandBy;
			Vtens{end} = cat(mc, Vtens{end}, zeros(d));
			para.d_opt(mc,sitej) = size(Vtens{mc},2);
		end
    end
end
d = size(mps{sitej});
dOBB1 = prod(para.d_opt(1:end-1,sitej)); dOBB2 = para.d_opt(end,sitej);
if (dOBB1-2 >= dOBB2) && para.tdvp.expandOBB
	expandBy = min([ceil(dOBB2*0.4),para.tdvp.maxOBBDim*2-dOBB2,dOBB1-dOBB2, d(1)*d(2)-dOBB2]);			% expand by at least double as above!
	if expandBy == 1, expandBy = 2; end
	if expandBy >= 1
		mps{sitej}(end,end,end+expandBy) = 0;
		d = size(Vtens{end}); d(end) = expandBy;
		Vtens{end} = cat(NC+1, Vtens{end}, zeros(d));
		para.d_opt(end,sitej) = size(Vtens{end},NC+1);
	end
end



%% 1. SVD to move focus from MPS onto V_S
% Structure of Vmat{sitej} : {V(1), V(2), ..., V(NC), VS}
[Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);		% left-normalize A, SVD in n.
para.d_opt(end,sitej) = size(Amat,3);						% could have changed -> save
Vtens{end} = contracttensors(Vtens{end}, NC+1, NC+1, V, 2, 2);
para.Vtens.focus  = NC+1;                                   % current focus on VS

op = H_Eff(Amat, [] , 'V' , op, para);						% create effective H terms for V ([H/Op][left/right]A)
op = H_Eff([], Vtens, 'MC-OBB', op, para);					% create chain-wise OBB terms

%% time-evolve all chain-V and CV matrices
for mc = 1:NC
	% enter each chain!
	para.currentChain = mc;
	Vtens = prepare_Vtens(Vtens, para);                    % move focus to Vtens{currentChain}
	if isempty(regexp(para.model,'SpinBoson\dCT','once')) || mod(mc,2)
		op = H_Eff([], Vtens, 'MC-V', op, para);               % create [H/Op][left/right]AS - terms and op.HnonInt
		[dk,dOBB] = size(Vtens{mc});

		%% Take matrix exponential V{i}
		% V{i}(t+dt) = exp(-i HAS_i dt)_(n',n~',n,n~) * V{i}(t)_(n,n~)
		[Vtens{mc}, err] = expvCustom(- 1i*t,'MC-HAS',Vtens{mc}, para, op);
		Vtens{mc} = reshape(Vtens{mc}, [dk,dOBB]);
		results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);
	end
	%% Move focus to VC{mc} & prepare H_Eff
	if para.tdvp.expandOBB
		minDim = max(para.d_opt_min, ceil(para.d_opt(end,sitej)^(1/NC)));
		[Vtens{mc}, CV, ~, para.d_opt(mc,sitej), sv] = prepare_onesiteVmat(Vtens{mc},para,[],[],minDim);         % SVD + truncation
	else
		[Vtens{mc}, CV, ~, para.d_opt(mc,sitej), sv] = prepare_onesiteVmat(Vtens{mc},para,results,sitej);		 % SVD
	end
	results.Vmat_sv{mc,sitej} = sv;
% TODO: delete empty dimensions! Only if expanded previously!
% 	% remove empty SV:
% 	keep = results.Vmat_sv{sitej} ~= 0;
% 	results.Vmat_sv{sitej} = results.Vmat_sv{sitej}(keep);
% 	Vmat{sitej} = Vmat{sitej}(:,keep); V = V(keep, :);
% 	[n1, n2] = size(V);
% 	OBBDimNew = n1;
	if isempty(regexp(para.model,'SpinBoson\dCT','once')) || mod(mc,2)
		op = H_Eff([], Vtens, 'MC-CV', op, para);                              % create h[1/2]jMCOBB
		[n1,n2] = size(CV);


		%% Take matrix exponential CV{i}
		% CV{i}(t) = exp(+i HAVS dt) * CV{i}(t+dt)
		[CV, err] = expvCustom(+ 1i*t,'MC-HASV',CV, para, op);
		CV = reshape(CV,[n1,n2]);
		results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);
	end
	% contract CV with VS and finish
	Vtens{end} = contracttensors(CV,2,2, Vtens{end},NC+1,mc);
	Vtens{end} = permute(Vtens{end}, [2:mc,1,(mc+1):(NC+1)]);
end

% op = H_Eff([],[],'MC-VS',op,para);                                         % create h12jAV for VS % not needed anymore!
n = size(Vtens{end});
%% Take matrix exponential VS
% VS(t+dt) = exp(-i HAV dt) * VS(t)
[Vtens{end}, err] = expvCustom(- 1i*t,'MC-HAV',Vtens{end}, para, op);
Vtens{end} = reshape(Vtens{end}, [prod(n(1:end-1)),n(end)]);                % 2D reshape to take focus away
results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

% Take focus to CVS
if para.tdvp.expandOBB
	[Vtens{end}, CVS, ~, n(end), sv] = prepare_onesiteVmat(Vtens{end},para,[],[],para.d_opt_min);
else
	[Vtens{end}, CVS, ~, n(end), sv] = prepare_onesiteVmat(Vtens{end},para,results,sitej);
end
results.Vmat_sv{NC+1,sitej} = sv;
results.Vmat_vNE(1,sitej)   = vonNeumannEntropy(diag(sv));
para.d_opt(end,sitej) = n(end);
Vtens{end} = reshape(Vtens{end}, n);                                        % now reshape completely
Vmat{sitej} = Vtens;
% Vtens = Vmat{sitej} is now completely evolved and normalised

% create effective OBB terms
op = H_Eff([],Vtens, 'MC-A',op,para);                                   % only need VS and all h[1/2]jMCOBB terms
n = size(CVS);

%% Take matrix exponential CVS
% CVS(t) = exp(+i HAVS dt) * CVS(t+dt)
[CVS, err] = expvCustom(+ 1i*t,'MC-HAVS',CVS, para, op);
CVS = reshape(CVS,n);
results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

mps{sitej} = contracttensors(Amat, 3, 3, CVS, 2, 2);
[BondDimLeft, BondDimRight, OBBDim]  = size(mps{sitej});

%% Take and apply Matrix exponential
% A(t+dt) = exp(-i Hn dt)_(l',r',n',l,r,n) * A(t)_(l,r,n)

[mpsNew,err] = expvCustom(- 1i*t, 'Hn',mps{sitej}, para,op);
Hn = [];		% dummy return value;
mps{sitej} = reshape(mpsNew,[BondDimLeft,BondDimRight,OBBDim]);
results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

% now: A and V are time-evolved, A is focused
% if sitej = L, then start lr sweep with decomposition of mps

% Only return current-site matrices!
mps = mps{sitej};
Vmat = Vmat{sitej};

end

%% NESTED FUNCTIONS
%   can be exported if needed elsewhere
function Vtens = prepare_Vtens(Vtens, para)
%% Input: focused VS = Vtens{end}
%  sets focus onto Vtens{mc}
%  no truncation for now!
%
mc = para.currentChain;                                % the chain to focus on
d = size(Vtens{end});

[VS]       = tensShape(Vtens{end}, 'unfold', mc, d);   % dOBB(mc) x rest

[VS, CV]   = prepare_onesiteVmat(VS.', para);          % need for transpose since m x n, m > n input needed

Vtens{end} = tensShape(VS.', 'fold', mc, d);
Vtens{mc}  = Vtens{mc} * CV.';

end

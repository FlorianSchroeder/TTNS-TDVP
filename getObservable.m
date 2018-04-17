function out = getObservable(type,mps,Vmat,para)
%% Calculates the following Observables:
%	SBM:    'spin', 'occupation', 'current', 'shift', 'rdm', 'staroccupation', 'starpolaron',
%			'energy', 'displacement', 'bath2correlators'
%   MLSBM:  'participation', 'tunnelenergy', 'staroccupation', 'starpolaron', 'energy'
%
%   Use as: first argument is cell, with descriptor as type{1}
%               all other type{n} are options defined in each function
%           getObservable({'spin'}					,mps,Vmat,para)
%           getObservable({'occupation'}			,mps,Vmat,para)
%           getObservable({'current',AnAm}			,mps,Vmat,para)
%           getObservable({'shift'}					,mps,Vmat,para)
%           getObservable({'rdm',2}					,mps,Vmat,para)
%           getObservable({'tunnelenergy',op}		,mps,Vmat,para)
%           getObservable({'displacement'}			,mps,Vmat,para)
%           getObservable({'bath2correlators'}		,mps,Vmat,para)
%           getObservable({'starpolaron'}			,mps,Vmat,para)
%           getObservable({'staroccupation',AnAm}	,mps,Vmat,para)
%           getObservable({'energy',op}				,mps,Vmat,para)
%           getObservable({'hshi',op}				,mps,Vmat,para)		% op is optional
%			getObservable({'stateproject',i,j}      ,mps,Vmat,para)		% i,j optional
%			getObservable({'sysheff',i}				,mps,Vmat,para)		% i optional; possibly only TreeMPS operational
%
% In 'current' and 'staroccupation': AnAm optional!
%
%   The output 'out' can be of different form:
%       'spin':             struct, fields: sx, sy, sz
%       'occupation':       array	(L x NC)
%       'shift':            array
%       'rdm':              matrix
%       'participation':    number
%       'tunnelenergy':     number
%		'displacement':		vector	(L x 1)
%		'bath2correlators':	matrix	(L x L x NC)
%		'staroccupation':	array	(3 x X x NC)
%		'current':			array	(NC x L-1)
%		'starpolaron':		array	(3 x X x NC)
%		'energy'			scalar
%		'hshi'				scalar
%		'coherence'			not implemented yet
%		'stateproject'		scalar
%		'state'				matrix  (dk x dk)
%
%   Created 03/06/2014 by Florian Schroeder @ University of Cambridge
%   TODO:   - Implement Boson Site Shift to export the routine from optimizesite.m. Only get a single shift value.
%
% Modified:
%	- 21/12/14 FS: replaced OBB contractions with faster matrix products.
%	- 21/08/15 FS: added Multi-Chain Vmat / Vtens capability
%	- 28/02/16 FS: added TreeMPS support
switch type{1}
    case 'spin'
        % applicable for spin-boson model and for folded SBM2
        out = calSpin(mps,Vmat,para);

    case 'occupation'
        % applicable to all single stranded chains. Extensible to 2-chains
		if para.useTreeMPS == 1
			out = calBosonOcc_Tree(mps,para);			% (L x NC)
			return;
		elseif para.foldedChain == 0
			if para.nChains == 1 && ~para.useStarMPS
				out = calBosonOcc(mps,Vmat,para);
			else										% Multi-chain with Vmat / V-tensor-network
				out = calBosonOcc_MC(mps,Vmat,para);	% (NC x L)
			end
		elseif para.foldedChain == 1
			out(1,:) = calBosonOcc(mps,Vmat,para,1);
			out(2,:) = calBosonOcc(mps,Vmat,para,2);
		end
		if size(out,1) > 1
			out = out.';								% (L x NC) for compatibility with tresults
		end

    case 'shift'
        % applicable to all single stranded chains. Extendable to 2-chains
        out = calBosonShift(mps,Vmat,para);

	case '1siteshift'
		% intended for use in optimizesite. Only for focused single-site
		% mps and Vmat are single-site matrices
		out = cal1siteShift(mps,Vmat,para);
		
    case 'rdm'
        % applicable to all systems.
		if type{2} <= para.L
			if para.useTreeMPS
				out = calRDM(mps.mps,mps.Vmat,para,1);	% only for site 1 for now.
			else
	            out = calRDM(mps,Vmat,para,type{2});
			end
		else
            out = [];
		end
		
	case 'rdm_adiabatic'
		% applicable to site 1 for now only!
		if type{2} ~= 1
			error('rdm_adiabatic only available for site 1');
		end
		if para.nChains > 1
			error('rdm_adiabatic only for single chains!');
		end
		assert( type{3} <= para.dk(1,1), 'Need valid number for state selection');
		k = type{3};		% the (adiabatic) bond state to construct the (adiabatic) RDM from
		
		% SVD of mps{1} and delete unwanted bond states!
		A = mps{1};
		[D1, D2, d] = size(A);				% d is site dimension
		A = permute(A, [3, 1, 2]);
		A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
		[B, S, U] = svd2(A);				% Could also use QR decomposition if nargin !=5
		DB = size(S,1);
		B = B* sparse(k,k,1,DB,DB) *U;
		B = reshape(B, [d, D1, D2]);
		mps{1} = permute(B, [2, 3, 1]);		% now only contains the desired bond state
		
		out = calRDM(mps,Vmat,para,type{2});
		
	case 'rdm_adiabatic2'
		% applicable to site 1 for now only!
		% A different way of selecting the bond-state, without SVD.
		% Seems to be equivalent to 'rdm_adiabatic'
		if type{2} ~= 1
			error('rdm_adiabatic only available for site 1');
		end
		if para.nChains > 1
			error('rdm_adiabatic only for single chains!');
		end
		assert( type{3} <= para.dk(1,1), 'Need valid number for state selection');
		k = type{3};		% the (adiabatic) bond state to construct the (adiabatic) RDM from
		
		[~, D2, ~] = size(mps{1});				% d is site dimension
		for ii = 1:D2
			if ii ~= k
				mps{1}(1,ii,:) = 0;				% delete 2 for k=1 or 1 for k=2
			end
		end
		mps{1} = mps{1}./norm(reshape(mps{1},[],1));		% normalise the state!
		
		out = calRDM(mps,Vmat,para,type{2});
		
	case 'state'
		% applicable to site 1 for now only!
		% saves the MPS of the site with diabatic & adiabatic information!
		% closely related to 'rdm_adiabatic'
		% {'state',sitej}
		if type{2} ~= 1
			error('VMPS:getObservable:state','state only available for site 1');
		end
		if para.nChains > 1 || para.useTreeMPS
			if para.useStarMPS
				d = size(mps{1});
				A = reshape(mps{1},[],d(end));		% ... x dk
				[U,S,~] = svd2(A.');				% transpose to benefit from speedup in svd2
			elseif para.useTreeMPS
				d = size(mps.mps{1});
				[~,V] = prepare_onesiteAmat(reshape(mps.mps{1},d(1),[],d(end)),para,1);		% A_(1,D(1)*..*D(NC),dOBB') * V_(dOBB',dk)
				if diff(size(V)) ~= 0														% if dOBB < dk then expand
					V(d(end),d(end)) = 0;
				end
				out = V';
				return;
			else
				error('VMPS:getObservable:NotImplemented','Needs to be implemented!')
			end
		else
			d = size(mps{1});
			A = reshape(mps{1},[],d(end));		% ... x dk
			[U,S,~] = svd2(A.');				% transpose to benefit from speedup in svd2
			if size(S,1) < d(end)
				% add zeros
				U(d(end),d(end)) = 0;
				S(d(end),d(end)) = 0;
			end
		end
		out = U*S;								% should be dk x D(1:dk) for D > dk
		
	case 'sys-env-state'
		% saves the MPS of site 1 & 2 containing all diabatic & adiabatic information!
		% decompose also 1st site to have adiabatic information available
		% returns out: struct
		% {'sys-env-state'}
		if para.nChains > 1
			if para.useStarMPS
				d = size(mps{1});
				A = reshape(mps{1},[],d(end));		% ... x dk
				[U,S,~] = svd2(A.');				% transpose to benefit from speedup in svd2
			else
				error('VMPS:getObservable:NotImplemented','Needs to be implemented!')
			end
		else
			d = size(mps{1});
			A = reshape(mps{1},[],d(end));		% ... x dk
			[U,S,~] = svd2(A.');				% transpose to benefit from speedup in svd2
			if size(S,1) < d(end)
				% add zeros
				U(d(end),d(end)) = 0;
				S(d(end),d(end)) = 0;
			end
		end
		Vmat{1} = U*S;							% should be dk x D(1:dk) for D > dk; sets focus on Vmat!! This is simpler for StarMPS
		% save the system and reaction coord state
		out.mps = mps([1,2]);					
		out.Vmat = Vmat([1,2]);
		
    case 'participation'
        % applicable to all systems, only for first site.
        out = calParticipation(calRDM(mps,Vmat,para,1));
		
    case 'tunnelenergy'
        % useful for MLSBM, applicable to all
        % needs type{2} = op
		if length(type) == 2
            out = calTunnelingEnergy(mps,Vmat,para,type{2});
		else
            out = calTunnelingEnergy(mps,Vmat,para);
		end

	case 'bath1correlators'
		% = [state projected] chain observable
		% returns a L [x dk(1,1)] x nChains Vector
		%
		%	obs: ['bp'|'x'|'bp^2'|'x^2'|'n']	selects observable to calculate
		%
		% by default:
		% {'bath1correlators',obs}            : no projection - TODO!
		% {'bath1correlators',obs,'diabatic'} : projects onto diabatic states
		% {'bath1correlators',obs,'adiabatic'}: projects onto adiabatic states
		
		NC = para.nChains;
		
		if para.useTreeMPS && isstruct(mps)
			if length(type) == 2
				out = calBath1SiteCorrelators_Tree(mps,para,type{2});		% L x NC
			else
				out = zeros(para.L, mps.dk(1), NC);
				% One SVD to dk to separate diabatic or adiabatic basis!
				% Assume, Vmat of node = eye; mps == treeMPS
				assert(all(all(mps.Vmat{1} == eye(mps.dk(1)))),'Please implement for Vmat ~= eye, if this happens!');
				d = size(mps.mps{1});
				[A,V] = prepare_onesiteAmat(reshape(mps.mps{1},d(1),[],d(end)),para,1);		% A_(1,D(1)*..*D(NC),dOBB') * V_(dOBB',dk)
				A = reshape(A,[],size(A,3));
				% Now dOBB is adiabatic basis, dk is diabatic. V carries probabilities for each state, which will be removed later!
				for kk = 1:mps.dk(1)
					newV = zeros(size(V));
					if strcmpi(type{3},'diabatic')
					% Project onto root nodes' diabatic states
						newV(:,kk) = V(:,kk);
					elseif strcmpi(type{3},'adiabatic')
					% Project onto root nodes' adiabatic states
						newV(kk,:) = V(kk,:);
					elseif strcmpi(type{3},'lettcoherence')		% get |TT><LE+| mode-mediated coherence terms
						newV = V;
						projOp = zeros(size(V,2));
						projOp(1,2) = 1;% projOp(2,1) = 1;
					end
					if norm(newV) ~= 0
						newV = newV./norm(newV);				% remove weight by /norm()
					end
					mps.mps{1} = reshape(A*newV,d);
					if strcmpi(type{3},'lettcoherence')
						out(:,kk,:) = calBath1SiteCorrelators_Tree(mps,para,type{2},projOp);		% projection needed for cross terms
						break;																		% only do once!
					else
						out(:,kk,:) = calBath1SiteCorrelators_Tree(mps,para,type{2});				% no projection needed anymore inside calBath1SiteCorrelators
					end
				end
			end
			return;
		end
		
		if length(type) == 3
			out = zeros(para.L, para.dk(1,1), NC);
		elseif length(type) == 2
			out = zeros(para.L,NC);
		else
			error('VMPS:getObservable:bath1correlators','Please call with correct input');
		end
		
		if para.useStarMPS
			for mc = 1:NC
				% extract MPS for each chain and calculate <anam> separately
				cL = para.chain{mc}.L;
				mpsC  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
				VmatC = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
				mpsC{1} = permute(mpsC{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
				mpsC{1} = reshape(mpsC{1}, [],para.D(mc,1),para.dk(1,1));

				paraC = para;
				paraC.nChains = 1; paraC.L = para.chain{mc}.L;
				paraC.dk = para.dk(mc,:); paraC.dk(1) = para.dk(1,1);	% need to preserve system dk
				paraC.shift = para.shift(mc,:);
				% calculate <anam>
				if length(type) == 2
					out(1:paraC.L,mc) = calBath1SiteCorrelatorsProject(mpsC,VmatC,paraC,type{2},[],[]);
				elseif strcmpi(type{3},'diabatic')
					% diabatic projections (H0 eigenstates)
					for kk = 1:paraC.dk(1,1)
						stateProj = zeros(paraC.dk(1,1)); stateProj(kk,kk) = 1;
						bondProj  = [];
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelatorsProject(mpsC,VmatC,paraC,type{2},stateProj,bondProj);
					end
				elseif strcmpi(type{3},'adiabatic')
					% adiabatic states
					for kk = 1:para.D(mc,1)
						stateProj = [];
						bondProj  = zeros(para.D(mc,1)); bondProj(kk,kk) = 1;
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelatorsProject(mpsC,VmatC,paraC,type{2},stateProj,bondProj);
					end
				end
			end
		else
			if length(type) == 2
				out(1:para.L,:) = calBath1SiteCorrelatorsProject(mps,Vmat,para,type{2},[],[]);
			elseif strcmpi(type{3},'diabatic')
				% diabatic projections (H0 eigenstates)
				for kk = 1:para.dk(1,1)
					stateProj = zeros(para.dk(1,1)); stateProj(kk,kk) = 1;
					bondProj  = [];
					out(1:para.L,kk,:) = calBath1SiteCorrelatorsProject(mps,Vmat,para,type{2},stateProj,bondProj);
				end
			elseif strcmpi(type{3},'adiabatic')
				% adiabatic states
				for kk = 1:para.dk(1,1)
					stateProj = [];
					bondProj  = zeros(para.dk(1,1)); bondProj(kk,kk) = 1;
					out(1:para.L,kk,:) = calBath1SiteCorrelatorsProject(mps,Vmat,para,type{2},stateProj,bondProj);
				end
			end
% 			out(:,1,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[1,0;0,0]);	% spin up
% 			out(:,2,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[0,0;0,1]);	% spin down
		end

	case 'bath2correlators'
		% needed for mapping from chain to star
		% returns triangular L x L Matrix
		if para.useStarMPS
			NC = para.nChains;
			out = zeros(para.L, para.L, NC);
			for mc = 1:NC
				% extract MPS for each chain and calculate <anam> separately
				cL = para.chain{mc}.L;
				mpsC  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
				VmatC = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
				mpsC{1} = permute(mpsC{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
				mpsC{1} = reshape(mpsC{1}, [],para.D(mc,1),para.dk(1,1));

				paraC = para;
				paraC.nChains = 1; paraC.L = para.chain{mc}.L;
				paraC.dk = para.dk(mc,:); paraC.dk(1) = para.dk(1,1);	% need to preserve system dk
				paraC.shift = para.shift(mc,:);
				% calculate <anam>
				out(1:paraC.L,1:paraC.L,mc) = calBath2SiteCorrelators_MC(mpsC,VmatC,paraC);
			end
		else
			out = calBath2SiteCorrelators_MC(mps,Vmat,para);
		end

	case 'staroccupation'
		%% does mapping chain -> star
		% uses bath2correlators
		if para.foldedChain == 1, error('Not yet implemented'); end
		% Get the correlator
		if length(type) == 2
			AmAn = real(type{2});
		else
			AmAn = real(getObservable({'bath2correlators'},mps,Vmat,para));		% L x L x nChains
		end
		% preallocate output array
		[x,~,~] = getStarMapping(para,1);
		out      = zeros(2,length(x),para.nChains);
		for mc = 1:para.nChains
			[x,hsquared,pxn] = getStarMapping(para,mc);
			if strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
				result = diag(pxn * (AmAn(:,:,mc) + triu(AmAn(:,:,mc),1).') * pxn.');	% pxn == U map -> no need for hsquared
				out(:,1:length(x),mc) = [x,result(2:end)].';
			else
% 				out(:,1:length(x),mc) = [x, hsquared.*((pxn.^2)*diag(AmAn(:,:,mc)) + 2.* diag(pxn*(AmAn(:,:,mc)-diag(diag(AmAn(:,:,mc))))*pxn.') )]';
				out(:,1:length(x),mc) = [x, hsquared.*diag(pxn*(AmAn(:,:,mc)+triu(AmAn(:,:,mc),1).')*pxn.')].';
			end
		end

	case 'current'
		%% gets current through each bond.
		% works with Single and Multi-Chain Vmat / Vtens
		L = para.L; NC = para.nChains;
		if length(type) == 2
			AmAn = imag(type{2});												% (L x L x NC)
		else
			AmAn = imag(getObservable({'bath2correlators'},mps,Vmat,para));		% tridiagonal
		end
		out = zeros(L-1,NC);
		for mc = 1:NC
			out(:,mc) = (para.chain{mc}.t.*diag(AmAn(1:end-1,2:end,mc)));		% t(n)*a(n)*a(n+1)^+
		end

	case 'displacement'
		% = state projected chain polaron
		% needed for mapping from chain to star, starpolaron
		% returns a L x dk(1,1) x nChains Vector
		%
		% by default:
		% {'displacement'}            : no projection - TODO!
		% {'displacement','diabatic'} : projects onto diabatic states
		% {'displacement','adiabatic'}: projects onto adiabatic states
		
		NC = para.nChains;
		out = zeros(para.L, para.dk(1,1), NC);
			
		if para.useStarMPS
			for mc = 1:NC
				% extract MPS for each chain and calculate <anam> separately
				cL = para.chain{mc}.L;
				mpsC  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
				VmatC = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
				mpsC{1} = permute(mpsC{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
				mpsC{1} = reshape(mpsC{1}, [],para.D(mc,1),para.dk(1,1));

				paraC = para;
				paraC.nChains = 1; paraC.L = para.chain{mc}.L;
				paraC.dk = para.dk(mc,:); paraC.dk(1) = para.dk(1,1);	% need to preserve system dk
				paraC.shift = para.shift(mc,:);
				% calculate <anam>
				if length(type) == 1 || strcmpi(type{2},'diabatic')
					% diabatic projections (H0 eigenstates)
					for kk = 1:paraC.dk(1,1)
						stateProj = zeros(paraC.dk(1,1)); stateProj(kk,kk) = 1;
						bondProj  = [];
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelatorsProject(mpsC,VmatC,paraC,'x',stateProj,bondProj);
					end
				elseif strcmpi(type{2},'adiabatic')
					% adiabatic states
					for kk = 1:para.D(mc,1)
						stateProj = [];
						bondProj  = zeros(para.D(mc,1)); bondProj(kk,kk) = 1;
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelatorsProject(mpsC,VmatC,paraC,'x',stateProj,bondProj);
					end
				end
			end
		else
			if length(type) == 1 || strcmpi(type{2},'diabatic')
				% diabatic projections (H0 eigenstates)
				for kk = 1:para.dk(1,1)
					stateProj = zeros(para.dk(1,1)); stateProj(kk,kk) = 1;
					bondProj  = [];
					out(1:para.L,kk,:) = calBath1SiteCorrelatorsProject(mps,Vmat,para,'x',stateProj,bondProj);
				end
			elseif strcmpi(type{2},'adiabatic')
				% adiabatic states
				for kk = 1:para.dk(1,1)
					stateProj = [];
					bondProj  = zeros(para.dk(1,1)); bondProj(kk,kk) = 1;
					out(1:para.L,kk,:) = calBath1SiteCorrelatorsProject(mps,Vmat,para,'x',stateProj,bondProj);
				end
			end
% 			out(:,1,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[1,0;0,0]);	% spin up
% 			out(:,2,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[0,0;0,1]);	% spin down
		end

	case 'starpolaron'
		%% does mapping chain -> star
		% uses displacement
		% x stands for continuous variable (momentum k)
		% Does state projection!
		% Works with Single and Multi Chain Vmat / Vtens
		%
		%
		if para.foldedChain == 1, error('FoldedChain not yet implemented for polaron'); end
		NC = para.nChains;
		% Get correlators
		if length(type) == 1 || strcmpi(type{2},'diabatic')
			An = real(getObservable({'displacement','diabatic'},mps,Vmat,para));		% L x states x nChains
		elseif strcmpi(type{2},'adiabatic')
			An = real(getObservable({'displacement','adiabatic'},mps,Vmat,para));		% L x states x nChains
		end
% 		AnUp   = real(calBath1SiteCorrelators_MC(mps,Vmat,para,1));			% L x nChains
% 		AnDown = real(calBath1SiteCorrelators_MC(mps,Vmat,para,-1));		% L x nChains

		[x,~,~] = getStarMapping(para,1);
		out = zeros(size(An,2)+1, length(x), NC);
		for mc = 1:NC % for each chain:
			[x,h,pxn] = getStarMapping(para,mc);
			if strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
				result = 2.* pxn * An(:,:,mc);	% pxn == U map -> no need for hsquared
				out(:,1:length(x),mc) = [x,result(2:end,:)].';
			else
% 				out(:,:,mc) = [x, h.*(pxn*An(:,1,mc)),h.*(pxn*An(:,2,mc))]';
				out(:,:,mc) = [x, (h*ones(1,size(An,2))).*(pxn*An(:,:,mc))]';
			end

		end

    case 'energy'
        % Calculates the entire energy of the chain
		% only usable for focus on spin site!
		% needs type{2} = op
		if length(type) == 2
            out = calEnergy(mps,Vmat,para,type{2});
		else
            error('Need op as 2nd argument');
		end

	case 'coherence'
        % only for spin in SBM

	case 'hshi'
        % Calculates the entire energy of the system hamiltonian and first interaction terms
		% needs type{2} = op
		if length(type) == 2
            out = calHsHi(mps,Vmat,para,type{2});
		else
            op.h1term = cell(1,2); % One body terms of the Hamiltonian
			op.h2term = cell(para.M, 2,2);			% default: nearest neighbor interaction

			for s=1:2
				op=genh1h2term_onesite(para,op,s);
			end
			out = calHsHi(mps,Vmat,para,op);
		end
		
	case 'stateproject'
		% Calculates the amplitude of the projection <P|Psi> of the MPS
		% Here |P> is a simple state like |1>|0000000> (System-Bath)
		%
		% by default:
		% {'stateproject'}          : projects onto |2>|0000000>
		% {'stateproject',i}        : projects onto |i>|0000000>
		% {'stateproject',i,j}      : projects onto |i>|jjjjjjj>
		
		% default settings
		systemState =  2;		% take by default state 2
		envState    = 1;		% the state on each boson to project on
		if length(type) >= 2
			systemState = type{2};
		end
		if length(type) == 3
			envState = type{3};
		end
		
		out = calStateProject(mps,Vmat,para,systemState,envState);
		
	case 'sysheff'
		% Calculates the effective Hamiltonian representing the adiabatic total energy surfaces (TES)
		% on which the system evolves.
		% For now: optimised for treeMPS
		%
		% {'sysheff'}				: returns (dk x dk)^2 operator and C
		% {'sysheff',i}				: returns (dk x i)^2, operator with the i strongest adiabats
		% out = {Heff,sysState}		returns the effective hamiltonian together with |psi> mapping from diabatic to the adiabatic basis
		out = cell(1,3);
		if para.useTreeMPS
			[out{1},out{2},out{3}] = calSysHeff(mps,para);
		else
			error('Needs to be implemented');
		end
		
	case 'syspes'
		% Calculates the effective potenital representing the adiabatic potential energy surfaces (PES)
		% on which the system evolves.
		out = cell(1,3);
		if para.useTreeMPS
			[out{1},out{2},out{3}] = calSysPES(mps,para,type{2});
		else
			error('Needs to be implemented');
		end
end

end

function spin = calSpin(mps,Vmat,para)
% Calculate the spin expectation value
% Has to be modified if site 1 changes dimension!!
% Modified:
%   FS 24/05/2014:  excluded MLSpinBoson. Use with try-catch block
%   FS 19/08/2015: - Can be used with Multi-Chain Vmat and Vtens code
%
%   TODO: Needs extension to use information about dimension of spin site.
%
if para.useVtens || para.useStarMPS               % only Quick Fix for Vtens code, Could be used as Standard code!
	McOp = cell(1,1,3);        % L x NC x N
	[McOp{1}, McOp{2}, McOp{3}] = spinop(para.spinbase);
	if para.useVtens
		spinVal = real(expectation_allsites_MC(McOp,mps,Vmat,para));
	else
		spinVal = real(expectation_allsites_StarMPS(McOp,mps,Vmat,para));
	end
	spin.sx = spinVal(1);
	spin.sy = spinVal(2);
	spin.sz = spinVal(3);
	return;
end
if para.useTreeMPS
	McOp = cell(1,1,3);        % L x NC x N
	[McOp{1}, McOp{2}, McOp{3}] = spinop(para.spinbase);
	error('VMPS:getObservable:NotImplementedYet','Please implement this feature!')
end
N=para.L;
assert(N==length(mps) && N==length(Vmat));
assert(~strcmp(para.model,'MLSpinBoson'),'not possible for MLSBM');

ndset=cell(1,N);

for j=1:N
    ndset{1,j}=eye(size(Vmat{j},1));
end
sx=ndset;sy=sx;sz=sy;

%debug:
% sx{1,1}

[sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
sx{para.spinposition(end)}=sigmaX;
sy{para.spinposition(end)}=sigmaY;
sz{para.spinposition(end)}=sigmaZ;
if strcmp(para.model,'2SpinPhononModel')
    sx{para.spinposition}=kron(sigmaZ,eye(2));  %measures excitation of site1
    sz{para.spinposition}=kron(eye(2),sigmaZ);  %measures excitation of site2
    %try to find a good way to measure this!
    sy{para.spinposition}=kron(sigmaY,eye(2));  %measures only sy of site1
end
spin.sx=expectationvalue(sx,mps,Vmat,mps,Vmat);
spin.sy=expectationvalue(sy,mps,Vmat,mps,Vmat);
spin.sz=expectationvalue(sz,mps,Vmat,mps,Vmat);
spin.sx=real(spin.sx);
spin.sy=real(spin.sy);
spin.sz=real(spin.sz);

end

function nx = calBosonOcc(mps,Vmat,para,varargin)
% Calculate the boson occupation on x chain
% The operator on the spin site is set to zero
% left or right chain occupation via varargin{1} = 1 or 2
%
% Modified:
%	FS 23/01/2014:	- Introduced '~' to ignore unused returned values
%					- support for folded Chain models
%   FS 10/03/2014:  - Using correlator_allsites which is more general
%   FS 04/06/2014:  - updated for spinposition array.
%	FS 17/07/2015:	- varargin = {1,2} to calculate nx/nz for folded chain

n_op = cell(1,para.L);

for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 0
            % 1-chain SBM
            [~,~,n_op{j}] = bosonop(para.dk(j),para.shift(j),para.parity);
        %Modification for 2chain model!! Not perfect or right yet!
        elseif para.foldedChain == 1
            % only kron(n,1) chain occupation calculated.
            if para.parity == 'n'
                [~,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
                idm = eye(size(n));
				if varargin{1} == 1
					nr = kron(n,idm);
				elseif varargin{1} == 2
					nr = kron(idm,n);
				end
            else
                [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);   % Why without sqrt??
                [~,~,nr,~,~,~]=paritykron(bp,para.bosonparity);
            end
            n_op{j} = nr;
        else
            disp('Not Implemented: getObservable, calBosonOcc, foldedchain >1');
        end
    else
        n_op{j}=zeros(para.dk(j));      % don't measure spin.
    end
end

%
nx = expectation_allsites(n_op,mps,Vmat);
nx = real(nx);			% imag(nx) = eps -> neglect

end

function n = calBosonOcc_MC(mps,Vtens,para)
%% calculate the boson occupation for Single / Multi-Chain + Vmat / V-tensor-network
% Zero operator for spin site
% calculate for all chains simultaneously
%
% Created 18/08/2015 by FS
assert(para.foldedChain == 0, 'not for folded chains');
L = para.L; NC = para.nChains;

% create Operator for expectationvalue: need nc^2 each (i,:,j) is one operator [] x [] x n x [] x []
McOp = cell(L, NC, NC);
for mc = 1:NC
	for j = 1:para.chain{mc}.L
		if j ~= para.spinposition                     % works with array
			[~,~,McOp{j,mc,mc}] = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
		else
			McOp{j,mc,mc} = zeros(para.dk(1,j));      % don't measure spin.
		end
	end
end

% if para.useStarMPS
% 	n = real(expectation_allsites_StarMPS(McOp,mps,Vtens,para));
% else
	n = real(expectation_allsites_MC(McOp,mps,Vtens,para));  % (NC x L)
% end

end

function n = calBosonOcc_Tree(treeMPS,para)
%% calculate the boson occupation for treeMPS
% Zero operator for spin site
% calculate for all chains simultaneously
%
% Created by FS 28/02/2016
%

% L = para.L; NC = para.nChains;				% total number of site and chains

% Implement recursively for now!
if treeMPS.height == 0
	% this is leaf / chain
	% create Operator for expectationvalue: need nc^2 each (i,:,j) is one operator [] x [] x n x [] x []
	Lc = treeMPS.L;							% length of chain
	Op = cell(1,Lc);
	for j = 1:Lc
		if isempty(treeMPS.spinposition) || j ~= treeMPS.spinposition					% works with array
			[~,~,Op{1,j}] = bosonop(treeMPS.dk(1,j),treeMPS.shift(1,j),para.parity);
		else
			Op{1,j} = zeros(treeMPS.dk(1,j));		% don't measure spin.
		end
	end

	n = real(expectation_allsites(Op,treeMPS.mps,treeMPS.Vmat,treeMPS.BondCenter));		% (NC x L)
	n = reshape(n,[],1);					% L x 1
else
	% this is node -> more recursive calls
	% assume, Focused on MPS!
	% save n into temp cell array, since total number of chains at this node is unknown
	nc    = treeMPS.degree;										% number of subchains at node
	nTemp = cell(1,nc);
	Atemp = contracttensors(treeMPS.BondCenter,2,2,treeMPS.mps{1},nc+2,1);
	for ii = 1:nc
		treeMPS.child(ii).BondCenter = contracttensors(conj(treeMPS.mps{1}),nc+2,[1:ii,ii+2:nc+2],Atemp,nc+2,[1:ii,ii+2:nc+2]);	% contract all except D_(ii+1)
		nTemp{ii}                    = calBosonOcc_Tree(treeMPS.child(ii),para);		% L x nc(child)
	end
	lengths = cellfun(@(x) size(x,1),nTemp);
	nChains = cellfun(@(x) size(x,2),nTemp);
	n  = zeros(max(lengths)+1,sum(nChains));											% L x nc_max
	mc = 0;
	for ii = 1:nc
		n(1+(1:lengths(ii)),mc+(1:nChains(ii))) = nTemp{ii};							% copy all together
		mc = mc + nChains(ii);
	end
end

end

function bosonshift = calBosonShift(mps,Vmat,para)
% Calculate boson shift x, x^2, var(x)
% The operator on the spin site is set to zero
% Modified:
%	FS 22/01/2014:	- changed to using para.foldedChain.
%   FS 04/06/2014:  - updated for spinposition array
%

L = para.L;
x_opx  = cell(1,L);
x2_opx = cell(1,L);
for j = 1:L
    if j ~= para.spinposition
        if para.foldedChain == 1
            % Constructs Supersite Operators
            % Only measures kron(x,1) chain part
            [bp,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
            idm = eye(size(n));
            bpx = kron(bp,idm); bmx = bpx';         %nx = kron(n,idm); unused
            x_opx{j}  = sqrt(2)/2.*(bpx+bmx);
            x2_opx{j} = x_opx{j}*x_opx{j};
        elseif para.foldedChain == 0
            [bp,bm,~] = bosonop(para.dk(1,j),para.shift(1,j),para.parity);
            x_opx{j}  = sqrt(2)/2.*(bp+bm);
            x2_opx{j} = x_opx{j}*x_opx{j};
        end
    else
        x_opx{j}  = zeros(para.dk(j));
        x2_opx{j} = zeros(para.dk(j));
    end
end
bosonshift.x        = expectation_allsites(x_opx,mps,Vmat);
bosonshift.xsquare  = expectation_allsites(x2_opx,mps,Vmat);
bosonshift.xvariant = sqrt(bosonshift.xsquare-bosonshift.x.^2);
bosonshift.xerror   = mean(abs(para.shift-bosonshift.x));
end

function shift = cal1siteShift(mps,Vmat,para)
%% Calculates shift for a single site, para.sitej!
% for use of shifting procedure in optimizesite.m
% mps and Vmat are single-site matrices
% mps has to be focused

NC = para.nChains; j = para.sitej;
xOp_OBB  = cell(1,NC);			% one x-OBB operator per chain
shift    = zeros(1,NC);			% one value per chain
% Generate operators and contract with OBB
if para.foldedChain == 1
	switch para.model
		case 'SpinDoubleBoson' % calculates shift only for x chain!
			bp = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
			if para.parity == 'n'
				idm = eye(size(bp));
				bpx = kron(bp,idm); bmx = bpx';
			else
				[bpx,bmx,~,~,~,~]=paritykron(bp,para.bosonparity);
			end
			x = sqrt(2)/2*(bpx+bmx);
		case '2SpinPhononModel'
			bp = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
			if para.parity=='n'
				idm=eye(size(bp));
				bpr=kron(bp,idm);   bmr=bpr';	% right chain
% 				bpl=kron(idm,bp);   bml=bpl';	% left chain
			else
				[bpr,bmr,~,~,~,~]=paritykron(bp,para.bosonparity);
			end
			x = sqrt(2)/2*(bpr+bmr);                                % why only evaluate it for right part?
	end
	if para.useVmat
		xOp_OBB{1} = Vmat' * x * Vmat;
	else
		xOp_OBB{1} = x;
	end
else
	for mc = 1:NC
		x     = cell(1,NC);
		[bp,bm,~]   = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
		x{mc}       = sqrt(2)/2*(bp+bm);
		xOp_OBB{mc} = contractMultiChainOBB(Vmat, x, para);						% _(n~',n~)
	end
end
temp = contracttensors(conj(mps),3,[1,2],mps,3,[1,2]);							% _(n~',n~)
for mc = 1:NC
	shift(mc) = real(sum(sum(temp.*xOp_OBB{mc})));								% fast tr(A * B.')
end

end

function reducedDensity = calRDM(mps,Vmat,para,k)
% calculates the reduced density matrix of any single site, mps{k} for size(mps{k})=a_{k-1} x a_k x n_k
%
%   created 24/05/2014 by Florian Schroeder @ University of Cambridge
%
%
% copied from prepare.m:
% does l -> r sweep to create state in local picture of k
para.sweepto = 'r';
if length(k) == 1		% single-site RDM
	for i = 1:k-1
		if para.useVmat==1
			[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
			mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
		end
		[mps{i}, U] = prepare_onesite(mps{i}, para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
		mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
		para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
		para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
	end

	% now in form: Al{1}...Al{k-1} M{k} Ar{k+1}...Ar{L}
	%   with Al = left-normalized, Ar: right-normalized.

	n = ndims(mps{k});
	reducedDensity = contracttensors(mps{k},n,1:n-1,conj(mps{k}),n,1:n-1);		% contract rD_nm = Mk_abn Mk*_abm

	% reducedDensity = Vmat{k} * (reducedDensity * Vmat{k}');
	reducedDensity = contracttensors(reducedDensity,2,2,conj(Vmat{k}),2,2);     % contract rD_nj = rD_nm Vmat*_jm
	reducedDensity = contracttensors(Vmat{k},2,2,reducedDensity,2,1);			% contract rD_ij = Vmat_in rd_nj
elseif length(k) == 2 && k(2) == k(1) + 1
	% 2-site RDM of nearest neighbours!
	% move to k(1)
	for i = 1:k(1)-1
		if para.useVmat==1
			[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
			mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
		end
		[mps{i}, U] = prepare_onesite(mps{i}, para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
		mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
		para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
		para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
	end
	% use ncon, since easy and not performance-critical!
	reducedDensity = ncon({mps{k(1)}, mps{k(2)}, conj(mps{k(1)}), conj(mps{k(2)}), Vmat{k(1)}, Vmat{k(2)}, conj(Vmat{k(1)}), conj(Vmat{k(2)})},...
						  {[1,2,3],   [2,4,5],   [1,6,7],         [6,4,8],         [-2,3],     [-1,5],     [-4,7],           [-3,8]});
end
end

function participation = calParticipation(rdm)
% rdm: a reduced density matrix of a single site
%   e.g. rdm = calRDM(mps,Vmat,para,k)
%
% takes a rdm and calculates the participation ratio

participation = 1/sum(diag(rdm).^2);

end

function tunnelE = calTunnelingEnergy(mps,Vmat,para,op)
%% calculates the tunneling energy in MLSBM of the PPC system
%   <Psi|H0-diag(H0)|Psi>
%
%   could be applied to any system.
%   Interacting System Hamiltonian: HI


HI = cell(1);
if ~exist('op','var') || ~isfield(op,'h1term')
    disp('guessing H0 from para for MLSBM');
    op.h1term{1,1} = Hamiltonian_PPC(para);
end
HI{1} = op.h1term{1,1}-diag(diag(op.h1term{1,1}));
tunnelE = expectationvalue(HI,mps,Vmat,mps,Vmat);

end

function An = calBath1SiteCorrelatorsProject(mps,Vmat,para,opSelect,stateProj,bondProj)
% calculates the (projected) single-site expectation value
% output is matrix containing all values for all chains (L x nChains)
% custom made solution for biggest speedup, exact solution!
%	opSelect: ['bp'|'x'|'bp^2'|'x^2'|'n']
% 	stateProj: equals state-projection operator, normalised!
%				[]: no selection, take all
%	bondProj:  equals bond-projection operator, normalised!
%				[]: no selection, take entire bond
%
% supports Multi-Chain with Vmat and Vtens, Works with single chain!
%  does NOT support StarMPS or TreeMPS input! This needs pre-processing

assert(para.foldedChain == 0, 'Please use single-chain code for folded Chains');

An     = zeros(para.L, para.nChains);									% initialize results array
op_OBB = cell( para.L, para.nChains);
%% generate all operators & contract with OBB Vmat / Vtens:
for j = 1:para.L
	para.sitej = j;
	for mc = 1:para.nChains
		op = cell(1, para.nChains);						            % containing single-site, single-chain operator
		if j ~= para.spinposition
			switch opSelect
				case 'bp'
					op{mc}        = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
				case 'x'
					bp            = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
					op{mc}        = (bp+bp')/2;
				case 'bp^2'
					bp            = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
					op{mc}        = bp^2;
				case 'x^2'
					bp            = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
					op{mc}        = (bp+bp')^2/4;
				case 'n'
					[~,~,op{mc}]  = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
			end
			op_OBB{j, mc} = contractMultiChainOBB(Vmat{j}, op, para);
		else	% Spin & Multi-Level System
			if ~isempty(stateProj)
				op_OBB{j,mc} = Vmat{j}' * stateProj * Vmat{j};					% Spin-site Vmat, project onto system state!
			else
				op_OBB{j,mc} = Vmat{j}' * Vmat{j};
			end
		end
	end
end

if isempty(stateProj) && ~isempty(bondProj)
	% do SVD from site 1 to 2 in order to get projections onto the dominant system states!
	%	otherwise have projection onto dominant Bath states (from previous sweep to left)
	A = mps{1};
	[D1, D2, d] = size(A);				% d is site dimension
	A = permute(A, [3, 1, 2]);
    A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
    [B, S, U] = svd2(A);				% Could also use QR decomposition if nargin !=5
	
	DB = size(S, 1);					% new a2 dimension of B should == d
	B = reshape(B, [d, D1, DB]);
	mps{1} = permute(B, [2, 3, 1]);		% not used anymore!
	% do not contract S into mps{2}, to avoid the displacement to be contaminated with the state amplitude!
	mps{2} = contracttensors(U,2,2, mps{2},3,1);
end

%% compute all partial contractions
% such that (1+sz/4) gets updated to the right in row 1;
% Trace(contraction for result) = Am(j)
opContract = cell(para.L, 1);										% same for all chains! since only 1 common MPS backbone
for j = 1:para.L
	for mc = 1:para.nChains
		if j ~= 1
			An(j,mc)    = trace(updateCleft(opContract{j-1,1},mps{j},[],op_OBB{j,mc},mps{j},[]));	% this is the result!
		end
	end

	if j ~=1
		opContract{j,1} = updateCleft(opContract{j-1,1},mps{j},[],[],mps{j},[]);				% Vmat is normalised -> leave out!
	else % j == 1
		if isempty(bondProj)
			opContract{1,1} = updateCleft([],mps{j},[],op_OBB{j,1},mps{j},[]);					% for state-projection; = (r',r)
		else
			opContract{1,1} = bondProj;
		end
	end
end

end

function An = calBath1SiteCorrelators_Tree(treeMPS,para,opSelect,projOp)
% calculates the single-site expectation value
% Zero operator for spin sites / nodes
% output is matrix containing all values for all chains (L x nChains)
% custom made solution for biggest speedup, exact solution!
%	opSelect: ['bp'|'x'|'bp^2'|'x^2'|'n']
%
% Projection operator only applied to dk of root node!
% supports only TreeMPS input!
%
% Created by FS 11/03/2016
%

% Implement recursively for now!
if treeMPS.height == 0
	% this is leaf / chain
	% create Operator for expectationvalue: need nc^2 each (i,:,j) is one operator [] x [] x n x [] x []
	Lc = treeMPS.L;							% length of chain
	Op = cell(1,Lc);
	for j = 1:Lc
		if isempty(treeMPS.spinposition) || j ~= treeMPS.spinposition					% works with array
			bp              = bosonop(treeMPS.dk(1,j),treeMPS.shift(1,j),para.parity);
			switch opSelect
				case 'bp'
					Op{1,j} = bp;
				case 'x'
					Op{1,j} = (bp+bp')/2;
				case 'bp^2'
					Op{1,j} = bp^2;
				case 'x^2'
					Op{1,j} = (bp+bp')^2/4;
				case 'n'
					Op{1,j} = bp*bp';
			end
		else
% 			if nargin > 3
% 				Op{1,j} = projOp;						% use Projection operator if specified
% 			else
				Op{1,j} = zeros(treeMPS.dk(1,j));		% no projection
% 			end
		end
	end

	An = expectation_allsites(Op,treeMPS.mps,treeMPS.Vmat,treeMPS.BondCenter);			% (NC x L)
	An = reshape(An,[],1);					% L x 1
else
	%% this is node -> more recursive calls
	% assume, Focused on MPS!
	% save n into temp cell array, since total number of chains at this node is unknown
	nc    = treeMPS.degree;										% number of subchains at node
	AnTemp = cell(1,nc);
	Atemp = contracttensors(treeMPS.BondCenter,2,2,treeMPS.mps{1},nc+2,1);
	if nargin > 3 && ~isempty(projOp)
		Atemp = contracttensors(Atemp,nc+2,nc+2,projOp.',2,1);							% Apply projection to first site
	end
		
	for ii = 1:nc
		treeMPS.child(ii).BondCenter = contracttensors(conj(treeMPS.mps{1}),nc+2,[1:ii,ii+2:nc+2],Atemp,nc+2,[1:ii,ii+2:nc+2]);	% contract all except D_(ii+1)
		AnTemp{ii}                   = calBath1SiteCorrelators_Tree(treeMPS.child(ii),para,opSelect);		% L x nc(child)
	end
	lengths = cellfun(@(x) size(x,1),AnTemp);
	nChains = cellfun(@(x) size(x,2),AnTemp);
	An  = zeros(max(lengths)+1,sum(nChains));											% L x nc_max
	mc = 0;
	for ii = 1:nc
		An(1+(1:lengths(ii)),mc+(1:nChains(ii))) = AnTemp{ii};							% copy all together
		mc = mc + nChains(ii);
	end
end

end

function An = calBath1SiteCorrelators_MC(mps,Vmat,para,stateProj,bondProj)
% calculates Re<(1+-sz)/4 * a_n^+>.
% output is matrix containing all values for all chains (L x nChains)
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!
% 	stateProj: equals state-projection operator, normalised!
%	bondProj:  selecting nth-dominant state of bond
%				[]: no selection, take entire bond
%				n : select nth-state
% supports Multi-Chain with Vmat and Vtens, Works with single chain!
%  does NOT support StarMPS input!

assert(para.foldedChain == 0, 'Please use single-chain code for folded Chains');

An     = zeros(para.L, para.nChains);									% initialize results array
bp_OBB = cell( para.L, para.nChains);
%% generate all operators & contract with OBB Vmat / Vtens:
for j = 1:para.L
	para.sitej = j;
	for mc = 1:para.nChains
		bp = cell(1, para.nChains);						            % containing single-site, single-chain operator
		if j ~= para.spinposition
			bp{mc}        = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
			bp_OBB{j, mc} = contractMultiChainOBB(Vmat{j}, bp, para);
		else	% Spin & Multi-Level System
			if ~isempty(stateProj)
				bp_OBB{j,mc} = Vmat{j}' * (stateProj./2) * Vmat{j};				% Spin-site Vmat, 1/2 since f = (a + a^+)/2
			else
				bp_OBB{j,mc} = (Vmat{j}' * Vmat{j})/2;							% = eye(d_opt)/2; remains unused, since using bondProj below!
			end
		end
	end
end

if isempty(stateProj) && ~isempty(bondProj)
	% do SVD from site 1 to 2 in order to get projections onto the dominant system states!
	%	otherwise have projection onto dominant Bath states (from previous sweep to left)
	A = mps{1};
	[D1, D2, d] = size(A);				% d is site dimension
	A = permute(A, [3, 1, 2]);
    A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
    [B, S, U] = svd2(A);				% Could also use QR decomposition if nargin !=5
	
	DB = size(S, 1);					% new a2 dimension of B should == d
	B = reshape(B, [d, D1, DB]);
	mps{1} = permute(B, [2, 3, 1]);		% not used anymore!
	% do not contract S into mps{2}, to avoid the displacement to be contaminated with the state amplitude!
	mps{2} = contracttensors(U,2,2, mps{2},3,1);
end

%% compute all partial contractions
% such that (1+sz/4) gets updated to the right in row 1;
% Trace(contraction for result) = Am(j)
bpContract = cell(para.L, 1);										% same for all chains! since only 1 common MPS backbone
for j = 1:para.L
	for mc = 1:para.nChains
		if j ~= 1
			An(j,mc)    = trace(updateCleft(bpContract{j-1,1},mps{j},[],bp_OBB{j,mc},mps{j},[]));	% this is the result!
		end
	end

	if j ~=1
		bpContract{j,1} = updateCleft(bpContract{j-1,1},mps{j},[],[],mps{j},[]);				% Vmat is normalised -> leave out!
	else % j == 1
		if isempty(bondProj)
			bpContract{1,1} = updateCleft([],mps{j},[],bp_OBB{j,1},mps{j},[]);					% for state-projection; = (r',r)
		else
			bpContract{1,1} = bondProj/2;														% 1/2 since f = (a + a^+)/2
		end
	end
end

end

function AmAn = calBath2SiteCorrelators_MC(mps,Vmat,para)
% calculates <a_m^+ a_n> for any combination, where n>m.
% output is matrix containing all values.
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!
%
% Only reduction for Multi-Chain: effective left density matrix
assert(para.foldedChain == 0, 'FoldedChain not yet supported');
AmAn = zeros(para.L,para.L,para.nChains);	% initialize results array
bpbm = cell(3,para.L,para.nChains);			% containing all necessary operators, contracted with OBB
leftDM  = cell(para.L,1);					% left-effective Density Matrix for each site; right-effective should always be 1
rightDM = cell(para.L,1);
%% generate all operators and contract OBB
%  also produce left-effective DM
for j = 1:para.L
	para.sitej = j;
	for mc = 1:para.nChains
		bp = cell(1,para.nChains); bm = cell(1,para.nChains); n = cell(1,para.nChains);
		if j ~= para.spinposition
			[bp{mc}, bm{mc}, n{mc}] = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
			bpbm{1,j,mc} = contractMultiChainOBB(Vmat{j}, bp, para);
			bpbm{2,j,mc} = contractMultiChainOBB(Vmat{j}, bm, para);
			bpbm{3,j,mc} = contractMultiChainOBB(Vmat{j},  n, para);
		else
			bpbm{1,j,mc} = zeros(para.dk(1,j));			% zero operators for spin sites
			bpbm{2,j,mc} = zeros(para.dk(1,j));			% dk = d_opt
			bpbm{3,j,mc} = zeros(para.dk(1,j));
		end
	end
	if j ~= 1
		leftDM{j} = updateCleft(leftDM{j-1}, mps{j}, [], [], mps{j}, []);
	else
		leftDM{j} = updateCleft(    []     , mps{j}, [], [], mps{j}, []);
	end
end

%% compute all partial contractions
% such that <m n> is meeting in the middle -> reduce overhead!
for mc = 1:para.nChains
	bpContract = cell(para.L);
	for j = 1:para.L-1
		if j == 1 || mod(j,10) == 0
			fprintf('%g-',j);
		end
		for i = 1:j+2
			if (i < j) && (j <= floor( (para.L+i)/2))
				% take a^+ to next effective basis
				bpContract{i,j} = updateCleft(bpContract{i,j-1},mps{j},[],[],mps{j},[]);
			elseif i == j % contruct <leftDM * a+
				if j ~=1
					bpContract{i,i} = updateCleft(leftDM{j-1},mps{j},[],bpbm{1,j,mc},mps{j},[]);
				else
					bpContract{i,i} = updateCleft(    []     ,mps{j},[],bpbm{1,j,mc},mps{j},[]);
				end
			elseif i == j+1 % construct left-effective DM
				% do nothing, leave empty
			elseif i == j+2 % a^+ a of same site
				if j ~= 1
					bpContract{i,j} = trace(updateCleft(leftDM{j-1},mps{j},[],bpbm{3,j,mc},mps{j},[]));
				else
					bpContract{i,j} = trace(updateCleft(     []    ,mps{j},[],bpbm{3,j,mc},mps{j},[]));
				end
			end
		end
	end
	L = para.L;
	% The last occupation:
	bpContract{L+2,L} = trace(updateCleft(leftDM{L-1},mps{L},[],bpbm{3,L,mc},mps{L},[]));
	fprintf('\n');

	bmContract = cell(para.L);
	for j = para.L:-1:2
		if j == 1 || mod(j,10) == 0
			fprintf('%g-',j);
		end
		for i = para.L:-1:j
			if (i > j) && (j > ceil( i/2 ) )
				bmContract{i,j} = updateCright(bmContract{i,j+1},mps{j},[],[],mps{j},[]);
			elseif i == j
				if j ~= para.L
					bmContract{i,i} = updateCright(rightDM{j+1},mps{j},[],bpbm{2,j,mc},mps{j},[]);
				else
					bmContract{i,i} = updateCright(     []     ,mps{j},[],bpbm{2,j,mc},mps{j},[]);
				end
			end
		end
		if mc == 1 % construct right effective DM only once
			if j ~= para.L
				rightDM{j} = updateCright(rightDM{j+1},mps{j},[],[],mps{j},[]);
			else
				rightDM{j} = updateCright(     []     ,mps{j},[],[],mps{j},[]);
			end
		end
	end

	for j = 1:para.L
		for i = 1:j
			midPoint = floor((j+i)/2);
			if  i ~= j
				AmAn(i,j,mc) = sum(sum(bpContract{i,midPoint} .* bmContract{j,midPoint+1}));		% faster tr(A * B.')
			else
				AmAn(i,j,mc) = bpContract{i+2,j};
			end
		end
	end
end
end

function E = calEnergy(mps,Vmat,para,op)
	%%
	sitej=1;
	if para.useStarMPS
		E = 0;			% nut supported for now!
		return;
	end
% 	[~,BondDimRight] = size(op.Hright);
% 	[~,BondDimLeft]  = size(op.Hleft);
% 	op = h1h2toOBB(Vmat{sitej},para,op);		% DEPRECATED
	op  = H_Eff([],Vmat{sitej},'A',op,para);
% 	[~,OBBDim]	 = size(op.h1j);
	[BondDimLeft, BondDimRight, OBBDim] = size(mps{sitej});
	M = size(op.h2j,1);
	A = reshape(mps{sitej},[numel(mps{sitej}),1]);
	E = A'*HmultA(A, op, BondDimLeft, BondDimRight, OBBDim, M,para.parity,[]);
	E = real(E);		% imag(E) = eps -> negligible!
end

function n = expectation_allsites_MC(McOp,mps,Vmat,para)
%% Multi Chain Expectation values for Vmat and Vtens
% copied from expectation_allsites since faster to write
%
% McOp : L x NC x N
% NC = para.nChains;
if para.useStarMPS		% shortcut for StarMPS
	n = expectation_allsites_StarMPS(McOp,mps,Vmat,para);
	return;
end
[L,~,N] = size(McOp);               % N = number of "observable chains"
% assert(L == length(mps));			% would also work for N < dim(mps)

% Contract McOp with OBB of Vmat or Vtens
McOp_OBB = cell(N,L);
for j = 1:L
	para.sitej = j;
	for n = 1:N
		if j ~= para.spinposition
			McOp_OBB{n,j} = contractMultiChainOBB(Vmat{j}, McOp(j,:,n), para);
		else
			Ind = ~cellfun('isempty', McOp(j,:,n));
			McOp_OBB{n,j} = Vmat{j}' * McOp{j,Ind,n} * Vmat{j};
		end
	end
end

% call with empty Vmat since OBB is already contracted
n = expectation_allsites(McOp_OBB, mps, cell(1,para.L));	% N x L

end

function n = expectation_allsites_StarMPS(McOp,mps,Vmat,para)
%% Multi Chain Expectation values for StarMPS code
% copied from expectation_allsites since faster to write
%
% McOp : L x NC x N
NC = para.nChains;
[L,M,N] = size(McOp);               % N = number of "observable chains";
% assert(L == length(mps));			% would also work for N < dim(mps)

n = zeros(N,L);
% Here: iterate over chains and determine whether McOp is empty
for mc = 1:M
	tempOp = reshape(McOp(:,mc,:),L,N);				% L x N
	empty  = sum(cellfun('isempty',tempOp),1);		% 0 where Op, L where empty
	if all(empty == L), continue; end
	% else: extract mpsChain, mpsVmat including system site
	cL = para.chain{mc}.L;
	mpsChain  = [mps(1), cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
	VmatChain = [Vmat(1), cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
	mpsChain{1} = permute(mpsChain{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
	mpsChain{1} = reshape(mpsChain{1}, [],para.D(mc,1),para.dk(1,1));
	for ii = find(empty ~= L)
		n(ii,1:min(L,cL)) = expectation_allsites(tempOp(1:min(L,cL),ii).', mpsChain, VmatChain);
	end
end

end

function [x,hsquared,pxn] = getStarMapping(para,mc)
	%% finds mapping chain -> star
		n = 0:(para.L-2);		% length(n) = L-1
		alphaN   = zeros(length(n),1);
		betaN    = zeros(length(n),1);			% this is sqrt(betaN); betaN(1) = beta_0 is meaningless!

		% find bath parameters:
		if strcmp(para.chain{mc}.mapping,'OrthogonalPolynomials')
			step = max(1/para.L/10,10^-3);
			x = (0:step:1)';		% resolution
			hsquared = ones(length(x),1);
			pxn = zeros(length(x),para.L);			% these are orthonormalized polynomials. examples:
					% pxn(:,1) = 0; pxn(:,2) = px(0) = 1/t(1); pxn(:,3) =
					% x-alpha(0)/norm();
			s = para.chain{mc}.s;
			if strcmp(para.chain{mc}.spectralDensity,'Leggett_Hard')
				hsquared = 2.*para.chain{mc}.alpha.*(x.^s);
				alphaN = (1+( (s^2)./((s+2.*n).*(2+s+2.*n)) ))./2;
				betaN  = ((1+n).*(1+s+n))./(s+2+2.*n)./(3+s+2.*n).*sqrt((3+s+2.*n)./(1+s+2.*n));
			elseif strcmp(para.chain{mc}.spectralDensity,'Leggett_Soft')
				x = x.*3;						% need some longer range since spectralDensity is not limited in bandwidth!
				hsquared = 2.*para.chain{mc}.alpha.*(x.^s).*exp(-x);
				alphaN = 2.*n+1+s;
				betaN  = sqrt((n+1).*(n+1+s));
			end
		elseif strcmp(para.chain{mc}.mapping, 'Stieltjes')
			x        = reshape(para.chain{mc}.xi, [],1);
			hsquared = reshape(para.chain{mc}.Gamma.^2,[],1);
			alphaN   = reshape(para.chain{mc}.epsilon, 1,[]);
			betaN    = reshape(para.chain{mc}.t(2:end),1,[]);
			pxn      = zeros(length(x),para.chain{mc}.L);
		elseif strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
			x        = reshape(para.chain{mc}.xi, [],1);			% frequency Axis
			hsquared = reshape(para.chain{mc}.Gamma.^2,[],1);		% not needed
			pxn      = para.chain{mc}.U;							% not poly, but map U
			return;
		end
		% find OrthPol:
		pxn(:,2) = ones(length(x),1).*1/para.chain{mc}.t(1);		% p(0) = constant normalized pol.
		for i = 1:para.L-2
			% p(i+1) = ((x - a(i))p(i)-b(i)p(i-1))/b(i+1)
			pxn(:,i+2) = ( (x - alphaN(i)).*pxn(:,i+1) - betaN(i).*pxn(:,i) )./betaN(i+1);
		end


end

function E = calHsHi(mps,Vmat,para,op)
%% function E = calHsHi(mps,Vmat,para,op)
% calculates the expectation value of <H_s + H_i(1,2)>
% of only System and first bond.
if ~exist('op','var') || ~isfield(op,'h1term')
	error('op is needed to use calHsHi');
end

% create n_op for expectation_general_StarMPS
% n_op = M+1 x L x NC
%	where H1term = n_op(1,1,1)
M  = para.M;
NC = para.nChains;
L  = para.L;

n_op = cell(M+1,L,NC);

n_op{1,1,1} = op.h1term{1,1};
for m = 1:M
	idx = squeeze(~cellfun('isempty',op.h2term(m,2,2,:)));
	n_op(m+1,1,idx) = op.h2term(m,1,1,1);
	n_op(m+1,2,idx) = op.h2term(m,2,2,idx);
end

E_i = expectation_general_StarMPS(n_op,mps,Vmat,para);
E = [real(E_i(1));2*real(E_i(2:2:end))]';				% find a better solution for this! But not now..

% E = sum(E_i);

end

function n = expectation_general(n_op, mps, Vmat, para)
%% function n = expectation_general(n_op, mps, Vmat, para)
% computes the expectation value of an arbitrary (general) multi-site operator defined as:
% O = 1 x L cell
%		where site-operators are stored in respective cells. Other cells: []
%	e.g:  < O_i O_j O_k >
%		O = {[],...,[],O_i,[],...,[],O_j,[],...,[],O_k,[],...}
%
% if multiple operators are desired on the same chain
% n_Op = M x L cell
%       where n_Op{1,:} = O{1}
%             n_Op{2,:} = O{2}  ...
%
% only for single-chain MPS
% output contains no site-dependence, such as expectation_allsites(). This has to be dealt with inside the calling function!
%
% written by FS 09/11/2015

[M,~] = size(n_op);						% M x L
empty = cellfun('isempty',n_op);

tempC = cell(M,1);						% will contain temporary right-effective operators.
lastnz = find(sum(~empty,1),1,'last');	% last non-zero

for ii = lastnz:-1:1
	% sweep from back
	for m = 1:M
		if ii == 1 && size(mps{1},1) > size(mps{1},2)
			n(m) = trace(updateCleft([],mps{ii},Vmat{ii},n_op{m,ii},mps{ii},Vmat{ii})*tempC{m}.');		% contract left first, since memory can blow up!
		else
			tempC{m} = updateCright(tempC{m},mps{ii},Vmat{ii},n_op{m,ii},mps{ii},Vmat{ii});
			if ii == 1
				n(m) = trace(tempC{m});
			end
		end
	end
end

end

function n = expectation_general_StarMPS(n_op, mps, Vmat, para)
%% function n = expectation_general(n_op, mps, Vmat, para)
% computes the expectation value of an arbitrary (general) multi-site operator defined as:
% O = 1 x L cell
%		where site-operators are stored in respective cells. Other cells: []
%	e.g:  < O_i O_j O_k >
%		O = {[],...,[],O_i,[],...,[],O_j,[],...,[],O_k,[],...}
%
% if multiple operators are desired on any chain mc
% n_Op = M x L x NC cell
%       where n_Op{1,:,mc} = O{1}
%             n_Op{2,:,mc} = O{2}  ...
%
%
% only for single-chain MPS
% output contains no site-dependence, such as expectation_allsites(). This has to be dealt with inside the calling function!
%
% written by FS 09/11/2015

[M,~,NC] = size(n_op);					% M x Lmax x NC
empty = cellfun('isempty',n_op);

n = zeros(M,1);

for mc = 1:NC
	% find which m can be computed
	allnz = find(sum(~empty(:,:,mc),2));	% all non-zero for mc

	% reshape MPS into single-chain
	cL = para.chain{mc}.L;
	mpsChain  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
	VmatChain = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
	mpsChain{1} = permute(mpsChain{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
	mpsChain{1} = reshape(mpsChain{1}, [],para.D(mc,1),para.dk(1,1));

	% calculate the observables for that chain:
	n(allnz) = expectation_general(n_op(allnz,:,mc),mpsChain,VmatChain,para);
end


end

function A = calStateProject(mps,Vmat,para,systemState,envState)
% Calculates the amplitude of the projection <P|Psi> of the MPS
% Here |P> is a simple state like |1>|0000000> (System-Bath)

Pmps = zeros(para.dk(1,1),1);		% projector for system state with all D=1
Pmps(systemState,1) = 1;
if para.useStarMPS
	A = contracttensors(mps{1},para.nChains+2,para.nChains+2,Pmps,2,1);	% A_(1,D1,D2,...,Dnc,1) = A_(1,D1,D2,...,Dnc,d) * P_(d,1)
	A = squeeze(A);					% A_(D1,D2,...,Dnc)
elseif para.useVtens
	% needs implementation
	error('VMPS:getObservable:stateproject:NotImplemented','Not yet implemented')
else
	A = contracttensors(mps{1},3,3,Pmps,2,1);							% A_(1,D,1) = A_(1,D1,d) * P_(d,1)
end

% Do following reverse for less code
for mc = para.nChains:-1:1
	Cright = [];
	% define Chain Projector as MPS/Vmat
	L = para.chain{mc}.L-1;					% allows different chain lengths, exclude system!
	dk = para.dk(mc,2:end);
	Pmps  = 1;					% 1x1x1 tensor
	for kk = L:-1:1
		PVmat = zeros(dk(kk),1);	% dk x dOBB
		PVmat(envState,1) = 1;
		if para.useStarMPS
			Cright = updateCright(Cright,Pmps,PVmat,[],mps{kk+1}{mc},Vmat{kk+1}{mc});		% 1 x D(kk)
		elseif para.useVtens
			% needs implementation
			error('VMPS:getObservable:stateproject:NotImplemented','Not yet implemented');
		else
			Cright = updateCright(Cright,Pmps,PVmat,[],mps{kk+1},Vmat{kk+1});
		end
	end
	if mc > 1
		A = squeeze(contracttensors(A,mc,mc,Cright,2,2));
	else
		A = reshape(A,1,[])*reshape(Cright,[],1);			% upper expression fails for mc = 1, since ndims(vec) = 2 in matlab! -> special treatment
	end
end

end

function [Heff, C, E] = calSysHeff(treeMPS,para,n)
	%% [Heff, C, E] = calSysHeff(treeMPS,para,n)
	%	calculates the effective Hamiltonian and corresponding system state to obtain the total energy surfaces
	%
	%	Output:
	%	 Heff: D' x dk' x D x dk is the effective Hamiltonian; D are adiabatic (bond) states, dk are the diabats
	%		C: dk(diab) x D(adiab) maps diabatic to adiabatic states;
	%		E: Energy when contracting Heff with C
	if nargin == 2
		n = [];
	end
	if strcmp(para.model,'SpinBoson')
		%% shortcut for simpler code
		OBBDim = treeMPS.d_opt;
		BondDimRight = treeMPS.D(2);	% works only for single chain models now! e.g. SBM
		
		% Heff_(r',n',r,n) = kron(OBB,BondRight); this follows convention in tdvp_1site_evolveHn
		temp = 0;
		temp = temp + kron(treeMPS.op.h1jOBB,eye(BondDimRight));
		temp = temp + kron(eye(OBBDim),treeMPS.child(1).op.Hlrstorage{1});
		for mm = 1:para.M
			temp = temp + kron(treeMPS.op.h2jOBB{mm,1,1},treeMPS.child(1).op.Opstorage{mm,2,1,1});
		end
		C = squeeze(treeMPS.mps{1}).';					% dk x D
	else
		%% try to write general approach for treeMPS; Copied from expvCustom/TreeMultA
		% 1. split off center to obtain isometry from chains to adiabatic system state
		d = size(treeMPS.mps{1});							% 1 x D1...Dn x dk
		OBBDim = d(end);
		BondDimRight = d(end);
		
		para.sweepto = 'l';
		
		Atens = reshape(treeMPS.mps{1},1,[],d(end));	%  1 x prod(D) x dk
		Atens = permute(Atens,[3,2,1]);					% dk x prod(D) x  1
		[A,C] = prepare_onesite(Atens,para);			% C: dk x dk maps diabatic to adiabatic states; A: dk x prod(D) x 1
		A = permute(A,[3,2,1]);							%  1 x prod(D) x dk
		A = reshape(A,d);								% 1 x D1...Dn x dk; now only isometry to map into systems diabatic states.
		
		NC = treeMPS.degree;
		
		temp = 0;				% takes kron(diabatic,adiabatic=BondDimRight)
		
		% 1. on-site H1
		temp = temp + kron(treeMPS.op.h1jOBB,eye(BondDimRight));
		
		NCoupOffset = 0;									% needed for child nodes being entanglement renormalisation tensors
		% 2. iterate through children
		for mc = 1:NC
			ordTo = [mc+1:NC+2,1:mc];
			A1 = permute(A,ordTo);
			A1 = reshape(A1,d(mc+1),[]);				% D_mc x ( ... dk x 1 x ... x D_{mc-1}); Useful for 3) as well
			
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = treeMPS.child(mc).op.Hlrstorage{1} * A1;
			OpTemp = reshape(OpTemp, d(ordTo));
			OpTemp = permute(OpTemp, [NC+3-mc:NC+2,1:NC+2-mc]);	% 1 x D1...Dn x dk
			OpTemp = contracttensors(conj(A),NC+2,1:NC+1, OpTemp, NC+2, 1:NC+1);
			temp = temp + kron(eye(OBBDim),OpTemp);
			
			% 3. all interacting parts
			NCoup = size(treeMPS.child(mc).op.Opstorage,4);
			for nn = 1:NCoup								% iterate through chains of child
				for mm = 1:para.M
					OpTemp = treeMPS.child(mc).op.Opstorage{mm,2,1,nn} * A1;
					OpTemp = reshape(OpTemp, d(ordTo));
					OpTemp = permute(OpTemp, [NC+3-mc:NC+2,1:NC+2-mc]);	% 1 x D1...Dn x dk
					OpTemp = contracttensors(conj(A),NC+2,1:NC+1, OpTemp, NC+2, 1:NC+1);
					
					temp = temp + kron(treeMPS.op.h2jOBB{mm,1,nn+NCoupOffset}, OpTemp);
				end
			end
			NCoupOffset = NCoupOffset + NCoup;
		end
	end
	Heff = reshape(temp,[BondDimRight,OBBDim,BondDimRight,OBBDim]);
	
	% also get energy
	E = contracttensors(Heff,4,[3,4], C,2,[2,1]);
	E = real(contracttensors(E,2,[1,2], conj(C),2,[2,1]));
	
	if ~isempty(n)
		Heff = Heff(1:n,:,1:n,:);
		C = C(:,1:n);
	end
end

function [Heff, C, E] = calSysPES(treeMPS,para,Dmax)
	% New approach, trying to remove kinetic energy contributions from the adiabatic potentials
	% To obtain the potential energy surfaces
	
	%% 1. Traverse tree and remove kinetic terms from h1term
	treeMPSPot = setPotOp(treeMPS);
	
	% update everything based on modified h1term,h2term
	treeMPSPot = initstorage(treeMPSPot,[],[],para);		% new call pattern for treeMPS
	[Heff, C, E] = calSysHeff(treeMPSPot,para,Dmax);
	
	
	function treeMPS = setPotOp(treeMPS)
		% set the potential operator by removing kinetic terms
		% now only works for harmonic modes with omega*n inside h1term!
		% h1term{2,2} contains omega for shift == 0
		
		idx      = num2cell(treeMPS.treeIdx+1);                       % index in para.*
		
		if treeMPS.height == 0
			for ii = 1:treeMPS.L
				treeMPS.op.h1j = []; treeMPS.op.h1jOBB = [];
% 				temp = genh1h2term_onesite(para,treeMPS.treeIdx,ii);		% would be better way
				[bp,bm,n]		  = bosonop(treeMPS.dk(ii),treeMPS.shift(ii),para.parity);
				% need: w/2*(x^2-1/2)
				treeMPS.op.h1term{ii}     = treeMPS.op.h1term{ii}(2,2)./2*((bp+bm)*(bp+bm)/2-eye(treeMPS.dk(ii))/2);				% Since have x = (bp+bm)/sqrt(2)
				
				% Also set correct interaction terms:  t * Q_k * Q_k+1
				treeMPS.op.h2term{1,1,ii} = (treeMPS.op.h2term{1,1,ii} + treeMPS.op.h2term{1,1,ii}')./2;
				treeMPS.op.h2term{2,1,ii} = (treeMPS.op.h2term{2,1,ii} + treeMPS.op.h2term{2,1,ii}')./2;
			end
		else
			%% for now nodes don't have to be changed since they don't carry Harmonic Modes
			treeMPS.op.h1j = []; treeMPS.op.h1jOBB = [];
			for mc = 1:treeMPS.degree
				treeMPS.child(mc) = setPotOp(treeMPS.child(mc));
			end
		end
		
	end
end
function [w, err, hump, E] = expvCustom( t, A, v, para, op)
%% A is string telling the operator to construct:
%    A = {'Hn',   'HAA',  'HAV',   'Kn'}
%  v contains the tensor to evolve in time:
%    v = {mps{j}, Vmat{j}, CV{j},   C{j}}
%  mps contains mps A-tensor if v ~= Amat
%    v2= {[],     mps{j},  mps{j},  mps{j}}
%  Vmat contans V-tensor if v ~=Vmat
%    v3= {[],      [],     Vmat{j},  [] }

dv = size(v);					% Dimensions of input tensor 
n = numel(v);
tol = para.tdvp.expvTol;
m   = min(n,para.tdvp.expvM);
M = para.M;
method = 1;						% 0: Arnoldi; 1: Lanczos with part. orth.

%% set Function handle
switch A
	case 'HAA'
		AFUN = @HAAmultV;		% also works with Multi-Chain
	case 'HAV'
		% copy OBB into h1j, h2j to match HmultVmat behaviour. This is only called for Bosons, so is fine!
		% needed to use HmultVmat for CV
		op.h1j = op.h1jOBB;
		op.h2j = op.h2jOBB;

		AFUN = @HAVmultCV;
	case 'Hn'
		AFUN = @HnmultA;
% 		AFUN = @(x) HmultA(x, op, dv(1), dv(2), dv(3), M,para.parity,[]);		% execution speed comparable to nested function above! has approx. 10% overhead
	case 'Kn'
		AFUN = @KnmultCA;
	case 'MC-HAS'               % multi-chain HOSVD exponentials
		AFUN = @MCmultV;
	case 'MC-HASV'
		AFUN = @MCmultCV;
	case 'MC-HAV'
		AFUN = @MCmultVS;
	case 'MC-HAVS'
		AFUN = @HAVmultCV;      % same as 'HAV'
	case 'STAR-Hn1'
		AFUN = @STARmultA;
	case 'STAR-Hn1Trotter'
		AFUN = @STARmultATrotter;
	case 'TREE-Hn1'
		if para.hasSite
			AFUN = @TREEmultA;
			assert( size(op.h2jOBB,3) == length([para.child.chainIdx]), 'VMPS:expvCustom:TREE-Hn1','There are not enough chain terms in h2jOBB for all child nodes!');
		else
			AFUN = @TREEmultA_NoSite;
			assert( size(op.Opstorage,4) == length([para.child.chainIdx]), 'VMPS:expvCustom:TREE-Hn1','There are not enough interaction terms from parent node for all child nodes!');
		end
	case 'TREE-Hn1Trotter'
		AFUN = @TREEmultATrotter;
	
end

%%
% anorm, rndoff are estimated using norm(H) later
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;

k1 = 2; xm = 1/m;
w = reshape(v,[n,1]);										% reshape since v = tensor
normv = norm(w); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = t_out-t_now;
% t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;			% properly approximated using norm(H)
% s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;
sgn = sign(t);

hump = normv;
while t_now < t_out
  nstep = nstep + 1;
  t_step = min( t_out-t_now,t_new );
  V = zeros(n,m+1);
  H = zeros(m+2,m+2);

% Start finding Krylov subspace using v_n = A v_n-1
  V(:,1) = (1/beta)*w;
  if method == 0						% Arnoldi Iteration for general Operator A
	  for j = 1:m
		 p = AFUN(V(:,j));				% AFUN is handle to specific function
		 for i = 1:j					% Modified Gram-Schmidt orthogonalisation
			H(i,j) = V(:,i)'*p;
			p = p-H(i,j)*V(:,i);
		 end;
		% Reorthogonalize to ensure orthogonal V-matrix. This is time consuming and might not be essential to the algorithm
		% Nevertheless, without it, V can deviate quite strongly!
	% 		nonOrthRest = V(:,1:j)' * p;	% vectorised Modified Gram-Schmidt orthogonalisation, based on orthogonal V(:,1:j), slower than for-loop!
	% 		while norm(nonOrthRest) > btol
	% 			p = p - V(:,1:j)*nonOrthRest;
	% 			nonOrthRest = V(:,1:j)' * p;
	% 		end
		 s = norm(p);
		 if s < btol,					% if residual orthogonal vector shorter than btol -> invariant subspace -> step out
			k1 = 0;
			mb = j;
			t_step = t_out-t_now;
			break;
		 end;
		 H(j+1,j) = s;
		 V(:,j+1) = (1/s)*p;
	  end;
  elseif method == 1					% Lanczos Iteration with partial orthogonalisation
	  W = zeros(m,m);					% measure of orthogonality
	  W(1,1) = 1;
	  theta = eps;						% TODO: replace anorm
	  psi   = eps;						% TODO: replace anorm
	  ALPHA = zeros(m+2,1);
	  BETA  = zeros(m+1,1);
	  BETA(1) = beta;
	  for j = 1:m
		 p = AFUN(V(:,j));				% r_(j+1)
		 if j > 1						% V(:,0) = 0
			 p = p - BETA(j)*V(:,j-1);  % r_(j+1) = r_(j+1) - beta_j * v_j
		 end
		 ALPHA(j) = V(:,j)'*p;			% alpha_j
		 p = p - ALPHA(j)*V(:,j);		% r_(j+1) = r_(j+1) - alpha_j * v_j
		 
		 s = norm(p);
		 if s < btol							% first short cut
			k1 = 0;
			mb = j;
			t_step = t_out-t_now;
			break;
		 end
		 
		 BETA(j+1) = s;
		 
		 % update measure of loss of orthogonality
		 theta = eps*(BETA(2:end) + BETA(j+1));		%.*randn(m,1)*0.3;			% randn takes too much time
		 if j == 1
			 W(j+1,1)     = BETA(2)*W(j,2);
		 else
			 W(j+1,1)     = BETA(2)*W(j,2) + (ALPHA(1)-ALPHA(j))*W(j,1) - BETA(j)*W(j-1,1);
		 end
		 W(j+1,1) = W(j+1,1)/BETA(j+1) + theta(1);
		 
		 if j > 2
			 idx = 2:j-1;
			 W(j+1,idx) = (BETA(idx+1).'.*W(j,idx+1) + (ALPHA(idx)-ALPHA(j)).'.*W(j,idx) - BETA(idx).'.*W(j,idx-1) - BETA(j)*W(j-1,idx))/BETA(j+1) + theta(idx)';
		 end
		 psi = eps*BETA(2)/BETA(j+1)*n;		%*randn(1)*0.6;							% psi_(j+1)
		 W(j+1,[j,j+1]) = [psi,1];
		 
		 % reorthogonalise if needed:
		 nonOrth = abs(W(j+1,1:j)) > 10^-6;
		 if any(nonOrth)
			 %nonOrth = logical(conv(single(nonOrth),[1,1,1],'same'));						% extend re-orth region
			 nonOrth = find(abs(W(j+1,1:j)) > 10^-10);
			 for i = nonOrth
				 p = p - V(:,i)*(V(:,i)' * p);
			 end
			 W(j+1,nonOrth) = eps;
		 end
		 
		 s = norm(p);
		 if s < btol 							% second short cut
			k1 = 0;
			mb = j;
			t_step = t_out-t_now;
			break;
		 end
		 
		 BETA(j+1) = s;					% beta_(j+1)
		 V(:,j+1) = p/s;
	  end
	  H = sparse(2:m+1, 1:m, BETA(2:end), m+2, m+2);
	  H = H + H.'+diag(ALPHA);
  end
  if nstep == 1				% was moved down from initial variable definitions
	  anorm = norm(H);
	  t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
	  s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;				% round to 2 significant digits!
	  t_step = min( t_out-t_now,t_new );		% re-estimate best t_step
	  rndoff= anorm*eps;
  end
%   fprintf('step: %s, anorm: %g, norm(H): %g, j: %d\n', A, anorm, norm(H), j);		% use to confirm: anorm == norm(H) for the corresponding subspace of v
  if k1 ~= 0,						% i.e. k1 = 2; = 0 if Krylov subspace < m
     H(m+2,m+1) = 1;
     avnorm = norm(AFUN(V(:,m+1)));		% was: norm(A*V(:,m+1));
  end;
  ireject = 0;
  while ireject <= mxrej,			% mxrej = 10 -> 10 tries for good t_step
     mx = mb + k1;							% mb: number of basis vectors -> dim of Krylov
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
		err_loc = btol;						% if dim(Krylov) < m -> precise result -> break while
        break;
	 else % estimate local error
        phi1 = abs( beta*F(m+1,1) );
        phi2 = abs( beta*F(m+2,1) * avnorm );
        if phi1 > 10*phi2,
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2,
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
           xm = 1/(m-1);
        end;
     end;
     if err_loc <= delta * t_step*tol,
        break;
	 elseif ireject == mxrej,
	   error('The requested tolerance is too high.');
	 else
		warning('expvCustom needs smaller t_step, error too large');
		t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        ireject = ireject + 1;
     end;
  end;
  % err_loc small enough @t_step -> apply the exponential
  mx = mb + max( 0,k1-1 );
  w = V(:,1:mx)*(beta*F(1:mx,1));					% calculate resulting w = exp(tA)v
  beta = norm( w );
  hump = max(hump,beta);

  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;	% estimate t_step for next loop
  s = 10^(floor(log10(t_new))-1);
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;

if nargout == 4
	% mainly useful for imaginary time evolution to monitor convergence.
	E = w'*AFUN(w)/(w'*w);
end
%%	Nested functions
%	used to have shared workspace for dv, op, para etc.
%	at least faster than anonymous functions
%
	function w = HAAmultV(V)
		%% HAAmultV(V)
		% Needs previous calculation of op.HrightA, op.HleftA, op.OprightA, op.OpleftA
		% works with multi-chain model!
		
		w = HmultVmat(V, op, dv(1),dv(2), M,para.parity);
	end

	function w = HAVmultCV(CV)
		%% For backward evolution of Center V
		% HAV_(n^',n~',n^,n~) = Vmat*_(n',n^')* HAA_(n',n~',n,n~) Vmat_(n,n^)
		%   w(n^',n~') = HAV_(n^',n~',n^,n~) * CV_(n^,n~)
		% expects h1j, h2j in OBB
		% Vmat = Vmat_(n,n^); since focus is taken to CV_(n^,n~)

		% CV = reshape(CV,dv);
		
		% reuse HmultVmat. Save all temp-results in w to save memory!
		w = HmultVmat(CV, op, dv(1),dv(2), M,para.parity);

	end

	function w = HnmultA(A)
		%% Multiplies H(n) with MPS-A
		%	uses old function HmultA.m since equivalent! Only now for
		%	parity = 'n'
		%	expect op.h1j, op.h2j transformed to OBB if para.useVmat = 1
		%   since para should not be passed to here.
		
		w = HmultA(A, op, dv(1), dv(2), dv(3), M,para.parity,[]);
		

	end

	function w = KnmultCA(CA)
		%%
		% CA_(r^,r) given as vector, r^ indicates now BondDim from SVD on A
		%	is sweep direction dependent! C_(l,l^) and C_(r^,r)
		% A_(l,r^,n) in mps as tensor.
		w = 0;
		NC = size(op.h2jAV,2);		% for multiple chain-channels through bond
		CA = reshape(CA, dv);
		switch para.sweepto
			case 'r'
				w = w + op.HleftAV * CA + CA * op.Hright.';
				for l = 1:NC
					for k = 1:M
						if ~isempty(op.Opright{k,l}) && ~isempty(op.h2jAV{k,l})
							w = w + op.h2jAV{k,l} * CA * op.Opright{k,l}.';
						end
					end
				end
			case 'l'
				w = w + op.Hleft * CA + CA * op.HrightAV.';
				for l = 1:NC
					for k = 1:M
						if ~isempty(op.Opleft{k,l}) && ~isempty(op.h2jAV{k,l})
							w = w + op.Opleft{k,l} * CA * op.h2jAV{k,l}.';
						end
					end
				end
		end

		% input V is always vectorized -> reshape
		w = reshape(w, [numel(w),1]);
	end

	function w = MCmultV(V)
		%% w = MCmultV(V)
		% called by 'MC-HAS'
		nc = para.currentChain;
		nTerms = para.M / para.nChains;
		iM      = 1:para.M;
		iM(ceil(iM/nTerms) ~= nc) = [];                       % could be moved to the top!

		V = reshape(V, dv);

		w = V * op.HnonInt.' + op.h1j{nc} * V;                % V_(n,n~) * HnonInt_(n~',n~)
		for k = iM
			w = w + op.h2j{k,2,nc} * (V * op.OpleftAS{k}.');
			w = w + op.h2j{k,1,nc} * (V * op.OprightAS{k}.');
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = MCmultCV(CV)
		%% w = MCmultCV(CV)
		% called by 'MC-HASV'
		% Similar to MCmultV only h[1/2]j <-> h[1/2]jMCOBB
		nc = para.currentChain;
		nTerms = para.M / para.nChains;
		iM      = 1:para.M;
		iM(ceil(iM/nTerms) ~= nc) = [];                       % could be moved to the top!

		CV = reshape(CV, dv);

		w = CV * op.HnonInt.' + op.h1jMCOBB{nc} * CV;         % CV_(n,n~) * HnonInt_(n~',n~)
		for k = iM
			w = w + op.h2jMCOBB{k,2,nc} * (CV * op.OpleftAS{k}.');
			w = w + op.h2jMCOBB{k,1,nc} * (CV * op.OprightAS{k}.');
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = MCmultVSOld(VS)
		%% w = MCmultVS(VS)
		% called by 'MC-HAV'
		d = dv;													% these are the dimensions of VS
		NC = para.nChains;
		vs = reshape(VS,[prod(d(1:end-1)), d(end)]);
		w = vs * (op.HleftA.' + op.HrightA.');
		w = reshape(w,d);										% store w as fully ordered tensor

		VS = reshape(VS,d);
		for k = 1:NC
			if isempty(op.h12jAV{k}), continue; end;
			order = [k, NC+1, 1:(k-1), (k+1):NC];
			dnow = d(order);
			vs = permute(VS, order);
			vs = reshape(vs, [prod(dnow([1,2])), prod(dnow(3:end))]);       % now ready for multiplications

			wk = op.h12jAV{k} * vs;
			wk = reshape(wk, dnow);
			w  = w + ipermute(wk,order);
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = MCmultVS(VS)
		%% w = MCmultVS(VS)
		% called by 'MC-HAV'
		d = dv;												% these are the dimensions of VS
		nTerms = para.M/para.nChains;
		NC = para.nChains;
		vs = reshape(VS,[prod(d(1:end-1)), d(end)]);
		w = vs * (op.HleftA.' + op.HrightA.');
		w = reshape(w,d);                                   % store w as fully ordered tensor

		VS = reshape(VS,d);
		for k = 1:NC
			vs = tensShape(VS, 'unfold', k, d);				% chain index to front!
			if ~isempty(op.h1jMCOBB{k})
				vsTemp = op.h1jMCOBB{k} * vs;
				w = w + tensShape(vsTemp, 'fold',k,d);		% reshape and add to result
			end
			Mc = 1:para.M;
			Mc(ceil(Mc/nTerms) ~=k) = [];					% find which m to operate on for this chain

			for mm = Mc										% rather slow loop since taking out of M
				if ~isempty(op.h2jMCOBB{mm,2,k})
					vsTemp = op.h2jMCOBB{mm,2,k} * vs;
					vsTemp = tensShape(vsTemp, 'fold',k,d);
					w = w + contracttensors(vsTemp, NC+1, NC+1, op.OpleftA{mm}.',2,1);
				end

				if ~isempty(op.h2jMCOBB{mm,1,k})
					vsTemp = op.h2jMCOBB{mm,1,k} * vs;
					vsTemp = tensShape(vsTemp, 'fold',k,d);
					w = w + contracttensors(vsTemp, NC+1, NC+1, op.OprightA{mm}.',2,1);
				end
			end
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = STARmultA(A)
		%% w = STARmultA(A)
		% called by 'STAR-Hn1'
		d = dv;												% size(A)
		NC = para.nChains;
		nTerms = para.M/NC;

		A = reshape(A,d);

		% 1. on-site H1
		w =	contracttensors(A, NC+2, NC+2, op.h1jOBB.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(A, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.chain(mc).Hlrstorage{1} * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:nTerms
				systemM = nTerms*(mc-1) + mm;				% position in op.h2j

				OpTemp = op.chain(mc).Opstorage{mm,2,1} * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w+ contracttensors(OpTemp, NC+2, NC+2, op.h2jOBB{systemM,1}.',2,1);
% 				w = w + OpTemp;

% 				OpTemp = contracttensors(A, NC+2, NC+2, op.h2jOBB{systemM,1}.',2,1);
% 				OpTemp = contracttensors(OpTemp, NC+2, mc+1, op.chain(mc).Opstorage{mm,2,1}.',2,1);
% 				OpTemp = permute(OpTemp,ord);

			end
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = STARmultATrotter(A)
		%% w = STARmultATrotter(A)
		% Trotter splitting in chains
		% called by 'STAR-Hn1Trotter'
		mc = para.currentChain;																% chain to be evolved
		dA = dv;																			% current shape of A
		NC = para.nChains;
		nTerms = para.M/NC;

		A = reshape(A,dA);

		% 1. on-site H1
		w =	contracttensors(A, 3, 3, (op.h1jOBB.')./NC, 2,1);								% Symmetrically divide h1jOBB onto each chain

		Atemp = tensShape(A, 'unfold', 2, dA);												% chain index to front
		% 2. non-interacting Hlrstorage (Hright)
		OpTemp = op.chain(mc).Hlrstorage{1} * Atemp;
		w = w + tensShape(OpTemp,'fold',2, dA);

		% 3. all interacting parts
		for mm = 1:nTerms
			systemM = nTerms*(mc-1) + mm;													% position in op.h2j

			OpTemp = op.chain(mc).Opstorage{mm,2,1} * Atemp;								% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
			OpTemp = tensShape(OpTemp,'fold',2, dA);
			w = w + contracttensors(OpTemp, 3, 3, op.h2jOBB{systemM,1}.',2,1);
		end

		w = reshape(w, [numel(w),1]);
	end

	function w = TREEmultA(A)
		%% w = TREEmultA(A)
		% called by 'TREE-Hn1'
		% for node which hasSite
		% para == treeMPS
		d = dv;												% size(A)
		NC = para.degree;
		
		A1 = reshape(A,[],d(NC+2));							% good shape for 3)
		A = reshape(A,d);
		
		% 1. on-site H1
% 		w =	contracttensors(A, NC+2, NC+2, op.h1jOBB.', 2,1);
		w = reshape(A1 * op.h1jOBB.',d);
		
		
		NCoupOffset = 0;									% needed for child nodes being entanglement renormalisation tensors
		% 2. iterate through children
		for mc = 1:NC
			ordTo = [mc+1:NC+2,1:mc];
			Atemp = permute(A,ordTo);
			Atemp = reshape(Atemp,d(mc+1),[]);
			% 2. non-interacting Hlrstorage (Hright)
			wTemp = para.child(mc).op.Hlrstorage{1} * Atemp;
			
			% 3. all interacting parts
			NCoup = size(para.child(mc).op.Opstorage,4);
			for nn = 1:NCoup								% iterate through chains of child
				for mm = 1:para.M
					OpTemp = A1 * op.h2jOBB{mm,1,nn+NCoupOffset}.';
					OpTemp = reshape(OpTemp,d);
					OpTemp = permute(OpTemp,ordTo);
					OpTemp = reshape(OpTemp,d(mc+1),[]);
					wTemp  = wTemp + para.child(mc).op.Opstorage{mm,2,1,nn} * OpTemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				end
			end
			
			wTemp = reshape(wTemp, d(ordTo));
			w = w + permute(wTemp, [NC+3-mc:NC+2,1:NC+2-mc]);									% see tensShape 'fold' for details. here: i = mc+1, r = NC+1
			
			NCoupOffset = NCoupOffset + NCoup;
		end
		if ~para.isRoot
			% also contract with parent operators
			% 1. Contract to non-interacting parts
			w = w + contracttensors(op.Hlrstorage{1},2,2, A, NC+2, 1);								% order unchanged, Focused on Node -> Hlrstorage{1} is Hleft of Node
			
			% 2. Contract the other chain-system interacting parts
			for mm = 1:para.M
				OpTemp = contracttensors(A, NC+2, NC+2, op.h2jOBB{mm,2,1}.',2,1);					% parent is 1st, node is 2nd position in H_int !!
				w      = w + contracttensors(op.Opstorage{mm,1,1},2,2, OpTemp, NC+2, 1);			% (m,1,1) should be the operator of parent in the effective left basis for node site 1
			end
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = TREEmultATrotter(A)
		%% w = TREEmultATrotter(A)
		% Trotter splitting in chains
		% called by 'TREE-Hn1Trotter'
		mc = para.currentChain;								% chain to be evolved; 0: parent
		dA = dv;											% current shape of A
		NC = para.degree;
		if ~para.isRoot
			NC = NC+1;		% divide by 1 more due to additional evolution step
		end
		
		A = reshape(A,dA);
		
		% 1. on-site H1
		w =	contracttensors(A, 3, 3, (op.h1jOBB.')./NC, 2,1);		% Symmetrically divide h1jOBB onto each chain
		Atemp = tensShape(A, 'unfold', 2, dA);		% chain index to front
		
		if mc ~= 0
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = para.child(mc).op.Hlrstorage{1} * Atemp;
			w = w + tensShape(OpTemp,'fold',2, dA);

			% 3. all interacting parts
			for mm = 1:para.M
				OpTemp = para.child(mc).op.Opstorage{mm,2,1} * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',2, dA);
				w = w + contracttensors(OpTemp, 3, 3, op.h2jOBB{mm,1,mc}.',2,1);
			end
		elseif mc == 0 && ~para.isRoot
			% 2. non-interacting Hlrstorage (Hleft)
			OpTemp = op.Hlrstorage{1} * Atemp;
			w = w + tensShape(OpTemp,'fold',2, dA);

			% 3. all interacting parts
			for mm = 1:para.M
				OpTemp = op.Opstorage{mm,1,1} * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',2, dA);
				w = w + contracttensors(OpTemp, 3, 3, op.h2jOBB{mm,2,1}.',2,1);
			end
		end
		w = reshape(w, [numel(w),1]);
	end

	function w = TREEmultA_NoSite(A)
		%% w = TREEmultA_NoSite(A)
		% called by 'TREE-Hn1'
		% for node which ~hasSite; allow only non-root nodes to have no sites
		% para == treeMPS
		d  = dv;							% size(A)
		NC = para.degree;
		
		A1 = reshape(A,d(1),[]);			% good shape for 3)
		A = reshape(A,d);
		
		
		% 1. initialise w
		w =	0;
		
		NCoupOffset = 0;									% needed for child nodes being entanglement renormalisation tensors
		% 2. iterate through children
		for mc = 1:NC
			% Order for permute to contract with child mc
			ordTo = [mc+1:NC+1,1:mc];
			Atemp = permute(A,ordTo);
			Atemp = reshape(Atemp,d(mc+1),[]);
			% 2. non-interacting Hlrstorage (Hright)
			wTemp = para.child(mc).op.Hlrstorage{1} * Atemp;
			
			% 3. all interacting parts
			NCoup = size(para.child(mc).op.Opstorage,4);
			for nn = 1:NCoup								% iterate through chains of child
				for mm = 1:para.M
					OpTemp = op.Opstorage{mm,1,1,nn+NCoupOffset} * A1;
					OpTemp = reshape(OpTemp,d);
					OpTemp = permute(OpTemp, ordTo);
					OpTemp = reshape(OpTemp, d(mc+1),[]);
					wTemp  = wTemp + para.child(mc).op.Opstorage{mm,2,1,nn} * OpTemp;			% (m,2,1,n) should be the operator of child site 1, chain n in the effective right basis for node
				end
			end
			% reshape wTemp back, bring into correct order and add to w
			wTemp = reshape(wTemp, d(ordTo));
			w = w + permute(wTemp, [NC+2-mc:NC+1,1:NC+1-mc]);									% see tensShape 'fold' for details. here: i = mc+1, r = NC+1
			
			NCoupOffset = NCoupOffset + NCoup;
		end
		
		% also contract with non-interacting parent operator
		% 1. Contract to non-interacting parts
		w = w + reshape(op.Hlrstorage{1} * A1, d);												% order unchanged, Focused on Node -> Hlrstorage{1} is Hleft of Node
		
		w = reshape(w, [numel(w),1]);
	end
end


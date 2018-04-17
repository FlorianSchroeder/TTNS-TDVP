function [chainpara,xi,Gamma]=SBM_genpara(chainpara)
%Calculate the wilson chain parameters for spin-boson model using logarithmic discretization.
%See Appendix A of Bulla 2005 PRB 71, 045122 for details.
%
% Modified:
%   FS 19/02/2014:  - added automatic L for L = 0;
%	FS 25/05/2014:	- added Spectral function of defined in Renger 2002,
%                     following Müh 2012
%	FS 27/05/2014:	- Changed to J(w) as in Renger 2002. This defines for
%                     B777, a monomer of B850, while Müh 2012 needs a
%                     prefactor to use it for pigments in PSII
%   FS 16/11/2014:  - Added analytical result for Orthogonal Polynomial
%                     chain mapping of J(w) of the SBM.
%	FS 12/02/2015:	- Adding Stieltje to OrthPol methods.
%					- Simplifying the conditionals for different methods
%					- Added more control over Analytic extraction and Numerical extraction of chain parameters.
%					  Allows for many more combinations!
%	FS 09/07/2015:	- Tidied up the necessary variables and shortened code
%					- Introduced init_() for cutoff initialization
%	FS 05/10/2015:	- Added J(w) interpolated from points: spectralDensity = 'PointsInterp'
%					  Needed definitions: w_cutoff, datapoints = [w,J(w)]
% 					- Added J(w) interp from broadened peaks: spectralDensity = 'CoupBroad'
% 					  Needed defs: w_cutoff, dataPoints = [w,lambda(w)], peakWidth
%					- Added spectralDensity = 'CoupDiscr'
% 					  direct use of given dataPoints = [w,lambda(w)]

if isfield(chainpara,'s')
	s = chainpara.s;
end
wc = 1;
% initialise fixed model parameters
eval(sprintf('init_%s()',chainpara.spectralDensity));

% first: OrthPol, analytic from recurrence relation.
if strcmp(chainpara.mapping,'OrthogonalPolynomials')
	assert(isfield(chainpara,'L') && chainpara.L>=1,'Please state required finite chain length!');
    assert(isfield(chainpara,'alpha'),'Please state para.alpha, the strength of coupling!');

	[chainpara.epsilon, chainpara.t] = chainParams_OrthogonalPolynomials();			% returns length L parameters to allow pure boson chains
	return;
end

%% Start Discretization
% Used for mappings
%	- Stieltjes
%	- LanzcosTriDiag
% either numerically using integrations (Renger)
% or from analytical result (Leggett)
% if para.Lambda > 1 -> Logarithmic, Zitko
% if para.Lambda = 1 -> Linear

if chainpara.Lambda > 1
	z = chainpara.z;
	bigL = ceil(floor(-1*log(realmin)/log(chainpara.Lambda))/2)
	%Use a large enough start H to make sure the accuracy after transformed to Wilson chain
elseif chainpara.Lambda == 1
	assert(chainpara.L > 1 || strcmp(chainpara.spectralDensity,'CoupDiscr'),'You need to define a chain length > 1!');
	% TODO: make sensible estimate!
	bigL = 10*chainpara.L		% arbitrary now, was 10*
	if strcmp(chainpara.discrMethod,'Direct') && strcmp(chainpara.mapping,'Stieltjes')
		bigL = 100*bigL;							% can deal with much more! was: 100
	end
else
	error('VMPS:SBM_genpara','Please define Lambda >= 1 for discretization!');
end

if strcmp(chainpara.discrMethod,'Analytic')
	assert(strcmp(chainpara.discretization,'LogZ'), 'Analytic discretization only available for LogZ!');
	assert(strcmp(chainpara.spectralDensity,'Leggett_Hard'), 'Analytic discretization only available for Leggett_Hard!');
	% If analytic discretization:
%% use J(w) from A. Leggett et al., Rev. Mod. Phys. 59, 1 (1987). DOI: 10.1103/RevModPhys.59.1
    %           and R. Bulla, N.-H. Tong, and M. Vojta, Phys. Rev. Lett. 91, 170601 (2003). DOI: 10.1103/PhysRevLett.91.170601
    %           and R. Žitko and T. Pruschke, Phys. Rev. B 79, 085106 (2009). DOI: 10.1103/PhysRevB.79.085106
    % Obtained by solving Eq. C4 in Žitko 2009, with f'(x) == 0;
    fac=2*pi*chainpara.alpha/(1+chainpara.s);
    tempfac=1/((1+chainpara.s)*log(chainpara.Lambda));
    tempLam=chainpara.Lambda^(1+chainpara.s);
    tempexp=1/(1+chainpara.s);

    xi=zeros(bigL,1);
    Gamma=zeros(bigL,1);

    % calculate discretized energy levels and couplings
	for j=2:bigL+1
		if j==2
            xi(j-1)=(tempfac*(1-tempLam^(-z))-z+1)^tempexp;
            Gamma(j-1)=sqrt(fac*(1-tempLam^(-z)));                      % modify this after the vector calculation.
        else
            xi(j-1)=(tempfac*tempLam^(2-j-z)*(tempLam-1))^tempexp;
            Gamma(j-1)=sqrt(fac*(tempLam-1)*tempLam^(2-j-z));           % replace this by a vector notation
		end
	end
elseif strcmp(chainpara.discrMethod,'None') && strcmp(chainpara.spectralDensity,'CoupDiscr')
    % if directly tridiagonalizing given couplings!
    % define values via the init_CoupDiscr function.
	bigL = length(xi);
elseif any(strcmp(chainpara.discrMethod,{'Direct','Numeric'}))
	%% from here: either discretize J(w) for Lanczos or h^2(x) for Stieltjes
	% Uses numerical evaluation of Integrals. Works for any spectral function!
	% numerical LogZ uses not the Zitko scheme, which is more difficult!!

	% assign spectral function to J:
	switch chainpara.spectralDensity
		case 'Leggett_Hard'
			J = @J_Leggett_Hard;
		case 'Leggett_Soft'
			J = @J_Leggett_Soft;
		case 'Renger'
			J = @J_Renger;
		case 'PointsInterp'
% 			J = @J_PointsInterp;   Do not do that!
            J = @(w,i) interp1(w_i, j_i, w,'spline')./(w.^i).*ceil(heaviside(w_cutoff-w));
        case 'CoupBroad'
%             J = @(w,i) interp1(w_i, j_i, w,'pchip')./(w.^i).*ceil(heaviside(w_cutoff-w));
			J = @(w,i) interp1(w_i, j_i, w,'pchip')./(w.^i).*(w_cutoff>=w);
		otherwise
			error('VMPS:SBM_genpara:UnknownSpectralDensity',...
				  'cannot not find J_%s(w,i); Please define para.spectralDensity properly!',chainpara.spectralDensity);
	end
	wmax = w_cutoff;

    xi	  = zeros(bigL,1);
    Gamma = zeros(bigL,1);			% not overloading gamma function!

	if strcmp(chainpara.discretization,'LogZ')
		% gives different result than 'analytic', 'Leggett_Hard'
		% Due to additional Discr. detail of differential equation!
		% TODO: Problem: 0 is not included as boundary! need: w_limits = [w_limits,0];
		w_limits = chainpara.Lambda.^(1-z-(0:bigL)).*wmax;      % Define limits of all the discretization intervals
		w_limits(1) = wmax;
	elseif strcmp(chainpara.discretization,'Linear')
		% Divide spectrum linearly
		w_limits = wmax.*(1:-1/bigL:0);
	end

	if strcmp(chainpara.discrMethod,'Numeric')
		for j=1:bigL
			intJ	  = integral(@(w) J(w,0)./pi,w_limits(j+1),w_limits(j));
			intJoverW = integral(@(w) J(w,1)./pi,w_limits(j+1),w_limits(j));
			xi(j)	  = intJ/intJoverW;							% bath energy levels
			Gamma(j)  = sqrt(intJ);								% bath couplings.
		end
% 		figure(1);plot(xi,Gamma);
	elseif strcmp(chainpara.discrMethod,'Direct')
		xi    = flipud(w_limits.');								% flip necessary to make vectorised code work
		difXi = diff(xi);										% width of each interval
		xi    = xi(1:end-1)+difXi./2;							% take middle of each interval
		Gamma = sqrt(J(xi,0).*difXi./pi);						% need to calc sqrt(area)
%         figure(1);plot(xi,Gamma);
	else
		error('VMPS:SBM_genpara:WrongDiscrMethod','Please choose discrMethod = [Direct | Numeric]');
	end
else
	error('VMPS:SBM_genpara','Unsupported combination of discrMethod & spectralDensity');
end

xi(isnan(xi))=0;                                    % if numbers < e-161 they become NaN. This is fix for it.
Gamma(isnan(Gamma))=0;
% Needed for mapping chain->star
chainpara.xi = xi;
chainpara.Gamma = Gamma;

%% Start mapping to Chain
% keep one more t to allow chain->star mapping: t = L+1 x 1
if strcmp(chainpara.mapping, 'LanczosTriDiag')
	%% Start Tridiagonalization
	[chainpara.epsilon, chainpara.t, chainpara.U] = chainParams_Lanczos(xi, Gamma);

	%chainpara.t=abs(chainpara.t); %I noticed that the r form hess() function sometimes changed sign.
elseif strcmp(chainpara.mapping, 'Stieltjes')
	%% Start Stieltjes
	% assume: xi = levels = x
	%		  gamma = coupling strengths = h^2(x)
    [chainpara.epsilon, chainpara.t] = chainParams_Stieltjes(xi, Gamma);

else
	error('VMPS:SBM_genpara:l179','This is really hard to reach!');
end

    function [w, t] = chainParams_OrthogonalPolynomials()
        %% returns w = epsilon and t for J(w) from A. Leggett
        % using orthogonal polynomial mapping from A. W. Chin et al. J. Math. Phys. 51, 092109 (2010). DOI: 10.1063/1.3490188
        % No need for tridiagonalization since analytical result.
        % t(1) == sqrt(eta_0/pi)/2
% 		wc = 1;     % cutoff is set to 1 always!
        n = (0:chainpara.L-1)';
		if ~isempty(strfind(chainpara.spectralDensity,'Leggett')) && strcmp(chainpara.discretization, 'None')
			if ~isempty(strfind(chainpara.spectralDensity,'Hard'))
				w = w_cutoff/2.*(1+ (s^2./((s+2.*n).*(2+s+2.*n))));                                     % w(1) = w(n=0)
				t = (w_cutoff.*(1+n).*(1+s+n))./(s+2+2.*n)./(3+s+2.*n).*sqrt((3+s+2.*n)./(1+s+2.*n));   % t(1) = t(n=0) != coupling to system
				t = [ w_cutoff*sqrt(2*pi*chainpara.alpha/(1+s))/(sqrt(pi)); t];				% t(1) = sqrt(eta_0/pi)/2, 1/2 from sigma_z; t(2) = t(n=0)
			elseif ~isempty(strfind(chainpara.spectralDensity,'Soft'))
				w = wc.*(2.*n+1+s);																		% w(1) = w(n=0)
				t = wc.*sqrt((n+1).*(n+s+1));															% t(1) = t(n=0) != coupling to system
				t = [ wc*sqrt(2*pi*chainpara.alpha*gamma(s+1))/(sqrt(pi)); t];				% t(1) = sqrt(eta_0/pi)/2, 1/2 from sigma_z; t(2) = t(n=0)
			else
				error('VMPS:SBM_genpara:chainParams_OrthogonalPolynomials','Choose spectralDensity = [Leggett_Soft|Leggett_Hard]');
			end
		else
			error('VMPS:SBM_genpara:chainParams_OrthogonalPolynomials','Only available for power-law spectral functions and discretization = ''None''.');
		end
    end

    function [w, t, U] = chainParams_Lanczos(xi, Gamma)
        % Tridiagonalize using Lanczos algorithm with complete reorthogonalization
		% U is transformation U_(k,n), useful for chain->star mapping
        indiag=zeros(length(Gamma)+1,1);
        inrow = indiag;

        inrow(2:end)  = Gamma;                         % factor of 1/2 only affected t(1) -> was moved into Hamiltonian; 1/sqrt(pi) changes only t(1)
        indiag(2:end) = xi;

        [epsilon,t,U]    = star2tridiag(indiag,inrow);

		if chainpara.L == 0
            chainpara.L = find(epsilon > chainpara.precision, 1,'Last') + 1;    % +1 as L includes spin site
            dispif(sprintf('Optimum chain length is: %u',chainpara.L),chainpara.logging)
		end

		if ~strcmp(chainpara.spectralDensity,'CoupDiscr')
	        [epsilon,t]			= extroplate(epsilon,t,chainpara.L);              % extrapolate very small levels for higher precision
		end
        w = epsilon(1:chainpara.L);
        t = t(1:chainpara.L);
		U = U(:,1:chainpara.L);
    end

    function [w, t] = chainParams_Stieltjes(xi, Gamma)
        ab = stieltjes(chainpara.L,[xi,Gamma.^2]);
        w  = ab(:,1);
		t  = ab(:,2);
		t  = [sqrt(sum(Gamma.^2)); t];                             % put sqrt(eta_0/pi) = t(1) by hand in!
    end

	function y = J_Renger(w,i)
		%% Spectral function for LH2 of Rhodoblastus acidophilus  from: F. Müh and T. Renger, Biochim. Biophys. Acta 1817, 1446 (2012). DOI: 10.1016/j.bbabio.2012.02.016
		% this code section can be used with any definition of J(w)
		% Units: eV
		% need for numerical integration.
		% multiplied by \pi * \omega^2: \pi to equal out division of \sqrt\pi later in the \gamma. Problem due to different
		%               definition of J(w) in Leggett et al.
		%               \omega^2 because Leggett defined J(w) as spectrum of the dimension-ful coupling constant! Renger
		%               took one \hbar\omega out.
		% This J(w) gives a Huang-Rhys factor of 1.3
		% i to get function J*w^(-i)

		%    tic; chainpara.Lambda = 2; bigL = 511; z=1;        % testing
% 		  w_c = 1000; accurate up to 10^-3 eV; w_c = 2000: 10^-6 eV; w_c = 3000: 10^-8 eV
% 			w_cutoff = cmToeV(3000);                            % upper cutoff frequency


% 		w=0:cmToeV(1):cmToeV(3000); i=0;                                                                % testing
		y = pi/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^(5-i).*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/(cmToeV(1.94).^4) .* w.^(5-i).*exp(-(w/cmToeV(1.94)).^(1/2))));
% 		plot(y);xlabel('$\omega / cm^{-1}$');ylabel('$J(\omega)/eV$'); formatPlot(1);                             % testing

		% original definition: J = @(w) 1/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^3.*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/cmToeV(1.94).^4 .* w.^3.*exp(-(w/cmToeV(1.94)).^(1/2))));
	end

	function init_Renger()
% 		wc = ?;
		w_cutoff = cmToeV(16000);                            % upper cutoff frequency
	end

	function y = J_Leggett_Hard(w,i)
		%% Section for Spectral function for Spin-Boson Model from A. Leggett
		% Hard cutoff at wc using stepfunction
		% i to get function J*w^(-i)
		y = 2*pi*chainpara.alpha*w_cutoff.^(1-chainpara.s).*w.^(chainpara.s-i).*(w_cutoff>=w);
	end

	function init_Leggett_Hard()
		wc		 = 1;
		if isfield(chainpara,'w_cutoff')
			w_cutoff = chainpara.w_cutoff;
		else
			w_cutoff = 1;
		end
	end

	function y = J_Leggett_Soft(w,i)
		%% Section for Spectral function for Spin-Boson Model from A. Leggett
		% Soft cutoff at wc using exponential function
		% additional hard cutoff for numeric integration
		% i to get function J*w^(-i)

% 		y = 2*pi*chainpara.alpha*wc^(1-chainpara.s).*w.^(chainpara.s-i).*exp(-w./wc).*~stepfun(w,w_cutoff);
		y = 2*pi*chainpara.alpha*wc^(1-chainpara.s).*w.^(chainpara.s-i).*exp(-w./wc).*ceil(heaviside(w_cutoff-w));	% non-zero for w = w_cutoff
	end

	function init_Leggett_Soft()
		wc		 = 1;				% characteristic frequency ~= w_cutoff!
		w_cutoff = chainpara.w_cutoff;
    end

%     function y = J_Coupling_Discrete(w,i)
%         %%
%         x_i = [0,0.1,0.3,0.5,0.8,1];
%         w_i = [0,0.1,0.5,0.8,0.2,0];
%         y = interp1(x_i,w_i, w,'pchip')./(w.^i).*ceil(heaviside(w_cutoff-w));
%     end

    function init_PointsInterp()
        %default values:
        w_cutoff = 1;
        w_i = [0,0.1,0.3,0.5,0.8,1];
        j_i = [0,0.1,0.3,0.5,0.8,1];
        % overrides:
        if isfield(chainpara,'w_cutoff')
			w_cutoff = chainpara.w_cutoff;
        end
        if isfield(chainpara,'dataPoints')
            w_i = chainpara.dataPoints(:,1);
            j_i = chainpara.dataPoints(:,2);
        end
    end

    function init_CoupDiscr()
		if isfield(chainpara,'dataPoints')
			chainpara.dataPoints = sortrows(chainpara.dataPoints,1);
            xi = chainpara.dataPoints(:,1);
            Gamma = chainpara.dataPoints(:,2);
		else
			error('VMPS:SBM_genpara:init_CoupDiscr','Please give dataPoints with CoupDiscr');
		end
		chainpara.discrMethod = 'None';			% nothing else possible!
	end

	function init_CoupBroad()
		if isfield(chainpara,'dataPoints')
            w_i = chainpara.dataPoints(:,1);			% omega
            j_i = pi.*chainpara.dataPoints(:,2).^2;			% lambda
		end
		w_cutoff = 1;
		if isfield(chainpara,'w_cutoff')
			w_cutoff = chainpara.w_cutoff;
		end
		sigma = 0.01;									% width of broadened peaks
		if isfield(chainpara,'peakWidth')
			sigma = chainpara.peakWidth;
		end
		minRes = min(diff(w_i),sigma)/10;				% the needed freq resolution
		omega = (0:minRes:w_cutoff).';				% omega-axis
		y = 0;									% the final broadened J(w)

		for ii = 1:length(w_i)
			A = 1/(integral(@(w) cauchypdf(w,w_i(ii),sigma)./w, 1e-5,w_cutoff)*w_i(ii));		% normalization constant to keep constant reorganisation Energy
%			A = 1/(integral(@(w) cauchypdf(w,w_i(ii),sigma)./w^2, 1e-5,w_cutoff)*w_i(ii)^2);	% normalization constant to keep constant Huang-Rhys
% 			A = 1;																				% normalization constant to keep constant coupling strength
			y = y + A .* j_i(ii) .*cauchypdf(omega,w_i(ii),sigma);
		end

		% overwrite interpolation points for use above
		w_i = omega;
		j_i = y;

% 		J_CoupBroad = @(w,i) interp1(omega,y,w,'pchip')./(w.^i).*ceil(heaviside(w_cutoff-w));		% Define function for use later
	end

    function y = normpdf(X,mu,sigma)
		% gaussian shape
        y = exp(-(X-mu).^2 ./2 ./sigma.^2)./sigma ./sqrt(2*pi);
	end

	function y = cauchypdf(X,mu,gamma)
		% lorentzian shape
		y = 1./(pi .* gamma .* (1+((X-mu)./gamma ).^2));
	end
end
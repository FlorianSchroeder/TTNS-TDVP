classdef TreeTDVP
	% The TreeTDVP class is a helper class to provide methods used in the TDVP and Ground State search of a TreeMPS
	% It is not meant to be instantiated, instead it provides static methods.
	
	methods(Static=true)
		function [bp,bm,n] = bosonop(dim, shift, parity, invertedOrder)
			%% bosonop - Generate Boson Operators
			%	[bp,bm,n] = bosonop(dim, shift, parity, invertedOrder)
			%
			%	Generate local boson operators
			%	shift sets only offset to diagonal entries
			%	Modified:
			%       FS 24/05/2014: replaced slow for-loop by gallery() to create bp
			%       FS 25/10/2014: replaced slow gallery by faster direct sparse
			
			if nargin > 3 && (isempty(invertedOrder) || ~invertedOrder)
				bp = sparse(2:dim, 1:dim-1, sqrt(1:dim-1), dim, dim);			% 1:n-1
			else
				bp = sparse(1:dim-1, 2:dim, sqrt(dim-1:-1:1), dim, dim);		% n-1:-1:1 inverted
			end
			
			if parity~='n'
				bp = TreeTDVP.parityorderOP(bp);
			end
			
			if nargout == 1 && shift == 0, return, end
			
			bm = bp';
			n  = bp*bm;
			
			if shift~=0
				bp   = full(bp); bm = full(bm); n = full(n);
				x    = sqrt(2)/2.*(bp+bm);
				iden = eye(dim);
				bp   = bp + sqrt(2)/2.*shift*iden;
				bm   = bm + sqrt(2)/2.*shift*iden;
				n    = n  + shift.*x + shift^2/2.*iden;
			end
		end
		
		function newmat = parityorderOP(mat)
			%% parityorderOP - Reorders Boson Operators
			%	B = parityorderOP(A)
			%
			%	For m:0 ordered Boson operators:
			%
			%	A: Boson op in local basis [m:0]
			%	B: Boson op in reordered parity basis [o,e]
			%
			%	Reordering the local boson operators into 2*2 block matrices according to
			%	parity odd and even (number). The shape of the block matrix is :
			%
			%	 oo oe
			%	 eo ee
			%
			%	 So operator
			%	 0 sqrt(3)   0       0
			%	 0   0     sqrt(2)   0
			%	 0   0       0     sqrt(1)
			%	 0   0       0       0
			%
			%	 will become:
			%	 0   0     sqrt(3)   0
			%	 0   0       0     sqrt(1)
			%	 0 sqrt(2)   0       0
			%	 0   0       0       0
			%
			%	For 0:m ordered Boson operators:
			%
			%	A: Boson op in local basis [0:m]
			%	B: Boson op in reordered parity basis [e,o]
			%
			%	 So operator
			%	 0 sqrt(1)   0       0
			%	 0   0     sqrt(2)   0
			%	 0   0       0     sqrt(3)
			%	 0   0       0       0
			%
			%	 will become:
			%	 0   0     sqrt(1)   0
			%	 0   0       0     sqrt(3)
			%	 0 sqrt(2)   0       0
			%	 0   0       0       0
			%
			% Modified:
			%		FS 07/07/15: 43% speedup using rotation operator

			% assert input
			[m,n] = size(mat);
			assert(mod(m,2)==0 && m==n, 'parityorderOP needs square matrix inputs of even dimensions');

			% create reorder operator, perhaps export into toolbox?
			U = sparse(1:m, [1:2:m,2:2:m],1,m,n);			% nicer
			% U has 1s in (row,column) = (1:m, [odd,even])

			% Apply reordering
			newmat = (U*mat*U');

		end
		
		function U = fockToX(dk, shift, parity, invertedOrder, omega, x)
			%% fockToX - Generate Map from Fock space to Position space
			%	U = fockToX(dk, shift, parity, invertedOrder)
			%	Generates the dim(x) x dk operator to map fock states into the position space.
			%	Can be used as: X = U*V
			%	Where V is the dk x dOBB isometry for the OBB
			%	Spatial resolution given in x 
			%
			%	Created by FS 28/04/2017
			
% 			U = arrayfun(@(n) (omega/pi)^(1/4)/(2^n * factorial(n)) * hermiteH(n,x) .* exp(-omega*x.^2/2), 0:dk-1,'UniformOutput',false); 
			x = reshape(x,[],1);									% make sure x is column vector!
			H = TreeTDVP.hermite(dk-1,sqrt(omega)*x);
			U = bsxfun(@times, H, (omega/pi)^(1/4).*exp(-omega/2 * (x.^2)));			% Gaussian envelope
% 			U = bsxfun(@rdivide, U, sqrt((2.^(0:dk-1) .* factorial(0:dk-1))));				% normalise
			U = bsxfun(@rdivide, U, sqrt(sum(U.^2,1)));				% normalise numerically instead of using prefactors.
		end
		
		function H = hermite(n,x)
			%% hermite - Generate Hermite Polynomials
			%	H = hermite(n,x)
			%	Generates all Hermite Polynomials up to order n
			%	Uses the recurrence relation H[n+1](x) = 2x H[n](x) - 2n H[n-1](x)
			%	
			%	H is x by n+1 matrix
			%
			%	Created by FS 28/04/2017
			
			x = reshape(x,[],1);		% make sure x is column vector!
			
			if n == 0
				H = ones(length(x),1);
				return;
			elseif n == 1
				H = 2*x;
				return;
			end
			
			H = ones(length(x), n+1);			% H(:,n+1) contains H[n]
			H(:,2) = 2*x;
			
			for ii = 2:n
				H(:,ii+1) = 2*x.*H(:,ii) - 2*(ii-1)*H(:,ii-1);			% H[n](x) = 2x H[n-1](x) - 2(n-1) H[n-2](x)
			end
		end
		
		function p = plotHermite(n)
			%% plotHermite - Plots the Hermite polynomials
			%	p = plotHermite(n)
			%	Plots the Hermite polynomials up to degree n and returns the plot handle
			%
			%	Created by FS 28/04/2017
			
			x = -3:0.1:3;
			p = plot(x,TreeTDVP.hermite(n,x)); 
			grid on;
		end
		
		function h = plotQHE(n,omega)
			%% plotQHE - Plots the Quantum Harmonic Oscillator wave functions
			%	p = plotQHE(n)
			%	Plots the Quantum Harmonic Oscillator wave functions up to energy level n
			%
			%	p = plotQHE(n,omega)
			%	Plots for defined energy omega
			%
			%	Created by FS 28/04/2017
			
			if nargin == 1
				omega = 1;
			end
			
			xmax = sqrt((2*n+1))/omega*1.4;					% range needed to display all levels + 10% more
			x = -xmax:0.1:xmax; x = x';
			
			% Get the wave functions
			U = TreeTDVP.fockToX(n+1,0,'n',0,omega,x);
			
			% Energy levels
			E = 1/2+(0:n);
			UE = bsxfun(@plus,U,E);			% shifted wave functions to plot ontop of energy levels
			
			h.f = figure; hold all;
			col = get(0,'defaultaxescolororder');
			
			% Plot Parabola
			h.pP(1) = plot(x, omega^2/2.*x.^2, 'k');

			% plot wavefunctions
			% h.pWF	= plot(x,UE);
			h.pWF	= patch(repmat([x;x(1)],1,n+1), [UE;UE(1,:)], (1:n+1),...
				'LineStyle','none');
			
			xlim([-1,1]*xmax*0.9)
			ylim(omega*[-0.5,n+2]);
		end
	end
	
end
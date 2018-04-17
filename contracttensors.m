function X = contracttensors(X, numindX, indX, Y, numindY, indY)
% Contracts tensors X and Y about the indices (dimensions) indX, indY
%	1. Permutes dimensions indX to the last and indY to the first place of X and Y respectively
%	2. Reshapes X and Y into 2-dim matrices s.t. all not contracting dimension are in one dimension. All to be contracted dimensions are reshaped into the other dimension.
% 	3. Matrixproduct X*Y
%	4. Reshape result into remainders: dimXl x dimYr x 1  ; here contracted dimension is at the end
%
% eg: Cleft = contracttensors(conj(B), 3, [1, 3], Cleft, 3, [3, 1]);      % Cleft_fb = B*_dfe Cleft_ebd


Xsize = ones(1, numindX); Xsize(1:ndims(X)) = size(X);
Ysize = ones(1, numindY); Ysize(1:ndims(Y)) = size(Y);

indXl = 1:numindX; indXl(indX) = []; 			% remove indX from index number array for proper permutation
indYr = 1:numindY; indYr(indY) = [];

sizeXl = Xsize(indXl);
sizeX = Xsize(indX); 						% should be just a number otherwise have: indX = [2,5] ; indY = [1,3]  to contract [2 with 1] and [5 with 3]
sizeYr = Ysize(indYr);
sizeY = Ysize(indY);

if prod(sizeX)~= prod(sizeY)
    error('indX and indY are not of same dimension.');
end

if isempty(indYr) 								% if Y will be fully contracted
    if isempty(indXl) 							% if X and Y are to be fully contracted
        X = permute(X, indX); 					% permutes array dimensions
        X = reshape(X, [1, prod(sizeX)]); 		% rearranges columnwise into new stated array shape
        Y = permute(Y, indY);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y; Xsize = 1;
        return;
    else
        X = permute(X, [indXl, indX]);
        X = reshape(X, [prod(sizeXl), prod(sizeX)]);
        Y = permute(Y, [indY]);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y;
        Xsize = Xsize(indXl);
        X = reshape(X, [Xsize, 1]);
        return
    end
end

% if parts from X and Y will be left over
		% conditionals may help for larger tensors!
% if ~all([indXl, indX] == 1:numindX)
	X = permute(X, [indXl, indX]);
% end
X = reshape(X, [prod(sizeXl), prod(sizeX)]);
% if ~all([indY, indYr] == 1:numindY)
	Y = permute(Y, [indY, indYr]);
% end
Y = reshape(Y, [prod(sizeY), prod(sizeYr)]);
% sprintf('cond(X) = %.10g, cond(Y) = %.10g', cond(X),cond(Y))
% a=norm(full(X));b=norm(full(Y));
% sprintf('norm(X) = %.10g, norm(Y) = %.10g', a,b)

X = X * Y;
% fprintf('||X|| ||Y|| /||XY|| = %.10g\n', a*b/norm(full(X)))

% Xsize = [Xsize(indXl), Ysize(indYr)]; 			% new array dimensions are from X to the left, and Y to the right
% numindX = length(Xsize); 					% new total dimension number
X = reshape(X, [Xsize(indXl), Ysize(indYr)]);					% why 1 at the end? To conserve the other dimensions?

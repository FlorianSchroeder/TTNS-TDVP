function [B, U, para, results, sv, vNE] = prepare_onesite_truncate(A,para,sitej,results)
% Truncate Bond dimension via SV of A
%
% Changed:
%   FS 20/10/2014: - removed explicit direction, using para.sweepto now!
%   FS 28/10/2014: - implemented sweep to right for TDVP
%                  - Added keeping 1 more dim if smallest SV too large,
%                    limited by para.tdvp.maxBondDim in 'r' sweep
%                  - major restructuring using truncation function
%	FS 09/04/2015: - in 'r' sweep: apply to site L to normalize MPS. No
%					 truncation or para.D saving!

[D1, D2, d] = size(A);
D_old = para.D;
err = 0;
switch para.sweepto
    case 'r'
        if para.parity=='n'
            %% normal Parity
            A = permute(A, [3, 1, 2]);
            A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
            [B, S, U] = svd2(A);				% Anew with orthonormal columns
			if norm(diag(S)) ~= 1
				S = S./norm(diag(S));					% normalisation necessary for thermal steps
			end
            sv = diag(S);
            if para.tdvp.truncateExpandBonds && sitej ~= para.L         % only used in TDVP. For VMPS add: && sitej <= para.trustsite(end)
                if sv(end) > para.svmintol
					%% Expand A dims
					if sv(end) > para.svmaxtol
						dimincratio = log10(sv(end)/para.svmaxtol)/2;						% aggressive expansion if min(SV) > SVmaxtol.
					else
						dimincratio = log10(sv(end)/para.svmintol)/(abs(log10(sv(end)))-1);		% new scheme, also weighted by current lowest exponent, less aggressive than above
					end
					adddim = ceil(D2*dimincratio);
					adddim = min(adddim, para.tdvp.maxBondDim(end) - para.D(sitej));

					A = cat(2,A,zeros(d * D1,adddim));
                    [B, S, U] = svd2(A);                % better since B contains no zeros but orthonormal vectors!
					if norm(diag(S)) ~= 1
						S = S./norm(diag(S));					% normalisation necessary for thermal steps
					end
                    U = U(:,1:D2);                      % form old Bond Dim in columns == U=cat(2,U,addmat);
				else
                    %% Truncate A dims
                    [B, S, U, err] = truncate(B,S,U);
                end
            end
            DB = size(S, 1);				% new a2 dimension of B
			if sitej ~= para.L
	            para.D(sitej)=DB;
			end
            B = reshape(B, [d, D1, DB]);
            B = permute(B, [2, 3, 1]);
            U = S * U;
        else
            error('Have not implemented different parity yet!');
        end
        resultsSite = sitej;

    case 'l'
        if para.parity=='n'
            A = permute(A, [1, 3, 2]);
            A = reshape(A, [D1, d * D2]);
            [U, S, B] = svd2(A);
			if norm(diag(S)) ~= 1
				S = S./norm(diag(S));					% normalisation necessary for thermal steps
			end
            sv = diag(S);
            %% Start Truncation / Expansion
            if sitej <= min(para.trustsite(end),para.L)
                if sv(end) > para.svmintol
                    %% Expand A dims
                    % Expand A then do SVD since this leaves B with orthonormal
                    % rows. Padding after SVD would only leave zeros.
					if sv(end) > para.svmaxtol
						dimincratio = log10(sv(end)/para.svmaxtol)/2;						% aggressive expansion if min(SV) > SVmaxtol.
					elseif abs(log10(sv(end))) > 2
						dimincratio = log10(sv(end)/para.svmintol)/(abs(log10(sv(end)))-1);	% new scheme, also weighted by current lowest exponent, less aggressive than above
					else
						dimincratio = log10(sv(end)/para.svmintol)/2;						% new scheme for big SV, even less aggressive than above
					end
                    adddim = ceil(D1*dimincratio);          % at least 1 for dimincratio != 0
					if isfield(para,'tdvp')
						if sitej ~= 1 && length(para.tdvp.maxBondDim) > 1
							adddim = min(adddim, para.tdvp.maxBondDim(2) - D1);
						else
							adddim = min(adddim, para.tdvp.maxBondDim(1) - D1);
						end
					end
                    A = cat(1,A,zeros(adddim, d * D2));
                    [U, S, B] = svd2(A);                % better since B contains no zeros but orthonormal vectors!
					if norm(diag(S)) ~= 1
						S = S./norm(diag(S));					% normalisation necessary for thermal steps
					end
                    U = U(1:D1,:);                      % form old Bond Dim in columns == U=cat(2,U,addmat);
                else
                    %% Truncate A dims
                    [U, S, B, err] = truncate(U,S,B);
                    % If Bond Dim truncated then sv(end) < para.svmaxtol
                end
            end

            DB = size(S, 1);
			if sitej ~= 1
	            para.D(sitej-1) = DB;
			end
            B = reshape(B, [DB, d, D2]); B = permute(B, [1, 3, 2]);
            U = U * S;

        else
            error('Have not implemented different parity yet!');
        end
        resultsSite = sitej-1;
end

sv  = diag(S);
vNE = vonNeumannEntropy(S);

if nargin == 4 && resultsSite >= 1
    results.Amat_vNE(resultsSite) = vNE;
    sv(sv < eps) = 0;					% need to set to zero, since otherwise might have artifacts after expansion!
	if length(sv)==para.D(resultsSite) || (length(sv)==para.D(resultsSite)/2 && (para.parity~='n'))
        results.Amat_sv{resultsSite} = sv;
	end
	results.Amat_truncErr(resultsSite) = err;
end

    function [U, S, V, err] = truncate(U,S,V)
    %% Truncates Bond Dimensions
    % keeps always one SV in range [para.svmaxtol, para.svmintol]
    %   if not possible -> Do nothing, expand later
    % keeps always at least para.Dmin many dimensions
    % by Florian Schroeder 29/10/2014
        sv = diag(S);
		err = 0;
		%% Truncate A dims
        keepdims = find(sv > para.svmintol);
        % If smallest SV too large, keep 1 more
		while (sv(keepdims(end))/norm(sv(keepdims)) > para.svmaxtol)
			if keepdims(end)+1 <= length(sv)
                keepdims = [keepdims;keepdims(end)+1];
			else
%                 newBondDim = length(sv);        % = old Dim
                % expand Bond dims! Happens automatically!
			end
		end
        if length(keepdims) < para.Dmin                     % keep at least Dmin bonds
            keepdims = 1:para.Dmin;
        end
        if length(keepdims) < length(sv)				% && sitej < para.trustsite(end)
            U = U(:,keepdims);                          % keep columns
            S = diag(sv(keepdims)./norm(sv(keepdims)));
            V = V(keepdims,:);                          % keep rows
%             newBondDim = length(keepdims);

			% truncation error estimation
			truncdims = 1:length(sv);
			truncdims(keepdims) = [];
			err = sum(sv(truncdims).^2);				% all disregarded SV^2
        end
    end

end

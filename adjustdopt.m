function [op,para,results,mps,Vmat] = adjustdopt(op,para,results,mps,Vmat)
% Adaptively change $d_{opt}$ so that the smallest sigular value s is below the thredhold svmintol
%   Either truncate or expand d_opt per sweep!
% Also change dk using schemes:
%   - para.useDkExpand1
%   - para.useDkExpand2
%
% Cheng Guo
% 17 Feb 2012
%
% Modified:
%	FS 03/02/2014:	- added abs() for calculation of d_opt_change, as this prevents elimination of positive and negative contributions.
%   FS 30/05/2014:  - added that always one sv < svmaxtol is kept, to prevent oscillating expansion / deletion
%                   - added para.d_opt_min as parameter.
%                   - added if-statements to check dimincratio
%                   - Expand dk if tail of Vmat indicates, that expansion is needed.
%   FS 25/10/2014:  - Dk expansion for para.useVmat=0 corrected (eye for
%                     Vmat)
%	FS 21/08/2015:	- Exported estimateDkExpand to use it also for TDVP
%					- added lines for inverse/forward Boson ordering
%
%   TODO: replace for by while loop and remove limitation by trustsite, as this adds artificial step.


if para.logging													% save d_opt history
		results.d_opt{para.loop-1} = para.d_opt;
end

d_opt_old = para.d_opt;

adjSites = 1:para.L;
adjSites(para.spinposition) = [];						% remove spins
adjSites(adjSites > para.trustsite(end)) = [];			% remove j > trustsite

for s = adjSites 			% Only modify boson sites
    %% Adjust d_opt OBB dimensions
    if para.useVmat
        para.increasedk = 0;
        discarddims     = find(results.Vmat_sv{s}<para.svmintol);       % lowest SV below certain threshold --> discard these d_opt dimensions

		dispif(sprintf('Site %d:',s),para.logging);

        if ~isempty(discarddims)
            %% Truncate d_opt
            if para.parity~='n'
                error('Haven not implemented parity code for adaptive d_opt!');
    %                 adddim=ceil(maxp*para.dk(s)/2); %This is according to experience,  if maxp=0.5 then this will add 50 new states.
    %                 if para.dk(s)+2*adddim>para.dkmax
    %                     adddim=(para.dkmax-para.dk(s))/2;
    %                 end
    %                 para.dk(s)=para.dk(s)+2*adddim;
    %                 addmat=zeros(adddim,d_opt);
    %                 Vmat{s}=cat(1,addmat,Vmat{s}(1:Db/2,:),addmat,Vmat{s}(Db/2+1:end,:));
    %                 Vmat{s}(1:adddim,1:d_opt/2)=1e-15*randn(adddim,d_opt/2); %add a lit perturbation in the new added dims
    %                 Vmat{s}(Db/2+adddim+1:Db/2+2*adddim,d_opt/2+1:end)=1e-15*randn(adddim,d_opt/2);
    %            else if para.d_opt(s)-length(discarddims) >= para.d_opt_min            %old 30/05/2014
            else
                % parity == 'n'
                if results.Vmat_sv{s}(discarddims(1)-1) > para.svmaxtol         % if next highest not-discarded element too large (causing expansion in next sweep)
                    discarddims = discarddims(2:end);                           % always keep one element < svmaxtol
                end
                % if Dif > 0: would remove too many dims; if Dif < 0, no
                % problem so set Dif = 0
                difference = para.d_opt_min + length(discarddims) - para.d_opt(s);
                if difference <= 0      % discarddims does not violate d_opt_min
                    difference = 0;
                end                     % else: discarddims would remove too much by amount = diff;
                dispif(sprintf('Truncate OBB by: %d', length(discarddims(difference+1:end))),para.logging)
%                 para.d_opt(s)=para.d_opt(s)-length(discarddims-difference);			% -length(discarddims)+difference
                Vmat{s}(:,discarddims(difference+1:end))=[];                          % only cut the end
                mps{s}(:,:,discarddims(difference+1:end))=[];
                results.Vmat_sv{s}(discarddims(difference+1:end))=[];
				para.d_opt(s) = size(Vmat{s},2);
    %                else                                                       % Don't need this case anymore.
    %                    dispif('Done nothing',para.logging)
    %                end
    %                 adddim=ceil(maxp*para.dk(s)); %This is according to experience,  if maxp=0.5 then this will add 50 new states.
    %                 if para.dk(s)+adddim>para.dkmax
    %                     adddim=para.dkmax-para.dk(s);
    %                 end
    %                 para.dk(s)=para.dk(s)+adddim;
    %                 addmat=zeros(adddim,d_opt);
    %                 Vmat{s}=cat(1,addmat,Vmat{s});
            end
        else
            %% Expand d_opt.
            % from here: no SV < para.svmintol
			[dk, dOBB] = size(Vmat{s});
            dimincratio=log10(results.Vmat_sv{s}(end)/para.svmaxtol)/2;						% Imperial relation
            if dimincratio > 0
                adddim=ceil(dOBB*dimincratio);												% increase dim (only if addim > 0)
                if dOBB + adddim >= dk
                    exceeded = dOBB + adddim - dk +1;										% By how much dk is too small. +1 as d_opt != dk
                    adddim = adddim - exceeded;												% Only increase to: d_opt = dk - 1
                    if exceeded > para.increasedk
                        para.increasedk = exceeded;											% log highest needed addition
                    end
                    dispif(sprintf('Expand OBB by: %d;  More needed: %d', adddim, exceeded),para.logging)
                else
                    dispif(sprintf('Expand OBB by: %d', adddim),para.logging)
                end
                if results.Vmat_sv{s}(end) > para.svmaxtol && dOBB + adddim < dk			% Expand d_opt; this if could be removed.
                    dispif('..Expanding OBB ',para.logging)
                    para.d_opt(s) = dOBB + adddim;
                    addmat        = zeros(dk,adddim);
                    Vmat{s}       = cat(2,Vmat{s},addmat);									% TODO: also expand SV vector for correctness?
                    [a,b,~]       = size(mps{s});
                    addmat        = zeros(a,b,adddim);
                    mps{s}        = cat(3,mps{s},addmat);
                end
            %else: svmaxtol > smallest SV > svmintol --> good, nothing to do
            end
        end
    end
    %% Expand dk
	% only works with OBB now!
    if para.useDkExpand == 1
        adddim = estimateDkExpand(Vmat{s}, results.Vmat_sv{1,s}, para);
		if adddim > 0
        %% Apply the expansion
			[dk,dOBB]   = size(Vmat{s});
			para.dk(s)  = dk+adddim;						% operators will be expanded in genh1h2term?
			if para.useVmat
				addmat  = zeros(adddim,dOBB);
% 				Vmat{s} = cat(1,addmat,Vmat{s});			% N:1 ordering
				Vmat{s} = cat(1,Vmat{s},addmat);			% 1:N ordering
			else
				Vmat{s} = eye(para.dk(s));
			end
			para.hasexpanded = 1;
			para.increasedk  = 0;							% reset this value
% 			dispif('Increased dk',para.logging)
		end
    end
end

% Expand operators
[op,para] = genh1h2term(para,op);
op.h1jOBB = []; op.h2jOBB = [];			% invalidate OBB trafos for proper updateop
op.h1j = []; op.h2j = [];				% also for old on-site operators

disp('(para.d_opt-d_opt_old)./para.d_opt=')
disp(mat2str((para.d_opt-d_opt_old)./para.d_opt,2));
para.d_opt_change=mean(abs(para.d_opt-d_opt_old)./para.d_opt);

if para.logging													% save d_opt history
		results.d_opt{para.loop} =para.d_opt;
		results.dk{para.loop}    = para.dk;
end

end

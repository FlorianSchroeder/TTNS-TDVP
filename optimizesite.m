function [Amat_focused,Vmat,results,para,op]=optimizesite(mps,Vmat,op,para,results,sitej)
% Uses Vmat only on not-spin-sites.
% Uses shift only on not-spin-sites.
% If at spinsite, sitej = 1: Amatlaststep = mps{sitej}; [Aj, E] = minimizeE_onesiteA(op, Vmat{sitej}, Amatlaststep,para,sitej);
% TODO: Improve speed for optV loop by using QR.
%		SVD result needed in lines: 77-79, so perhaps do it before that in a single run only if useshift == 1?
%		    will produce overhead if shift not used as (QR+SVD), but great reduction with loops >=2
%		QR usable in Lines: 13, 22
%
% Commented by Florian Schroeder 13/01/2014
% Modified:
%   FS 26/05/2014:  Shift only if sitej not in para.spinposition --> bosonic site.

optV=1;
while optV
    % Use Vmat only for bosonic sites
    if  prod(sitej ~= para.spinposition) && para.useVmat==1                 % Only use Vmat{j} and optimize for the boson sites. Old 05/05/14: (para.dk(sitej)>2 && para.useVmat == 1); Now: ready for array in spinposition
        [Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);				% left-normalize A, as r -> l sweep did right normalize.

		Vmat_focused = Vmat{sitej} * V.';
        [Vmat_focused,E] = minimizeE_onesiteVmat(op, Amat, Vmat_focused,para);			% first Energy Optimization
%         if sitej>=3 && para.rescaling==1
%             results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.Lambda;
%             E=(E+results.geoffset(sitej))./(para.Lambda.^(sitej-2));
%         end
%         %results.leftge(sitej)
%         fprintf('\nE = %.10g\n', E);
        [Vmat{sitej}, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);
        Amat_focused = contracttensors(Amat, 3, 3, V, 2, 2);
        clear('Vmat_focused');
    else
        Amat_focused = mps{sitej}; optV=0;
    end
    [Amat_focused, E, op] = minimizeE_onesiteA(op, Vmat{sitej}, Amat_focused,para);

% Shifting basis of bosonic sites only. Different criteria explained in declaration file.
% prod(sitej ~= para.spinposition) excludes all spin-sites
    if para.loop > 1 && para.useshift == 1 && prod(sitej ~= para.spinposition) && para.nChains == 1          % old: para.dk(sitej)>2
		if para.useFloShift == 1 && ~(para.trustsite(end)>0)						% stop loops if criterium not fulfilled
			optV = 0;
%		elseif para.useFloShift2 == 1 && ~(min(cellfun(@(x) x(1,1), results.Vmat_sv(2:end))) < para.FloShift2minSV)	% cellfun takes all maximum SV of each site into an array
		elseif para.useFloShift2 == 1 && ~(min(results.Vmat_svLog{end}) < para.FloShift2minSV)	% results.Vmat_svLog{end} contains max SV from last sweep
			optV = 0;
%		elseif para.useFloShift3 == 1 && (~(mod(para.loop,para.FloShift3loops)==0) || ~(results.Vmat_svLog{end}(sitej-1) < (1.+para.FloShift3minMaxSV)/2.))		%FloShift3minMaxSV from last sweep. (sitej-1) as spinsite not included
		elseif para.useFloShift3 == 1 && (~(mod(para.loop,para.FloShift3loops)==0) || ~((results.Vmat_svLog{end}(sitej-1) < (2.+para.FloShift3minMaxSV)/3.) || (sitej-1 < find(results.Vmat_svLog{end} <= para.FloShift3minMaxSV,1) && results.Vmat_svLog{end}(sitej-1) ~= 1)))
		% Apply shift if: every nth loop; (site has had a smaller SV than criterium) or (site is below worst site and SV ~= 1)
			optV = 0;
		elseif para.useChengShift == 1 && ~(sitej<=para.trustsite(end))		% -- original condition from Cheng.
			optV = 0;
		elseif ~para.useEveryShift && ~para.useChengShift && ~para.useFloShift && ~para.useFloShift2 && ~para.useFloShift3	% wrong definition
			optV = 0;
			disp('Define a shift method!');
	% && results.Eerror < sqrt(para.precision)  ---  inserted energy criterium myself!
		else
			shift = getObservable({'1siteshift'}, Amat_focused, Vmat{sitej}, para); % calculate applicable shift

			para.relativeshift(sitej) = max(abs(shift-para.shift(:,sitej)))/para.maxshift(sitej);
			if  para.relativeshift(sitej)>para.relativeshiftprecision                       % if relevant shift, update and rerun optimization
				para.shift(:,sitej) = shift;
				para.shifted        = 1;
				op                  = update_sitej_h1h2(para,op,sitej);						% shift boson operators
				mps{sitej}          = Amat_focused;											% store Amat in mps to use in next loop
			else
				optV=0;
                % Amat will be returned and normalised + saved later!
			end
		end
    else
        optV=0;
    end
end

if sitej>=3 && para.rescaling==1
    results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.chain{1}.Lambda;
    E=(E+results.geoffset(sitej))./(para.chain{1}.Lambda.^(sitej-2));
end
%results.leftge(sitej)
%fprintf('E = %.10g\n', E);
results.E=E;

end

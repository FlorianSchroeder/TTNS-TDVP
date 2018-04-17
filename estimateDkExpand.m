function NDimExpand = estimateDkExpand(V, sv, para)
%% NDimExpand = estimateDkExpand(Vmat, sv, para)
%	estimates by which amount Dk should be expanded
%
%	para.useDkExpand1: expand by 40% if largest SV < para.expandBelowSV
%	para.useDkExpand2: Tail analysis, expand if Vmat tail indicates need for more Dk
%					   More advanced technique than DkExpand1!
%
%	Vmat: single-site V-OBB matrix
%	sv  : corresponding SV-vector
%
%	Created by FS 21/08/2015
%

[dk,dOBB] = size(V);
dOBB = min(dOBB,length(sv));		% in case V was recently expanded/truncated

if para.useDkExpand1 == 1 && para.useVmat && max(sv) < para.expandBelowSV
	% Expand dk if highest SV < para.expandBelowSV
	dimincratio = 0.4;							% find a better ratio from SV
elseif para.useDkExpand2 == 1 && para.useVmat
	% Expand dk according to Vmat SV-tail analysis: Much better than DkExpand1
	% vector of length dk:
	% calculate the sum over the weighted Vmat contributions. This tells which states are occupied / needed

% 	occupation = sum(abs(real(V(:,1:dOBB)*diag(sv(1:dOBB)))),2);		% Old! real should be wrong!
	occupation = sum(abs( V(:,1:dOBB)*diag(sv(1:dOBB)) ),2);			% sum over OBB

% 	occTail = occupation(1:ceil(para.dkEx2_tail*length(occupation)));									% high energy tail N:1 ordering
	occTail = occupation(ceil((1-para.dkEx2_tail)*length(occupation)):end);								% high energy tail 1:N ordering
	if std(log10(occTail)) > para.dkEx2_maxDev && abs(mean(log10(occTail))) < para.svmintol^2			% measure order of fluctuations on tail. should be small for flat tail.
		dimincratio = 1/mean(abs(diff(log10(occTail))));												% the steeper the less increase needed.
		%dimincratio = log10(std(occTail)/para.dkEx2_maxDev)/3;											% similar to imperial relation, slightly modified (estimated). always positive.
	elseif mean(occTail) > para.svmintol																% if absolute value is too high
		dimincratio = 1/mean(abs(diff(log10(occTail))));												% some increase
	else
		dimincratio = 0;                        % no increase
	end
else
	dimincratio = 0;
end
dimincratio = max(dimincratio,para.increasedk/dk);      % if d_opt needs more dk, this assures it.

if dimincratio < 0
	NDimExpand = 0;
	return;
end;

NDimExpand = min(ceil(dimincratio*dk),para.dkmax-dk);
if para.foldedChain
	% ensure that sqrt(dk) = integer
	NDimExpand = ceil(sqrt(dk+NDimExpand))^2 - dk;
end


end
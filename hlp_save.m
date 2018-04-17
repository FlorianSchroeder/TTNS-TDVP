function hlp_save(fname,para,op,results,tresults,mps,Vmat,treeMPS)
% saves above variables serialized into a file.

para     = hlp_serialize(para);
save(fname,'para','-v7.3');
fprintf('.');

if ~isempty(op)
	op       = hlp_serialize(op);
	save(fname,'op','-append');
	clear('op');
	fprintf('.');
end

results  = hlp_serialize(results);
save(fname,'results','-append');
clear('results');
fprintf('.');

tresults = hlp_serialize(tresults);
save(fname,'tresults','-append');
clear('tresults');
fprintf('.');

if ~isempty(mps)
	mps      = hlp_serialize(mps);
	save(fname,'mps','-append');
	clear('mps');
	fprintf('.');
end

if ~isempty(Vmat)
	Vmat     = hlp_serialize(Vmat);
	save(fname,'Vmat','-append');
	clear('Vmat');
	fprintf('.\n');
end

if nargin > 7 && ~isempty(treeMPS)
	treeMPS  = hlp_serialize(treeMPS);
	save(fname,'treeMPS','-append');
	clear('treeMPS');
	fprintf('.');
end
end
function r = initresults(para)
%
% Modified:
%   FS 25/10/2014: - added line for ~para.useVmat to prepare Vmat_sv == 1
%	FS 12/10/2015: - better multi chain support: Vtens and StarMPS
L     = para.L;
d_opt = para.d_opt;
if L ~= 1
	D     = para.D(end);
end
NC    = para.nChains;

r.Vmat_sv  = cell(1,L);
r.Amat_sv  = cell(1,L-1);
r.Vmat_vNE = zeros(1,L); %phonon basis von Neumann Amat_vNEj
r.Amat_vNE = zeros(1,L-1);
r.Vmat_truncErr = zeros(1,L);
r.Amat_truncErr = zeros(1,L-1);

if para.parity~='n'
    % old code:
    %r.Vmat_sv=zeros(L,max(d_opt)/2);
    %r.Amat_sv=zeros(L,D/2);
else
	if para.useVtens
		r.Vmat_sv = cell(NC+1,L);
	elseif para.useStarMPS
		r.Vmat_sv  = cell(NC,L);
		r.Vmat_vNE = zeros(NC,L);
		r.Amat_sv  = cell(NC,L-1);
		r.Amat_vNE = zeros(NC,L-1);
	else
		r.Vmat_sv = cell(1,L);
	end
    %r.Vmat_sv=zeros(L,max(d_opt));
    %r.Amat_sv=zeros(L,D);
end

if ~para.useVmat
    % initialise with any dummy values to prevent error
    for i = 1:L
        r.Vmat_sv{i} = eye(para.dk(i));
    end
end

r.leftge = zeros(L,1);
r.geoffset=zeros(L,1);
r.Eerror(1)=1;
r.lastVmat_vNE=zeros(1,L);
if para.logging
    r.d_opt = cell(1);
    r.d_opt{1} = para.d_opt;
    r.D = cell(1);
    r.D{1} = para.D;
	r.dk = cell(1);
    r.dk{1} = para.dk;
    r.shift = cell(1);
    r.shift{1} = para.shift;
    r.Vmat_svLog = cell(1);
	r.EvaluesLog = cell(1);
	r.flowdiag = cell(1);
end
function out = performanceTests(variant, varargin)
%% function performanceTests(variant, varargin)
%	runs certain benchmarks which are useful to assess performance for TDVP-VMPS
out = [];
switch lower(variant)
	case 'matfile'
		testMatFile();
	case 'ncon'
		testNCON();
	case 'cellfun'
		testCellfun();
	case 'matcontract1'
		testGPUMatContract1();
	case 'tenscontract1'
		testGPUTensContract1();
	case 'tenscontract2'
		out = testGPUTensContract2();
	case 'rrqr'
		out = testRandomizedQR();
	case 'rrqrsweep'
		out = testRandomizedQRsweep();
	case 'cellstructmps'
		out = testCellStructMPS();
	case 'copyonwrite'
		out = testCopyOnWrite();
end

end

function testMatFile()

%% Test performance of matfile command
timings = zeros(1,3);
for n = 10.^(1:3)
	system(['del "test.mat"'] );
	out = matfile('test.mat','Writable',true);
	start = tic;
	out.mps(n,n) = {[]};
	timings(log10(n)) = toc(start);
end
%%
	L = 100;
	N = [5,5,20];
	T = 1000;
    filespec = 'matfile_test.mat';
	%%
    mat = rand( N, 'double' );
	mps = cell(T,L);
	for i = 1:T
		for j = 1:L
			mps{i,j} = mat;
		end
	end
    save( filespec, 'mps','N','-v7.3' )
	clear mat mps
    %%
    obj = matfile( filespec ,'Writable', true );
    %%
    tic, mfm = obj.mps; toc
    tic, h5m = h5read( filespec, '/mps' ); toc
    %%
    dfm  = mfm-mat;
    d5m  = h5m-mat;
    max(abs(dfm(:)))
    max(abs(d5m(:)))
    %% column wise load
    tic, mfm = obj.mps( :, 1 ); toc
    tic, h5m = h5read( filespec, '/mps', [1,1], [N,1] ); toc
    %%
    dfm  = mfm-mat( :, 1 );
    d5m  = h5m-mat( :, 1 );
    max(abs(dfm(:)))
    max(abs(d5m(:)))
    %% row wise load
    tic, mfm = obj.mps( 2, : ); toc
    tic, h5m = h5read( filespec, '/mps', [2,1], [1,N] ); toc
	max(abs(mfm-h5m))
    %%
    dfm  = mfm-mat( 1, : );
    d5m  = h5m-mat( 1, : );
    max(abs(dfm(:)))
    max(abs(d5m(:)))
	%% writing test

	tic, obj.mat( 1, : ) = mfm; toc
%     tic, h5write( filespec, '/mat', h5m, [2,1], [1,N] ); toc

end

function testNCON()

%%  Benchmark updateCright vs ncon
tNCON = []; tUCR = []; step = 10;
for j = 4
	D = j*step; dOBB = 50; dk = 1000;
	mps = randn(D,D,dOBB);
	Vmat = randn(dk,dOBB);
	H = randn(dk,dk);
	C = randn(D,D);
	rep = 100;

	t1 = tic;
	for i = 1:rep
		cNCON =		ncon({mps,		conj(mps), Vmat,  conj(Vmat), H,	 C},...
			 {[-2,6,4], [-1,5,2],  [3,4], [1,2],	  [1,3], [5,6]});
	end
	tNCON(j) = toc(t1)

	t1 = tic;
	for i = 1:rep
		cUCR = updateCright(C, mps, Vmat, H, mps, Vmat);
	end

	tUCR(j) = toc(t1)
end

%% Plot results:
plot((1:5) .*step,[tUCR;tNCON]);

end

function testCellfun()
%% find(cellfun('isempty',mcOp)) --> fast enough!
rep = 1000;
mcOp = cell(10,1);
mcOp{4} = eye(10);
tic
for i = 1:rep
	ind = zeros(1,length(mcOp));
	ind = find(cellfun('isempty',mcOp));
end
toc

tic
for i = 1:rep
	ind = [];
	for j = 1:length(mcOp)
		if isempty(mcOp{j})
			ind = [ind;j];
		end
	end
end
toc

end

function testGPUMatContract1()
%% Test CPU vs GPU performance (matrix)
Nruns = 10;
dims = 10:20:5000;
tCPU = zeros(length(dims),Nruns);
tGPU = zeros(length(dims),Nruns);
gd = gpuDevice;

for ii = 1:length(dims)
	A = randn(dims(ii));
	B = randn(dims(ii));
	C = gpuArray(A);
	D = gpuArray(B);

	for jj = 1:Nruns
		t = tic;
		AB = A*B;
		tCPU(ii,jj) = toc(t);

		t = tic;
		CD = C*D;
		wait(gd);
		tGPU(ii,jj) = toc(t);
	end

	fprintf('D = %d, Average CPU time: %g, STD: %g\n', dims(ii), mean(tCPU(ii,:)), std(tCPU(ii,:)));
	fprintf('D = %d, Average GPU time: %g, STD: %g\n', dims(ii), mean(tGPU(ii,:)), std(tGPU(ii,:)));

end
end

function testGPUTensContract1()
%% Test CPU vs GPU performance (tensor)
Nruns = 10;
M = 6000;
N = 100;
P = 20;
K = 100;
dims = 10:20:K;
tCPU = zeros(length(dims),Nruns,4);
tGPU = zeros(length(dims),Nruns,4);
gd = gpuDevice;


for ii = 1:length(dims)
	A = rand(M,dims(ii),P);		% M K P
	B = rand(dims(ii),N);		% K N
	C = gpuArray(A);
	D = gpuArray(B);
	CD2 = gpuArray.zeros(M,N,P);

	for jj = 1:Nruns
		t = tic;
		AB = permute(reshape(reshape(permute(A,[1,3,2]),[],dims(ii))*B, M,P,N),[1,3,2]);
		tCPU(ii,jj,1) = toc(t);

		t = tic;
		ABcell = arrayfun(@(kk) A(:,:,kk) * B, 1:P,'UniformOutput',false);
		tCPU(ii,jj,2) = toc(t);

		t = tic;
		for kk = 1:P
			AB(:,:,kk) = A(:,:,kk) * B;
		end
		tCPU(ii,jj,3) = toc(t);

		% CONTRACTTENSORS
		t = tic;
		AB = permute(contracttensors(A,3,2,B,2,1),[1,3,2]);
		tCPU(ii,jj,4) = toc(t);

		% PAGEFUN
		t = tic;
		CD = pagefun(@mtimes, C,D);
		wait(gd);
		tGPU(ii,jj,1) = toc(t);

		% PERMUTE
		t = tic;
		CD = permute(reshape(reshape(permute(C,[1,3,2]),[],dims(ii))*D, M,P,N),[1,3,2]);
		wait(gd);
		tGPU(ii,jj,2) = toc(t);

		% FOR-LOOP
		t = tic;
		for kk = 1:P
			CD2(:,:,kk) = C(:,:,kk) * D;
		end
		wait(gd);
		tGPU(ii,jj,3) = toc(t);
		%parfor was pretty bad!

		% ARRAYFUN
		t = tic;
		CDcell = arrayfun(@(kk) C(:,:,kk) * D, 1:size(C,3),'UniformOutput',false);

		wait(gd);
		tGPU(ii,jj,4) = toc(t);

		% CONTRACTTENSORS
		t = tic;
		CD = permute(contracttensors(C,3,2,D,2,1),[1,3,2]);
		wait(gd);
		tGPU(ii,jj,5) = toc(t);
	end

	fprintf('D = %d, Average CPU time (permute): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,1)), std(tCPU(ii,:,1)));
	fprintf('D = %d, Average CPU time (arrayfu): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,2)), std(tCPU(ii,:,2)));
	fprintf('D = %d, Average CPU time (forloop): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,3)), std(tCPU(ii,:,3)));
	fprintf('D = %d, Average CPU time (contens): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,4)), std(tCPU(ii,:,4)));
	fprintf('D = %d, Average GPU time (pagefun): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,1)), std(tGPU(ii,:,1)));
	fprintf('D = %d, Average GPU time (permute): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,2)), std(tGPU(ii,:,2)));
	fprintf('D = %d, Average GPU time (forloop): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,3)), std(tGPU(ii,:,3)));
	fprintf('D = %d, Average GPU time (arrayfu): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,4)), std(tGPU(ii,:,4)));
	fprintf('D = %d, Average GPU time (contens): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,5)), std(tGPU(ii,:,5)));
	fprintf('\n');
end
	%%	Plot results
	f = figure(); f.Name = 'GPUbench'; ax=gca; hold all
	pl = {};
	pl{1} = plot(dims,mean(tCPU(:,:,1),2),'DisplayName','CPU, permute');
	pl{2} = plot(dims,mean(tCPU(:,:,2),2),'DisplayName','CPU, arrayfun');
	pl{3} = plot(dims,mean(tCPU(:,:,3),2),'DisplayName','CPU, for-loop');
	pl{4} = plot(dims,mean(tCPU(:,:,4),2),'DisplayName','CPU, con-tens');
	ax.ColorOrderIndex = 1;
	pl{5} = plot(dims,mean(tGPU(:,:,1),2),'-.','DisplayName','GPU, pagefun');
	pl{6} = plot(dims,mean(tGPU(:,:,2),2),'-.','DisplayName','GPU, permute');
	pl{7} = plot(dims,mean(tGPU(:,:,3),2),'-.','DisplayName','GPU, for-loop');
	pl{8} = plot(dims,mean(tGPU(:,:,4),2),'-.','DisplayName','GPU, arrayfun');
	pl{9} = plot(dims,mean(tGPU(:,:,5),2),'-.','DisplayName','GPU, con-tens');
	legend show
	xlabel('Array Dim');
	ylabel('Time in s');
	ax.Color = 'None';
	M*N*K*P

end

function out = testGPUTensContract2()
%% Test CPU vs GPU performance (tensor)
% test real StarmultA testcase
Nruns = 1000;
D = 20;		% max Dimension
dk = 4;		% local dimension
NC = 5;		% number of chains
d = [1, ones(1,NC)*D, dk];		% dimensions of A
dims = 1:D;
tCPU = nan(length(dims),Nruns);
tGPU = nan(length(dims),Nruns);
gd = gpuDevice;

for ii = 1:length(dims)
	d = [1, ones(1,NC)*dims(ii), dk];		% dimensions of A
	A = rand(d);
	op.h1jOBB = randn(dk);					% local single-site operator
	op.Hright = randn(dims(ii));			% right chain op
	op.Opright = randn(dims(ii));			% right chain op

% 	Agpu = gpuArray(A);

% 	CD2 = gpuArray.zeros(M,N,P);
	% CPU tensShape code
	jj = 0;
% 	while sum(tCPU(ii,:,1),'omitnan') < 20
	while (jj > 1 && (std(tCPU(ii,1:jj,1))/mean(tCPU(ii,1:jj,1)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% 1. on-site H1
		w =	contracttensors(A, NC+2, NC+2, op.h1jOBB.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(A, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hright * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Opright * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBB.',2,1);
			end
		end
		tCPU(ii,jj,1) = toc(t);
		if jj == Nruns
			break;
		end
	end

	fprintf('D = %d, Average CPU time (tensShape): %g\n', dims(ii), mean(tCPU(ii,1:jj,1),'omitnan'));

	% GPU tensShape code (excluding transfer time)
	% copy all onto GPU
	Ag = gpuArray(A);
	op.h1jOBBg = gpuArray(op.h1jOBB);					% local single-site operator
	op.Hrightg = gpuArray(op.Hright);			% right chain op
	op.Oprightg = gpuArray(op.Opright);

	jj = 0;
% 	while sum(tCPU(ii,:,1),'omitnan') < 20
	while (jj > 1 && (std(tGPU(ii,1:jj,1))/mean(tGPU(ii,1:jj,1)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% 1. on-site H1
		w =	contracttensors(Ag, NC+2, NC+2, op.h1jOBBg.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(Ag, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hrightg * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Oprightg * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBBg.',2,1);
			end
		end
		wait(gd);
		tGPU(ii,jj,1) = toc(t);
		if jj == Nruns
			break;
		end
	end

	fprintf('D = %d, Average GPU time (tShape-ex): %g\n', dims(ii), mean(tGPU(ii,1:jj,1),'omitnan'));

	% GPU tensShape code (including transfer time)
	jj = 0;
% 	while sum(tGPU(ii,:,2),'omitnan') < 20
	while (jj > 1 && (std(tGPU(ii,1:jj,2))/mean(tGPU(ii,1:jj,2)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% copy all onto GPU

		Ag = gpuArray(A);
		op.h1jOBBg = gpuArray(op.h1jOBB);					% local single-site operator
		op.Hrightg = gpuArray(op.Hright);			% right chain op
		op.Oprightg = gpuArray(op.Opright);

		% 1. on-site H1
		w =	contracttensors(Ag, NC+2, NC+2, op.h1jOBBg.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(Ag, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hrightg * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Oprightg * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBBg.',2,1);
			end
		end
		wait(gd);
		tGPU(ii,jj,2) = toc(t);
		if jj == Nruns
			break;
		end
		reset(gd);
	end

	fprintf('D = %d, Average GPU time (tShape-in): %g\n', dims(ii), mean(tGPU(ii,1:jj,2),'omitnan'));
	fprintf('\n');
end

	%%	Plot results
	f = figure(); f.Name = 'GPUbench'; ax=gca; hold all
	pl = {};
	pl{1} = plot(dims,mean(tCPU(:,:,1),2,'omitnan'),'DisplayName','CPU, tensShape');
% 	pl{2} = plot(dims,mean(tCPU(:,:,2),2),'DisplayName','CPU, arrayfun');
% 	pl{3} = plot(dims,mean(tCPU(:,:,3),2),'DisplayName','CPU, for-loop');
% 	pl{4} = plot(dims,mean(tCPU(:,:,4),2),'DisplayName','CPU, con-tens');
% 	ax.ColorOrderIndex = 1;
	pl{2} = plot(dims,mean(tGPU(:,:,1),2,'omitnan'),'-.','DisplayName','GPU, tensShape');
	pl{3} = plot(dims,mean(tGPU(:,:,2),2,'omitnan'),'-.','DisplayName','GPU, tensShape+transfer');
% 	pl{7} = plot(dims,mean(tGPU(:,:,3),2),'-.','DisplayName','GPU, for-loop');
% 	pl{8} = plot(dims,mean(tGPU(:,:,4),2),'-.','DisplayName','GPU, arrayfun');
% 	pl{9} = plot(dims,mean(tGPU(:,:,5),2),'-.','DisplayName','GPU, con-tens');
% 	legend show
	xlabel('Array Dim');
	ylabel('Time in s');
	ax.Color = 'None';
	set(gca,'YScale','log');

	out.tCPU = tCPU;
	out.tGPU = tGPU;
end

function out = testRandomizedQR(N,r)
% function testRandomizedQR()
%   call with performanceTests('rrQR')
%	tests the randomized reduced rank QR rrQR()

out = zeros(2);

if nargin == 0
	N = 1000;	% dimension of A
	r = 200;		% rank of A accurate up to eps. Only used to construct A
end
Nruns = 10;

% 1. create a randomized A
A = randn(N);
[U,S,V] = svd(A);
S=diag(exp(-30/r*(1:N)));		% make S such that for i>r S_i < eps. log(eps) = -36, but 30 gives the right rank(A) = r

A = U*S*V';		% now A is low rank!
rank(A)
t = zeros(Nruns,1);			% keep the timings for each run
for ii = 1:Nruns
	tic
	[Q,R] = qr(A,0);
	t(ii) = toc;
end
out(1,1) = mean(t);			% mean time per calculation
out(1,2) = norm(A-Q*R);		% error estimate

t = zeros(Nruns,1);
for ii = 1:Nruns
	tic
	[Q,R] = rrQR(A,floor(r*0.4),0);
	t(ii) = toc;
end
out(2,1) = mean(t);			% mean time per calculation
out(2,2) = norm(A-Q*R);		% error estimate

% compare SV:
B = Q'*A;
[u,s,v] = svd(B,0);
plot(log10(abs(diag(S(1:size(s,1),1:size(s,1)))-diag(s))));

end

function out = testRandomizedQRsweep()
%  function out = testRandomizedQRsweep()
% sweeps through parameter space!
% different N and different r
step = 10;
N = 10:step:1000;

tQR  = zeros(length(N));	% timing of QR
eQR  = tQR;					% error
trrQR = tQR;				% timing of rrQR
errQR = tQR;				% error of rrQR

for ii = 1:length(N)
	r = 10:step:N(ii);
	for jj = 1:length(r)
		A = testRandomizedQR(N(ii),r(jj));
		tQR(ii,jj) = A(1,1);
		eQR(ii,jj) = A(1,2);
		trrQR(ii,jj) = A(2,1);
		errQR(ii,jj) = A(2,2);
		fprintf('finished N = %d, r = %d\n', N(ii),r(jj));
	end
	figure(1);
	surf(tQR)
	figure(2);
	surf(eQR);

	figure(3);
	surf(trrQR)
	figure(4);
	surf(errQR);
	if mod(ii,10)
		drawnow
	end
end
% plot all results
figure(1);
surf(tQR)
figure(2);
surf(eQR);

figure(3);
surf(trrQR)
figure(4);
surf(errQR);

end

function out = testCellStructMPS()
%% function out = testCellStructMPS()
%	tests whether a struct MPS has any speed disadvantages compared to a cell-array MPS
%	useful to assess suitability of struct MPS for TreeMPS implementation.

% define test case as StarMPS
% do star(center + 2 chains + star(center + 2 chains))
Lchains = 5;
D = 30;			% one BondDim for all edges
Nruns = 300;
tCell = zeros(Nruns,1);
tStruct = zeros(Nruns,1);

% Define ChainMPS
chainMPS = cell(1,Lchains);
for ii = 1:Lchains
	if ii == Lchains
		chainMPS{ii} = randn(D,1,D);			% end-of-chain -> Dr = 1
	else
		chainMPS{ii} = randn(D,D,D);
	end
end
% Define inner StarMPS
starMPS = cell(3,1);	% {center, chain1, chain2}
starMPS{1} = randn(D,D,D,D);	% Dl, Dc1, Dc2, dk
starMPS{2} = chainMPS;			% copy chainMPS
starMPS{3} = chainMPS;			% copy chainMPS

% Build Cell MPS
cellMPS = [starMPS;{starMPS}];
cellMPS{1} = reshape(cellMPS{1},[1,size(cellMPS{1})]);		% include Dl = 1

% Build Struct MPS
structChainMPS = struct();		% smallest subunit
structChainMPS.mps = chainMPS;

structStarMPS = struct();
structStarMPS.mps = {starMPS{1};structChainMPS;structChainMPS};

structMPS = struct();
structMPS.mps = {cellMPS{1}; structChainMPS; structChainMPS; structStarMPS};

%% Get Timings for contraction calls
for ii = 1:Nruns
	tic
	calcCrightCell(cellMPS);
	timing = toc;
	tCell(ii,1) = timing;
	
	tic
	calcCrightStruct(structMPS);
	timing = toc;
	tStruct(ii,1) = timing;
end

out = [tCell,tStruct];

plot(out);
legend('Cell','Struct')
end

function out = testCopyOnWrite()
%% function out = testCopyOnWrite()
%	tests whether an assignment followed by read only impairs speed

% define test case
% do star(center + 2 chains + star(center + 2 chains))
NC = 3;
D = 200;			% one BondDim for all edges
Nruns = 30;
tDirect = zeros(Nruns,1);
tRead = zeros(Nruns,1);
tWrite = zeros(Nruns,1);

MPS = randn(ones(1,NC)*D);
op = randn(D,D);

%% Get Timings for contraction calls
for ii = 1:Nruns
	% Direct access
	tic
	B = contracttensors(op,2,2, MPS,NC,1);
	timing = toc;
	tDirect(ii,1) = timing;
	
	% Read only access
	tic
	Amat = MPS;				% pass on handle
	B = contracttensors(op,2,2, MPS,NC,1);
	timing = toc;
	tRead(ii,1) = timing;
	
	% write inbetween
	tic
	Amat = MPS;				% pass on handle
	Amat(1) = 2;			% arbitrary write access
	B = contracttensors(op,2,2, Amat,NC,1);
	timing = toc;
	tWrite(ii,1) = timing;
end

out = [tDirect,tRead,tWrite];

plot(out);
legend('Direct','Read','Write')
end

function out = calcCright(MPS)
%% calculates <MPS|MPS> into the left effective basis
if isstruct(MPS)
	% code for Struct MPS
	MPS = MPS.mps;
else
	% code for Cell MPS
end

if size(MPS,1) ~= 1
	%% get the Cright for each chain and contract with center
	NC = size(MPS,1)-1;			% first element is Tree root
% 	Cright = cell(NC,1);
	Atens = shiftdim(MPS{1},1);		% already shift leading singleton to the end.
	nd = ndims(MPS{1});
	for ii = 2:NC+1
		Cright = calcCright(MPS{ii});
		Atens = contracttensors(Atens,nd,1,Cright.',2,1);
	end
	Atens = shiftdim(Atens,1);		% shift dk to the end
	
	out = contracttensors(conj(MPS{1}),nd,2:nd,Atens,nd,2:nd);
else
	%% iterate through chain and updateCright()
	C = [];
	L = size(MPS,2);
	for ii = 1:L
		C = updateCright(C,MPS{L-ii+1},[],[],MPS{L-ii+1},[]);
	end
	out = C;
end

end

function out = calcCrightCell(MPS)
%% calculates <MPS|MPS> into the left effective basis
%	for Cell MPS
if size(MPS,1) ~= 1
	%% get the Cright for each chain and contract with center
	NC = size(MPS,1)-1;			% first element is Tree root
% 	Cright = cell(NC,1);
	Atens = shiftdim(MPS{1},1);		% already shift leading singleton to the end.
	nd = ndims(MPS{1});
	for ii = 2:NC+1
		Cright = calcCrightCell(MPS{ii});
		Atens = contracttensors(Atens,nd,1,Cright.',2,1);
	end
	Atens = shiftdim(Atens,1);		% shift dk to the end
	
	out = contracttensors(conj(MPS{1}),nd,2:nd,Atens,nd,2:nd);
else
	%% iterate through chain and updateCright()
	C = [];
	L = size(MPS,2);
	for ii = 1:L
		C = updateCright(C,MPS{L-ii+1},[],[],MPS{L-ii+1},[]);
	end
	out = C;
end

end

function out = calcCrightStruct(MPS)
%% calculates <MPS|MPS> into the left effective basis
%	for Struct MPS
if size(MPS.mps,1) ~= 1
	%% get the Cright for each chain and contract with center
	NC = size(MPS.mps,1)-1;			% first element is Tree root
% 	Cright = cell(NC,1);
	Atens = shiftdim(MPS.mps{1},1);		% already shift leading singleton to the end.
	nd = ndims(MPS.mps{1});
	for ii = 2:NC+1
		Cright = calcCrightStruct(MPS.mps{ii});
		Atens = contracttensors(Atens,nd,1,Cright.',2,1);
	end
	Atens = shiftdim(Atens,1);		% shift dk to the end
	
	out = contracttensors(conj(MPS.mps{1}),nd,2:nd,Atens,nd,2:nd);
else
	%% iterate through chain and updateCright()
	C = [];
	L = size(MPS.mps,2);
	for ii = 1:L
		C = updateCright(C,MPS.mps{L-ii+1},[],[],MPS.mps{L-ii+1},[]);
	end
	out = C;
end

end

function [outMPS] = moveFocus(mps,i,j)
	% Moves the focus from site i to site j
	% i needs to indicate the current state
	if i > j
		error('VMPS:tdvp_1site:NotYetImplemented','Moving Focus from right to left still needs to be implemented'),
	end
	outMPS = mps;
	for k = i:j-1
		d = size(outMPS{k});
		A = permute(outMPS{k}, [3, 1, 2]);
		A = reshape(A, d(3)*d(1), []);		% reshapes to (a1 d),a2
		[A, S, U] = svd2(A);
		A = reshape(A, d(3), d(1), []);
		outMPS{k} = permute(A, [2, 3, 1]);
		U = S * U;
		outMPS{k+1} = contracttensors(U,2,2, outMPS{k+1},3,1);
	end
end
%% How to replace tensor contraction with bsxfun

%%
para.chain{1}.mapping			= 'Stieltjes';
% para.chain{1}.mapping			= 'LanczosTriDiag';
% J given though points
para.chain{1}.spectralDensity	= 'CoupBroad'; para.chain{1}.w_cutoff = 1;
para.chain{1}.dataPoints        = [0.2:0.1:1;abs(sin(2*(0.2:0.1:1)))]';          % [w, J(w)]
para.chain{1}.peakWidth			= 0.005;

% para.chain{1}.spectralDensity	= 'PointsInterp'; para.chain{1}.w_cutoff = 1;
% para.chain{1}.dataPoints        = [0:0.1:1;abs(sin(6*(0:0.1:1)))]';          % [w, J(w)]
% % para.chain{1}.peakWidth			= 0.004;

% lambda given through points
para.chain{1}.discrMethod		= 'Direct';
para.chain{1}.discretization	= 'Linear';

para.chain{1}.s		 = 1;                            % SBM spectral function power law behaviour
para.chain{1}.alpha  = 0.2;                       % SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
para.chain{1}.Lambda = 1;
para.chain{1}.z      = 1;
para.chain{1}.L		 = 40;
%%
a=[0:0.1:1;abs(sin(6*(0:0.1:1)))]';          % [w,
b = [0.25, 2; 0.3, 3];
w = 0:0.0001:1;
sigma = 0.010;
figure(3); clf; hold all;
stem(b(:,1),b(:,2));
% normpdf = @(X,mu,sigma) exp(-(X-mu).^2 ./2 ./sigma.^2)./sigma ./sqrt(2*pi);			% normal dist
normpdf = @(X,mu,gamma) 1./(pi .* gamma .* (1+((X-mu)./gamma ).^2));				% cauchy -> not really good
y = 0;
for ii = 1:size(b,1)
	A = 1/(b(ii,1) * integral(@(w) normpdf(w,b(ii,1),sigma)./w, -inf,inf));
	y = y + A .* b(ii,2) .* normpdf(w,b(ii,1),sigma);			% good way!
% 	y = y + sqrt(A .* b(ii,2) .* normpdf(w,b(ii,1),sigma));		% bad way!
end
% A = arrayfun(@(mu) 1/(mu*integral(@(w) normpdf(w,mu,0.5)./w, -inf,inf)), b(:,1));
% y = A(1).*b(1,2).*normpdf(w,b(1,1),0.5)+A(2).*b(2,2).* normpdf(w,b(2,1),0.5);
% plot(w, y)
plot(w,sqrt(y))
%%
para.chain{1} = SBM_genpara(para.chain{1});
%%
para.chain{2} = para.chain{1};
% para.chain{2}.mapping = 'Stieltjes';
para.chain{2}.mapping = 'LanczosTriDiag';
% para.chain{2}.discretization = 'Linear'; para.chain{2}.Lambda = 1;
para.chain{2} = SBM_genpara(para.chain{2});
%%
f = figure(2); clf; hold on;
col = {'k','r','b','g'};
for i = 1:length(para.chain)
	plot(para.chain{i}.epsilon,col{i});
	plot(para.chain{i}.t,col{i});
end
% set(gca,'yscale','log');

%% test accuracy of Direct Stieltjes vs bigL
% need to adjust bigL for each chain by hand!
para.chain{1}.mapping			= 'Stieltjes';
para.chain{1}.spectralDensity	= 'Leggett_Soft';
para.chain{1}.discrMethod		= 'Direct';
para.chain{1}.discretization	= 'Linear';

para.chain{1}.s		 = 1;                            % SBM spectral function power law behaviour
para.chain{1}.alpha  = 0.3;                       % SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
para.chain{1}.Lambda = 1;
para.chain{1}.z      = 1;
para.chain{1}.L		 = 100;

%%
para.chain{1} = SBM_genpara(para.chain{1});			% bigL = 10*L
%%
para.chain{2} = para.chain{1};
para.chain{2} = SBM_genpara(para.chain{2});			% bigL = 100*L

%%
para.chain{3} = para.chain{1};
para.chain{3} = SBM_genpara(para.chain{3});			% bigL = 1000*L

%%
para.chain{4} = para.chain{1};
para.chain{4} = SBM_genpara(para.chain{4});			% bigL = 10000*L

%%
para.chain{5} = para.chain{1};
para.chain{5}.discretization = 'None';
para.chain{5}.mapping = 'OrthogonalPolynomials';
para.chain{5} = SBM_genpara(para.chain{5});			% bigL = 100000*L

%%
f = figure(2); clf; hold on;
col = {'k','r','b','g','k-'};
for i = 1:length(para.chain)-1
	p(i) = plot(abs(para.chain{end}.epsilon-para.chain{i}.epsilon),col{i});
	plot(abs(para.chain{end}.t-para.chain{i}.t),col{i});
% 	plot(para.chain{i}.epsilon,col{i});
% 	plot(para.chain{i}.t,col{i});
end
set(gca,'yscale','log');
legend(p,'10','100','1000','10000');

figure(3); clf; hold on;
errorBounds = zeros(length(para.chain),2);
for i = 1:length(para.chain)
	errorBounds(i,1) = std(para.chain{end}.epsilon-para.chain{i}.epsilon);
	errorBounds(i,2) = std(para.chain{end}.t-para.chain{i}.t);
end
plot(10.^(1:4),errorBounds);
set(gca,'xscale','log');set(gca,'yscale','log');

%% Find relation of w to DiscrMode input
% para.chain{1}.mapping			= 'Stieltjes';
para.chain{1}.mapping			= 'LanczosTriDiag';
% J given though points
para.chain{1}.spectralDensity	= 'CoupDiscr'; para.chain{1}.w_cutoff = 1;
para.chain{1}.dataPoints        = [0.2,0.4,0.6;0.2 0.5 1]';          % [w, lambda]
para.chain{1}.peakWidth			= 0;

% para.chain{1}.spectralDensity	= 'PointsInterp'; para.chain{1}.w_cutoff = 1;
% para.chain{1}.dataPoints        = [0:0.1:1;abs(sin(6*(0:0.1:1)))]';          % [w, J(w)]
% % para.chain{1}.peakWidth			= 0.004;

% lambda given through points
para.chain{1}.discrMethod		= 'Direct';
para.chain{1}.discretization	= 'Linear';

para.chain{1}.s		 = 1;                            % SBM spectral function power law behaviour
para.chain{1}.alpha  = 0.2;                       % SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
para.chain{1}.Lambda = 1;
para.chain{1}.z      = 1;
para.chain{1}.L		 = 3;

para.chain{1} = SBM_genpara(para.chain{1});
%% ncon test
	tNCON = []; tUCR = [];
for j = 50:50:1000
	D = 50; dOBB = 50; dk = j;
	mps = randn(D,D,dOBB);
	Vmat = randn(dk,dOBB);
	H = randn(dk,dk);
	C = randn(D,D);
	L = 10;

	t1 = tic;
	for i = 1:L
	% 	cNCON =		ncon({mps,		conj(mps), Vmat,  conj(Vmat), H,	 C},...
	% 					 {[-2,6,4], [-1,5,2],  [3,4], [1,2],	  [1,3], [5,6]});
		newH = (Vmat' * H) * Vmat;
		cNCON =		ncon({mps,		conj(mps), newH,  C},...
						 {[-2,4,2], [-1,3,1],  [1,2], [3,4]});

	end
	tNCON(j/50) = toc(t1)

	t1 = tic;
	for i = 1:L
		cUCR = updateCright(C, mps, Vmat, H, mps, Vmat);
	end
	tUCR(j/50) = toc(t1)
end
%
figure(1); clf; hold on;
plot(50:50:1000,tUCR);
plot(50:50:1000,tNCON);
% set(gca,'yscale','log');
xlabel('d_{k}'); ylabel('t in s');
legend('VMPS','NCON');

%% OBB multi-chain contraction test
for j = 1:5
	dOBB = 10; dk = j*2;
	nChains = 4;
	H = cell(1,nChains);
	H{3} = randn(dk,dk);
	Vmat = randn(dk^nChains,dOBB);
	L = 10;

	t1 = tic;
	for i = 1:L
		newH = 1;
		for k = 1:nChains
			if isempty(H{k})
				newH = kron(eye(dk),newH);
			else
				newH = kron(H{k},newH);
			end
		end
		cKRON = Vmat'*newH*Vmat;
	end
	tKRON(j) = toc(t1)

	t1 = tic;
	for i = 1:L
% 		Vmat = reshape(Vmat,[dk,dk*dOBB]);
% 		Vmat = Vmat'*Vmat;
		cHand = reshape(Vmat,[ones(1,nChains).*dk,dOBB]);
		%contract all empty parts
		ind = find(cellfun('isempty',H));		% finds all empty indices
		cHand = contracttensors(cHand,nChains+1,ind,conj(cHand),nChains+1,ind);
		ind = find(~cellfun('isempty',H));		% finds position to contract
		cHand = contracttensors(cHand,4,[1,3],H{ind},2,[1,2]);
	end
	tHand(j) = toc(t1)
end
%
figure(1); clf; hold on;
plot((1:5).*2,tHand);
plot((1:5).*2,tKRON);
set(gca,'yscale','log');
xlabel('d_{k}'); ylabel('t in s');
legend('contracttensors','kron');

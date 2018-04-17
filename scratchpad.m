%% Path where to save figures
path = pwd;
saveto = '..\Presentations\20140206 - Summary\';
%saveto = '..\..\Presentations\20140520 - MLSBM\20140525-CoupFunc\';
wantSave = 0;
%% Plot Vmat contributions for different Sites i (normalized)
i=2;
plot(abs((Vmat{1,i}(:,:))))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
if wantSave
    export_fig(sprintf('%sVmatNormalized%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end

%% Plot Vmat contributions for different Sites i (SV-weighted)
i = 5;
plot(real(Vmat{1,i}*diag(results.Vmat_sv{1,i})))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Sum over Vmat
i = 2;
a=sum(abs((Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2);
plot(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2))
set(gca,'YScale','log')
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Sum over Vtens
i = 2;
for j = 1:para.nChains
	a(1:size(Vmat{i}{j},1),j) = sum(abs((Vmat{i}{j}*diag(results.Vmat_sv{j,i}))),2)'./para.d_opt(j,i);
end
plot(a);
set(gca,'YScale','log')
title(['k = ',num2str(i)])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Sum over Vmat in 3D
plotMat = [];
mc = 1;						% choose the chain!
rotate3d off
zeroVals = -50;             % value of zeros for padding and -inf replacement
for i = 2:length(Vmat)
    if size(Vmat{1,i}{mc},2) < size(results.Vmat_sv{mc,i},1)
        a = log10(sum(abs((Vmat{1,i}{mc}*diag(results.Vmat_sv{mc,i}(1:size(Vmat{1,i}{mc},2),:)))),2));
    else
        a = log10(sum(abs((Vmat{1,i}{mc}(:,1:size(results.Vmat_sv{mc,i},1))*diag(results.Vmat_sv{mc,i}))),2));
    end
    a(a==-inf)=zeroVals;
    dim = max(length(a),size(plotMat,1));
    if length(a)< dim
        a = padarray(a,dim-length(a),zeroVals,'pre');
    elseif size(plotMat,1) < dim
        plotMat = padarray(plotMat,dim-size(plotMat,1),zeroVals,'pre');
    end
    plotMat = [plotMat,a];
end
surf(plotMat)
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{mc,i}(1,1))])
ylabel('$d_k$')
xlabel('Site $k$')
set(gca,'View',[9.5 40]);
formatPlot(1)
rotate3d on
axis tight
shading interp

%% Plot Results
f1=figure(1);
for i = 2:1:size(results.D,2)
    if isempty(results.D{i})
        results.D{i} = results.D{i-1};
    end
    if isempty(results.d_opt{i})
        results.d_opt{i} = results.d_opt{i-1};
    end
    if isempty(results.dk{i})
        results.dk{i} = results.dk{i-1};
    end
end
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
end
subplot(2,2,1);
    plot(real(results.nx));
    title('$$<n_x(k)>$$');
subplot(2,2,2);
    plot(para.trustsite);
    title('Trustsite')
% 3D-version elucidating change:
if para.useshift
    subplot(2,2,3);
        surf(cell2mat(results.shift'))
        set(gca,'View',[-25 10]);
        shading interp
        title('Bosonic shift');
end
subplot(2,2,4);
    surf(cell2mat(results.d_opt'));
    shading interp
    set(gca,'View',[0 90]);
    title('OBB dim')
% 2D-version showing final
% subplot(2,2,3);
%     plot(results.shift{end});
%     title('Bosonic shift');
% subplot(2,2,4);
%     plot(results.d_opt{end});
%     title('OBB dim')
%     text(-80,-30,sprintf(para.filename(1:38)))
if wantSave
    export_fig(sprintf('%sResultsSummary%s',saveto,para.filename(1:13)),'-transparent','-png','-painters')
end

%% Plot d_opt, D adjustments
% fill cell arrays
for i = 2:1:size(results.D,2)
    if isempty(results.D{i})
        results.D{i} = results.D{i-1};
    end
    if isempty(results.d_opt{i})
        results.d_opt{i} = results.d_opt{i-1};
    end
    if isempty(results.dk{i})
        results.dk{i} = results.dk{i-1};
    end
end
f2 = figure(2);
subplot(1,3,1);
surf(cell2mat(results.D'))
set(gca,'View',[0 90]);
shading interp
title('Change in bond dimension D')
subplot(1,3,2);
surf(cell2mat(results.d_opt'))
set(gca,'View',[0 90]);
shading interp
title('Change in $d_{opt}$')
subplot(1,3,3);
surf(cell2mat(results.dk'))
set(gca,'View',[0 90]);
shading interp
title('Change in $d_{k}$')
rotate3d on
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
%export_fig(['png/',para.folder,'-D-dopt.png'],'-transparent',f2)
if wantSave
    export_fig(sprintf('%sChangeOfdoptDdk%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% Plot shift adjustment history of results.shift
saveto = '..\Presentations\20140206 - Summary\';
% fill cell array
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
end
f3 = figure(3);
%subplot(1,2,1);
%surf(cell2mat(results.D'))
%shading interp
%title('Change in bond dimension D')
%subplot(1,2,2);
surf(cell2mat(results.shift'))
set(gca,'View',[-25 10]);
shading interp
title('Change of shift $\delta_k$')
xlabel('site k')
ylabel('sweep')
zlabel('$\delta_k$')
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
if wantSave
    export_fig(sprintf('%sChangeOfShift%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% Simulate shifting procedure
dim = 200; [bp,bm,n]=bosonop(dim,0,'n');
v1 = zeros(dim,1); v1(2:4) = 1/sqrt(3);
% 1. iteration
x = (bp+bm)./sqrt(2);
shift1 = v1'*x*v1;
[bp1n,bm1n,n1n]=bosonop(dim,shift1,'n');
bm1t = expm((bp-bm).*shift1./sqrt(2))'*bm*expm((bp-bm).*shift1./sqrt(2));
v1n = expm((bp-bm).*shift1./sqrt(2))'*v1;
% 2. iteration
x1n = (bp1n+bm1n)./sqrt(2);
shift2 = v1n'*x1n*v1n;
[bp2n,bm2n,n2n]=bosonop(dim,shift2,'n');
v2n = expm((bm1n-bp1n).*shift2./sqrt(2))*v1n;

plot([v1,v1n,v2n])
%% Try to analyse Differential Shift
diffShift = cell(1);
diffShift{1} = zeros(1,para.L);
diffShift{2} = zeros(1,para.L);
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
    if i~=2
        diffShift{i} = results.shift{i}-results.shift{i-1};
    end
end
f3 = figure(3);
%subplot(1,2,1);
%surf(cell2mat(results.D'))
%shading interp
%title('Change in bond dimension D')
%subplot(1,2,2);
surf(abs(cell2mat(diffShift')))
%set(gca,'View',[-25 10]);
set(gca,'View',[0 80]);
shading interp
title('Differential shift $\Delta (\delta_k)$')
xlabel('site k')
ylabel('sweep')
zlabel('$\delta_k$')
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
if wantSave
    export_fig(sprintf('%sDifferentialShift%s',saveto,para.filename(1:13)),'-transparent','-png')
end
%% plot results.Vmat_svLog
f4 = figure(4);
surf(cell2mat(results.Vmat_svLog'))
shading interp
title('Maximum SV of V')
xlabel('site k')
ylabel('sweep')
zlabel('SV')
set(gca,'View',[0 90]); %top
%set(gca,'View',[90 0]); %side
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
rotate3d on
if wantSave
    export_fig(sprintf('%sChangeOfVmatSV%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% Plot results.Amat_sv - logPlot
figure;
svCell = results.Amat_sv;
maxN = max(cellfun(@(x) length(x),svCell));
pTemp = zeros(maxN,length(svCell));
for ii = 1:length(svCell)
	pTemp(1:length(svCell{ii}),ii) = svCell{ii};
end
surf(log10(pTemp));
% shading interp;
axis tight; rotate3d on;
set(gca,'View',[0,0]);

%% plotting first rows of uneven cell arrays
cellfun(@(x) x(1,1), results.Vmat_sv(2:end))
%% try to calculate Shift analytically:
t = para.chain{1}.t; e = para.chain{1}.epsilon;
A = gallery('tridiag',t(2:end-1),e(1:end),t(2:end-1));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
B = zeros(para.L-1,1);
B(1) = -t(1)*1;%(results.spin.sz);  % -t1*sigmaZ
shift = A\B.*sqrt(2);
%% Plot my Shift versus VMPS Shift
colors={[0 0 1];[0 1 0];[1 0 0];[1 1 0];[1 0 1];[0 1 1];[0 0 0]};
hold all
pl(1) = plot([0;shift]);
set(pl(1), 'Marker', 'none','Color',colors{1}); % blue my
pl(2) = plot(para.shift);
set(pl(2), 'Marker', 'none','Color',colors{2}); % green VMPS
if wantSave
    export_fig(sprintf('%sAnalyticShift%s',saveto,para.filename(1:13)),'-transparent','-png','-eps')
end
%% Plot < n > of chain
figure(1); hold on;
results.nx = real(getObservable({'bath1correlators','n'}, mps,Vmat,para));
pl = plot(real(results.nx));
%set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<n_{k,VMPS}>$')
% set(gca,'yscale','log')
formatPlot(1)
if wantSave
    export_fig(sprintf('%s%s-Occupation',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot < x > of chain
figure(1); hold on;
results.xx = real(getObservable({'bath1correlators','x'}, mps,Vmat,para));
pl = plot(real(results.xx));
%set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<x_{k,VMPS}>$')
% set(gca,'yscale','log')
formatPlot(1)
if wantSave
    export_fig(sprintf('%s%s-Occupation',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot Chain epsilon-t
f = figure(3)
f.Name = 'Chain hopping and site energies';
subplot(1,2,1); hold all
% pl = plot(cell2mat(cellfun(@(x) x.t, para.chain, 'UniformOutput',false)));
cellfun(@(x) plot(x.t), para.chain, 'UniformOutput',false);
% set(gca,'YScale','log');
xlabel('Site k')
ylabel('$t_k$')
axis tight
legend(strsplit(sprintf('%d,', 1:length(para.chain)),','))
subplot(1,2,2); hold all;
cellfun(@(x) plot(x.epsilon), para.chain, 'UniformOutput',false);
% set(gca,'YScale','log');
xlabel('Site k')
ylabel('$\epsilon_k$')

%% Plot RC coupling over energy
% to estimate reorg. Energy
f = figure(4);
f.Name = 'RC reorganisation energies';
subplot(1,2,1);
pl = plot(cellfun(@(x) x.t(1)/x.epsilon(1), para.chain));
xlabel('Chain k');
ylabel('$\lambda/\omega$');

subplot(1,2,2);
pl = plot(cellfun(@(x) x.t(1)^2/4/x.epsilon(1), para.chain));
xlabel('Chain k');
ylabel('$\lambda^2/4\omega$');

% formatPlot(2)

%% Plot relative deviation from calculated < n >
%  Much more important as here the Wavefunction corrects also for different shift.
nx = [0 (shift.*shift./2)'];
figure(1)
plot((nx-real(results.nx))./nx,'LineStyle','none','Marker','*');
%set(gca,'YScale','log');
%title('Relative Deviation of $<n_k>$');
xlabel('Site k');
ylabel('$\frac{<n_k>-<n_{k,VMPS}>}{<n_k>}$')
formatPlot(1)
if wantSave
    export_fig(sprintf('%s%s-RelativeDeviationN',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot relative deviation of shift
relShift = ((shift-para.shift(2:end)')./shift);
figure(1)
pl(1) = plot([0; relShift],'LineStyle','none','Marker','*');
set(gca,'YScale','log');
%title('Relative deviation of shift from calculation')
xlabel('Site k')
ylabel('$\frac{\delta_k-\delta_{k,VMPS}}{\delta_k}$')
formatPlot(1)
if wantSave
    export_fig(sprintf('%sRelativeDeviationShift%s',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot ShiftUp vs ShiftDown
figure(5);
wantSave=0;
hold on
plot(results.bosonshift.x,'b')
plot(results.bosonshiftPerSpin.xUp,'r')
plot(results.bosonshiftPerSpin.xDown,'g')
title('Bosonic Shift vs Up or Down')
xlabel('site k')
ylabel('$\delta_k$')
legend('$x$','$x_\uparrow$','$x_\downarrow$');
ylim=get(gca,'YLim');
text(3,ylim(2)*2/3,...
    ['$x-(x_{\uparrow} +x_\downarrow) = $', sprintf('%.4g; \n', norm(results.bosonshift.x-results.bosonshiftPerSpin.xUp-results.bosonshiftPerSpin.xDown)),...
    '$|x|-|x_{\uparrow}|-|x_{\downarrow}| = $',sprintf('%.4g; ', norm(abs(results.bosonshift.x)-abs(results.bosonshiftPerSpin.xUp)-abs(results.bosonshiftPerSpin.xDown)))]);
formatPlot(5)
if wantSave
    export_fig(sprintf('%s%s-SBM-subOhmic-a0036622-shiftUpvsDown',saveto,para.filename(1:13)),'-transparent','-png','-painters')
end
%% Plot s-alpha relation
figure(5)
x = 0:0.01:1;
y = 0.005;
z = y.^(x./(1-x));
plot(x,z);

PlotData(:,PlotData(:,7)==0.001)

%% For PPC MLSBM:

%% Plot Energy convergence
figure(2);
plot(cell2mat(results.EvaluesLog)-min(cell2mat(results.EvaluesLog)));
% disp(sprintf('%.15e',results.E))
set(gca,'YScale','log');
% try
% title(sprintf('$E_0 = %.10g, \\Lambda =  %.2g, z =  %.2g$',results.E, para.Lambda, para.z));catch end
xlabel('Site$\cdot$Loop');
ylabel('$E-E_0$');
formatPlot(2)
yLim = get(gca,'YLim');
for i = 1:para.loop
%     line([para.L*i para.L*i],yLim,'LineWidth',1,'Color','black');
end
if wantSave
    export_fig(sprintf('%s%s-MLSBM-Econvergence-Lambda%.2gz%.2gp16',saveto,para.filename(1:13),para.Lambda,para.z),'-transparent','-png','-painters')
end

%% Plot Energy Error convergence
figure(3);
plot(results.Eerror);
set(gca,'YScale','log');
xlabel('Sweep');
ylabel('E error');
formatPlot(3);
%% Get System-Bath coupling
figure(2);
plot(real(diag(op.h2term{1,1,1})./para.t(1)));
xlabel('Site k');
ylabel('$\hat\eta$');
formatPlot(2);

%% Plot System Wavefunction
figure(3);
plot(diag(calReducedDensity(mps,Vmat,para,1)))
if para.MLSB_staticDisorder
    hold all
    plot(para.MLSB_disDiag/500);        % /500 scales to see both
end
%%
%% Plot Flowdiagram
loop = para.loop;
plot(results.flowdiag{1,loop});

%% Plot the CoupDiscr dataPoints
figure(1); clf; hold all;
pl = {};
for ii = 1:length(para.chain)
	pl{ii} = stem(para.chain{ii}.dataPoints(:,1),para.chain{ii}.dataPoints(:,2))

end

%% Get the RC parameters in cm^-1
x = res(1);
[cellfun(@(x) x.epsilon(1), x.para.chain);cellfun(@(x) x.t(1), x.para.chain)]./1.23984e-4

%% For TDVP analysis:

%% TDVP (1) SBM: Plot evolution of the spin
figure(11);clf;
hold all
sphereon = true;
if sphereon
    sphere
    daspect([1 1 1])
    alpha(0.2)
	set(get(gca,'children'),'linestyle',':')
end
col = parula(size(tresults.spin.sx,1));
scatter3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz,20,col,'filled');
plot3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz);
set(gca,'xlim',[-1,1]);
set(gca,'ylim',[-1,1]);
set(gca,'zlim',[-1,1]);
set(gca,'view',[-29,16]);
rotate3d on

%% TDVP (2) SBM: Plot <sz> / Visibility / Coherence
figure(20); hold all;
% tresults = res{3}.tresults;	% for res-extraction
% plotOpt = {'LineWidth',1.5};
% plotOpt = {'black','LineWidth',1};
plotOpt = {};
if isfield(tresults,'t')
	t=tresults.t;		% for the new convention when extracting in intervals >= rev42
else
	t=para.tdvp.t;		% for the old files
end
% n = find(tresults.n(:,3),1,'last');
n = find(t,1,'last');
plot(t(1:n), tresults.spin.sx(1:n), plotOpt{:});
plot(t(1:n), tresults.spin.sy(1:n), plotOpt{:});
plot(t(1:n), tresults.spin.sz(1:n), plotOpt{:});
plot(t(1:n), tresults.spin.visibility(1:n), plotOpt{:});
set(gca,'ylim',[-1,1]);
% set(gca,'xscale','log');
xlabel('t');
ylabel('$\left<\sigma_z\right>$');
l=legend('$\left< \sigma_x \right>$','$\left< \sigma_y \right>$','$\left< \sigma_z \right>$');
% l=legend('$\left< \sigma_x \right>$','$\left< \sigma_y \right>$','$\left| D(t) \right|$');
l.Interpreter = 'latex';

%% TDVP (2.1) SBM: Plot spin for adiabatic RDM:
figure(21); clf; hold all;

n = size(tresults.rho,1);
rho = squeeze(tresults.rho(:,:,:,3));		% 1: full rdm; 2: RDM of dominant state
% s2 = arrayfun(@(x) trace(squeeze(tresults.rho(x,:,:,1))*rho(x,:,:)),1:n);   % SV^2 corresponding to the bond state
s2 = diag(reshape(tresults.rho(:,:,:,1),n,[])*reshape(tresults.rho(:,:,:,2),n,[])');	% achieves the same
rho = reshape(rho,n,[]);			% reshape into t x 4
[sigmaX,sigmaY,sigmaZ] = spinop('Z');
sigmaX = reshape(sigmaX.',[],1);
sigmaY = reshape(sigmaY.',[],1);
sigmaZ = reshape(sigmaZ.',[],1);

plot(tresults.t,rho*sigmaX);
plot(tresults.t,rho*sigmaY);
plot(tresults.t,real(rho*sigmaZ));

set(gca,'ylim',[-1,1]);

%% TDVP (2) SBM: Plot Bloch length
figure(2); hold all;
plot(para.tdvp.t(1:length(tresults.spin.sz)), sqrt(tresults.spin.sz.^2+tresults.spin.sx.^2+tresults.spin.sy.^2));
plot(para.tdvp.t, tresults.spin.visibility);
set(gca,'ylim',[0,1]);
% set(gca,'xscale','log');
xlabel('t');
ylabel('$\sqrt{<s_x>^2+<s_y>^2+<s_z>^2}$');
legend('Bloch length','Visibility');

%% TDVP (2.3) RDM analysis separate plots
plotOpt = {'-fsev'};
f = figure(1); clf; hold all;
x.plot('rhoii',plotOpt{:});
formatPlot(f,'twocolumn-single')
grid on;

f = figure(2); clf; hold all;
x.plot('rhoij-imag',plotOpt{:});
formatPlot(f,'twocolumn-single')
grid on;

f = figure(3); clf; hold all;
x.plot('rhoij-real',plotOpt{:});
formatPlot(f,'twocolumn-single')
grid on;

%% TDVP (2.3) RDM analysis DPMES joint plots
plotOpt = {'-fsev'};
x(2).plot('rho-dpmes',plotOpt{:});
f = gcf; f.Name = 'Rho of DPMES Overview';
formatPlot(f,'twocolumn-single');
% diff(rhoii) overlay: select upper plot first!
% plot(x.t(1:x.lastIdx-1)*0.658,-25*diff(real(x.rho(1:x.lastIdx,1,1))),'-.','Color',[0.5,0.5,0.5],'LineWidth',1.5); % for LE+

%% TDVP Get TreeMPS vNE across all ER-node
% open a full .mat file and deserialise treeMPS
TDVPData.printSvNETreeMPS(treeMPS)

%% TDVP (3) Environment Plots
%% TDVP (3.1): Plot <n> CHAIN
mode = 0;		% 0: lin, 1: log
f=figure(312); clf; f.Name = 'Chain Occupation';
% x = res{9,1}; tresults = x.tresults; para = x.para;
x = res(9); tresults = x.tresults; para = x.para;
mc = 1;							% choose chain for display!
if str2double(para.tdvp.version(2:end)) < 50
	tresults.n = tresults.nx;
end
n = tresults.n(:,:,mc);
l = find(tresults.n(:,3),1,'last');
if isfield(tresults,'t')
	t=tresults.t;		% for the new convention when extracting in intervals >= rev42
else
	t=para.tdvp.t;		% for the old files
end
if mode
	surf(1:size(n,2),t(1:l),log10(abs(n(1:l,:))));
% 	surf(1:size(n,2),t(1:l),log10(abs(real(n(1:l,:))-ones(l,1)*real(n(1,:)))));		% subtract intitial population
	zlabel('$\log_{10}\left<n_k\right>$');
else
	surf(1:size(n,2),t(1:l),real(n(1:l,:)));
% 	surf(1:size(n,2),t(1:l),real(n(1:l,:))-ones(l,1)*real(n(1,:)));			% subtract initial population
	zlabel('$\left<n_k\right>$');
end
ax = gca;
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = ax.ZLabel.String;
xlabel('Site $k$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');se
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	ax.ZLim = [-30, max(max(ax.Children.ZData))];
% 	ax.ZLim = [-6, 0];
else
% 	ax.ZLim = [0.1,1].*10^-26;
end
ax.CLim = ax.ZLim;
% ax.YLim = [0,100];
%% TDVP (3.1): 2D Plot <n> CHAIN - L x t - TDVPData
mode = 0;		% 0: lin, 1: log
f=figure(313); clf; f.Name = 'Chain Occupation';
del = findobj(f,'Type','hgjavacomponent'); del.delete;	% get rid of old sliders
hold all;

% x = res(35);

if mode
	h = x.plotSld2D('chain-n','-log');
else
	h = x.plotSld2D('chain-n','-fsev');
end
%% TDVP (3.1): 1D Plot <n> Chain coherence - TDVPData
f=figure(315);  f.Name = ''; clf; hold all; ax = gca;
% x = res(6);
pl = x.plotSld1D('chain-n-c-t','-fsev',f);
leg = legend('$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$');
ylabel('$|TT\rangle\langle LE^+|\hat n_k$');
grid on
axis tight
formatPlot(f,'twocolumn-single')

%% TDVP (3.1): 1D Plot <n> Chain DFT - norm dist - TDVPData
f=figure(315);  f.Name = ''; clf; hold all; ax = gca;
x = res(4);
h = x.plotSld1DFT('chain-n','-cmev','-norm','-dist','-smoothres',f);
set(h.sld{2},'Value',20);	% start with 20 zero padding
leg = legend('$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$');
ylabel('$|FT (n_2)|$');
grid on
axis tight
[h.ax.XLim] = deal([0,2000]);
formatPlot(f,'twocolumn-single')
%% TDVP (3.1): 1D Plot <n> Chain DFT - dist - TDVPData
f=figure(316);  f.Name = ''; clf; hold all; ax = gca;
% x = res(4);
h = x.plotSld1DFT('chain-n','-cmev','-dist','-smoothres',f);
set(h.sld{2},'Value',20);	% start with 20 zero padding
set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
ylabel('$|FT (n_{RC})|$');
[h.ax.XLim] = deal([0,2000]);
formatPlot(f,'twocolumn-single')
drawnow; h.controlPanel.delete;
%% TDVP (3.1): 1D Plot <n> Chain DFT - dist - TDVPData
f=figure(316);  f.Name = ''; clf; hold all; ax = gca;
x = res(1);
h = x.plot('chain-n-rc-ft','-cmev','-resetColorOrder','-norm','-dist',f)
% set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
set(h.pl,{'Displayname'},{'$A_{2}$','$A_{1,1}$','$A_{1,2}$','$B_{1,1}$','$B_{1,2}$','$B_{2,2}$','$B_{2,3}$'}');  %for DPMES-Tree3
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
ylabel('$|FT (n_{RC})|$');
[h.ax.XLim] = deal([0,2000]);
formatPlot(f,'twocolumn-single')
drawnow;
%% TDVP (3.1): 1D Plot <n> RC - diabatic all chains, grid - TDVPData
% x = y(end);
h = x.plot('chain-n-d-rc','-fsev');
h.f.Units = 'pixels'; h.f.Position = get(0,'ScreenSize'); h.f.Name = 'Diabatic Occupation RC';
chainlbl = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
for ii = 1:length(chainlbl)
	text(h.f.Children(ii),0.9,0.9,chainlbl{end-ii+1},'sc');
	text(h.f.Children(ii),0.8,0.85,sprintf('$\\langle n \\rangle_{max} \\approx %g$',h.f.Children(ii).YLim(end)),'sc');
end
ax = [h.f.Children]; [ax.XLim] = deal([0,199]); %[ax.YLim] = deal([-0.07,0.3]);
legend(popName{:});
ax(end).Visible = 'off';
t = title(x.folder); t.Units = 'norm'; t.Position = [0.5,-0.13,0];
%% TDVP (3.1): 1D Plot <n> RC - adiabatic all chains, grid - TDVPData
% x = y(1);
h = x.plot('chain-n-a-rc','-fsev');
h.f.Units = 'pixels'; h.f.Position = get(0,'ScreenSize'); h.f.Name = 'Adiabatic Occupation RC';
chainlbl = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};
for ii = 1:length(chainlbl)
	text(h.f.Children(ii),0.9,0.9,chainlbl{end-ii+1},'sc');
	text(h.f.Children(ii),0.8,0.85,sprintf('$\\langle n \\rangle_{max} \\approx %g$',h.f.Children(ii).YLim(end)),'sc');
end
ax = [h.f.Children]; [ax.XLim] = deal([0,199]); %[ax.YLim] = deal([-0.07,0.3]);
legend toggle;
ax(end).Visible = 'off';
t = title(x.folder); t.Units = 'norm'; t.Position = [0.5,-0.13,0];
%% TDVP (3.1): 1D Plot <n> RC - adiabatic & diabatic, grid - TDVPData
h = x.plot('chain-n-a-rc','-fsev','-xlim',[0,200],'-rowheight',10,'-rowwidth',40);
h.f.OuterPosition = [0 0.05 1 0.95];
h.f.Name = 'Bond State Occupation RC';
chainlbl = x.chainLabel;
for ii = 1:length(chainlbl)
	text(h.ax(ii),0.9,0.9,chainlbl{ii},'sc');
	text(h.ax(ii),0.05,0.95,sprintf('$\\langle n \\rangle \\in [%g,%g]$',h.ax(ii).YLim(1),h.ax(ii).YLim(2)),'sc');
end
ax = h.ax; %[ax.YLim] = deal([-0.07,0.3]);
legend toggle;
t = title(x.folder);t.Units = 'normalized'; t.Position = [0.5,-0.13,0];
isLine = arrayfun(@(x) isa(x,'matlab.graphics.chart.primitive.Line'), h.pl);
[h.pl(isLine).LineStyle] = deal(':');
[h.pl(isLine).LineWidth] = deal(2);

cpStates = [1,2];		% copy States TT and LE+ for now!
h2 = x.plot('chain-n-d-rc','-fsev','-xlim',[0,200]);
h2.f.Visible = 'off';
pl = reshape(h2.pl,[],size(h2.pl,3));
for ii = 1:length(ax)
	copyobj(pl(ii,cpStates),ax(ii));
end
close(h2.f);
%% TDVP (3.1): 2D Plot <n> RC tFFT
nc = 2;
f = figure(310+nc+4); clf;f.Name = sprintf('Chain %d',nc);
x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;			% convert from simulation timescale to cm
h.data = real(squeeze(real(x.occC(1:h.m,2,nc))));		% t x L x nChain
h.data = TDVPData.movAvgRes(h.data,350/x.dt);
h.tdata = x.t;
h.Fs = 1/x.dt;
h.nwind = 1024*6;
%spectrogram(x     ,window        ,Nolap,NFFT/F,h.Fs,'yaxis','reassigned','power','MinThreshold',-40); 
[h.S,h.F,h.T,h.P] = spectrogram(h.data,kaiser(h.nwind,5),h.nwind-2,[0:1600]*h.tTocm,h.Fs,'yaxis','reassigned','power','MinThreshold',-50);	% only calculate for specified points!
%  spectrogram(h.data,kaiser(1024,5),1000,0:0.0001:0.04,h.Fs,'yaxis','reassigned','power','MinThreshold',-40);	% only calculate for specified points!
%  spectrogram(h.data,kaiser(1024,5),1000,(h.m*7),h.Fs,'yaxis','reassigned','power','MinThreshold',-40);		% to zero padding across entire range!
% ax=gca; ax.Children.YData = ax.Children.YData/0.658*4.135*8065.73*1e-3; ylim([0,2000]);
% surf(h.T,h.F,10*log10(abs(h.P).^2))
surf(h.T*0.658,h.F,mag2db(h.P))
axis tight; shading interp; view([0,90]);
ylabel('$E/cm^{-1}$');xlabel('$t/fs$');

ax=gca; ax.Children.YData = ax.Children.YData/h.tTocm; ylim([0,1600]); xlim([0,960]);
%% TDVP (3.1): 2D Plot <n> RC WSST
%	Wavelet Synchrosqueezed Transform
% not optimal, since needs better scale resolution! padding is not controllable!
nc = 3; fwind = 1;
popName = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};
for nc = [1:7]
% for fwind = [1,2,4,6,8]
f = figure(300+nc*10+fwind); clf;f.Name = sprintf('Chain WSST - %s',popName{nc});
x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;					% convert from simulation timescale to cm
h.data = real(squeeze(real(x.occC(1:h.m,2,nc))));		% t x L x nChain
h.data = TDVPData.movAvgRes(h.data,760/x.dt);
h.tdata = x.t*h.tTocm; h.dt = x.dt*h.tTocm;
h.f0 = 6/(2*pi);				% 5 for bump, 6 for morl
% h.scales = helperCWTTimeFreqVector(30,1600,h.f0,h.dt,128)/h.dt;	%	(fmin, fmax, f0, dt, nVoices), returns s*dt and not pure scale
h.scales = fliplr(h.f0/h.dt./(30:1600));
[h.sst,h.f]=mywsst(double(h.data),1/h.dt,'VoicesPerOctave',128,'amor'); ylim([0,1.6]);
% [h.sst,h.f] = mywsst(double(h.data),1/h.dt,'scales',h.scales,'amor'); h.f = fliplr(h.f0./h.scales./h.dt);% does not work! frequencies are stretched
% [h.sst,h.f] = sswt(double(h.data),1/h.dt,'fmin',50,'fmax',1600,'plot','amp'); 
%
h.sstDB = mag2db(abs(h.sst)/max(max(abs(h.sst))));
h.sstDB(h.sstDB<-90) = -inf;
surf(x.t*0.658,h.f,h.sstDB);
% surf(x.t*0.658,Freq,abs(h.sst));
axis tight; shading interp; view([0,90]); ylim([0,10000]);xlim([0,987]);
ylabel('$E/cm^{-1}$');xlabel('$t/fs$');
title(sprintf('Chain WSST - %s',popName{nc}));

export_fig(sprintf('img/%d',get(gcf,'Number')),'-transparent','-png','-m3');
close(f);
% end
end
%% TDVP (3.1): 2D Plot rhoii tFFT
pop = 1; fwind = 3;
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
for pop = 1:5
for fwind = [1,2,4,6,8]
f = figure(300+pop*10+fwind); clf;f.Name = sprintf('Population STFT - %s',popName{pop});
x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;					% convert from simulation timescale to cm
h.data = abs(x.getData('rhoii')); h.data = h.data(:,pop);		% t x L x nChain
h.data = TDVPData.movAvgRes(h.data,350/x.dt);
h.tdata = x.t;
h.Fs = 1/x.dt;
h.nwind = 1024*fwind;
%spectrogram(x     ,window        ,Nolap,NFFT/F,h.Fs,'yaxis','reassigned','power','MinThreshold',-40); 
[h.S,h.F,h.T,h.P] = spectrogram(h.data,kaiser(h.nwind,5),h.nwind-2,[0:1600]*h.tTocm,h.Fs,'yaxis','reassigned','power','MinThreshold',-50);	% only calculate for specified points!
%  spectrogram(h.data,kaiser(1024,5),1000,0:0.0001:0.04,h.Fs,'yaxis','reassigned','power','MinThreshold',-40);	% only calculate for specified points!
%  spectrogram(h.data,kaiser(1024,5),1000,(h.m*7),h.Fs,'yaxis','reassigned','power','MinThreshold',-40);		% to zero padding across entire range!
% ax=gca; ax.Children.YData = ax.Children.YData/0.658*4.135*8065.73*1e-3; ylim([0,2000]);
% surf(h.T,h.F,10*log10(abs(h.P).^2))
surf(h.T*0.658,h.F,mag2db(h.P))
axis tight; shading interp; view([0,90]);
ylabel('$E/cm^{-1}$');xlabel('$t/fs$');
title(sprintf('Population STFT - %s - Nwind %d',popName{pop}, h.nwind));

ax=gca; ax.Children.YData = ax.Children.YData/h.tTocm; ylim([0,1600]); xlim([0,960]);
% export_fig(sprintf('img/%d',get(gcf,'Number')),'-transparent','-png','-m3');
close(f);
end
end

%% TDVP (3.1): 2D Plot rhoii CWTft
% not optimal since needs reassignment!
pop = 3; fwind = 2;
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
for pop = 1:5
% for fwind = [1,2,4,6,8]
f = figure(400+pop*10+fwind); clf;f.Name = sprintf('Population CWTft - %s',popName{pop});
% x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;					% convert from simulation timescale to cm
h.data = abs(x.getData('rhoii')); h.data = h.data(:,pop);		% t x L x nChain
h.data = TDVPData.movAvgRes(h.data,760/x.dt);
h.tdata = x.t*h.tTocm; h.dt = x.dt*h.tTocm;
h.s0 = 2*h.dt;
h.f0 = 5/(2*pi);				% 5 for bump, 6 for morl
h.scales = helperCWTTimeFreqVector(30,1600,h.f0,h.dt,128);	%	(fmin, fmax, f0, dt, nVoices)
% h.scales = fliplr(h.f0/h.dt./(30:1600));
h.cwtdata = cwtft({h.data,h.dt},'wavelet','bump','scales',h.scales,'padmode','symw');
helperCWTTimeFreqPlot((h.cwtdata.cfs),h.tdata/h.tTocm*0.658,h.cwtdata.frequencies,...
    'surf',sprintf('Population CWTFT - %s',popName{pop}),'$t/fs$','$E/cm^{-1}$')
% export_fig(sprintf('img/%d',get(gcf,'Number')),'-transparent','-png','-m3');
% close(f);
end
%% TDVP (3.1): 2D Plot rhoii WSST
%	Wavelet Synchrosqueezed Transform
% not optimal, since needs better scale resolution! padding is not controllable!
pop = 3; fwind = 1;
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
for pop = 1
% for fwind = [1,2,4,6,8]
f = figure(300+pop*10+fwind); clf;f.Name = sprintf('Population WSST - %s',popName{pop});
% x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;					% convert from simulation timescale to cm
h.data = abs(x.getData('rhoii')); h.data = h.data(:,pop);		% t x L x nChain
h.data = TDVPData.movAvgRes(h.data,760/x.dt);
h.tdata = x.t*h.tTocm; h.dt = max(diff(x.t))*h.tTocm;
h.f0 = 6/(2*pi);				% 5 for bump, 6 for morl
% h.scales = helperCWTTimeFreqVector(30,1600,h.f0,h.dt,128)/h.dt;	%	(fmin, fmax, f0, dt, nVoices), returns s*dt and not pure scale
% h.scales = fliplr(h.f0/h.dt./(30:1600));
% [h.sst,h.f]=mywsst(double(h.data),1/h.dt,'VoicesPerOctave',128,'amor'); %ylim([0,1.6]);
% [h.sst,h.f] = mywsst(double(h.data),1/h.dt,'scales',h.scales,'amor'); h.f = fliplr(h.f0./h.scales./h.dt);% does not work! frequencies are stretched
[h.sst,h.f,h.wopt,h.wt,h.wtf,h.ifr] = sswt(double(h.data),1/h.dt,'fmin',50,'fmax',1600,'plot','amp++-wr','Padding',0,'nv',256,'Wavelet','Morlet'); 

% separate plot, if needed
h.sstDB = mag2db(abs(h.sst)/max(max(abs(h.sst))));
h.sstDB(h.sstDB<-90) = -inf;
f = figure(300+pop*10+fwind); 
surf(x.t(1:h.m)*0.658,h.f,h.sstDB); shading interp, view([0,90]), axis tight
% surf(x.t*0.658,Freq,abs(h.sst));
end
%%
for pop = 1:5
figure(300+pop*10+fwind);
axis tight; shading interp; view([0,90]); ylim([0,1600]);xlim([0,987]);
ylabel('$E/cm^{-1}$');xlabel('$t/fs$');
title(sprintf('Population WSST - %s',popName{pop}));

export_fig(sprintf('img/%d',get(gcf,'Number')),'-transparent','-png','-m3');
% close(f);
% end
end
%% TDVP (3.1): 2D Plot rhoij-imag CWTft
% not optimal since needs reassignment!
pop = 4; fwind = 2;
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
% for pop = 1:5
% for fwind = [1,2,4,6,8]
f = figure(400+pop*10+fwind); clf;%f.Name = sprintf('Coherence CWTft - %s',popName{pop});
% x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73;					% convert from simulation timescale to cm
out = x.getData('rhoij'); h.data = imag(out{1}(:,pop));         % out: { t x coherences, state numbers}
h.data = TDVPData.movAvgRes(h.data,760/x.dt);
h.tdata = x.t(1:x.lastIdx)*h.tTocm; h.dt = x.dt*h.tTocm;
h.s0 = 2*h.dt;
h.f0 = 5/(2*pi);				% 5 for bump, 6 for morl
h.scales = helperCWTTimeFreqVector(30,8000,h.f0,h.dt,128);	%	(fmin, fmax, f0, dt, nVoices)
% h.scales = fliplr(h.f0/h.dt./(30:1600));
h.cwtdata = cwtft({h.data,h.dt},'wavelet','bump','scales',h.scales,'padmode','symw');
helperCWTTimeFreqPlot((h.cwtdata.cfs),h.tdata/h.tTocm*0.658,h.cwtdata.frequencies,...
    'surf',sprintf('Coherence CWTFT - %s-%s',popName{out{2}(1,pop)},popName{out{2}(2,pop)}),'$t/fs$','$E/cm^{-1}$')
% export_fig(sprintf('img/%d',get(gcf,'Number')),'-transparent','-png','-m3');
% close(f);
% end
%% TDVP (3.1) Wavelet Coherence between states and <n> RC
x = res(1); 
h.m = x.lastIdx; h.tTocm = 0.658/4.135/8065.73; h.dt = max(diff(x.t))*h.tTocm;
h.tdata = x.t(1:h.m);
h.Odata = real(squeeze(real(x.occC(1:h.m,2,:)))); h.Odata = TDVPData.movAvgRes(h.Odata,760/x.dt);
h.pdata = abs(x.getData('rhoii'));                h.pdata = TDVPData.movAvgRes(h.pdata,760/x.dt);
%% state-RC coherences
popName = {'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'}; chainName = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};
[~,m] = size(h.pdata); [~,n] = size(h.Odata);
hPl = TDVPData.plotGrid(m,n); hPl.f.Units = 'pixels';
hPl.f.Position = [1,1,1920,1080]; hPl.f.Name = 'Wavelet Coherence state-RC';
for kk = 1:m*n
	[ii,jj] = ind2sub([m,n],kk);
	axes(hPl.ax(kk));
	wcoherence(h.pdata(:,ii),h.Odata(:,jj),1/h.dt,'VoicesPerOctave',32,'PhaseDisplayThreshold',0.5,'NumScalesToSmooth',10);
	hPl.ax(kk).YLim(end) = 4;		% 2^4 = 16 kcm^-1
	hPl.ax(kk).YTick = hPl.ax(kk).YLim(1):(hPl.ax(kk).YLim(end)-1);
	title('');
	if ii < m 
		xlabel(''); 
		hPl.ax(kk).XTickLabel = '';
		if ii == 1
			title(chainName{jj});
		end
	else
		hPl.ax(kk).XTickLabel = strsplit(sprintf('%d ',round(hPl.ax(kk).XTick*10^-3/h.tTocm*0.658)));
		xlabel('$t/fs$');
	end
	if jj > 1
		ylabel(''); 
		hPl.ax(kk).YTickLabel = '';
		if jj == n
			hPl.ax(kk).YAxisLocation = 'right';
			ylabel(popName{ii});
		end
	else
		ylabel('$E/(10^3 cm^{-1})$');
	end
	colorbar off
	drawnow
end

%% inter-state coherences
[m,n] = size(h.pdata);
hPl = TDVPData.plotGrid(n-1,n-1); hPl.f.Units = 'pixels';
hPl.f.Position = [1,1,1920,1080]; hPl.f.Name = 'Wavelet Coherence inter-state';
hPl.ax = reshape(hPl.ax,n-1,n-1);
for kk = 1:n*n
	[ii,jj] = ind2sub([n,n],kk);
	if ii >= jj
		if all([ii,jj-1] > 0) && ii < n
			hPl.ax(ii,jj-1).Visible = 'off';
		end
		continue;
	end
	axes(hPl.ax(ii,jj-1));
	wcoherence(h.pdata(:,ii),h.pdata(:,jj),1/h.dt,'VoicesPerOctave',32,'PhaseDisplayThreshold',0.5,'NumScalesToSmooth',10);
	hPl.ax(ii,jj-1).YLim(end) = 4;		% 2^4 = 16 kcm^-1
	hPl.ax(ii,jj-1).YTick = hPl.ax(ii,jj-1).YLim(1):(hPl.ax(ii,jj-1).YLim(end)-1);
	title('');
	if ii < n 
		if ii == 1
			title(popName{jj});
		end
		xlabel(''); 
		hPl.ax(ii,jj-1).XTickLabel = '';
	end
	if jj > 2
		if jj ~= ii+1
			ylabel(''); 
			hPl.ax(ii,jj-1).YTickLabel = '';
		end
		if jj == n
			hPl.ax(ii,jj-1).YAxisLocation = 'right';
			ylabel(popName{ii});
		end
	end
	if ii == jj-1 && ii ~= n-1
		ylabel('$E/(10^3 cm^{-1})$');
		hPl.ax(ii,jj-1).XTickLabel = strsplit(sprintf('%d ',round(hPl.ax(ii,jj-1).XTick*10^-3/h.tTocm*0.658)));
		xlabel('$t/fs$');
	end
	if ii == n-1
		hPl.ax(ii,jj-1).YTickLabel = '';
		hPl.ax(ii,jj-1).XTickLabel = strsplit(sprintf('%d ',round(hPl.ax(ii,jj-1).XTick*10^-3/h.tTocm*0.658)));
		xlabel('$t/fs$');
	end
	colorbar off
	drawnow
end

%% TDVP (3.2): Plot <n> STAR
mode = 0;		% 0: lin, 1: log
f=figure(320); clf; f.Name = 'Star Occupation';
l = find(tresults.star.t,1,'last');
mc = 2;
n = tresults.star.n(:,:,mc);
if numel(tresults.star.omega) == length(tresults.star.omega)
	omega = tresults.star.omega;
else
	omega = tresults.star.omega(:,mc);
end
if mode
	surf(omega,tresults.star.t(1:l),log10(n(1:l,:)));
% 	surf(tresults.star.omega,tresults.star.t(1:n),log10(tresults.star.n(1:n,:))-log10(ones(n,1)*tresults.star.n(1,:)));
	zlabel('$\log_{10}\left<n_k\right>$');
else
	surf(omega,tresults.star.t(1:l),n(1:l,:));
	zlabel('$\left<n_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-3,1]);
% 	set(gca,'clim',get(gca,'zlim'));
end

	%% TDVPData

%% TDVP (3.3): Animate <n> propagation
figure(3); clf;
ax = axes('units','pixels');
pl = plot(1:para.L,log10(real(tresults.n(1,:,1))));
set(gca,'ylimmode','manual');
set(gca,'ylim',[-5,1]);
set(gca,'xlimmode','manual','xlim',[1,para.L]);
xlabel('Site $k$');
ylabel('$\left<n_k\right>$');
% shading interp
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(tresults.n,1),'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'ydata',log10(real(tresults.n(round(source.Value),:,1)))));
%% TDVP (3.4): Animate <n> star formation
mode = 0;		% 0: lin, 1: log
figure(3); clf;
ax = axes('units','pixels');
n = find(tresults.star.t,1,'last');
if mode
	pl = plot(tresults.star.omega,log10(tresults.star.n(1,:)));
	ylabel('$log_{10}\left<n_k\right>$');
	set(gca,'ylim',[-5,log10(max(max(tresults.star.n)))]);
else
	pl = plot(tresults.star.omega,tresults.star.n(1,:));
	ylabel('$\left<n_k\right>$');
	set(gca,'ylim',[0,max(max(tresults.star.n))]);
end
set(gca,'ylimmode','manual');
set(gca,'xlimmode','manual','xlim',[0,tresults.star.omega(end)]);
xlabel('Site $k$');
% shading interp
if mode
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',n,'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',log10(tresults.star.n(round(source.Value),:))));
else
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',n,'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',tresults.star.n(round(source.Value),:)));
end
%% TDVP (3.5): Animate <n> star single mode kinetics
mode = 0;		% 0: lin, 1: log
figure(3); clf;
n = find(tresults.star.t,1,'last');
if mode
	pl = plot(tresults.star.t(1:n),log10(tresults.star.n(1:n,1)));
	ylabel('$log_{10}\left<n_k\right>$');
	set(gca,'ylim',[-5,log10(max(max(tresults.star.n(1:n))))]);
else
	pl = plot(tresults.star.t(1:n),tresults.star.n(1:n,1));
	ylabel('$\left<n_k\right>$');
% 	set(gca,'ylim',[0,max(max(tresults.star.n))]);
end
% set(gca,'ylimmode','manual');
% set(gca,'xlimmode','manual','xlim',[0,tresults.star.t(end)]);
xlabel('Time $\omega_c t$');
ax = gca; ax.UserData = 1;
% OnChange actions:
hPl = handle(ax); hProp = findprop(hPl,'UserData');
if mode
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl,'ydata',log10(tresults.star.n(1:n,ax.UserData))));
else
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl,'ydata',tresults.star.n(1:n,ax.UserData)));
end
	hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$ \\omega_k = %g $',tresults.star.omega(ax.UserData))));

% set(gca,'ylimmode','manual');
% set(gca,'xlimmode','manual','xlim',[0,tresults.star.t(end)]);
xlabel('Time $\omega_c t$');

% % slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmax = size(tresults.star.n,2);
sld = javax.swing.JScrollBar(0,1,1,1,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax,'UserData',round(source.Value)));

%%			 : Capture Movie from existing slider
fignum = 4;
f = figure(fignum);
nFrames = sld.getMaximum;
F(nFrames) = struct('cdata',[],'colormap',[]);
for i = 1:nFrames
	ax.UserData = i;
	F(i) = getframe(f);
end
%% Playback Movie
fig = figure;
movie(fig,F,1,200)
%% Save Movie
writerObj = VideoWriter('img/Polaron001-dt01','MPEG-4');
writerObj.FrameRate = 60;
open(writerObj);
writeVideo(writerObj,F);
close(writerObj);

%% TDVP (3.6): Plot FT(<n>) star
% possible to restrict in time domain!
mode = 1;		% 0: lin, 1: log
f=figure(7); clf; f.Name = 'Star Occupation Fourier';
maxN = floor(find(tresults.star.t,1,'last')/2)*2;
% maxN = 500;
freq = 2*pi/tresults.star.t(2)/2 * linspace(0,1,maxN/2+1);		% fs * k/N  where k=0... N/2
FT = fft(tresults.star.n(1:maxN,:),maxN,1)/maxN;
if mode
	surf(tresults.star.omega,freq,log10(2*abs(FT(1:maxN/2+1,:))));
	zlabel('$\log_{10} \mathcal{F} \left\{ \left<n_k\right>(t)\right\} $');
else
	surf(tresults.star.omega,freq,2*abs(FT(1:(size(FT,1)/2+1),:)));
	zlabel('$\mathcal{F} \left\{ \left<n_k\right>(t)\right\}$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Frequency $\omega$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1.2]);
% 	set(gca,'clim',[-4,1.2]);
end

%% TDVP (3.7): Extract resonance of renormalized splitting
mode = 1;		% 0: lin, 1: log
figure(6); clf;
% tresults = res{42,1}.tresults;
if mode
	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n));
% 	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n)-log10(ones(size(tresults.star.n,1),1)*tresults.star.n(1,:)));
	zlabel('$\log_{10}\left<n_k\right>$');
else
	surf(tresults.star.omega,tresults.star.t,tresults.star.n);
	zlabel('$\left<n_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1]);
% 	set(gca,'clim',[-4,1]);
end
if exist('axIn','var')
	axIn.delete;
end
ax=gca; oldPos = ax.Position;
axIn = axes('Position',[oldPos(1)+0.32, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);		% use this if main with colorbar
% axIn = axes('Position',[oldPos(1)+0.36, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);	% use this if main without colorbar
axIn.TickDir = 'out';
axIn.FontSize=14; axIn.XColor=[0.7,0.7,0.7];axIn.YColor=axIn.XColor;

hold all;
surf(tresults.star.omega,tresults.star.t, ax.Children.ZData);
set(gca,'View',[0 90]);
shading interp
axis tight
axIn.XLim = [0,0.2]; %axIn.YLim = [0,800];
[M,I] = max(ax.Children.ZData(:,1:300),[],2);
plot3(tresults.star.omega(I(2:end)),tresults.star.t(2:end),ones(1,length(tresults.star.t)-1).*1000,':r')
plot3([0.1,0.1],axIn.YLim,[1000,1000],':black');

tresults.star.omega(I(2:end))

%% TDVP: Draw <n> propagation fit
figure(3); clf; hold on
n = size(tresults.nx,1); nCont = 30; picks = 8;
contourf(1:para.L,para.tdvp.t(1:n),real(tresults.nx),nCont)
%%
coneParam = zeros(2,picks);
hold on
for m=1:picks
	[coneParam(1,m),coneParam(2,m)] = ginput(1);
	plot(coneParam(1,m),coneParam(2,m),'ko','markerSize',10,'lineWidth',1,'markerFaceColor','w');
end
%%
% pf = polyfit(coneParam(1,:),coneParam(2,:),3);
% coneFunct = polyval(pf,para.tdvp.t(1:n));
pf = fnxtr(csaps(coneParam(2,:),coneParam(1,:)));
L = 800;
coneFunct = fnval(pf,1:L);
plot(coneFunct,1:L,'w','linewidth',4)
plot(coneFunct,1:L,'k','linewidth',2)

%% TDVP (3.8.1): 2D Plot <x> CHAIN - TDVPData
% x = res(35);
h = x.plotSld2D('chain-x','-fsev');
formatPlot(gcf,'twocolumn-single')

%% TDVP (3.8.2): 2D Plot <x> CHAIN - State avg - TDVPData
% x = res(35);
h = x.plotSld2D('chain-x-avg','-fsev');
formatPlot(gcf,'twocolumn-single')

%% TDVP (3.8.3): 1D Plot <x> CHAIN - TDVPData
% x = res(35);
h = x.plotSld1D('chain-x-t','-fsev');
formatPlot(gcf,'twocolumn-single')

%% TDVP (3.1): 1D Plot <x> RC - diabatic all chains, grid - TDVPData
%	displays <x> range and annotates Chain Labels
%	
% x = y(1);
h = x.plot('chain-x-d-rc','-fsev','-xlim',[0,110],'-rowwidth',40,'-rowheight',10);
h.f.Units = 'pixels'; h.f.Position = get(0,'ScreenSize'); h.f.Name = 'Diabatic Displacement RC';
% chainlbl = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};
chainlbl = x.chainLabel;
% popName = {'$\langle x \rangle=0$','TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'};
for ii = 1:length(chainlbl)
	text(h.ax(ii),0.9,0.9,chainlbl{ii},'sc');
	text(h.ax(ii),0.05,0.95,sprintf('$\\langle x \\rangle \\in [%g,%g]$',h.ax(ii).YLim(1),h.ax(ii).YLim(2)),'sc');
end
ax = [h.ax]; %[ax.YLim] = deal([-0.07,0.3]);
% legend(popName{:});
legend toggle;
t = title(x.folder);t.Units = 'normalized'; t.Position = [0.5,-0.13,0];
%% TDVP (3.1): 1D Plot <x> RC - adiabatic all chains, grid - TDVPData
% x = y(1);
h = x.plot('chain-x-a-rc','-fsev','-xlim',[0,110],'-rowheight',10,'-rowwidth',40);
h.f.OuterPosition = [0 0.05 1 0.95];
h.f.Name = 'Adiabatic Displacement RC';
% chainlbl = {'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'};  % for old DPMES
% chainlbl = {'$B_{1,1}$','$B_{1,2}$','$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{2,1}$','$B_{2,2}$'};		% for DPMES-Tree4
chainlbl = x.chainLabel;
for ii = 1:length(chainlbl)
	text(h.ax(ii),0.9,0.9,chainlbl{ii},'sc');
	text(h.ax(ii),0.05,0.95,sprintf('$\\langle x \\rangle \\in [%g,%g]$',h.ax(ii).YLim(1),h.ax(ii).YLim(2)),'sc');
end
ax = [h.ax]; %[ax.YLim] = deal([-0.07,0.3]);
legend toggle;
t = title(x.folder);t.Units = 'normalized'; t.Position = [0.5,-0.13,0];

%% TDVP (3.1.1): 1D Plot <x> RC - adiabatic & diabatic, grid - TDVPData
h = x.plot('chain-x-a-rc','-fsev','-xlim',[0,200],'-rowheight',10,'-rowwidth',40);
h.f.OuterPosition = [0 0.05 1 0.95];
h.f.Name = 'Bond State Displacement RC';
chainlbl = x.chainLabel;
for ii = 1:length(chainlbl)
	text(h.ax(ii),0.9,0.9,chainlbl{ii},'sc');
	text(h.ax(ii),0.05,0.95,sprintf('$\\langle x \\rangle \\in [%g,%g]$',h.ax(ii).YLim(1),h.ax(ii).YLim(2)),'sc');
end
ax = h.ax; %[ax.YLim] = deal([-0.07,0.3]);
legend toggle;
t = title(x.folder);t.Units = 'normalized'; t.Position = [0.5,-0.13,0];
isLine = arrayfun(@(x) isa(x,'matlab.graphics.chart.primitive.Line'), h.pl);
[h.pl(isLine).LineStyle] = deal(':');
[h.pl(isLine).LineWidth] = deal(2);

cpStates = [1,2];		% copy States TT and LE+ for now!
h2 = x.plot('chain-x-d-rc','-fsev','-xlim',[0,110]);
h2.f.Visible = 'off';
for ii = 1:length(ax)
	copyobj(h2.pl(cpStates,ii),ax(ii));
end
close(h2.f);
%% TDVP (3.8.3): 1D Plot <x> CHAIN - Diabatic - TDVPData
% x = res(35);
h = x.plotSld1D('chain-x-d-t','-fsev');
% formatPlot(gcf,'twocolumn-single')
formatPlot(gcf);

%% TDVP (3.8.4): 1D Plot <x> CHAIN - Coherence - TDVPData
f=figure(381);  f.Name = ''; clf; hold all; ax = gca;
% x = res(6);
pl = x.plotSld1D('chain-x-c-t','-fsev',f);
leg = legend('$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$');
ylabel('$|TT\rangle\langle LE^+|\hat x_k$');
grid on
axis tight
formatPlot(f,'twocolumn-single')

%% TDVP (3.8.1): Plot FT(<x>) CHAIN - TDVPData
% x = res(35);
% h = x.plotSld2D('chain-x-ft','-fsev'); h.ax.YLim = [0,0.6];
h = x.plotSld2D('chain-x-ft','-fsev','-cm'); h.ax.YLim = [0,3e3];
% 

%% TDVP (3.8): Plot polaron STAR
f=figure(9); clf; f.Name = 'Star Polaron';
% x = res{31,1};
% tresults = x.tresults; para = x.para;
nStat = para.dk(1,1);
mc = 1;							% choose which chain!
n = find(tresults.star.t,1,'last');
m = find(tresults.star.omega(:,mc),1,'last');
x = tresults.star.x(1:n,1:m,:,mc);
omega = tresults.star.omega(1:m,mc);
ax = cell(2,1);
ax{1} = axes();
surf(omega,tresults.star.t(1:n),x(:,:,1));
ax{1}.Position(1) = ax{1}.Position(1)-0.02; ax{1}.Position(3) = ax{1}.Position(3)/nStat;
zlabel('$f_k^\uparrow$');

ax{2} = axes('Position',[ax{1}.Position(1)+ax{1}.Position(3)+0.02, ax{1}.Position(2:4)]);
surf(omega,tresults.star.t(1:n),x(1:n,:,2));
zlabel('$f_k^\downarrow$');
ax{2}.YAxisLocation = 'right';
if nStat > 2
	ax{2}.YTickLabel = [];
	ax{3} = axes('Position',[ax{2}.Position(1)+ax{2}.Position(3)+0.02, ax{2}.Position(2:4)]);
	surf(omega,tresults.star.t(1:n),x(1:n,:,3));
	zlabel('$f_k^\downarrow$');
	ax{3}.YAxisLocation = 'right';
end
for i = 1:length(ax)
	axes(ax{i});
	cb = colorbar;cb.Title.Interpreter = 'latex';cb.Location = 'North';
	cb.Title.String = ax{i}.ZLabel.String;
	cb.Position(4) = cb.Position(4)/2;
	xlabel('Mode $\omega_k / \omega_c$');
	ylabel('Time $\omega_c t$');
% 	set(gca,'View',[0 42]);
	set(gca,'View',[0 90]);
	shading interp
	rotate3d on
	axis tight
end
%% TDVP (3.9): Plot polaron STAR spin-corrected
mode = 0;		% 0: lin, 1: log
f=figure(9); clf; f.Name = 'Star Polaron corrected';
% tresults = res{14,1}.tresults;
n = find(tresults.star.t,1,'last');
normFact = tresults.spin.sz(ismember(tresults.t, tresults.star.t(1:n)),1)./2+1/2;
normFact = normFact(1:n)*ones(1,size(tresults.star.n,2));
if mode
	surf(tresults.star.omega,tresults.star.t(1:n),log10(abs(tresults.star.x(1:n,:)./normFact)).*sign(tresults.star.x(1:n,:)));
	zlabel('$\log_{10}\left<x_k\right>$');
else
	surf(tresults.star.omega,tresults.star.t(1:n),(-tresults.star.x(1:n,:)./normFact));
	zlabel('$\left<x_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'xscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1]);
% 	set(gca,'clim',[-4,1]);
end

%% TDVP (3.10): Animate STAR polaron kinetics up/down
guides = 1;		% plot parabola mid-points
guide1 = 'off'; guide2 = 'off'; guide3 = 'on';		% 'on' / 'off', 1: wl, 2: max Displacement, 3: Silbey-Harris
fignum = 4; f = figure(fignum); clf;
ax = axes('box','on'); hold all; ax.UserData = 1;
hPl = handle(ax); hProp = findprop(hPl,'UserData');
n = find(tresults.star.t,1,'last');
ax.XLim = tresults.star.t([1,n]);
x = tresults.star.x(:,:,:,2);		% Choose Chain!
% Sim Parameters
a = para.chain{1}.alpha; s = para.chain{1}.s;
% Plots & OnChange actions on ax.UserData
pl1 = plot(tresults.star.t(1:n),(x(1:n,1,1)));
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl1,'ydata',x(1:n,ax.UserData,1)));
if size(x,3) == 2
	pl2 = plot(tresults.star.t(1:n),(x(1:n,1,2)));
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl2,'ydata',x(1:n,ax.UserData,2)));
end
ylabel('$\left<f_k\right>$');
set(gca,'ylim',[min(min(min(tresults.star.x(1:n,:,:))));max(max(max(x(1:n,:,:))))]);
if s == 1
	deltaR = abs(para.hx)^(1/(1-a));
else
	deltaR = para.hx;		% Need good formula for renormalized amplitude!
end
if guides
	% expected Position of Parabolas:
	% complete displacement
	pl3 = plot(tresults.star.t(1:n),-(tresults.spin.sz(1:n)+1)./4.*a./tresults.star.omega(1),'black--');
	pl4 = plot(tresults.star.t(1:n),(1-tresults.spin.sz(1:n))./4.*a./tresults.star.omega(1),'black--');
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl3,'ydata',-(tresults.spin.sz(1:n)+1)./4.*sqrt(2*a/tresults.star.omega(ax.UserData))));
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl4,'ydata',(1-tresults.spin.sz(1:n))./4.*sqrt(2*a/tresults.star.omega(ax.UserData))));
	pl3.Visible = guide2; pl4.Visible = guide2;
	% one wavelength
	pl5 = plot([1 1]*2*pi/tresults.star.omega(1),get(gca,'ylim'),'black--');
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl5,'xdata',[1 1]*2*pi/tresults.star.omega(ax.UserData)));
	pl5.Visible = guide1;
	% Silbey-Harris
	pl6 = plot(tresults.star.t(1:n),-(tresults.spin.sz(1:n)+1)./4.*a./tresults.star.omega(1),'--', 'color',[0.5,0.5,0.5]);
	pl7 = plot(tresults.star.t(1:n),(1-tresults.spin.sz(1:n))./4.*a./tresults.star.omega(1),'--','color',pl6.Color);
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl6,'ydata',-(tresults.spin.sz(1:n)+1)./4.*sqrt(2*a*tresults.star.omega(ax.UserData))/(tresults.star.omega(ax.UserData)+deltaR)));
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl7,'ydata',(1-tresults.spin.sz(1:n))./4.*sqrt(2*a*tresults.star.omega(ax.UserData))/(tresults.star.omega(ax.UserData)+deltaR)));
	pl6.Visible = guide3; pl7.Visible = guide3;
end
% labels and comments
% set(gca,'ylimmode','manual');
% set(gca,'xlimmode','manual','xlim',[0,tresults.star.t(end)]);
xlabel('Time $\omega_c t$');
l=legend([pl1,pl2,pl3,pl6],'$f_k^{\uparrow}$','$f_k^{\downarrow}$','$\frac{\pm 1-\sigma_z}{2} \frac{g_k}{2\omega_k}$','Silbey-Harris', 'Location','NorthEastOutside');
l.Interpreter = 'latex';
f.Position(3:4) = [720,416];	% correction needed for legend
% Text axes
axT = axes(	'Position',[l.Position(1),ax.Position(1),l.Position(3),0.6],...
			'box','on', 'Visible','off');
axT.YLim = [0,0.6*f.Position(4)];
t1 = text(0.1,210,{sprintf('$s=%g$',s);...
				   sprintf('$\\alpha=%g$',a);...
				   sprintf('$\\Delta t=%g$',para.tdvp.deltaT);...
				   sprintf('$ \\Delta_r = %.3g $',deltaR)});
t2 = text(0.1,t1.Extent(2)-9,sprintf('$ \\omega_k = %g$',tresults.star.omega(ax.UserData)));
hPl.addlistener(hProp,'PostSet',@(src,event) set(t2, 'String',sprintf('$ \\omega_k = %g $',tresults.star.omega(ax.UserData))));

% slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmax = size(tresults.star.x,2);
sld = javax.swing.JScrollBar(0,1,1,1,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax,'UserData',round(source.Value)));
%% TDVP (3.10): Animate STAR polaron kinetics corrected
mode = 0;		% 0: lin, 1: log
figure(3); clf;
ax = axes('units','pixels');
n = find(tresults.star.t,1,'last');
normFact = tresults.spin.sz(ismember(tresults.t, tresults.star.t(1:n)),1)./2+1/2;
normFact = normFact(1:n)*ones(1,size(tresults.star.n,2));
tresults.star.x(1:n,:)= tresults.star.x(1:n,:)./normFact;
if mode
	pl = plot(tresults.star.t(1:n),log10(tresults.star.x(1:n,1)));
	ylabel('$log_{10}\left<x_k\right>$');
	set(gca,'ylim',[-5,log10(max(max(tresults.star.x(1:n,:))))]);
else
	pl = plot(tresults.star.t(1:n),tresults.star.x(1:n,1));
	ylabel('$\left<x_k\right>$');
	set(gca,'ylim',[min(min(tresults.star.x(1:n,:))),max(max(tresults.star.x(1:n,:)))]);
end
% set(gca,'ylimmode','manual');
% set(gca,'xlimmode','manual','xlim',[0,tresults.star.t(end)]);
xlabel('Time $\omega_c t$');
% shading interp
if mode
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',size(tresults.star.x,2),'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',log10(tresults.star.x(1:n,round(source.Value)))));
else
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',size(tresults.star.x,2),'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',tresults.star.x(1:n,round(source.Value))));
end
ax.Units = 'norm';
%% TDVP (3.10): Animate STAR <n> <x> kinetics side-by-side
fignum = 5; figure(fignum); clf;
% x = res{9,1}; tresults = x.tresults; para = x.para;
nc = 1;		% Choose chain!
guideSH = 0;
width = 0.8; height = 0.375; posx = 0.1; posy = 0.13;
ax = axes(	'Position',[posx,posy,width,height],...
			'box', 'on');
ax2 = axes(	'Position',[posx,posy+height,width,height],...
			'box','on');
% Sim Parameters
if isfield(para,'alpha')
	a = para.alpha; s = para.s;
else
% 	a = para.chain{nc}.alpha; s = para.chain{nc}.s;
end
% polaron plot
axes(ax); ax.UserData = 1; hold all;
hPl = handle(ax); hProp = findprop(hPl,'UserData');
n = find(tresults.star.t,1,'last');
pl1 = plot(tresults.star.t(1:n),(tresults.star.x(1:n,1,1,nc)));
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl1,'ydata',(tresults.star.x(1:n,ax.UserData,1,nc))));
if size(tresults.star.x,3) == 2
	pl2 = plot(tresults.star.t(1:n),(tresults.star.x(1:n,1,2,nc)));
	hPl.addlistener(hProp,'PostSet',@(src,event) set(pl2,'ydata',tresults.star.x(1:n,ax.UserData,2,nc)));
	if guideSH
		if s == 1
			deltaR = abs(para.hx)^(1/(1-a));
		else
			deltaR = para.hx;		% Need good formula for renormalized amplitude!
		end
		pl3 = plot(tresults.star.t(1:n),-(tresults.spin.sz(1:n)+1)./4.*s./tresults.star.omega(1),'red--');
		pl4 = plot(tresults.star.t(1:n),(1-tresults.spin.sz(1:n))./4.*s./tresults.star.omega(1),'red--');
		hPl.addlistener(hProp,'PostSet',@(src,event) set(pl3,'ydata',-(tresults.spin.sz(1:n)+1)./4.*sqrt(2*a*tresults.star.omega(ax.UserData))/(tresults.star.omega(ax.UserData)+deltaR)));
		hPl.addlistener(hProp,'PostSet',@(src,event) set(pl4,'ydata',(1-tresults.spin.sz(1:n))./4.*sqrt(2*a*tresults.star.omega(ax.UserData))/(tresults.star.omega(ax.UserData)+deltaR)));
	end
end
xlabel('Time $\omega_c t$');
ylabel('$\left<f_k\right>$');
set(gca,'ylim',[min(min(min(tresults.star.x(1:n,:,:,nc)))),max(max(max(tresults.star.x(1:n,:,:,nc))))]);
plot(ax.XLim,[0,0],'black');
% OnChange actions:
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$ \\omega_k = %g $',tresults.star.omega(ax.UserData))));
% set(gca,'ylimmode','manual');
% set(gca,'xlimmode','manual','xlim',[0,tresults.star.t(end)]);

axes(ax2); hold all;
pl3 = plot(tresults.star.t(1:n), tresults.star.n(1:n,ax.UserData,nc));
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl3,'ydata',tresults.star.n(1:n,ax.UserData,nc)));
% xlabel('Time $\omega_c t$');
ylabel('$\left<n_k\right>$');
% ax2.YAxisLocation = 'right';
ax2.XAxisLocation = 'top';

% slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmax = size(tresults.star.x,2);
sld = javax.swing.JScrollBar(0,1,1,1,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax,'UserData',round(source.Value)));

%% TDVP (3.11): Plot FT(<x>) star
% possible to restrict in time domain!
mode = 1;		% 0: lin, 1: log
f=figure(12); clf; f.Name = 'Star Polaron Fourier'; ax = gca;
maxN = floor(find(tresults.star.t,1,'last')/2)*2;
% maxN = 1000;
freq = 2*pi/tresults.star.t(2)/2 * linspace(0,1,maxN/2+1);		% fs * k/N  where k=0... N/2
FT = fft(tresults.star.x(1:maxN,:,1),maxN,1)/maxN;
if mode
	surf(tresults.star.omega,freq,log10(2*abs((FT(1:maxN/2+1,:)))));
	zlabel('$\log_{10} |\mathcal{F} \{ f_k^\uparrow \}| $');
else
	surf(tresults.star.omega,freq,2*abs(FT(1:(maxN/2+1),:)));
	zlabel('$\mathcal{F} \left\{ f_k^\uparrow(t)\right\}$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = ax.ZLabel.String;
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Frequency $\omega/ \omega_c$');
shading interp
rotate3d on
axis tight
ax.View = [0,90];
ax.YLim = [0,1.2];
if mode
% 	ax.ZLim = [-1,1];
	ax.CLim = ax.ZLim;
end
%%      (3.11a): Sin-Cos contrib
fignum = 12; fh = figure(fignum);clf; ax=gca; hold all
mode = 1;		% 0: lin, 1: log
fh.UserData = FT;
ax.UserData = 2000;
pl1 = plot(freq,(real(fh.UserData(1:maxN/2+1,ax.UserData))));
pl2 = plot(freq,(imag(fh.UserData(1:maxN/2+1,ax.UserData))));
ax.XLim = [0,1.2];

%
% pl = surf(log10(abs(x{1}.tresults.nx)));
% if mode
% 	zlabel('$log_{10}\left<n_k\right>$');
% 	pl.ZDataSource = 'log10(abs(pl.UserData{1}.tresults.nx));';
% else
% 	zlabel('$\left<n_k\right>$');
% 	pl.ZDataSource = 'abs(pl.UserData{1}.tresults.nx);';
% end
% pl.XDataSource = '1:pl.UserData{1}.para.L';
% % pl.YDataSource = 'pl.UserData{1}.para.tdvp.t';
% pl.YDataSource = 'pl.UserData{1}.tresults.t';
% xlabel('Mode $\omega_k / \omega_c$');
% ylabel('Time $\omega_c t$');
% axis tight;
% if mode
% 	set(gca,'zlim',[-4,0]);
% end
% com = ', ';
% % slider definition:
% sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
% javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
% sld.setUnitIncrement(1); sld.setBlockIncrement(1);
% hsld = handle(sld,'CallbackProperties');
% set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'userdata',res(round(source.Value),:)));
% hPl = handle(pl); hProp = findprop(hPl,'UserData');
% hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(gcf));
% hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$s = %g, \\alpha = %g$',pl.UserData{1}.para.s, pl.UserData{1}.para.alpha)));
% pl.UserData=x;

%% TDVP (3.12): Plot FT(<x>) star corrected
% possible to restrict in time domain!
mode = 1;		% 0: lin, 1: log
f=figure(12); clf; f.Name = 'Star Polaron Fourier corrected';
maxN = floor(find(tresults.star.t,1,'last')/2)*2;
% maxN = 1000;
normFact = tresults.spin.sz(ismember(tresults.t, tresults.star.t(1:maxN)),1)./2+1/2;
normFact = normFact(1:maxN)*ones(1,size(tresults.star.n,2));

freq = 2*pi/tresults.star.t(2)/2 * linspace(0,1,maxN/2+1);		% fs * k/N  where k=0... N/2
FT = fft(tresults.star.x(1:maxN,:)./normFact,maxN,1)/maxN;
if mode
	surf(tresults.star.omega,freq,log10(2*abs(FT(1:maxN/2+1,:))));
	zlabel('$sgn(\left<x_k\right>)\log_{10} \mathcal{F} \left\{ |\left<x_k\right>(t) |\right\} $');
else
	surf(tresults.star.omega,freq,2*abs(FT(1:(size(FT,1)/2+1),:)));
	zlabel('$\mathcal{F} \left\{ \left<x_k\right>(t)\right\}$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Frequency $\omega$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-1,1]);
	set(gca,'clim',get(gca,'zlim'));
end

%% TDVP (3.13): Plot Polaron-Antipolaron proportion
mode = 1;		% 0: lin, 1: log
f=figure(9); clf; f.Name = 'Star Polaron corrected';
% tresults = res{14,1}.tresults;
n = find(tresults.star.t,1,'last');
normFact = tresults.spin.sz(ismember(tresults.t, tresults.star.t(1:n)),1)./2+1/2;
normFact = normFact(1:n)*ones(1,size(tresults.star.n,2));
if mode
	surf(tresults.star.omega,tresults.star.t(1:n),log10((tresults.star.n(1:n,:)./((tresults.star.x(1:n,:)).^2).*(normFact.^2))));
% 	zlabel('$\log_{10}\left<x_k\right>$');
else
	surf(tresults.star.omega,tresults.star.t(1:n),(tresults.star.n(1:n,:)./((tresults.star.x(1:n,:)).^2).*(normFact.^2)));
% 	zlabel('$\left<x_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'xscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1]);
% 	set(gca,'clim',[-4,1]);
end

%% TDVP (3.14): Plot Current <j> CHAIN
mode = 0;		% 0: lin, 1: log
f=figure(14); clf; f.Name = 'Chain Current';
% tresults = res{1,1}.tresults;
n = tresults.lastIdx;
j = tresults.j(:,:,4);
t=tresults.t;		% for the new convention when extracting in intervals >= rev42
if mode
	surf(1:size(tresults.j,2),t(1:n),log10(abs(j(1:n,:))));
	zlabel('$\log_{10}\left<j_k\right>$');
else
	surf(1:size(tresults.j,2),t(1:n),-real(j(1:n,:)));
	zlabel('$\left<j_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = get(get(gca,'zlabel'),'String');
xlabel('Site $k$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');se
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
	set(gca,'zlim',[-4,-1.5]);
	set(gca,'clim',get(gca,'zlim'));
end
% set(gca,'Xlim',[1,15]);

%% TDVP (4) Bond Dimension Plots
%% TDVP (4.1): Plot d_opt - per chain
f=figure(411); clf; f.Name = 'OBB Dimension';
mc = 2;				% choose Chain!
if issparse(results.tdvp.d_opt)							% convert into full array & reshape
	n = find(sum(abs(results.tdvp.d_opt),2),1,'last');	%tresults.lastIdx;
	results.tdvp.d_opt = uint8(full(cumsum(results.tdvp.d_opt(1:n,:))));
	results.tdvp.d_opt = reshape(results.tdvp.d_opt,n,para.L,[]);
else
	n = size(results.tdvp.d_opt,1);
end
surf(1:para.L,para.tdvp.t(1:n),results.tdvp.d_opt(:,:,mc))
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$d_{opt}$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
%% TDVP (4.1): Plot d_opt - all chains
f=figure(411); clf; f.Name = 'OBB Dimension';
if issparse(results.tdvp.d_opt)							% convert into full array & reshape
	n = find(sum(abs(results.tdvp.d_opt),2),1,'last');	%tresults.lastIdx;
else
	n = size(results.tdvp.d_opt,1);
end
surf(1:size(results.tdvp.d_opt,2),para.tdvp.t(1:n),cumsum(results.tdvp.d_opt(1:n,:)))
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$d_{opt}$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight

%% TDVP (4.2): Plot D - per chain
f=figure(421); clf; f.Name = 'Bond Dimension';
mc = 3;				% choose Chain!
if issparse(results.tdvp.D)						% convert into full array & reshape
	n = find(sum(abs(results.tdvp.D),2),1,'last'); %tresults.lastIdx;
	results.tdvp.D = uint8(full(cumsum(results.tdvp.D(1:n,:))));
	results.tdvp.D = reshape(results.tdvp.D,n,para.L-1,[]);
else
	n = size(results.tdvp.D,1);
end
surf(1:para.L-1,para.tdvp.t(1:n),results.tdvp.D(:,:,mc))
colorbar
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$D$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
%% TDVP (4.2): Plot D - all chains
f=figure(421); clf; f.Name = 'Bond Dimension';
if issparse(results.tdvp.D)
	n = find(sum(abs(results.tdvp.D),2),1,'last'); %tresults.lastIdx;
else
	n = size(results.tdvp.D,1);
end
surf(1:size(results.tdvp.D,2),para.tdvp.t(1:n),cumsum(results.tdvp.D(1:n,:)))
colorbar
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$D$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight

%% TDVP (4.3): Plot d_k
f=figure(431); clf; f.Name = 'Local Dimension';
mc = 1;				% choose Chain!
if issparse(results.tdvp.dk)						% convert into full array & reshape
	n = find(sum(abs(results.tdvp.dk),2),1,'last'); %tresults.lastIdx;
	results.tdvp.dk = uint8(full(cumsum(results.tdvp.dk(1:n,:))));
	results.tdvp.dk = reshape(results.tdvp.dk,n,para.L-1,[]);
else
	n = size(results.tdvp.dk,1);
end
surf(1:para.L*para.nChains,para.tdvp.t(1:n),results.tdvp.dk(:,:,mc))

xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$D$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
%% TDVP (6): Plot vNE of A / V
ii = 4;
figure(600+ii); clf;
% x = res{1+ii,1}; results = x.results;
mc = 4;
if para.tdvp.logSV
	results.tdvp.Amat_vNE = cell2mat(cellfun(@(x) sum(-x.^2.*log(x.^2)), results.tdvp.Amat_sv, 'UniformOutput',false));
	results.tdvp.Vmat_vNE = cell2mat(cellfun(@(x) sum(-x.^2.*log(x.^2)), results.tdvp.Vmat_sv, 'UniformOutput',false));
end
lastIdx = find(results.tdvp.Amat_vNE(:,2,mc),1,'last');
subplot(1,2,1);
surf(1:size(results.tdvp.Amat_vNE,2),para.tdvp.t(1:lastIdx),results.tdvp.Amat_vNE(1:lastIdx,:,mc));
% contour3(1:size(results.tdvp.Amat_vNE,2),para.tdvp.t(1:size(results.tdvp.Amat_vNE,1)),results.tdvp.Amat_vNE,40,'k');
xlabel('Bond $k$'); ylabel('Time $\omega_c t$'); zlabel('$S_{vNE}(A)$');
set(gca,'View',[90 0]);  %top: 0 90
shading interp
rotate3d on
axis tight
subplot(1,2,2);
surf(1:size(results.tdvp.Vmat_vNE,2),para.tdvp.t(1:lastIdx),results.tdvp.Vmat_vNE(1:lastIdx,:,mc));
xlabel('Bond $k$'); ylabel('Time $\omega_c t$'); zlabel('$S_{vNE}(V)$');
set(gca,'View',[90 0]);
shading interp
rotate3d on
axis tight
% formatPlot(7)

%% TDVP (6.1): Vmat, use of dk Slider in L for Multi-Chain!
nc = 5;		% which chain?
ii = 0;
% x = res{10+ii,1}; Vmat = x.Vmat; results = x.results;
figure(610+ii); clf; ax = {};
ax{1} = gca;
hPl = handle(ax{1}); hProp = findprop(hPl,'UserData');
pl1 = surf(abs(Vmat{2}{nc}));	% plot resets UserData!
ax{1}.UserData = 2;
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl1,'zdata',abs(Vmat{ax{1}.UserData}{nc})));
xlabel('$d_{OBB}$');
ylabel('$d_k$');
ax{1}.View = [0,90];
shading interp; rotate3d on; axis tight;
origHeight = ax{1}.Position(4);
ax{1}.Position(4) = origHeight*0.8;
ax{2} = axes('Position',[ax{1}.Position(1:2)+[0,ax{1}.Position(4)], abs(ax{1}.Position(3:4)-[0,origHeight])]);
pl2 = plot(log10(results.Vmat_sv{nc,ax{1}.UserData})); axis tight
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl2,'ydata',log10(results.Vmat_sv{nc,ax{1}.UserData})));
ax{2}.XTickLabel = [];
% ax{2}.YLim=ax{2}.YLim;
ax{2}.XLim = ax{1}.XLim;
hPl.addlistener(hProp,'PostSet',@(src,event) set(ax{2},'XLim',ax{1}.XLim));
ylabel('$log_{10}(SV(V))$');
t = text(mean(ax{2}.XLim),ax{2}.YLim(1)/2, 'k = 2');
hPl.addlistener(hProp,'PostSet',@(src,event) set(t,'String',sprintf('k = %d',ax{1}.UserData)));

% slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmin = 2; sldmax = para.L+1;
sld = javax.swing.JScrollBar(0,sldmin,1,sldmin,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax{1},'UserData',round(source.Value)));

%% TDVP (6.2): Vmat, use of dk Slider in L for Single-Chain!
nc = 1;
ii = 0;
% x = res{10+ii,1}; Vmat = x.Vmat; results = x.results;
figure(620+ii); clf; ax = {};
ax{1} = gca;
hPl = handle(ax{1}); hProp = findprop(hPl,'UserData');
pl1 = surf(abs(Vmat{2}));	% plot resets UserData!
ax{1}.UserData = 2;
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl1,'zdata',abs(Vmat{ax{1}.UserData})));
xlabel('$d_{OBB}$');
ylabel('$d_k$');
ax{1}.View = [0,90];
shading interp; rotate3d on; axis tight;
origHeight = ax{1}.Position(4);
ax{1}.Position(4) = origHeight*0.8;
ax{2} = axes('Position',[ax{1}.Position(1:2)+[0,ax{1}.Position(4)], abs(ax{1}.Position(3:4)-[0,origHeight])]);
pl2 = plot(log10(results.Vmat_sv{nc,ax{1}.UserData})); axis tight
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl2,'ydata',log10(results.Vmat_sv{nc,ax{1}.UserData})));
% ax{2}.XTick = [];
ax{2}.YLim=ax{2}.YLim; ax{2}.XLim = ax{1}.XLim;
hPl.addlistener(hProp,'PostSet',@(src,event) set(ax{2},'XLim',ax{1}.XLim));
ylabel('$log_{10}(SV(V))$');
t = text(mean(ax{2}.XLim),ax{2}.YLim(1)/2, 'k = 2');
hPl.addlistener(hProp,'PostSet',@(src,event) set(t,'String',sprintf('k = %d',ax{1}.UserData)));

% slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmin = 2; sldmax = para.L;
sld = javax.swing.JScrollBar(0,sldmin,1,sldmin,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax{1},'UserData',round(source.Value)));

%% TDVP (6.3): Amat, use of D/dOBB Slider in L for Single-Chain! Vtens!
nc = 1;
ii = 0;
% x = res{10+ii,1}; Vmat = x.Vmat; results = x.results;
figure(630+ii); clf; ax = {};
% surface plot
ax{1} = gca;
hPl = handle(ax{1}); hProp = findprop(hPl,'UserData');
pl1 = surf(abs(reshape(mps{3},para.D(2),[]).'));	% plot resets UserData!
ax{1}.UserData = 1;
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl1,'zdata',abs(reshape(mps{ax{1}.UserData},para.D(ax{1}.UserData-1),[]).')));
xlabel('$D_{l}$');
ylabel('$D_r \times d_{OBB}$');
ax{1}.View = [0,90];
shading interp; rotate3d on; axis tight;
origHeight = ax{1}.Position(4);
ax{1}.Position(4) = origHeight*0.8;
% SV plot
ax{2} = axes('Position',[ax{1}.Position(1:2)+[0,ax{1}.Position(4)], abs(ax{1}.Position(3:4)-[0,origHeight])]);
pl2 = plot(log10(results.Amat_sv{1,ax{1}.UserData-1}));	axis tight
hPl.addlistener(hProp,'PostSet',@(src,event) set(pl2,'ydata',log10(results.Amat_sv{1,ax{1}.UserData-1})));
ax{2}.XAxisLocation = 'top';
ax{2}.XLim = ax{1}.XLim;
hPl.addlistener(hProp,'PostSet',@(src,event) set(ax{2},'XLim',ax{1}.XLim));
ylabel('$log_{10}(SV(V))$');
t = text(mean(ax{2}.XLim),ax{2}.YLim(1)/2, sprintf('$k = 3; D_l = %d, D_r = %d, d_{OBB} = %d$', para.D(2),para.D(3),para.d_opt(3)));
hPl.addlistener(hProp,'PostSet',@(src,event) set(t,'String', ...
	sprintf('$k = %d; D_l = %d, D_r = %d, d_{OBB} = %d$',ax{1}.UserData, para.D(ax{1}.UserData-1),para.D(ax{1}.UserData),para.d_opt(ax{1}.UserData))));

% slider definition and UserData setting:
f = gcf; pos = f.Position;
sldmin = 2; sldmax = para.L;
sld = javax.swing.JScrollBar(0,sldmin,1,sldmin,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax{1},'UserData',round(source.Value)));

%% TDVP: Plot temporal change in Vmat SV
figure(4); clf;
ax = axes('units','pixels');
pl = surf(cell2mat(results.tVmat_sv(1,:)));
set(gca,'View',[150 0]);
set(gca,'zscale','log');
set(gca,'zlimmode','manual');
set(gca,'zlim',[1e-15,1]);
xlabel('Site $k$');
ylabel('OBB Dimension');
% shading interp
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(tmps,1),'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'zdata',cell2mat(results.tVmat_sv(round(source.Value),:))));

%% TDVP: Plot Min of Amat SV over time
figure(5);hold on;
minSV = cell2mat(cellfun(@min, results.tdvp.Amat_sv, 'UniformOutput',false));
surf(1:(para.L-1),para.tdvp.t(1:size(minSV,1)),log10(minSV))
% minSV = cell2mat(cellfun(@(x) log10(x(end-2:end)), results.Amat_sv(:,2:end),'UniformOutput',false));
% plot(minSV')
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
shading interp
rotate3d on
set(gca,'zlim',[-10,0]);
%% TDVP: Plot Min of Vmat SV over time
figure(5);clf;
% minSV = cell2mat(cellfun(@(x) log10(min(x)), results.tdvp.Vmat_sv(:,2:end),'UniformOutput',false));
minSV = cell2mat(cellfun(@(x) log10(x(1:3)), results.Vmat_sv(:,2:end),'UniformOutput',false));
plot(minSV')
% surf(1:(para.L-1),para.tdvp.t(1:size(minSV,1)),minSV)
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
shading interp
rotate3d on
% set(gca,'zlim',[-10,0]);

%% TDVP: Plot Sum over Vmat in 3D - analyse dk expand
figure(6);
rotate3d off
zeroVals = -50;             % value of zeros for padding and -inf replacement
plotMatTime = cell(size(tmps,1),1);
for j = 1:size(tmps,1)      % timeslices
    Vmatj = tVmat(j,:);
    Vmatj_sv = results.tdvp.Vmat_sv(j,:);
    plotMat = [];
for i = 2:length(Vmatj)
    if size(Vmatj{1,i},2) < size(Vmatj_sv{1,i},1)
        a = log10(sum(abs(real(Vmatj{1,i}*diag(Vmatj_sv{1,i}(1:size(Vmatj{1,i},2),:)))),2));
    else
        a = log10(sum(abs(real(Vmatj{1,i}(:,1:size(Vmatj_sv{1,i},1))*diag(Vmatj_sv{1,i}))),2));
    end
    a(a==-inf)=zeroVals;
    %a = log10(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2));
    dim = max(length(a),size(plotMat,1));
    if length(a)< dim
        a = padarray(a,dim-length(a),zeroVals,'pre');
    elseif size(plotMat,1) < dim
        plotMat = padarray(plotMat,dim-size(plotMat,1),zeroVals,'pre');
    end
    plotMat = [plotMat,a];
end
plotMatTime{j} = plotMat;
end
%%
pl = surf(plotMatTime{2,1});
title('Temporal change of OBB')
ylabel('$d_k$')
xlabel('Site $k$')
set(gca,'View',[9.5 40]);
formatPlot(6)
axis tight
rotate3d on
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(plotMatTime,1),'Value',1,...
        'Position', [1100 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'zdata',plotMatTime{round(source.Value)}));

%% TDVP (7): System Density Matrix

%% TDVP (7.01): RDM Diagonal - TDVPData
f=figure(701); clf; f.Name = 'RDM diagonal'; ax = gca;
% x = res{15,1};
% ph = x.plot('rhoii',f);
ph = x.plot('rhoii','-fsev',f);
leg = legend('TT','LE+','LE-','CT+','CT-');
grid on; axis tight;
formatPlot(f,'twocolumn-single')
%% TDVP (7.02): RDM Off-Diagonal imag - TDVPData
f=figure(702); clf; f.Name = 'RDM off-diagonal imag'; ax = gca;
% x = res{15,1}; tresults = x.tresults; para = x.para;
ph = x.plot('rhoij-imag',f);
legend toggle;
% ax.YScale = 'log';
%% TDVP (7.03): RDM Off-Diagonal real - TDVPData
f=figure(703); clf; f.Name = 'RDM off-diagonal real'; ax = gca;
% x = res{15,1}; tresults = x.tresults; para = x.para;
ph = x.plot('rhoij-real',f);
legend toggle;
% ax.YScale = 'log';
%% TDVP (7.04): RDM Off-Diagonal abs  - TDVPData
f=figure(704); clf; f.Name = 'RDM off-diagonal abs'; ax = gca;
% x = res{15,1}; tresults = x.tresults; para = x.para;
ph = x.plot('rhoij-abs',f);
legend toggle;
% ax.YScale = 'log';
%% TDVP (7.05): RDM DPMES Overview - TDVPData
f=figure(705); clf; f.Name = 'RDM DPMES all'; ax = gca;
% x = res{15,1}; tresults = x.tresults; para = x.para;
ph = x.plot('rho-dpmes','-fsev',f);
formatPlot(f,'twocolumn-single');

%% TDVP (7.1): RDM Fourier Trafos
	%% TDVPData: DFT
f=figure(706);  f.Name = 'MLSBM Wave-Occupation DFT'; clf; hold all; ax = gca;
% x = res(6);
% pl = x.plot('rhoii'); xlabel('$t$');
% pl = x.plot('rhoii-ft','-cmev',f); xlabel('$E/cm^{-1}$');
x.plotSld1DFT('rhoii-ft','-cmev',f);
% leg = legend('TT','LE+','CT+','CT-');
leg = legend('TT','LE+','LE-','CT+','CT-');
ylabel('$|FT(\rho_{ii})|^2$');
grid on
axis tight
% ax.XLim = [0,9000];ax.YLim = [0,2];
ax.XLim = [0,2000];ax.YLim = [0,200];
formatPlot(f,'twocolumn-single')
	%% TDVPData: DFT Fit residuals
f=figure(707);  f.Name = 'Rhoii residual DFT'; clf; hold all; ax = gca;
% x.plotSld1DFT('rhoii-osc-res','-cmev',f); ax=gca; f=gcf;
ph=x.plot('rhoii-osc-res','-cmev',f); ax=gca; f=gcf;
leg = legend('TT','LE+','LE-','CT+','CT-');
ylabel('$|FT(\rho_{ii})|^2$');
grid on
axis tight
% ax.XLim = [0,9000];ax.YLim = [0,2];
ax.XLim = [0,2000];ax.YLim = [0,200];
formatPlot(f,'twocolumn-single')

for kk = 1:numel(ph)
	ph(kk).YData = ph(kk).YData./max(ph(kk).YData)+5-kk;
	ax.YLim = [0,5];
end

	%% TDVPData: DFT Smooth residuals norm dist SLD
f=figure(713);  f.Name = 'Rhoii residual DFT'; clf; hold all; ax = gca;
x = res(2);
x.plotSld1DFT('rhoii','-cmev','-smoothres','-norm','-dist',f);
leg = legend('TT','LE+','LE-','CT+','CT-');
ylabel('$|FT(\rho_{ii})|^2$');
grid on
axis tight
% ax.XLim = [0,9000];ax.YLim = [0,2];
ax.XLim = [0,2000];
formatPlot(f,'twocolumn-single')

	%% TDVPData: DFT Smooth residuals dist
f=figure(714);  f.Name = 'Rhoii residual DFT'; clf; hold all; ax = gca;
x = res(1);
h = x.plot('rhoii-osc-res-med','-cmev','-resetColorOrder','-dist',f); ax=gca; f=gcf;
set(h.pl,{'Displayname'},{'TT','LE+','LE-','CT+','CT-'}')
arrayfun(@(x) legend(x,'show'),h.ax);
ylabel('$|FT(\rho_{ii})|^2$');
% ax.XLim = [0,9000];ax.YLim = [0,2];

[h.ax.XLim] = deal([0,2000]);
formatPlot(f,'twocolumn-single')

	%% TDVPData: DFT Smooth residuals norm dist
f=figure(715);  f.Name = 'Rhoii residual DFT'; clf; hold all; ax = gca;
% x = res(2);
ph = x.plot('rhoii-osc-res-med','-cmev','-resetColorOrder','-norm','-dist',f); ax=gca; f=gcf;
leg = legend('TT','LE+','LE-','CT+','CT-');
ylabel('$|FT(\rho_{ii})|^2$');
grid on
axis tight
% ax.XLim = [0,9000];ax.YLim = [0,2];
ax.XLim = [0,2000];
formatPlot(f,'twocolumn-single')

	%% TDVPData: DFT Smooth residuals norm 1272
f=figure(711);  f.Name = 'Rhoii residual DFT'; clf;  hold all; ax = gca;
% 	a = x.getData('rhoii-osc-res');
% x.plotSld1DFT('rhoii-osc-res','-cmev',f);
h = x.plot('rhoii-osc-res-med','-cmev','-resetColorOrder',f);
ph = h.pl;
if numel(ph) ~= length(ph)
	set(ph(2,:),'LineStyle','-.');
end
% peak range:
rng = find(ph(1).XData > 1260 &ph(1).XData < 1290);
% for kk = 1:numel(ph)
% 	ph(kk).YData = ph(kk).YData./max(ph(kk).YData(rng));
% end
% leg = legend('TT','LE+','LE-','CT+','CT-');
ylabel('$|FT(\rho_{ii})|^2$');
grid on
axis tight
% ax.XLim = [0,9000];ax.YLim = [0,2];
ax.XLim = [0,2000];
formatPlot(f,'twocolumn-single')
	
%% TDVP (7.3) MLSBM HsHi
f=figure(732);  f.Name = 'DPMES HsHi'; clf; hold all; ax = gca;
% x = res(6);
% pl = x.plot('hshi');xlabel('$t$')
pl = x.plot('hshi','-fsev');xlabel('$t/fs$')
% leg = legend('$H_S$','$H_I(1)$','$H_I(2)$','$H_I(3)$','$H_I(4)$','$H_I(5)$','$H_I(6)$','$H_I(7)$','$H_S+\sum H_I$');
leg = legend('$H_S$','$H_I(A_{1,0})$','$H_I(A_{1,2})$','$H_I(A_2)$','$H_I(B_1)$','$H_I(B_{1,0})$','$H_I(B_{1,x})$','$H_I(B_{1,-x})$','$H_S+\sum H_I$');
ylabel('$\left<H\right>/eV$');
grid on
axis tight
% set(gca,'yscale','log')
formatPlot(f,'twocolumn-single')

%% TDVP (7.4) StateProjection - TDVPData
f=figure(741);  f.Name = 'State Projection Amplitude'; clf; hold all; ax = gca;
% x = res(6);
% pl = x.plot('stateproj'); xlabel('$t$');
pl = x.plot('stateproj','-fsev'); xlabel('$t/fs$');
leg = legend('abs','real','imag');
ylabel('$\langle LE^+|\left< 0|\Psi\right>$');
grid on
axis tight
formatPlot(f,'twocolumn-single')

%% TDVP (7.5) Linear Absorption - TDVPData
f=figure(753);  f.Name = 'Linear Absorption of LE+';hold all; ax = gca;
% x = res(6);
% pl = x.plot('linabs','-fsev'); xlabel('$E/eV$'); ax.XLim = [1.5,2.5];
% pl = x.plot('linabs','-nmev'); xlabel('$\lambda/nm$'); ax.XLim = [400, 1000];
pl = x.plotSld1DFT('linabs','-nmev'); xlabel('$E/eV$'); ax = gca; f = gcf; ax.XLim = [400,1000];
ylabel('$FT(\langle TT|\left< 0|\Psi_{TT}\right>)$');
grid on
% ax.XLim = [1000,100000];
formatPlot(f,'twocolumn-single')

%% TDVP (7.6) Adiabatic States - 2D TDVPData
% x = res(6);
h = x.plotSld2D('state-adiab');
% h.ax.ZLim = [0,1]; h.ax.CLim = [0,1];
set(h.ax,'TickLabelInterpreter', 'latex');
h.ax.XTick = 1:5;
h.ax.XTickLabel = {'$TT$','$LE^+$','$LE^-$','$CT^+$','$CT^-$'};

%% TDVP (7.6.1) Adiabatic States rearranged - 2D TDVPData
% especially done for: 20160208-1904-27-DPMES5-7C-Star-v72-dk20D5dopt5L8Delta0State2
% could work in others as well!
% x = res(6);
x1 = x;
r = [0,37.5,894.5,x.lastIdx-1]+1; from = [3,2,1];								% adiab state 1
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,1) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,894.5,2773.5,x.lastIdx-1]+1; from = [1,2,3];								% adiab state 2
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,2) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,122.5,294.5,327.5,356.5,2773.5,x.lastIdx-1]+1; from = [5,4,3,4,3,2];		% adiab state 3
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,3) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,37.5,294.5,327.5,356.5,807.5,1046.5,4661.5,x.lastIdx-1]+1; from = [2,3,4,3,4,5,4,5];		% adiab state 4
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,4) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,122.5,807.5,1046.5,4661.5,x.lastIdx-1]+1; from = [4,5,4,5,4];		% adiab state 4
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,5) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end

h = x1.plotSld2D('state-adiab');
% h.ax.ZLim = [0,1]; h.ax.CLim = [0,1];
% h.ax.ZLim = [-5,0]; h.ax.CLim = h.ax.ZLim;
set(h.ax,'TickLabelInterpreter', 'latex');
h.ax.XTick = 1:5;
h.ax.XTickLabel = {'$TT$','$LE^+$','$LE^-$','$CT^+$','$CT^-$'};

%% TDVP (7.6.2) Adiabatic States rearranged - 1D TDVPData
% especially done for: 20160208-1904-27-DPMES5-7C-Star-v72-dk20D5dopt5L8Delta0State2
% could work in others as well!
% x = res(6);
x1 = x;
r = [0,37.5,894.5,x.lastIdx-1]+1; from = [3,2,1];								% adiab state 1
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,1) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,894.5,2773.5,x.lastIdx-1]+1; from = [1,2,3];								% adiab state 2
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,2) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,122.5,294.5,327.5,356.5,2773.5,x.lastIdx-1]+1; from = [5,4,3,4,3,2];		% adiab state 3
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,3) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,37.5,294.5,327.5,356.5,807.5,1046.5,4661.5,x.lastIdx-1]+1; from = [2,3,4,3,4,5,4,5];		% adiab state 4
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,4) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end
r = [0,122.5,807.5,1046.5,4661.5,x.lastIdx-1]+1; from = [4,5,4,5,4];		% adiab state 4
for ii = 1:length(from)
	x1.sysState(ceil(r(ii)):floor(r(ii+1)),:,5) = x.sysState(ceil(r(ii)):floor(r(ii+1)),:,from(ii));
end

h = x1.plotSld1D('state-adiab');
leg = legend('TT','LE+','LE-','CT+','CT-');
grid on
axis tight
formatPlot(gcf,'twocolumn-single')


%% TDVP (8.1) TES/PES Energy Surfaces
x = res(1);
h = x.plot('heff-full-diab-v2','-fsev',figure(20)); h.pl(end).SizeData = 10; resizePlot(h.f,2.2,2);
% h = x.setHeffTo('pes').plot('heff-full-diab-v2','-fsev',figure(21),'-ylim',[1.3,3.4]); h.pl(end).SizeData = 10; resizePlot(h.f,2.2,2);
% h = x.plot('heff-full-pop-diab','-xlim',[0,40]); ylabel('$E$');h.pl(end).SizeData = 50;
% [h.pl(1:end-1).FaceAlpha] = deal(0.5);
% h.ax.Color = 'none'
% box on;
% xlim([0,420]);
% h = x.setHeffTo('pes').plot('heff-full-pop-diab','-xlim',[0,40]); ylabel('$E$');h.pl(end).SizeData = 50

%% TDVP (8.1) TES/PES Energy Surfaces SBM
% x = res(8);
% xlimits = [0,40]; arrayfun(@(x) clf(x),1:8);
ylimits = [-0.5,0.5];
h = x.setHeffTo('pes').plot('heff-pop-diab','-state',1,'-xlim',xlimits, figure(1)); title('PES (1)'); xlabel('$t$'); ylabel('$E$'); resizePlot(h.f,1.3);
h = x.setHeffTo('pes').plot('heff-pop-diab','-state',2,'-xlim',xlimits, figure(2)); title('PES (2)'); xlabel('$t$'); ylabel('$E$'); resizePlot(h.f,1.3);
h = x.setHeffTo('tes').plot('heff-pop-diab','-state',1,'-xlim',xlimits, figure(3)); title('TES (1)'); xlabel('$t$'); ylabel('$E$'); resizePlot(h.f,1.3);
h = x.setHeffTo('tes').plot('heff-pop-diab','-state',2,'-xlim',xlimits, figure(4)); title('TES (2)'); xlabel('$t$'); ylabel('$E$'); resizePlot(h.f,1.3);

% h = x.setHeffTo('tes').plot('heff-full-pop-diab','-xlim',xlimits, figure(5)); xlabel('$t$'); ylabel('$E$'); resizePlot(h.f,1.3);
% h = x.setHeffTo('pes').plot('heff-full-pop-diab','-xlim',xlimits, figure(6)); ylabel('$E$'); ylabel('$E$'); resizePlot(h.f,1.3);
h = x.setHeffTo('tes').plot('heff-full-pop-diab-v2','-xlim',xlimits, figure(5)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);
h = x.setHeffTo('pes').plot('heff-full-pop-diab-v2','-xlim',xlimits, figure(6)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);
% h = x.setHeffTo('tes').plot('heff-full-pop-diab-v2','-xlim',xlimits,'-ylim',ylimits, figure(5)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,2.6,2);
% h = x.setHeffTo('pes').plot('heff-full-pop-diab-v2','-xlim',xlimits,'-ylim',ylimits, figure(6)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,2.6,2);

% h = x.plot('tes-from-pes-pop-diab-v2','-xlim',xlimits, figure(7)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);
% h = x.plot('pes-from-tes-pop-diab-v2','-xlim',xlimits, figure(8)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);
h = x.plot('tes-from-pes-pop-diab-v2','-xlim',xlimits,'-ylim',ylimits, figure(7)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);
h = x.plot('pes-from-tes-pop-diab-v2','-xlim',xlimits,'-ylim',ylimits, figure(8)); ylabel('$E$');h.pl(end).SizeData = 50; resizePlot(h.f,1.3);

%% TDVP (8.2) j Currents
x = res(1);
h = x.plot('heff-current-nonzero','-fsev');
legend toggle;
xlim([0,420]);
TDVPData.resetColorOrder(h.ax);
%% TDVP (8.2.1) sum(j) Currents cumsum
x = res(1);
h = x.plot('heff-current-cumsum','-fsev');
delete(h.pl(~sum(abs(h.ydata),1)));
legend toggle;
xlim([0,420]);
TDVPData.resetColorOrder(h.ax);
%% export
export_fig(['img/test'],'-transparent','-painters','-png','-m4')

%% TDVP (8.3) System States: PES to environment <n> mapping
% try to find the environmental state corresponding to a specific surface in the TES/PES
x = res(1);
tes = x.getData('heff-full-pop');
popSurf = tes{2};
% tes{1} = E;																	% t x D*dk_eig			Energies
% tes{2} = pop;																	% t x D*dk_eig			population on surface
% tes{3} = popDiab;																% t x dk x D*dk_eig
% tes{4} = Vtemp;																% t x D*dk x D*dk_eig	map to surface
tIdx = find(x.t*0.658>20,1);		% want to look at 80fs.
eigIdx = find(popSurf(tIdx,:)== max(popSurf(tIdx,:)),1);						% largest population on surface

d = size(x.sysState);															% t x dk x D
stateMix = reshape(tes{4}(tIdx,:,eigIdx), d(3), d(2));							% D x dk
stateMixPop = stateMix.*conj(stateMix)
% TODO: need to continue!

%% TDVP (8.3) System States: ES breakdown
% try to plot the diabatic states corresponding to a specific surface in the TES/PES
x = res(1);
es = x.setHeffTo('tes').getData('heff-full-pop');
popSurf = es{2};
Vtemp = es{4};
% tes{1} = E;																	% t x D*dk_eig			Energies
% tes{2} = pop;																	% t x D*dk_eig			population on surface
% tes{3} = popDiab;																% t x dk x D*dk_eig
% tes{4} = Vtemp;																% t x D*dk x D*dk_eig	map to surface
t_fs = 80;								% want to look at 80fs.
tIdx = find(x.t*0.658>t_fs,1);
eigIdx = find(popSurf(tIdx,:)== max(popSurf(tIdx,:)),1);						% largest population on surface

% imagesc(abs(squeeze(Vtemp(tIdx
% TODO not fully done yet
%% TDVP z-averaging in one file
% naming scheme to find files:
%   take series filename and replace z-value by *
clear;
folder = '20150225-1826-SpinBoson-LogZ-v40TCMde9-alpha1.5delta0.1Lambda1.1dk20D5dopt5L100-artificial';
figure(7);clf;
% folder = '20141117-0531-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L84';
% folder = '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49';
filescheme = 'results-Till1000Step2v40-OBBExpand-noBondExpand-expvCustom800-1core-z*-small.mat';
files = ls(sprintf('%s/%s',folder,filescheme));
PlotData.spin.sz = [];PlotData.spin.sx = []; PlotData.spin.sy = []; PlotData.nx = [];
PlotData.z = [];
for k = 1:size(files,1)
    load([folder,'/',files(k,:)],'para','tresults');
	PlotData.spin.sx(k,:) = tresults.spin.sx;
	PlotData.spin.sy(k,:) = tresults.spin.sy;
    PlotData.spin.sz(k,:) = tresults.spin.sz;
	PlotData.nx(:,:,k)    = tresults.nx;
    PlotData.z(k) = para.z;
end
PlotData.t = para.tdvp.t;
plot(PlotData.t,[PlotData.spin.sz;mean(PlotData.spin.sz)]);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData.z)); legLabels{end} = 'z-Ave';
legend(legLabels);
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
para.tdvp.filename = [para.tdvp.filename(1:end-5),'Avg.mat'];
tresults.spin.sx = mean(PlotData.spin.sx,1);
tresults.spin.sy = mean(PlotData.spin.sy,1);
tresults.spin.sz = mean(PlotData.spin.sz,1);
tresults.nx = mean(PlotData.nx,3);
if ~exist(para.tdvp.filename,'file')
	save(para.tdvp.filename, 'para','tresults');
end

%% Plot only z-averaged sz
plot(PlotData.t,mean(PlotData.spin.sz));
ylim([-1,1]);set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');

%% TDVP (8) z-averaging create Orth 2010
% naming scheme to find files:
%   take series filename and replace z-value by *
figure(8);clf;
folder = {'20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49',...
          '20141114-2019-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0405-SpinBoson-alpha0.1delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0406-SpinBoson-alpha0.15delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49'};
filescheme = 'results-Till325Step4*-OBBExpand-noBondExpand*.mat';

if ~exist('PlotData8','var')
    PlotData8.spin.meanSz = [];
    PlotData8.alpha = [];
end
if isempty(PlotData8.spin.meanSz)
    for k = 1:length(folder)
        PlotData8.spin.sz = [];
        files = dir(sprintf('%s/%s',folder{k},filescheme));
        for l = 1:length(files)
            load([folder{k},'/',files(l).name]);
            PlotData8.spin.sz(l,:) = tresults.spin.sz;
        end
        PlotData8.alpha(k) = para.alpha;
        PlotData8.spin.meanSz(k,:) = mean(PlotData8.spin.sz);
    end
end

PlotData8.t = para.tdvp.t;
plot(PlotData8.t,PlotData8.spin.meanSz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData8.alpha)); legLabels{end} = 'z-Ave';
legend(legLabels);
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');

%% TDVP (9) Orthogonal Polynomials L = 50: Orth 2010
figN = 10;
figure(figN); clf;
folder50 = {'20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50',...
            '20141114-1617-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L50',...
            '20141117-0641-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L50',...
            '20141117-0642-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L50',...
            '20141116-0229-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L50'};
if ~exist('PlotData9','var')
    PlotData9.sz = [];
    PlotData9.alpha = [];
end
if isempty(PlotData9.sz)
    for k = 1:length(folder50)
		% No OBB No Bond Expand:
%         filename = dir(sprintf('%s/results-Till325Step4-noOBBExpand-noBondExpand.mat',folder50{k}));
		% OBB and Bond Expand, max Bond 20;
		filename = dir(sprintf('%s/results-Till325Step4-OBBandBondExpand20*.mat',folder50{k}));
		try
	        load(sprintf('%s/%s',folder50{k},filename.name),'tresults','para');
		catch
			continue
		end
        PlotData9.sz(k,1:length(tresults.spin.sz)) = tresults.spin.sz;
        PlotData9.alpha(k) = para.alpha;
    end
    PlotData9.t = para.tdvp.t;
end
plot(PlotData9.t,PlotData9.sz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData9.alpha));
legend(legLabels(1:end-1));
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
formatPlot(figN);

%% TDVP (10) Orthogonal Polynomials L = 200: Orth 2010
figN = 10;
figure(figN); clf;
folder50 = {'20141221-0148-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0148-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0151-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0151-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0153-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L200'};
if ~exist('PlotData9','var')
    PlotData9.sz = [];
    PlotData9.alpha = [];
end
if isempty(PlotData9.sz)
    for k = 1:length(folder50)
		% OBB and Bond Expand, max Bond 20;
		filename = dir(sprintf('%s/results-Till325Step4-OBB*Expand20*.mat',folder50{k}))
		try
	        load(sprintf('%s/%s',folder50{k},filename.name),'tresults','para');
		catch
			continue
		end
        PlotData9.sz(k,1:length(tresults.spin.sz)) = tresults.spin.sz;
        PlotData9.alpha(k) = para.alpha;
    end
    PlotData9.t = para.tdvp.t;
end
plot(PlotData9.t,PlotData9.sz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData9.alpha));
legend(legLabels(1:end-1));
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
formatPlot(figN);

%% TDVP (11.1) expvCustom Benchmarking: evolve V
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(9);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),results.tdvp.expvTime(:,1)+results.tdvp.expvTime(:,4),'+','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),results.tdvp.expvTime(:,2)+results.tdvp.expvTime(:,4),'*','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),results.tdvp.expvTime(:,3),'o','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(9,'twocolumn-single')
%% TDVP (11.1.1) expvCustom Benchmarking 3D: evolve V
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(9);clf;
hold all
scatter3(results.tdvp.expvTime(:,19),results.tdvp.expvTime(:,20),results.tdvp.expvTime(:,1)+results.tdvp.expvTime(:,4),'+','DisplayName','expM');
scatter3(results.tdvp.expvTime(:,19),results.tdvp.expvTime(:,20),results.tdvp.expvTime(:,2)+results.tdvp.expvTime(:,4),'*','DisplayName','expV');
scatter3(results.tdvp.expvTime(:,19),results.tdvp.expvTime(:,20),results.tdvp.expvTime(:,3),'o','DisplayName','expvCustom');
set(gca,'zscale','log')
% set(gca,'xscale','log')
% set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
axis tight
formatPlot(9,'twocolumn-single')

%% TDVP (11.2) expvCustom Benchmarking: evolve CV
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(10);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,[19,19]),2),results.tdvp.expvTime(:,5)+results.tdvp.expvTime(:,8),'+','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,[19,19]),2),results.tdvp.expvTime(:,6)+results.tdvp.expvTime(:,8),'*','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,[19,19]),2),results.tdvp.expvTime(:,7),'o','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(10,'twocolumn-single')
%% TDVP (11.3) expvCustom Benchmarking: evolve A
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(11);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,[17,18,19]),2),results.tdvp.expvTime(:, 9)+results.tdvp.expvTime(:,12),'+','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,[17,18,19]),2),results.tdvp.expvTime(:,10)+results.tdvp.expvTime(:,12),'*','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,[17,18,19]),2),results.tdvp.expvTime(:,11),'o','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(11,'twocolumn-single')

%% TDVP (11.4) expvCustom Benchmarking: evolve CA
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(12);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,[17,18]),2),results.tdvp.expvTime(:,13)+results.tdvp.expvTime(:,16),'+','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,[17,18]),2),results.tdvp.expvTime(:,14)+results.tdvp.expvTime(:,16),'*','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,[17,18]),2),results.tdvp.expvTime(:,15),'o','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
% set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(12,'twocolumn-single')

%% TDVP (11.5) expvCustom Benchmarking: evolve V+CV
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(15);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),sum(results.tdvp.expvTime(:,[1,4,5,8]),2),'.','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),sum(results.tdvp.expvTime(:,[2,4,6,8]),2),'.','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,[19,20]),2),sum(results.tdvp.expvTime(:,[3,7]),2),'.','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
% set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(15,'twocolumn-single')
%% TDVP (11.6) expvCustom Benchmarking: evolve A+CA
% [ expM, expV, exvCustom, Hn building, numel(Hn)]
figure(16);clf;
hold all
scatter(prod(results.tdvp.expvTime(:,17:19),2),sum(results.tdvp.expvTime(:,[1,4,5,8]+8),2),'.','DisplayName','expM');
scatter(prod(results.tdvp.expvTime(:,17:19),2),sum(results.tdvp.expvTime(:,[2,4,6,8]+8),2),'.','DisplayName','expV');
scatter(prod(results.tdvp.expvTime(:,17:19),2),sum(results.tdvp.expvTime(:,[3,7]+8),2),'.','DisplayName','expvCustom');
set(gca,'yscale','log')
set(gca,'xscale','log')
% set(gca,'xlim',[3,1e4]); set(gca,'xtick',[1,10,100,1e3,1e4]);
legend('Location','best')
xlabel('Matrix dimension n')
ylabel('Time/s')
% formatPlot(16,'twocolumn-single')

%% TDVP (12) expError analysis

figure(12); clf;
[n,m] = size(results.tdvp.expError);
plot(para.tdvp.t(2:end),log10(results.tdvp.expError))

xlabel('Time $t$');
ylabel('$log_{10}$(Exponential error)');
set(gca,'yscale','log');
% set(gca,'zscale','log');
set(gca,'View',[35 42]);
shading interp
rotate3d on
axis tight

%% TDVP (13) prepare artificial GroundState
% clear Workspace, load old result.mat and execute. Saves into new folder
% for fresh start! Initializes Spin in +Sz
para.folder = [para.folder,'-artificial'];
para.filename = sprintf('%s/results.mat',para.folder);
mps{1} = reshape([1,zeros(1,numel(mps{1})-1)],[1,para.D(1),para.d_opt(1)]);
Vmat{1} = eye(para.dk(1));
for j=2:para.L-1
	mps{j} = reshape([1, zeros(1,numel(mps{j})-1)],para.D(j-1),para.D(j),para.d_opt(j));
	Vmat{j} = [zeros(para.dk(j)-para.d_opt(j),para.d_opt(j));...
		fliplr(eye(para.d_opt(j)))];
end
mps{para.L} = reshape([1, zeros(1,numel(mps{para.L})-1)],para.D(para.L-1),1,para.d_opt(para.L));
Vmat{para.L} = [zeros(para.dk(para.L)-para.d_opt(para.L),para.d_opt(para.L));...
		fliplr(eye(para.d_opt(para.L)))];
% Vmat{para.L}(para.dk(para.L):-1:para.dk(para.L)-para.d_opt(para.L)+1,1:para.d_opt(para.L)) = -eye(para.d_opt(para.L));

% results.nx         = getObservable({'occupation'},mps,Vmat,para);
% results.bosonshift = getObservable({'shift'},mps,Vmat,para);

if strcmp(para.model,'SpinBoson')
%     results.spin   = getObservable({'spin'},mps,Vmat,para);
end
% mkdir(para.folder);
% save(para.filename,'mps','Vmat','para','results','op');

%% TDVP SBM multi load files: OrthPol rev22 threading vs perfect L=200
cd('./../cacheComputations/');
res = {};
res{1,1} = load('20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{1,2} = 'v22 multi-core, 1pass, rescaling=0';
res{2,1} = load('20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{2,2} = 'v22 multi-core, restarted, rescaling=0';
res{3,1} = load('20150115-1504-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{3,2} = 'v22 single-core, 1pass, rescaling=1';
res{4,1} = load('20150115-1504-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{4,2} = 'v22 single-core, restarted, rescaling=1';
res{5,1} = load('20150115-1510-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{5,2} = 'v22 single-core, 1pass, rescaling=0';
res{6,1} = load('20150115-1510-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{6,2} = 'v22 single-core, restarted, rescaling=0';
res{7,1} = load('20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand25-expvCustom800-small.mat');
res{7,2} = 'v22 multi-core, restarted, rescaling=0, Bond25';
res{8,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{8,2} = 'perfect';
res{9,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{9,2} = 'v22,expvCustom800,Bond15';
res{10,1}= load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpand-BondExpand15-small.mat');
res{10,2}= 'v22,noExvpVCustom,Bond15';
res{11,1}= load('20150123-1237-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{11,2}= 'v35, rescaling = 1, Lappy';
res{12,1}= load('20150123-1326-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{12,2}= 'v35, rescaling = 0, Lappy';
res{13,1}= load('20150123-1431-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{13,2}= 'v35, rescaling = 1, TCM';
res{14,1}= load('20150123-1432-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{14,2}= 'v35, rescaling = 0, TCM';
res{15,1}= load('20150123-1509-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{15,2}= 'v20, rescaling = 1, TCM';
res{16,1}= load('20150123-1510-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{16,2}= 'v20, rescaling = 0, TCM';

cd('./../TDVP/');
% cell2mat(cellfun(@(x) [x.results.time,x.results.tdvp.time], res(:,1), 'UniformOutput', false))

%% TDVP SBM multi load files: OrthPol L=200 rev20 vs rev35 GS
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{1,2} = 'perfect v19';
res{2,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1509-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{2,2}= 'v20, rescaling = 1, TCM';
res{3,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1510-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{3,2}= 'v20, rescaling = 0, TCM';
res{4,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1602-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{4,2}= 'v26, rescaling = 1, TCM';
res{5,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1602-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{5,2}= 'v26, rescaling = 0, TCM';
res{6,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1554-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{6,2}= 'v27, rescaling = 1, TCM';
res{7,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1555-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{7,2}= 'v27, rescaling = 0, TCM';
res{8,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1607-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{8,2}= 'v28, rescaling = 1, TCM';
res{9,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1609-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{9,2}= 'v28, rescaling = 0, TCM';
res{10,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1614-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{10,2}= 'v30, rescaling = 1, TCM';
res{11,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{11,2}= 'v30, rescaling = 0, TCM';
res{12,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1237-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{12,2}= 'v35, rescaling = 1, Lappy';
res{13,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1326-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{13,2}= 'v35, rescaling = 0, Lappy';
res{14,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1431-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{14,2}= 'v35, rescaling = 1, TCM';
res{15,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1432-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{15,2}= 'v35, rescaling = 0, TCM';

%% TDVP SBM multi load files: OrthPol compare L = 50 rev20, rev22 vs
res = {};
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4-noOBBExpand-noBondExpand.mat');
res{1,2} = '14/11/2014, noExpands, Exp.v19?';
res{2,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150120-1952-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{2,2} = 'GS Exp.v25, rescaling = 1';
res{3,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2047-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{3,2} = 'GS Exp.v23, rescaling = 1';
res{4,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2034-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{4,2} = 'GS Exp.v22, rescaling = 1';
res{5,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2059-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{5,2} = 'GS Exp.v21, rescaling = 1';
res{6,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2118-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{6,2} = 'GS Exp.v20, rescaling = 1';
res{7,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150120-2255-SpinBoson-OrthPol-rev27-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{7,2} = 'GS to22.v27, rescaling = 1';
res{8,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150120-2310-SpinBoson-OrthPol-rev28-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{8,2} = 'GS to22.v28, rescaling = 1';
res{9,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150120-2321-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{9,2} = 'GS to22.v29, rescaling = 1';
res{10,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150120-2346-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{10,2} = 'GS to22.v29, rescaling = 1, higher precision';
res{11,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-0138-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{11,2} = 'GS to22.v29, rescaling = 0';
res{12,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-0156-SpinBoson-OrthPol-rev28-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{12,2} = 'GS to22.v28, rescaling = 0';
res{13,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-1602-SpinBoson-OrthPol-rev30-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{13,2} = 'GS to22.v30, rescaling = 1';
res{14,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-1634-SpinBoson-OrthPol-rev30-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{14,2} = 'GS to22.v30, rescaling = 0';
res{15,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-1817-SpinBoson-OrthPol-rev31-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{15,2} = 'GS to22.v31, rescaling = 1';
res{16,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-1821-SpinBoson-OrthPol-rev31-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{16,2} = 'GS to22.v31, rescaling = 0';
res{17,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-2152-SpinBoson-OrthPol-rev32-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{17,2} = 'GS to22.v32, rescaling = 1';
res{18,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150121-2203-SpinBoson-OrthPol-rev32-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{18,2} = 'GS to22.v32, rescaling = 1, no braces';
res{19,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150122-0006-SpinBoson-OrthPol-rev33-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{19,2} = 'GS to22.v33, rescaling = 1, HmultA sum';
res{20,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150122-0018-SpinBoson-OrthPol-rev34-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{20,2} = 'GS to22.v34, rescaling = 1';
res{21,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1444-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{21,2} = 'GS Exp.v35, rescaling = 1';
res{22,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1447-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{22,2} = 'GS Exp.v35, rescaling = 0';
res{23,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1436-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{23,2} = 'GS Exp.v35, rescaling = 1, TCM';
res{24,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1438-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{24,2} = 'GS Exp.v35, rescaling = 0, TCM';
res{25,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1507-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{25,2} = 'GS Exp.v20, rescaling = 1, TCM';
res{26,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150123-1508-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{26,2} = 'GS Exp.v20, rescaling = 0, TCM';
res{27,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1600-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{27,2} = 'GS Exp.v26, rescaling = 1, TCM';
res{28,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1600-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{28,2} = 'GS Exp.v26, rescaling = 0, TCM';
res{29,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1552-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{29,2} = 'GS Exp.v27, rescaling = 1, TCM';
res{30,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1553-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{30,2} = 'GS Exp.v27, rescaling = 0, TCM';
res{31,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1605-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{31,2} = 'GS Exp.v28, rescaling = 1, TCM';
res{32,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1608-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{32,2} = 'GS Exp.v28, rescaling = 0, TCM';
res{33,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{33,2} = 'GS Exp.v30, rescaling = 1, TCM';
res{34,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150124-1611-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{34,2} = 'GS Exp.v30, rescaling = 0, TCM';
res{35,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\cacheComputations\20150126-1511-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{35,2} = 'GS Exp.v35, rescaling = 1, TCM,prec';

%% TDVP SBM multi load files: OrthPol L= 50 rev28, rev30 architecture sweep
cd('./../cacheComputations/');
res{1,1} = load('20150124-1605-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{1,2} = 'GS Exp.v28, rescaling = 1, Haswell';
res{2,1} = load('20150124-1608-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{2,2} = 'GS Exp.v28, rescaling = 0, Haswell';
res{3,1} = load('20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{3,2} = 'GS Exp.v30, rescaling = 1, Haswell';
res{4,1} = load('20150124-1611-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{4,2} = 'GS Exp.v30, rescaling = 0, Haswell';
res{5,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L50resc1'),'/results.mat']);
res{5,2} = 'GS Exp.v28, rescaling = 1, Sandy';
res{6,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L50resc0'),'/results.mat']);
res{6,2} = 'GS Exp.v28, rescaling = 0, Sandy';
res{7,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L50resc1'),'/results.mat']);
res{7,2} = 'GS Exp.v30, rescaling = 1, Sandy';
res{8,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L50resc0'),'/results.mat']);
res{8,2} = 'GS Exp.v30, rescaling = 0, Sandy';
res{9,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L50resc1'),'/results.mat']);
res{9,2} = 'GS Exp.v28, rescaling = 1, Nehalem';
res{10,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L50resc0'),'/results.mat']);
res{10,2} = 'GS Exp.v28, rescaling = 0, Nehalem';
res{11,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L50resc1'),'/results.mat']);
res{11,2} = 'GS Exp.v30, rescaling = 1, Nehalem';
res{12,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L50resc0'),'/results.mat']);
res{12,2} = 'GS Exp.v30, rescaling = 0, Nehalem';
res{13,1} = load([ls('201501*OrthPol*v28TCM*Core2*L50resc1'),'/results.mat']);
res{13,2} = 'GS Exp.v28, rescaling = 1, Core2';
res{14,1} = load([ls('201501*OrthPol*v28TCM*Core2*L50resc0'),'/results.mat']);
res{14,2} = 'GS Exp.v28, rescaling = 0, Core2';
res{15,1} = load([ls('201501*OrthPol*v30TCM*Core2*L50resc1'),'/results.mat']);
res{15,2} = 'GS Exp.v30, rescaling = 1, Core2';
res{16,1} = load([ls('201501*OrthPol*v30TCM*Core2*L50resc0'),'/results.mat']);
res{16,2} = 'GS Exp.v30, rescaling = 0, Core2';
cd('./../TDVP/');
%% TDVP SBM multi load files: OrthPol L=200 rev28, rev30 architecture sweep
cd('./../cacheComputations/');
res{1,1} = load('20150124-1607-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1/results.mat');
res{1,2} = 'GS Exp.v28, rescaling = 1, Haswell';
res{2,1} = load('20150124-1609-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0/results.mat');
res{2,2} = 'GS Exp.v28, rescaling = 0, Haswell';
res{3,1} = load('20150124-1614-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1/results.mat');
res{3,2} = 'GS Exp.v30, rescaling = 1, Haswell';
res{4,1} = load('20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0/results.mat');
res{4,2} = 'GS Exp.v30, rescaling = 0, Haswell';
res{5,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L200resc1'),'/results.mat']);
res{5,2} = 'GS Exp.v28, rescaling = 1, Sandy';
res{6,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L200resc0'),'/results.mat']);
res{6,2} = 'GS Exp.v28, rescaling = 0, Sandy';
res{7,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L200resc1'),'/results.mat']);
res{7,2} = 'GS Exp.v30, rescaling = 1, Sandy';
res{8,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L200resc0'),'/results.mat']);
res{8,2} = 'GS Exp.v30, rescaling = 0, Sandy';
res{9,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L200resc1'),'/results.mat']);
res{9,2} = 'GS Exp.v28, rescaling = 1, Nehalem';
res{10,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L200resc0'),'/results.mat']);
res{10,2} = 'GS Exp.v28, rescaling = 0, Nehalem';
res{11,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L200resc1'),'/results.mat']);
res{11,2} = 'GS Exp.v30, rescaling = 1, Nehalem';
res{12,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L200resc0'),'/results.mat']);
res{12,2} = 'GS Exp.v30, rescaling = 0, Nehalem';
res{13,1} = load([ls('201501*OrthPol*v28TCM*Core2*L200resc1'),'/results.mat']);
res{13,2} = 'GS Exp.v28, rescaling = 1, Core2';
res{14,1} = load([ls('201501*OrthPol*v28TCM*Core2*L200resc0'),'/results.mat']);
res{14,2} = 'GS Exp.v28, rescaling = 0, Core2';
res{15,1} = load([ls('201501*OrthPol*v30TCM*Core2*L200resc1'),'/results.mat']);
res{15,2} = 'GS Exp.v30, rescaling = 1, Core2';
res{16,1} = load([ls('201501*OrthPol*v30TCM*Core2*L200resc0'),'/results.mat']);
res{16,2} = 'GS Exp.v30, rescaling = 0, Core2';
cd('./../TDVP/');
%% TDVP SBM multi load files: OrthPol/LogZ v37: L=50, L=200; a = 0.01;  VMPS vs artificial		LabBook: 27/01/2015
% works only in cacheComputations!
defPlot(1,:) = {'OrthPol-TDVP-Benchmark-Sz-L50-v37',				[1,3:6]};
defPlot(2,:) = {'OrthPol-TDVP-Benchmark-Sz-L200-v37',				[1,7:9]};
defPlot(3,:) = {'OrthPol-TDVP-Benchmark-Sz-L50-artificial-v37',		[1,10:14]};
defPlot(4,:) = {'OrthPol-TDVP-Benchmark-Sz-L200-artificial-v37',	[1,15:18]};
defPlot(5,:) = {'OrthPol-TDVP-Benchmark-Sz-LogZ-Best-v37',			[1,2,11,16]};
% Legend font size: 16 or 18
cd('./../cacheComputations/');
res{1,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat','para','tresults');
res{1,2} = 'v19, Perfect';
res{2,1} = load('20150126-1800-SpinBoson-LogZ-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{2,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{3,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{3,2} = 'OrthPol L=50, no Expand';
res{4,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{4,2} = 'OrthPol L=50, OBB no Bond';
res{5,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{5,2} = 'OrthPol L=50, Bond no OBB';
res{6,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{6,2} = 'OrthPol L=50, All Expand';
res{7,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{7,2} = 'OrthPol L=200, no Expand';
res{8,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{8,2} = 'OrthPol L=200, OBB no Bond';
res{9,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{9,2} = 'OrthPol L=200, All Expand';
res{10,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{10,2} = 'OrthPol L=50, Art, no Expand';
res{11,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{11,2} = 'OrthPol L=50, Art, OBB no Bond';
res{12,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{12,2} = 'OrthPol L=50, Art, Bond no OBB';
res{13,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{13,2} = 'OrthPol L=50, Art, All, expvCustom800';
res{14,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom1-1core-small.mat','para','tresults');
res{14,2} = 'OrthPol L=50, Art, All, expvCustom1';
res{15,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{15,2} = 'OrthPol L=200,Art, no Expand';
res{16,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{16,2} = 'OrthPol L=200,Art, OBB no Bond';
res{17,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{17,2} = 'OrthPol L=200,Art, Bond no OBB';
res{18,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{18,2} = 'OrthPol L=200,Art, All';
res{19,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{19,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{20,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{20,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond, Art';

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	% ph{8}.LineStyle = ':';				%temp
	set(gca,'ylim',[-1,1]);
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);		% comment me
end
cd('./../TDVP/');
%% TDVP SBM multi load files: Architecture sweep, OrthPol, artificial, L=50						LabBook: 01/02/2015
col = get(groot,'defaultAxesColorOrder');
defPlot(1,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L50-architecture-artificial-v37',	[1:15]};
cd('./../cacheComputations/');
res{ 1,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 1,2} = '\alpha = 0.01, Haswell'; res{1,3} = 'Haswell';
res{ 2,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 2,2} = '\alpha = 0.05, Haswell';
res{ 3,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 3,2} = '\alpha = 0.1, Haswell';
res{ 4,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 4,2} = '\alpha = 0.15, Haswell';
res{ 5,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 5,2} = '\alpha = 0.2, Haswell';
res{ 6,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 6,2} = '\alpha = 0.01, Sandy'; res{6,3} = 'Sandy Bridge';
res{ 7,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 7,2} = '\alpha = 0.05, Sandy';
res{ 8,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 8,2} = '\alpha = 0.1, Sandy';
res{ 9,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 9,2} = '\alpha = 0.15, Sandy';
res{10,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.2, Sandy';
res{11,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{11,2} = '\alpha = 0.01, Nehalem'; res{11,3} = 'Nehalem';
res{12,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{12,2} = '\alpha = 0.05, Nehalem';
res{13,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{13,2} = '\alpha = 0.1, Nehalem';
res{14,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{14,2} = '\alpha = 0.15, Nehalem';
res{15,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{15,2} = '\alpha = 0.2, Nehalem';
res{16,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{16,2} = 'v19, Perfect';
cd('./../TDVP-Git/');
%%
for fignum = 1:size(defPlot,1)
	f=figure(fignum); clf; hold all; f.Name = 'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L50-architecture-artificial-v37-new';
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');

	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; box on;
	cellfun(@(x) set(x,'Color',col(1,:)), ph(1:5));
	cellfun(@(x) set(x,'Color',col(2,:)), ph(6:10));
	cellfun(@(x) set(x,'Color',col(3,:)), ph(11:15));
	set(gca,'ylim',[-1,1]);
	xlabel('$t$');
	ylabel('$\langle\sigma_z\rangle$');
	leg = legend([ph{[1,6,11]}],res{[1,6,11],3},'location','best');
	resizePlot(f,2.3,2,true);
% 	leg.FontSize = 18;
% 	formatPlot(fignum);
% 	title(defPlot{fignum,1},'fontsize',15);
end
%% TDVP SBM Ohmic  s1  a<0.5  Orth2010, OrthPol, art, L=50, L=200; weak -0.75 coupling		LabBook: 30/01/2015, 27/02/2015, POSTER April 2015
cd('./../cacheComputations/');
defPlot(1,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L50-artificial-v37',					[1:5],			{'ylim',[-1,1]}};
defPlot(2,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L200-artificial-v37',				[6:10],			{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-OrthPol-TDVP-s1-med-AllExpand-L200-artificial-v41',				[11:12, 20:21],	{'ylim',[-0.2,1],'xlim',[0,320]}};		% POSTER, Incomplete, PC67
defPlot(4,:) = {'Orth2010-13a-OrthPol-TDVP-OBBExpand-L200-DeltaT1-artificial-v40',			[15:19],		{'ylim',[-1,1],'xlim',[0,320]}};		% POSTER, This is the PERFECT one!

% 01-05: Problem, L = 50
foldPattern = '20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'OBBExpand, L=50, art.';

% 06-10: Problem, dt = 4
foldPattern = '20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art.';


% 11-14		Medium			Finer REDO!! look at PC67
foldPattern = '20150218-1340-SpinBoson-OrthPol-v39TCM1-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art.';

% 15-19: Pretty good weak ohmic!
foldPattern = '20150227-231*-SpinBoson-OrthPol-v40TCMde9-alpha*delta0.1Lambda2dk20D5dopt5L200-artificial';
filePattern = 'results-Till1000Step1v40-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art., T=1000step1';

% 20-21
foldPattern = '20150305-0328-SpinBoson-OrthPol-v41TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L400-artificial';
filePattern = 'results-Till2000Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=400, t=0.5';
cd('./../TDVP/');

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=1$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM Ohmic  s1  0.5<a  Orth2010, OrthPol, artificial, L=200; strong coupling, all FAILS
% clear;
defPlot(1,:) = {'Orth2010-13b-OrthPol-TDVP-noExpand-L200-artificial-v40',							[1:5],	{'ylim',[-1,1], 'xlim',[0,1000]}};	% no good result
defPlot(2,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-artificial-v40-expTime',				[6:10],	{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-DeltaT2-artificial-v40',				[11:16],{'ylim',[-1,1]}};		% Incomplete?
defPlot(4,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-0.5-0.75-DeltaT0.5-artificial-v40',   [17:18,21:22],{'ylim',[0,1], 'xscale','log'}};		% Incomplete pc67
defPlot(5,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-2-artificial-v40',			[17:18,14:16],{'ylim',[0,1], 'xscale','lin','xlim',[1,1000]}};		% Incomplete pc67

n = max(cell2mat(defPlot(:,2)'));
while true
	cd('../cacheComputations/');
% 1-5
foldPattern = '20150225-1301*OrthPol*'; filePattern = '*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'noExpand';
if size(res,1) >= n, break; end;

i = i+1;
res{i,1} = load('20150219-0043-SpinBoson-OrthPol-v39TCM67-alpha0.5delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till100000ExpStep0.1v39-OBBExpand-noBondExpand-expvCustom800-1core-small.mat','para','tresults');
res{i,2} = '\alpha = 0.5'; res{i,3} = 0.5;

% 7-10
foldPattern = '20150218-1609*OrthPol*'; filePattern = '*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-ExpTime';
if size(res,1) >= n, break; end;

% 11-16
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step2*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
if size(res,1) >= n, break; end;

% 17-18
foldPattern = '20150305-0328-SpinBoson-OrthPol-v41TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L400-artificial';
filePattern = 'results-Till2000Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-0.5new';
if size(res,1) >= n, break; end;

% 19-20
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step1*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-1';
if size(res,1) >= n, break; end;

% 21-22
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step0.5*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-0.5';
if size(res,1) >= n, break; end;
end
cd('../TDVP/');


for fignum = 4:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	% ph{8}.LineStyle = ':';				%temp
	set(gca,defPlot{fignum,3}{:});
	xlabel('t');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM Ohmic  s1  0.4<a  Orth2010, OrthPol, coup + art, L=500, NEW							LabBook: 08/05/15, Paper
clear;
defPlot(1,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L500-DeltaT0.01-coupled-v42',			[1:6],{'ylim',[-0.1,1],    'xscale','lin','xlim',[1,1e3]}};
defPlot(2,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L500-DeltaT0.10-artificial-v43',		[7:10],{'ylim',[-0.1,1],   'xscale','lin','xlim',[1,1e3]}};
defPlot(3,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L500-DeltaT0.01-coupled-v42',			[1:6],{'ylim',[-0.1,1.05], 'xscale','log','xlim',[1,2e3],'YTick',[0,0.25,0.5,0.75,1]}};

res = cell(0,4);
i = 0;

% 1-6
foldPattern = '20150409-1304-SpinBoson-OrthPol-v42TCMde9-alpha*delta0.1epsilon0dk30D5dopt5L500';
filePattern = 'results-Till2000Step0.01v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.alpha;
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-2000/0.01';

cd('./../cacheComputations/');
% 7-10
foldPattern = '20150519-1402-SpinBoson-OrthPol-v43TCMde11-alpha*delta0.1epsilon0dk30D5dopt5L500-artificial';
filePattern = 'results-Till1000Step0.1v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.alpha;
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-1000/0.1';
cd('./../TDVP/');

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(sum(x.tresults.t~=0)), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:sum(x.tresults.t~=0)), x.tresults.spin.sz(1:sum(x.tresults.t~=0))), res(pick,1), 'UniformOutput', false);
	axis tight, ax = gca;
	% ph{8}.LineStyle = ':';				%temp
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=1$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM subOhmic s0.5 SBM: Orth2010, OrthPol, L=200:										LabBook: 06/03/2015, Poster!
clear
defPlot( 1,:) = {'Orth2010-14ab-OrthPol-TDVP-noExpand-L200-DeltaT4-artificial-v40',		[1:9],	 {'ylim',[-1,1]}};
defPlot( 2,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L200-DeltaT4-artificial-v40',	[10:18], {'ylim',[-1,1]}};	% not good enough!
defPlot( 3,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[19:27], {'ylim',[-1,1]}};
defPlot( 4,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[28:36], {'ylim',[-1,1]}};
defPlot( 5,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-dk50-coupled-v41',[37:45], {'ylim',[-1,1]}};
defPlot( 6,:) = {'Orth2010-14ab-OrthPol-TDVP-AllExpand-L100-DeltaT0.5-dk50-coupled-v41',[46:54], {'ylim',[-1,1]}};
defPlot( 7,:) = {'Orth2010-14a-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[19:22], {'ylim',[-1,1]}};						% Paper
defPlot( 8,:) = {'Orth2010-14b-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[23:27], {'ylim',[-0.5,1]}};					% Paper
defPlot( 9,:) = {'Orth2010-14a-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[28:31], {'ylim',[-1,1],'xlim',[0,320]}};		% Paper
defPlot(10,:) = {'Orth2010-14b-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[32:36], {'ylim',[0,1],'xlim',[0,320]}};		% Paper

% 01-09
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till324ExpStep4v40-noOBBExpand-noBondExpand*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'noExpand';

% 10-18
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till324ExpStep4v40-OBBExpand-noBondExpand*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand';

% 19-27
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till500Step0.5v*-OBBExpand-noBondExpand*small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Delta t=0.5';

% 28-36
foldPattern = '20150311-172*-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Delta t=0.5';

% 37-45
foldPattern = '20150312-1059-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk50D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Delta t=0.5, dk50';

% 46-54
foldPattern = '20150312-1059-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk50D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-BondExpand10-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'AllExpand, \Delta t=0.5';

for fignum = 5:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM subOhmic s0.25, s0.75 SBM: Kast2013, OrthPol, L=50:									LabBook, Poster!
% use: export_fig(defPlot{get(gcf,'number'),1},'-transparent','-png')
clear
defPlot(1,:) = {'Kast2013-2-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		[1:8],	 {'ylim',[-0.5,1]}};
defPlot(2,:) = {'Kast2013-2-OrthPol-TDVP-AllExpand-L50-DeltaT0.5-artificial-v41',		[9:16],  {'ylim',[-0.5,1]}};
defPlot(3,:) = {'Kast2013-2-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-coupled-v41',			[17:24], {'ylim',[-0.5,1]}};
defPlot(4,:) = {'Kast2013-3-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		[25:31], {'ylim',[-0.5,1]}};
defPlot(5,:) = {'Kast2013-3-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-coupled-v41',			[32:38], {'ylim',[-0.5,1]}};
defPlot(6,:) = {'Kast2013-OrthPol-TDVP-OBBExpand-L50-200-coupled-v41',					[39:41], {'ylim',[-0.1,1]}};
defPlot(7,:) = {'Kast2013-3-OrthPol-TDVP-OBBExpand-L50-DeltaT0.2-5-artificial-v41',		[25:31,42:46], {'ylim',[-0.5,1]}};

cd('../cacheComputations/');
foldPattern = '20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),3);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'artificial';

foldPattern = '20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-BondExpand15-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'AllExpand';

foldPattern = '20150307-1406-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'coupled';

foldPattern = '20150307-1324-SpinBoson-OrthPol-v41TCMde9-s0.25-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 's0.25 artificial';

foldPattern = '20150307-1350-SpinBoson-OrthPol-v41TCMde9-s0.25-alpha*delta0.1epsilon0dk20D5dopt5L50';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 's0.25 coupled';

foldPattern = '20150307-154*-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha0.14delta0.1epsilon0dk20D5dopt5L*';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'opt5L')+5:strfind(file,'\')-1));
	res{i,2} = sprintf('L = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 's0.25 Longer';

% 42-46
foldPattern = '20150522-2155-SpinBoson-OrthPol-v43TCM51-s0.25-alpha*delta0.1epsilon0dk30D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.2v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
% res = [res;cell(size(folds,1),3)];
offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 's0.25 finer, Crossover';

cd('../TDVP/');

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,defPlot{fignum,3}{:});
	% set(gca,'xlim',[1,1e8]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	set(gca,'color','none');
	formatPlot(fignum);
% 	title(defPlot{fignum,1},'fontsize',15);
	if fignum == 5
		ax = axes('Position',[0.17,0.225,0.4,0.3275]);hold on;
		ax.ColorOrderIndex=4; ax.FontSize = 14;
		ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz,...
			'LineWidth',2), res(pick(4:end),1), 'UniformOutput', false);
		box(gca,'on');
		set(gca,'ylim',[0.93,1],'xlim',[0,30]);

	end
end

% figure(fignum+1);clf; hold all;
% ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM compare ML-MCTDH
clear
defPlot(1,:) = {'Wang2010-4-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		1:3,			{'ylim',[-1,1], 'xlim',[0,120]}};
defPlot(2,:) = {'Wang2010-7b-OrthPol-TDVP-AllExpand-L200-Trial-artificial-v42',			[2,4,6,8,11],		{'ylim',[0.6,1],  'xlim',[0,40]}};
defPlot(3,:) = {'Wang2010-7c-OrthPol-TDVP-AllExpand-L200-Trial-artificial-v42',			[3,5,7,9],		{'ylim',[0.6,1], 'xlim',[0,20]}};
defPlot(4,:) = {'Wang2010-7a-OrthPol-TDVP-AllExpand-L600-Trial-artificial-v42',			[10],			{'ylim',[0.4,1], 'xlim',[0,20]}};

cd('./../cacheComputations/');
res = {};
res{1,1} = load('20150316-1346-SpinBoson-OrthPol-v41TCM51-s0.5-alpha0.005delta0.2epsilon0dk20D5dopt5L300-artificial/results-Till1200Step1v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{1,2} = '\alpha = 0.005, Fig4';				res{1,3} = '\Deltat=1,OBB';
res{2,1} = load('20150316-1347-SpinBoson-OrthPol-v41TCM51-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L100-artificial/results-Till400Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat');
res{2,2} = '\alpha = 0.5, Fig7b, \Deltat=0.5';	res{2,3} = '\Deltat=0.5,All7';
res{3,1} = load('20150316-1348-SpinBoson-OrthPol-v41TCM51-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L100-artificial/results-Till200Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat');
res{3,2} = '\alpha = 1, Fig7c, \Deltat=0.5';	res{3,3} = '\Deltat=0.5,All7';
res{4,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till20Step0.1v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{4,2} = '\alpha = 0.5, Fig7b, \Deltat=0.1';	res{4,3} = '\Deltat=0.1,All20';
res{5,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till20Step0.1v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{5,2} = '\alpha = 1, Fig7c, \Deltat=0.1';	res{5,3} = '\Deltat=0.1,All20';
res{6,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{6,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02'; res{6,3} = '\Deltat=0.02,OBB';
res{7,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{7,2} = '\alpha = 1, Fig7c, \Deltat=0.02';	res{7,3} = '\Deltat=0.02,OBB';
res{8,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{8,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02';	res{8,3} = '\Deltat=0.02,All20';
res{9,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{9,2} = '\alpha = 1, Fig7c, \Deltat=0.02';	res{9,3} = '\Deltat=0.02,All20';
res{10,1} = load('20150325-1532-SpinBoson-OrthPol-v42TCM59-s0.5-alpha0.3delta0.2epsilon0dk20D5dopt5L600-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.3, Fig7a, \Deltat=0.02';	res{10,3} = '\Deltat=0.02,OBB';
res{11,1} = load('20150325-1532-SpinBoson-OrthPol-v42TCM59-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L600-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{11,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02';	res{11,3} = '\Deltat=0.02,OBB';
cd('./../TDVP/');

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)).*(-x.para.hx)./2, x.tresults.spin.sz,'*'), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,defPlot{fignum,3}{:});
	xlabel('$t\Delta$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,3},'location','best');
	leg.FontSize = 18;
	set(gca,'color','none');
	formatPlot(fignum);
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi files: Ohmic LogZ																	LabBook: 04/02/2015
clear
defPlot(1,:) = {'Orth2010-13b-LogZ1.1-TDVP-OBBExpand-L100-DeltaT2-artificial-v40',		[1:6]};

foldPattern = '20150225*LogZ*'; filePattern = '*Till1000Step2v40-OBBExpand-noBondExpand-expvCustom800-1core*-zAvg.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),3);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3);

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,'ylim',[-1,1]);
% 	set(gca,'xlim',[1,1e3]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM multi files: Star Bath occupation
clear
defPlot(1,:) = {'Orth2010-13a-s1-OrthPol-TDVP-OBB10-L200-325T1-art-v42',		[1:5],	{'ylim',[-1,1]}};
defPlot(2,:) = {'Orth2010-14a-s0.5-OrthPol-TDVP-OBB10-L100-325T1-art-v42',		[6:9],	{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-14b-s0.5-OrthPol-TDVP-OBB10-L100-325T1-art-v42',		[10:14],{'ylim',[-0.5,1]}};
defPlot(4,:) = {'Orth2010-13a-s1-OrthPol-TDVP-OBB10-L200-325T1-coup-v42',		[15:19],{'ylim',[-1,1]}};
defPlot(5,:) = {'Orth2010-14a-s0.5-OrthPol-TDVP-OBB10-L100-325T1-coup-v42',		[20:23],{'ylim',[-1,1]}};
defPlot(6,:) = {'Orth2010-14b-s0.5-OrthPol-TDVP-OBB10-L100-325T1-coup-v42',		[24:28],{'ylim',[0,1]}};

foldPattern = '20150327-1327-SpinBoson-OrthPol-v42TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res = sortrows(res,3); res(1:i,4) = {'$s=1$, vac'};

foldPattern = '20150327-1330-SpinBoson-OrthPol-v42TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100-artificial';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=0.5$, vac'};

foldPattern = '20150327-1430-SpinBoson-OrthPol-v42TCMde10-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=1$, coup'};

foldPattern = '20150327-1434-SpinBoson-OrthPol-v42TCMde10-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=0.5$, coup'};

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi files: Star Bath and POLARON NEW L=500										LabBook: 26/05/2015
clear
defPlot(1,:) = {'AlltheRest artificial',													[1:8],{'ylim',[-1,1],'xlim',[0,500]}};
defPlot(2,:) = {'AlltheRest coupled',														[9:16],{'ylim',[-1,1],'xlim',[0,500]}};

i=0;

%1-8
foldPattern = '20150529-100*-SpinBoson-OrthPol-v43TCM*-*alpha*delta0.1epsilon0dk30D5dopt5L500-artificial';
filePattern = 'results-Till500Step0.25v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),4} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.alpha;
	res{i,4} = res{i,1}.para.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),4);
% res(offset+1:i,4) = {'L=500,art,Polaron'};

%9-16
foldPattern = '20150531-1124-SpinBoson-OrthPol-v43TCM*-*alpha*delta0.1epsilon0dk30D5dopt5L500';
filePattern = 'results-Till500Step0.25v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),4} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.alpha;
	res{i,4} = res{i,1}.para.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),4);

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(find(x.tresults.t,1,'last')), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:find(x.tresults.t,1,'last')), x.tresults.spin.sz(1:find(x.tresults.t,1,'last'))), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM super-Ohmic s3  0.01<a<1  OrthPol, art, L=300											LabBook: 13/05/15,
clear;
defPlot(1,:) = {'OrthPol-TDVP-SBM-s3-OBBExpand-L300-DeltaT0.2-art-v43',						[1:7],{'ylim',[-1,1], 'xscale','lin','xlim',[1,500]}};		% Incomplete node9
defPlot(2,:) = {'OrthPol-TDVP-SBM-s3-OBBExpand-L500-DeltaT0.1-art-v43',						[8:14],{'ylim',[-1,1], 'xscale','lin','xlim',[1,500]}};		% Incomplete node9

res = cell(0,4);
i = 0;

foldPattern = '20150512-1658-SpinBoson-OrthPol-v43TCMde10-s3-alpha*delta0.1epsilon0dk30D5dopt5L300-artificial';
filePattern = 'results-Till500Step0.2v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-500/0.2';

foldPattern = '20150519-1350-SpinBoson-OrthPol-v43TCMde10-s3-alpha*delta0.1epsilon0dk30D5dopt5L500-artificial';
filePattern = 'results-Till1000Step0.1v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-1000/0.1';


for fignum = 1:size(defPlot,1)
	%%
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(sum(x.tresults.t~=0)), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:sum(x.tresults.t~=0)), x.tresults.spin.sz(1:sum(x.tresults.t~=0))), res(pick,1), 'UniformOutput', false);
	axis tight, ax = gca;
	% ph{8}.LineStyle = ':';				%temp
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=1$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);
% ph2 = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)-1), diff(x.para.tdvp.calcTime)), res(pick,1), 'UniformOutput', false);

%% TDVP SBM Ohmic: v52n-v62 TDVP Benchmark s=1 0.01 < a < 0.5									% LabBook 06/08/2015 TOOODOOOO!!!
% deleted v52 results, since no real extra information gain. If needed, rerun!
% includes: j, sx, sn
clear
defPlot(1,:) = {'20150805-Benchmark-v52n-dt01-Using-ExpvCustom-only',		[1:6],  {'ylim',[-1,1],'xlim',[0,500]}};
defPlot(2,:) = {'20150928-Benchmark-v62-dt01-Bond5',						[7:14], {'ylim',[-1,1],'xlim',[0,300]}};
defPlot(3,:) = {'20150928-Benchmark-v62-dt01-Bond20',						[15:22],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(4,:) = {'20150928-Benchmark-v62-dt01-Bond150',						[23:30],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(5,:) = {'20150928-Benchmark-v62-dt1-Bond5',							[31:38],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(6,:) = {'20150928-Benchmark-v62-dt1-Bond20',						[39:46],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(7,:) = {'20150928-Benchmark-v62-dt1-Bond150',						[47:54],{'ylim',[-1,1],'xlim',[0,300]}};


i=0; cols = 5;
n = max(cell2mat(defPlot(:,2)'));
while true
%1-6
foldPattern = '20150805-2056-SpinBoson-OrthPol-v52TCMde10-alpha*delta0.1epsilon0dk30D5dopt5L200-art-sz';
filePattern = 'results-Till500Step0.1v52n-OBBExpand-noBondExpand-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'expvCustom, v52n, dt0.1'};

%7-14
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step0.1v62-alpha*-OBBmax40-Dmax5-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB40E, D5E, v62, dt0.1'};

%15-22
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step0.1v62-alpha*-OBBmax40-Dmax20-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB40E, D20E, v62, dt0.1'};

%23-30
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step0.1v62-alpha*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB40E, D150E, v62, dt0.1'};

%31-38
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step1v62-alpha*-OBBmax15-Dmax5-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB15E, D5E, v62, dt1'};

%39-46
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step1v62-alpha*-OBBmax15-Dmax20-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB15E, D20E, v62, dt1'};

%47-52
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step1v62-alpha*-OBBmax15-Dmax150-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB15E, D150E, v62, dt1'};

if size(res,1) >= n, break; end;
end
%% TDVP SBM Ohmic: v52n-v62 TDVP Benchmark s=1 0.01 < a < 0.5		TDVPData					% LabBook 06/08/2015 TOOODOOOO!!!
% deleted v52 results, since no real extra information gain. If needed, rerun!
% includes: j, sx, sn
clear
defPlot(1,:) = {'20150805-Benchmark-v52n-dt01-Using-ExpvCustom-only',		[1:6],  {'ylim',[-1,1],'xlim',[0,500]}};
defPlot(2,:) = {'20150928-Benchmark-v62-dt01-Bond5',						[7:14], {'ylim',[-1,1],'xlim',[0,300]}};
defPlot(3,:) = {'20150928-Benchmark-v62-dt01-Bond20',						[15:22],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(4,:) = {'20150928-Benchmark-v62-dt01-Bond150',						[23:30],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(5,:) = {'20150928-Benchmark-v62-dt1-Bond5',							[31:38],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(6,:) = {'20150928-Benchmark-v62-dt1-Bond20',						[39:46],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(7,:) = {'20150928-Benchmark-v62-dt1-Bond150',						[47:54],{'ylim',[-1,1],'xlim',[0,300]}};


i=0; cols = 5;
n = max(cell2mat(defPlot(:,2)'));

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'SpinBoson')),TDVPfolds));
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'OrthPol')),TDVPfolds));
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'2015')),TDVPfolds));
matches = [];
legLabels = [];

%1-6
dirPat = '20150805-2056-SpinBoson-OrthPol-v52TCMde10-alpha.*delta0.1epsilon0dk30D5dopt5L200-art-sz';
filPat = 'results-Till500Step0.1v52n-OBBExpand-noBondExpand-expvCustom0-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'expvCustom, v52n, dt0.1'};

%7-14
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step0.1v62-alpha.*-OBBmax40-Dmax5-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB40E, D5E, v62, dt0.1'};

%15-22
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step0.1v62-alpha.*-OBBmax40-Dmax20-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB40E, D20E, v62, dt0.1'};

%23-30
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step0.1v62-alpha.*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB40E, D150E, v62, dt0.1'};

%31-38
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step1v62-alpha.*-OBBmax15-Dmax5-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB15E, D5E, v62, dt1'};

%39-46
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step1v62-alpha.*-OBBmax15-Dmax20-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB15E, D20E, v62, dt1'};

%47-52
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step1v62-alpha.*-OBBmax15-Dmax150-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v.*-Dmax(?<D>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$D = %s, dt = %s$',tokens{ii}.D,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$D = %s$',tokens{ii}.D), I, 'UniformOutput', false);
legLabels = [legLabels, y];
% res(offset+1:i,5) = {'OBB15E, D150E, v62, dt1'};

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(arrayfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.sz(1:x.tresults.lastIdx), 'DisplayName',num2str(x.alpha)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	set(gca,'color','none');
	if fignum ~= 4
		ax.ColorOrderIndex = 1;
		ph2 = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.sz(1:x.tresults.lastIdx),'.black'), res(defPlot{4,2},1), 'UniformOutput', false);
	end
	leg = legend([ph{:}]);
end

% plot time taken
figure(fignum+1);clf; hold all;
ph = arrayfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);
%% plot error wrt D150 run															% LabBook 25/09/2015

for jj = 2:3
	figure(fignum+jj); clf; hold all;
	pick = defPlot{jj,2};
	for ii = 1:8
		plot( res(pick(ii)).tresults.t, res(ii+22).tresults.spin.sz - res(pick(ii)).tresults.spin.sz);
	% 	plot( res{ii,1}.tresults.t, res{ii+16,1}.tresults.spin.sz -
	% 	res{ii+8}.tresults.spin.sz);	% D20
% 		plot( res{ii,1}.tresults.t, res{ii+16,1}.tresults.spin.sz - res{ii}.tresults.spin.sz);
	end
% 	leg = legend(cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
% 	formatPlot(fignum+jj,'twocolumn-single');
end

for jj = 5:6
	figure(fignum+jj); clf; hold all;
	pick = defPlot{jj,2};
	for ii = 1:8
		plot( res(pick(ii)).tresults.t, res(ii+46).tresults.spin.sz - res(pick(ii)).tresults.spin.sz);
	% 	plot( res{ii,1}.tresults.t, res{ii+16,1}.tresults.spin.sz -
	% 	res{ii+8}.tresults.spin.sz);	% D20
% 		plot( res{ii,1}.tresults.t, res{ii+16,1}.tresults.spin.sz - res{ii}.tresults.spin.sz);
	end
% 	leg = legend(cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
% 	formatPlot(fignum+jj,'twocolumn-single');
end

%% plot error wrt D150 run as mean error			% Dissertation Plot!			% LabBook 25/09/2015
h = struct();
h.f = figure(120); clf;
h.f.Name = 'SBM-Convergence';
h.ax = gca;
set(h.f,'DefaultAxesFontSize', 8,...
		'DefaultLineLineWidth',1);
resizePlot(h.f,2.3,2,true);
col = get(0,'defaultaxescolororder');		% not necessary?
col(8,:) = [0 0 0];
hold all; box on;
for ii = 1:8
pick = [7 15 31 39 47 23]+ii-1;
ttPop = cell2mat(arrayfun(@(x) TDVPData.subsample(real(x.getData('sz')),x.dt,1),res(pick)','UniformOutput',false));
	% estimate mean error
ttErr = bsxfun(@minus,-ttPop,-ttPop(:,end));
ttErr = bsxfun(@rdivide,ttErr,ttPop(:,end)); ttErr(isnan(ttErr)) = 0;		% makes error relative
ttMeanErr = abs(mean(ttErr,1));
	% get Runtimes
tSweep = arrayfun(@(x) 60*mean(diff(nonzeros(x.para.tdvp.calcTime))),res(pick));
h.pl = scatter(tSweep(1:2), ttMeanErr(1:2),60,[1;1]*col(ii,:),'filled');
h.pl = scatter(tSweep(3:5), ttMeanErr(3:5),60,[1;1;1]*col(ii,:),'x');
% set(h.ax(ii), defPlot{ii,3}{:});
h.ax.YScale = 'log';
xlabel('$t_{CPU}$/Sweep/min');
ylabel('$\varepsilon(\langle\sigma_z\rangle)$');

	% shift text labels
offX = +ones(1,length(pick))*0.015; offY = ones(1,length(pick))*0;
% offX(5) = 0.2; offY(5) = 0;
% offX(7) = -1.0; offY(7) = -0.005;
% offX(8) = -0.1; offY(8) = 0.008;
% offX(9) = -0.6; offY(9) = 0.003;
offY = offY;
for kk = 1:length(pick)
% 	t = text(tSweep(kk)+offX(kk),double(ttMeanErr(kk))+offY(kk),strrep(res(pick(kk)).LegLabel,', ','\\\\'),'FontSize',8);
	t = text(tSweep(kk)+offX(kk),double(ttMeanErr(kk))+offY(kk),['$',res(pick(kk)).LegLabel(5:end)],'FontSize',8);
end
end
% axis tight

%% TDVP SBM multi files: v58 SBM2C Benchmark s=1 0.01 < a < 1						% LabBook ??/08/2015
clear
defPlot(1,:) = {'20150812-SBM2C-dt1-FirstTry',								[1:8], {'ylim',[-1,1],'xlim',[0,100]}};

i=0; cols = 5;

%1-6
foldPattern = '20150812-1731-SpinBoson2C-OrthPol-v58TCMde10-alpha*delta0.1epsilon0dk20D5dopt5L200-art-sz';
filePattern = 'results-Till350Step1v58-OBBExpand-BondExpand5-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'expvCustom, dt0.1'};


for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.sz(1:x.tresults.lastIdx)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi load : v59 TTM Benchmark, s=1 0.01 < a < 1						% LabBook 16/08/2015
clear
defPlot(1,:) = {'20150816-SBMTTM-dt002-FirstTry',							[1:8], {'ylim',[1e-8,1],'xlim',[0,100],'yscale','log'}};
defPlot(2,:) = {'20150805-Benchmark-v52-dt01-Using-ExpvCustom-only',		[9:14,17,18], {'ylim',[-1,1],'xlim',[0,500]}};
defPlot(3,:) = {'20150929-SBMTTM-dt01-szcoup',								[19:21],{'ylim',[-1,1],'xlim',[0,400]}};
defPlot(4,:) = {'20150928-Benchmark-v62-dt01-Bond150',						[22:29],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(5,:) = {'20150928-SBMTTM-dt01-szcoup-SBM-compare',					[22,27,28],{'ylim',[-1,1],'xlim',[0,1000]}};
i=0; cols = 5;

%1-8: TTM Data
foldPattern = '20150816-0156-SpinBosonTTM-OrthPol-v59TCMde9-alpha*delta0.1epsilon0dk30D5dopt5L100';
filePattern = 'results-Till150Step0.02v59-OBBExpand-BondExpand10-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'TTM, dt0.02'};

%9-14: TDVP Data weak
foldPattern = '20150805-2056-SpinBoson-OrthPol-v52TCMde10-alpha*delta0.1epsilon0dk30D5dopt5L200-art-sz';
filePattern = 'results-Till500Step0.1v52n-OBBExpand-noBondExpand-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'TDVP, dt0.1'};

% 15-18: TDVP Data strong
cd('../cacheComputations/');
foldPattern = '20150519-1402-SpinBoson-OrthPol-v43TCMde11-alpha*delta0.1epsilon0dk30D5dopt5L500-artificial';
filePattern = 'results-Till1000Step0.1v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.alpha;
	res{i,4} = res{i,1}.para.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'TDVP, dt0.1'};
cd('./../TDVP/');

%19-21: TTM Data
foldPattern = '20150929-1353-SpinBosonTTM-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till100Step0.1v62-alpha*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'TTM, dt0.02'};

%22-29
foldPattern = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha*delta0.1epsilon0dk40D5dopt5L200';
filePattern = 'results-Till300Step0.1v62-alpha*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.chain{1}.s;
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,5) = {'OBB40E, D150E, v62, dt0.1'};

	%% plot TTM norm
for fignum = [1,3]%:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	plot([0,xmax],[0,0],'black');
% 	ph = cellfun(@(x) plot(x.tresults.t, x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	ph = cellfun(@(x) plot(x.tresults.t(2:x.tresults.lastIdx), x.tresults.TTM.Tnorm(1:x.tresults.lastIdx-1)./(x.para.tdvp.deltaT.^2)), res(pick,1), 'UniformOutput', false); % plot TTM norm
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
% 	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	ylabel('$|T|/\Delta t^2$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);
	%% reconstruct Dynamics using the TTM & plot
finalT = 1000;
[sx,sy,sz] = spinop('Z');
tic;
for k = 19:21%6:size(res,1)
	n = round(finalT/res{k,1}.para.tdvp.deltaT)+1;
	rhoT = zeros(length(res{k,1}.tresults.TTM.T)*4,1);
	Esigma = zeros(n,3);
	T = reshape(res{k,1}.tresults.TTM.T, 4,[]);			% creates [T(1) T(2) T(3) ...]
	for i = 1:n
		if i == 1
			rho = [1,0,0,0]';
		else
			rho = T*rhoT;
		end
		rhoT = [rho; rhoT(1:end-4)];					% prepend new vector rho(i)
		rho = reshape(rho,[2,2]);						% reshape rho(i) for observables
		Esigma(i,1) = trace(sx*rho);
		Esigma(i,2) = trace(sy*rho);
		Esigma(i,3) = trace(sz*rho);
	end
	res{k,1}.tresults.spin.sx = real(Esigma(:,1));
	res{k,1}.tresults.spin.sy = real(Esigma(:,2));
	res{k,1}.tresults.spin.sz = real(Esigma(:,3));
	res{k,1}.tresults.t       = 0:res{k,1}.para.tdvp.deltaT:finalT;
% 	plot(res{k,1}.tresults.TTM.t,res{k,1}.tresults.TTM.Esigma(:,3));
end
toc
	%% Plot Spin Dynamics for szcoup
	set(0,'defaulttextinterpreter','latex');
for fignum = 5 %:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph = cellfun(@(x) plot(x.tresults.t, x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(2:x.tresults.lastIdx), x.tresults.TTM.Tnorm(1:x.tresults.lastIdx-1)./(x.para.tdvp.deltaT.^2)), res(pick,1), 'UniformOutput', false); % plot TTM norm
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	plot(ax.XLim,[0,0],'k');
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
% 	ylabel('$|T|/\Delta t^2$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
	ph = cellfun(@(x) plot(x.tresults.t, x.tresults.spin.sz,'k.'), res(defPlot{3,2},1), 'UniformOutput', false);
end

%% TDVP SBM multi load : v59 TTM Benchmark, s=1 0.01 < a < 1	TDVPData			% LabBook 16/08/2015
clear
defPlot(1,:) = {'20150816-SBMTTM-dt0.02-FirstTry',							[1:8], {'ylim',[1e-5,1],'xlim',[0,150],'yscale','log'}};
defPlot(2,:) = {'20150805-Benchmark-v52-dt01-Using-ExpvCustom-only',		[9:14,17,18], {'ylim',[-1,1],'xlim',[0,500]}};
defPlot(3,:) = {'20150929-SBMTTM-dt01-szcoup',								[19:21],{'ylim',[-1,1],'xlim',[0,100],'yscale','log'}};
defPlot(4,:) = {'20150928-Benchmark-v62-dt01-Bond150',						[22:29],{'ylim',[-1,1],'xlim',[0,300]}};
defPlot(5,:) = {'20150928-SBMTTM-dt01-szcoup-SBM-compare',					[22,27,28],{'ylim',[-1,1],'xlim',[0,2300]}};
i=0; cols = 5;

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'SpinBoson')),TDVPfolds));
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'2015')),TDVPfolds));
matches = [];
legLabels = [];

% 1-8: TTM Data
dirPat = '20150816-0156-SpinBosonTTM-OrthPol-v59TCMde9-alpha.*delta0.1epsilon0dk30D5dopt5L100';
filPat = 'results-Till150Step0.02v59-OBBExpand-BondExpand10-expvCustom800-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
legLabels = [legLabels, y];

% 9-14: TDVP Data weak
dirPat = '20150805-2056-SpinBoson-OrthPol-v52TCMde10-alpha.*delta0.1epsilon0dk30D5dopt5L200-art-sz';
filPat = 'results-Till500Step0.1v52n-OBBExpand-noBondExpand-expvCustom0-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
legLabels = [legLabels, y];

% 15-18: TDVP Data strong
dirPat = '20150519-1402-SpinBoson-OrthPol-v43TCMde11-alpha.*delta0.1epsilon0dk30D5dopt5L500-artificial';
filPat = 'results-Till1000Step0.1v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
legLabels = [legLabels, y];

% 19-21: TTM Data
dirPat = '20150929-1353-SpinBosonTTM-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till100Step0.1v62-alpha.*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$\\alpha = %s$',tokens{ii}.alpha), I, 'UniformOutput', false);
legLabels = [legLabels, y];

% 22-29
dirPat = '20150925-1538-SpinBoson-OrthPol-v62TCMde9-alpha.*delta0.1epsilon0dk40D5dopt5L200';
filPat = 'results-Till300Step0.1v62-alpha.*-OBBmax40-Dmax150-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'alpha(?<alpha>[0-9\.]*)delta.*Step(?<dt>[0-9\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.alpha),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
% y = arrayfun(@(ii) sprintf( '$\\alpha = %s, dt = %s$',tokens{ii}.alpha,tokens{ii}.dt), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$%s$',tokens{ii}.alpha), I, 'UniformOutput', false);
legLabels = [legLabels, y];

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot TTM norm
for fignum = [1,3]%:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	resizePlot(f,2,2,true);
	ph = res(pick).plot('ttm-norm',f);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	legend toggle
% 	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	resizePlot(f,2.5,2,true);
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
	title('');
end

figure(fignum+1);clf; hold all;
ph = arrayfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);
	%% reconstruct Dynamics using the TTM & plot
finalT = 2300;
[sx,sy,sz] = spinop('Z');
tic;
for k = 19:21%6:size(res,1)
	n = round(finalT/res(k).para.tdvp.deltaT)+1;
	rhoT = zeros(length(res(k).tresults.TTM.T)*4,1);
	Esigma = zeros(n,3);
	T = reshape(res(k).tresults.TTM.T, 4,[]);			% creates [T(1) T(2) T(3) ...]
	for i = 1:n
		if i == 1
			rho = [1,0,0,0]';
		else
			rho = T*rhoT;
		end
		rhoT = [rho; rhoT(1:end-4)];					% prepend new vector rho(i)
		rho = reshape(rho,[2,2]);						% reshape rho(i) for observables
		Esigma(i,1) = trace(sx*rho);
		Esigma(i,2) = trace(sy*rho);
		Esigma(i,3) = trace(sz*rho);
	end
	res(k).tresults.spin.sx = real(Esigma(:,1));
	res(k).tresults.spin.sy = real(Esigma(:,2));
	res(k).tresults.spin.sz = real(Esigma(:,3));
	res(k).tresults.t       = 0:res(k).para.tdvp.deltaT:finalT;
% 	plot(res{k,1}.tresults.TTM.t,res{k,1}.tresults.TTM.Esigma(:,3));
end
toc
	%% Plot Spin Dynamics for szcoup
	set(0,'defaulttextinterpreter','latex');
for fignum = 5 %:size(defPlot,1)
	f = figure(fignum); clf; hold all; resizePlot(f,2,2,true);
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(arrayfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph2 = arrayfun(@(x) plot(x.tresults.t, x.tresults.spin.sz,'k'), res(defPlot{3,2},1), 'UniformOutput', false);
	h = res(pick).plot('sz',f);
	title('');
	set(h(1).ax,defPlot{fignum,3}{:});
	plot(h(1).ax.XLim,[0,0],'Color',[1 1 1]*0,'LineStyle',':');
	xlabel('$\omega_ct$');
	set(gca,'color','none');
	leg = legend([h.pl],cellfun(@(x) sprintf('%s\n',x),legLabels(pick),'UniformOutput',false),'location','southeast');
	t1 = text(leg.Position(1)+h(1).ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$','Units','norm','VerticalAlignment','bottom');
	legend boxoff
	resizePlot(f,2,2,true);
end

	%% plot TTM norm
for fignum = [1]%:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
% 	pick = defPlot{fignum,2};			% plot all
	pick = [1,2,5];			% plot all
	resizePlot(f,2,2,true);
	h = res(pick).plot('ttm-norm-hull',f);
	axis tight; ax = gca;
	ax.YTick=10.^[-4:0];
	set(ax,defPlot{fignum,3}{:});
	legend([h.pl],{'0.01','0.05','0.2'});
	legend boxoff
	resizePlot(f,2.5,2,true);
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
	title('');
end

	%% reconstruct Dynamics using the TTM & plot
finalT = 1000;
[sx,sy,sz] = spinop('Z');
tic;
for k = 1:8%6:size(res,1)
	n = round(finalT/res(k).para.tdvp.deltaT)+1;
	rhoT = zeros(length(res(k).tresults.TTM.T)*4,1);
	Esigma = zeros(n,3);
	T = reshape(res(k).tresults.TTM.T, 4,[]);			% creates [T(1) T(2) T(3) ...]
	for i = 1:n
		if i == 1
			rho = [1,0,0,0]';
		else
			rho = T*rhoT;
		end
		rhoT = [rho; rhoT(1:end-4)];					% prepend new vector rho(i)
		rho = reshape(rho,[2,2]);						% reshape rho(i) for observables
		Esigma(i,1) = trace(sx*rho);
		Esigma(i,2) = trace(sy*rho);
		Esigma(i,3) = trace(sz*rho);
	end
	res(k).tresults.spin.sx = real(Esigma(:,1));
	res(k).tresults.spin.sy = real(Esigma(:,2));
	res(k).tresults.spin.sz = real(Esigma(:,3));
	res(k).tresults.t       = 0:res(k).para.tdvp.deltaT:finalT;
% 	plot(res{k,1}.tresults.TTM.t,res{k,1}.tresults.TTM.Esigma(:,3));
end
toc
	%% Plot Spin Dynamics for szcoup
	set(0,'defaulttextinterpreter','latex');
for fignum = 2 %:size(defPlot,1)
	f = figure(fignum); clf; hold all; resizePlot(f,2,2,true);
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(arrayfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph2 = arrayfun(@(x) plot(x.tresults.t/2, x.tresults.spin.sz,'k'), res(defPlot{1,2},1), 'UniformOutput', false);
	h = res(pick).plot('sz',f);
	title('');
	set(h(1).ax,defPlot{fignum,3}{:});
	plot(h(1).ax.XLim,[0,0],'k');
	xlabel('$\omega_ct$');
	set(gca,'color','none');
	leg = legend([h.pl],cellfun(@(x) sprintf('%s\n',x),legLabels(pick),'UniformOutput',false),'location','southeast');
	t1 = text(leg.Position(1)+h(1).ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$','Units','norm','VerticalAlignment','bottom');
	legend boxoff
	resizePlot(f,2,2,true);
end


%% THERM SBM multi load : v61 T=300K Benchmark, s=1 a=0.01							% LabBook 28/08/2015
% Only contains the first THERM trials!
clear
defPlot(1,:) = {'20150828-SBM2CT-dt01-FirstTry',							[ 1:5], {'ylim',[0,0.8],'xlim',[0,200],'yscale','lin'}};
defPlot(2,:) = {'20150828-SBM2CT-dt01-sz',									[ 6:9], {'ylim',[0,0.8],'xlim',[0,200],'yscale','lin'}};
defPlot(3,:) = {'20150904-SBM2CT-dt01-sz-DimSweep',							[10:16], {'ylim',[0,4],'xlim',[1,20],'yscale','lin'}};
defPlot(4,:) = {'20150905-SBM2CT-L20-sz-DimExpandSweep',					[17:32], {'ylim',[0,0.4],'xlim',[1,20],'yscale','lin'}};
i=0; cols = 5;

%1-5: Thermal Evolution from -sx
foldPattern = '20150827-1936-SpinBoson2CT-VT-OrthPol-v61TCMde10-alpha0.01delta0epsilon0.1dk*D5dopt5L*-art--sx';
filePattern = 'results-Till20Step0.1v61-OBBExpand-BondExpand20-expvCustom0-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults','results');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('d_k = %g', res{i,1}.para.dk(1,2));
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%6-9: Thermal Evolution from -sz
foldPattern = '20150828-1406-SpinBoson2CT-VT-OrthPol-v61TCMde9-alpha0.01delta0epsilon0.1dk*D5dopt5L*-art--sz';
filePattern = 'results-Till20Step0.1v61-OBBExpand-BondExpand20-expvCustom0-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults','results');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('d_k = %g', res{i,1}.para.dk(1,2));
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%10-16: Thermal Evolution from -sz
foldPattern = '20150904-0006-SpinBoson2CT-VT-OrthPol-v61TCMde11-alpha0.01delta0epsilon0.1dk60D*dopt*L20-art--sz';
filePattern = 'results-Till20Step0.1v61-noOBBExpand-noBondExpand-expvCustom0-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	load(file,'para','tresults','results');			% comment first!
	Vars = whos;
	for ii = 1:size(Vars)
		if strcmp(Vars(ii).class,'uint8')
			eval(sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
		end
	end
 	res{i,1}.para = para; res{i,1}.results = results; res{i,1}.tresults = tresults;
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('dOBB=%g, D=%g', res{i,1}.para.d_opt(3,3),res{i,1}.para.D(7));
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%17-27: Thermal Evolution from -sz
foldPattern = '20150905-0139-SpinBoson2CT-VT-OrthPol-v61TCM33-alpha0.01delta0epsilon0.1dk60D10dopt10L20-art--sz';
filePattern = 'results-Till20*small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
% 	res{i,2} = sprintf('dOBB=%g, D=%g', res{i,1}.para.d_opt(3,3),res{i,1}.para.D(7));
	[~,res{i,2}] = fileparts(res{i,1}.para.tdvp.filename);
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%%
for fignum = 4:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
% 	plot([0,xmax],[0,0],'black');
% 	ph = cellfun(@(x) plot(x.tresults.TTM.t, x.tresults.TTM.Esigma(3,:)), res(pick,1), 'UniformOutput', false);
	ph = cellfun(@(x) plot(1:size(x.tresults.n,2), x.tresults.n(x.tresults.lastIdx,:,1)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('Site $k$');
	ylabel('$\left<n_k\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%s\n',x),res(pick,2),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi load  : v61 T=300K Benchmark, s=1 a=0.01							% LabBook 30/08/2015
% TDVP on THERM evolution
clear
defPlot(1,:) = {'20150830-SBM2CT-dt1-dk10-TH-TDVP',							[1:4], {'ylim',[0,1] ,'xlim',[0,95]}};
defPlot(2,:) = {'20150831-SBM2CT-dt1-dk20-TH-TDVP',							[5:8], {'ylim',[0,1] ,'xlim',[0,95]}};

i=0; cols = 5;

%1-4: TDVP Data
foldPattern = '20150827-1936-SpinBoson2CT-VT-OrthPol-v61TCMde10-alpha0.01delta0epsilon0.1dk10D5dopt5L50-art--sx';
filePattern = 'results-Till350Step1v61-alpha*-OBBExpand-noBondExpand-expvCustom0-1core.mat-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('d_k = %g', res{i,1}.para.dk(1,2));
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt1,dk10,L50'};

%5-8: Thermal Evolution from -sz
foldPattern = '20150828-1406-SpinBoson2CT-VT-OrthPol-v61TCMde9-alpha0.01delta0epsilon0.1dk20D5dopt5L200-art--sz';
filePattern = 'results-Till350Step1v61-alpha*-OBBExpand-BondExpand20-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('d_k = %g', res{i,1}.para.dk(1,2));
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt1,dk20,L200'};
%%
for fignum = 4:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.visibility(1:x.tresults.lastIdx)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('Site $k$');
	ylabel('$\left|D(t)\right|$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%s\n',x),res(pick,2),'UniformOutput',false),'location','best');
	leg.Interpreter = 'latex';
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TH-TDVP: v61 Evolution of thermal Bath alone										% LabBook 03/09/2015
clear;
defPlot(1,:) = {'20150902-SBM2CT-dt1-TH-TDVP-noExcitation',							[1:5], {'ylim',[0,0.7] ,'xlim',[1,50]}};

res{1,1} = load('20150827-1936-SpinBoson2CT-VT-OrthPol-v61TCMde10-alpha0.01delta0epsilon0.1dk10D5dopt5L50-art--sx\results-Till350Step1v61-alpha0.01-noOBBExpand-noBondExpand-expvCustom0-1core-small.mat');
res{2,1} = load('20150827-1936-SpinBoson2CT-VT-OrthPol-v61TCMde10-alpha0.01delta0epsilon0.1dk20D5dopt5L50-art--sx\results-Till350Step1v61-alpha0.01-noOBBExpand-noBondExpand-expvCustom0-1core-small.mat');
res{3,1} = load('20150827-1936-SpinBoson2CT-VT-OrthPol-v61TCMde10-alpha0.01delta0epsilon0.1dk40D5dopt5L50-art--sx\results-Till350Step1v61-alpha0.01-noOBBExpand-noBondExpand-expvCustom0-1core-small.mat');
res{4,1} = load('20150828-1406-SpinBoson2CT-VT-OrthPol-v61TCMde9-alpha0.01delta0epsilon0.1dk10D5dopt5L50-art--sz\results-Till350Step1v61-alpha0.01-noOBBExpand-noBondExpand-expvCustom0-1core-small.mat');
res{5,1} = load('20150828-1406-SpinBoson2CT-VT-OrthPol-v61TCMde9-alpha0.01delta0epsilon0.1dk20D5dopt5L50-art--sz\results-Till350Step1v61-alpha0.01-noOBBExpand-noBondExpand-expvCustom0-1core-small.mat');
for i = 1:length(res)
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('d_k = %g', res{i,1}.para.dk(1,2));
end

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph = cellfun(@(x) plot(1:size(x.tresults.n,2), x.tresults.n(x.tresults.lastIdx,:,1)), res(pick,1), 'UniformOutput', false);
	ph = cellfun(@(x) plot(1:size(x.tresults.n,2), x.tresults.n(1,:,1)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('Site $k$');
	ylabel('$\left<n_k\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%s\n',x),res(pick,2),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
end

%% THERM SBM multi load : v62 T=300K Benchmark, s=1 a=0.01							% LabBook 07/09/2015

clear
defPlot(1,:) = {'20150907-SBM2CT-dt01-NewPreExpansion-dk70',		[ 1: 7], {'ylim',[0,0.3],'xlim',[0,20],'yscale','lin'}};
defPlot(2,:) = {'20150907-SBM2CT-dt01-NewPreExpansion-dk40',		[ 8: 13], {'ylim',[0,0.3],'xlim',[0,20],'yscale','lin'}};
defPlot(3,:) = {'20150907-SBM2CT-dt01-NewPreExpansion-D60-80',		[4:7,11:14], {'ylim',[0,0.3],'xlim',[0,20],'yscale','lin'}};
defPlot(4,:) = {'20150912-SBM2CT-dt01-NewPreExpansion-THTDVP',		[15: 21], {'ylim',[0.84,1],'xlim',[0,70],'yscale','lin'}};
i=0; cols = 5;

%1-7: Thermal Evolution from -sx
foldPattern = '20150907-0014-SpinBoson2CT-VT-OrthPol-v62TCMde10-alpha0.01delta0epsilon0.1dk70D5dopt5L20-art--sz';
filePattern = 'results-Till20Step0.1v62-OBBmax*-Dmax*-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('$D_{max} = %g$', res{i,1}.para.tdvp.maxBondDim);
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%8-13: Thermal Evolution from -sx
foldPattern = '20150907-1114-SpinBoson2CT-VT-OrthPol-v62TCMde11-alpha0.01delta0epsilon0.1dk40D5dopt5L20-art--sz';
filePattern = 'results-Till20Step0.1v62-OBBmax*-Dmax*-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('$D_{max} = %g$', res{i,1}.para.tdvp.maxBondDim);
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%15-21: Thermal Evolution from -sx
foldPattern = '20150907-0014-SpinBoson2CT-VT-OrthPol-v62TCMde10-alpha0.01delta0epsilon0.1dk70D5dopt5L20-art--sz';
filePattern = 'results-Till350Step1v62-alpha0.01-OBBmax*-Dmax*-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = sprintf('$D_{max} = %g$', res{i,1}.para.tdvp.maxBondDim);
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph = cellfun(@(x) plot(1:size(x.tresults.n,2), x.tresults.n(x.tresults.lastIdx,:,1)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('Site $k$');
	ylabel('$\left<n_k\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%s\n',x),res(pick,2),'UniformOutput',false),'location','best');
	leg.Interpreter = 'latex';
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	% plot analytic prediction
	para = res{pick(1),1}.para;
	n = para.L; beta = 40;
	sz= -1;
	A = full(sparse(1:n,1:n,[-para.hz*sz/2; para.chain{1}.epsilon],n,n)+...
		sparse(2:n,1:n-1,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n)+...
		sparse(1:n-1,2:n,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n));
	[V,D] = eig(A);
	occDiag = 1./(exp(beta.*diag(D))-1); occDiag(1) = 0;

	Aocc = V*diag(occDiag)*V.'; Aocc(1,1) = 0;				% 1-site is meaningless!
	occN = diag(Aocc);
	plot(occN);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% THERM SBM multi load : v61 T=300K, s=1 a=0.01, delta ~= 0							% LabBook 07/09/2015
% Only contains THERM with tunneling
clear
defPlot(1,:) = {'20150905-SBM2CT-dt01-WithTunneling',					[ 1: 4], {'xlim',[0,0.3],'yscale','lin'}};
i=0; cols = 5;

%1-4: Thermal Evolution from -sz
foldPattern = '20150905-1707-SpinBoson2CT-VT-OrthPol-v61TCM51-alpha0.01delta*epsilon0.1dk40D10dopt10L20-art--sz';
filePattern = 'results-Till20Step0.05v61-OBBmax20-Dmax20-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = res{i,1}.para.chain{1}.alpha;
	res{i,4} = res{i,1}.para.L;
	res{i,2} = -res{i,1}.para.hx;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'dt0.1'};

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	ph = plot(cell2mat(res(pick,2)),...
			  cell2mat(cellfun(@(x) x.tresults.spin.sz(end), res(pick,1),'UniformOutput',false)));
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\Delta$');
	ylabel('$\left<\sigma_z\right>$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
end

%% DPMES v63/v64, Vtens vs StarMPS													% LabBook 07/10/2015
% Purpose to Benchmark StarMPS
clear
defPlot(1,:) = {'20151010-DPMES-v63-Vtens-dk40-300',			[ 1: 2], {'xlim',[0,100],'yscale','lin'}};
defPlot(2,:) = {'20151010-DPMES-v64-StarMPS',					[ 3: 6], {'xlim',[0,100],'yscale','lin'}};
defPlot(3,:) = {'20151006-DPMES-v63-Vtens-dk100',				[ 7:11], {'xlim',[0,100],'yscale','lin'}};
defPlot(4,:) = {'20151012-DPMES-v64-StarMPS-Dsweep',			[12:15], {'xlim',[0,330],'yscale','lin'}};
% defPlot(4,:) = {'20151012-DPMES-v64-StarMPS-Dsweep',			[12:16], {'xlim',[0,100],'yscale','lin'}};
defPlot(5,:) = {'20151016-DPMES-v64-StarMPS-DsweepChain-D10',	[16:19], {'xlim',[0,1e3],'yscale','lin'}};
defPlot(6,:) = {'20151019-DPMES-v64-StarMPS-DsweepChain-D30',	[20:23], {'xlim',[0,330],'yscale','lin'}};
defPlot(7,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D30',		[24:27], {'xlim',[0,330],'yscale','lin'}};
defPlot(8,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D20',		[28:31], {'xlim',[0,330],'yscale','lin'}};
defPlot(9,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D10',		[32:35], {'xlim',[0,330],'yscale','lin'}};
i=0; cols = 5;

%1-2: Old Vtens code
foldPattern = '20151007-1429-DPMES3-4C-VT-v63TCMde9-dk40D10dopt5L11';
filePattern = 'results-Till150Step*v63-OBBmax20-Dmax60-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'Vtens, OBB20,D60'};

%3-6: New StarMPS code
foldPattern = '20151010-2042-DPMES3-4C-Star-v64TCMde9-dk40D10dopt5L11';
filePattern = 'results-Till150Step*v64-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, OBB20,D60'};

%7-11: Vtens
foldPattern = '20151006-1741-DPMES3-4C-VT-v63TCMde9-dk100D20dopt5L11';
filePattern = 'results-Till150Step*v63*-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'Vtens, dk100'};

%12-15: Star
foldPattern = '20151012-0208-DPMES3-4C-Star-v64TCMde9-dk60D5dopt5L11';
filePattern = 'results-Till500Step0.1v64-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, dk60'};

%16-19: Star
foldPattern = '20151016-1518-DPMES3-4C-Star-v64TCMde9-dk100D10dopt5L11';
filePattern = 'results-Till1500Step0.1v64-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim(end), res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, dk100, D1-10'};

%20-23: Star
foldPattern = '20151019-1216-DPMES3-4C-Star-v64TCMde9-dk100D10dopt5L11';
filePattern = 'results-Till500Step0.1v64-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim(end), res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, dk100,D1-30'};

%24-27: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax*0-Dmax*0-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim(end), res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, L3, D1-30'};

%28-31: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax*0-Dmax*0-expvCustom699-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim(end), res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, L3, D1-20'};

%32-35: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax*0-Dmax*0-expvCustom698-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, D%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,1}.para.tdvp.maxBondDim(end), res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, L3, D1-10'};


for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	col = ax.ColorOrder;
	ax.ColorOrder = kron(col,[1;1;1]);
	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = [ph{:}];
	leg = legend(phArr(1,:),res{pick,3},'location','Northwest');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
end
%% DPMES v63/v64, Vtens vs StarMPS - TDVPData										% LabBook 07/10/2015
% Purpose to Benchmark StarMPS
clear
defPlot(1,:) = {'20151010-DPMES-v63-Vtens-dk40-300',			[ 1: 2], {'xlim',[0,100],'yscale','lin'}};
defPlot(2,:) = {'20151010-DPMES-v64-StarMPS',					[ 3: 6], {'xlim',[0,100],'yscale','lin'}};
defPlot(3,:) = {'20151006-DPMES-v63-Vtens-dk100',				[ 7:11], {'xlim',[0,100],'yscale','lin'}};
defPlot(4,:) = {'20151012-DPMES-v64-StarMPS-Dsweep',			[12:15], {'xlim',[0,330],'yscale','lin'}};
% defPlot(4,:) = {'20151012-DPMES-v64-StarMPS-Dsweep',			[12:16], {'xlim',[0,100],'yscale','lin'}};
defPlot(5,:) = {'20151016-DPMES-v64-StarMPS-DsweepChain-D10',	[16:19], {'xlim',[0,1e3],'yscale','lin'}};
defPlot(6,:) = {'20151019-DPMES-v64-StarMPS-DsweepChain-D30',	[20:23], {'xlim',[0,330],'yscale','lin'}};
defPlot(7,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D30',		[24:27], {'xlim',[0,330],'yscale','lin'}};
defPlot(8,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D20',		[28:31], {'xlim',[0,330],'yscale','lin'}};
defPlot(9,:) = {'20151019-DPMES-v64-StarMPS-L3DSweep-D10',		[32:35], {'xlim',[0,330],'yscale','lin'}};
i=0; cols = 5;
%
TDVPfolds = TDVPData.getTDVPLib();
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

%1-2: Old Vtens code
foldPattern = '20151007-1429-DPMES3-4C-VT-v63TCMde9-dk40D10dopt5L11';
filePattern = 'results-Till150Step.*v63-OBBmax20-Dmax60-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim, res(i).dt));
	res(i) = res(i).setComment('Vtens, OBB20,D60');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%3-6: New StarMPS code
foldPattern = '20151010-2042-DPMES3-4C-Star-v64TCMde9-dk40D10dopt5L11';
filePattern = 'results-Till150Step.*v64-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, OBB20,D60');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%7-11: Vtens
foldPattern = '20151006-1741-DPMES3-4C-VT-v63TCMde9-dk100D20dopt5L11';
filePattern = 'results-Till150Step.*v63.*-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim, res(i).dt));
	res(i) = res(i).setComment('Vtens, dk100');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);


%12-15: Star
foldPattern = '20151012-0208-DPMES3-4C-Star-v64TCMde9-dk60D5dopt5L11';
filePattern = 'results-Till500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, dk60');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%16-19: Star
foldPattern = '20151016-1518-DPMES3-4C-Star-v64TCMde9-dk100D10dopt5L11';
filePattern = 'results-Till1500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, dk100, D1-10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);


%20-23: Star
foldPattern = '20151019-1216-DPMES3-4C-Star-v64TCMde9-dk100D10dopt5L11';
filePattern = 'results-Till500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, dk100, D1-30');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%24-27: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, L3, D1-30');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%28-31: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom699-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, L3, D1-20');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%32-35: Star
foldPattern = '20151019-1318-DPMES3-4C-Star-v64TCMde10-dk100D10dopt5L3';
filePattern = 'results-Till500Step0.1v64-OBBmax.*0-Dmax.*0-expvCustom698-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('OBB%d, D%d, dt%g',res(i).para.tdvp.maxOBBDim,res(i).para.tdvp.maxBondDim(end), res(i).dt));
	res(i) = res(i).setComment('StarMPS, L3, D1-10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(arrayfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	col = ax.ColorOrder;
	ax.ColorOrder = kron(col,[1;1;1]);
	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = [ph{:}];
	leg = legend(phArr(1,:),res(pick).LegLabel,'location','Northwest');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
end

%% DPMES4-5C v64  StarMPS - TDVPData												% LabBook 09/11/2015
% Purpose to Benchmark StarMPS
clear
defPlot(1,:) = {'20151109-DPMES-v64-L7DSweep-D10-2ndParams',				[ 1: 2], {'xlim',[0,1.5e3],'yscale','lin'}};
defPlot(2,:) = {'20151109-DPMES4-5C-v64-L7DSweep-D10-3rdParams',			[ 3: 7], {'xlim',[0,1.5e3],'yscale','lin'}};
defPlot(3,:) = {'20151110-DPMES4-5C-v65-L7DSweep-D10-3rdParams-high',		[ 8:13], {'xlim',[0,1.5e3],'yscale','lin'}};
defPlot(4,:) = {'20151113-DPMES4-5C-v65-L7DSweep-D10-3rdParams-high-dk200',	[14:16], {'xlim',[0,1.5e3],'yscale','lin'}};
% defPlot(5,:) = {'20151123-DPMES4-5C-v66-L7DSweep-D10-4thParams-high',		[17:20], {'xlim',[0,1.5e3],'yscale','lin'}};
% defPlot(6,:) = {'20151123-DPMES4-5C-v66-L7DSweep-D10-4thParams-high-v2',	[21:24], {'xlim',[0,1.5e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

%1: 3-4C simulation
foldPattern = '20151103-1203-DPMES3-4C-Star-v64-dk100D10dopt5L7';
filePattern = 'results-Till1500Step0.1v64-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%d, 3-4C',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('dk100, OBB60, D10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%2: 4-5C simulation
foldPattern = '20151104-1552-DPMES4-5C-Star-v64-dk100D10dopt5L7';
filePattern = 'results-Till1500Step0.1v64-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%d, 4-5C',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('dk100, OBB60, D10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 3-7: Params v3 simulation
foldPattern = '20151109-1356-DPMES4-5C-Star-v64TCMde10-dk100D10dopt5L7';
filePattern = 'results-Till1500Step0.1v64-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('dk100, OBB60, D10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 8-13: Params v3-high simulation + HsHi
foldPattern = '20151110-1125-DPMES4-5C-Star-v65TCMde10-dk100D10dopt5L7';
filePattern = 'results-Till1500Step0.1v65-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('dk100, OBB60, D10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 14-16: Params v3-high simulation + HsHi + dk200/1000
foldPattern = '20151113-1357-DPMES4-5C-Star-v65TCMde10-dk200D10dopt5L7';
filePattern = 'results-Till1500Step0.1v65-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('dk200, OBB60, D10');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% % 17-20: Params v4-high simulation + HsHi
% foldPattern = '20151123-1154-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
% filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
% matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
% res(i+size(matches,1),1) = TDVPData(); offset = i;
% for file = {matches.name}
% 	file = file{1};
% 	i = i+1;
% 	res(i) = TDVPData(file);			% comment first!
% 	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
% 	res(i) = res(i).setComment('Params-v4-high');
% end
% [y,I] = sort([res((offset+1):end).dt]);
% res((offset+1):end,1) = res(offset+I,1);
%
% % 21-24: Params v4-high-v2 simulation + HsHi
% foldPattern = '20151123-1201-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
% filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
% matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
% res(i+size(matches,1),1) = TDVPData(); offset = i;
% for file = {matches.name}
% 	file = file{1};
% 	i = i+1;
% 	res(i) = TDVPData(file);			% comment first!
% 	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
% 	res(i) = res(i).setComment('Params-v4-high-v2');
% end
% [y,I] = sort([res((offset+1):end).dt]);
% res((offset+1):end,1) = res(offset+I,1);
%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = cellfun(@(x) x(1),ph,'UniformOutput',false);
	leg = legend([phArr{:}],res(pick).LegLabel,'location','Northwest');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');ax.XLim(2) = 0.658*ax.XLim(2);
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
end

%% DPMES4-5C v66  StarMPS - TDVPData												% LabBook 23/11/2015
% Purpose to Benchmark StarMPS
clear
defPlot(1,:) = {'20151123-DPMES4-5C-v66-L7DSweep-D10-4thParams-high',		[ 2: 5], {'xlim',[0,1.5e3],'yscale','lin'}};
defPlot(2,:) = {'20151123-DPMES4-5C-v66-L7DSweep-D10-4thParams-high-v2',	[ 6:10], {'xlim',[0,1.5e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v3-high simulation + HsHi: for reference
foldPattern = '20151110-1125-DPMES4-5C-Star-v65TCMde10-dk100D10dopt5L7';
filePattern = 'results-Till1500Step0.1v65-OBBmax60-Dmax.*0-expvCustom698-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v3-high');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-5: Params v4-high simulation + HsHi
foldPattern = '20151123-1154-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax.*0-expvCustom.*-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v4-high');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 6-10: Params v4-high-v2 simulation + HsHi
foldPattern = '20151123-1201-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax.*0-expvCustom.*-1.*core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v4-high-v2');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);
%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = cellfun(@(x) x(1),ph,'UniformOutput',false);
	leg = legend([phArr{:}],res(pick).LegLabel,'location','Northwest');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	res(1).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
end

%% DPMES4-5C v66  StarMPS Disorder- TDVPData										% LabBook 04/12/2015
% Purpose to Benchmark StarMPS
clear
defPlot(1,:) = {'20151204-DPMES4-5C-v66-disorder10-D10-4thParams-highv2',		[ 2: 33], {'xlim',[0,1e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v4-high-v2 non-disorder
foldPattern = '20151123-1201-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v4-high-v2');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-17:
foldPattern = '20151204-1732-.*-DPMES4-5C-Star-v66TCMde10-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.2v66-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v4-highv2');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 18-33:
foldPattern = '20151204-1734-.*-DPMES4-5C-Star-v66TCMde11-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.2v66-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D%g',res(i).para.tdvp.maxBondDim(end)));
	res(i) = res(i).setComment('Params-v4-highv2');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-resetColorOrder'), res(pick), 'UniformOutput', false);legend('TT','LE+','CT+','CT-')
% 	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = cellfun(@(x) x(1),ph,'UniformOutput',false);
% 	leg = legend([phArr{:}],res(pick).LegLabel,'location','Northwest');
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
% 	xlabel('$t/fs$');
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	res(1).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
end

%% DPMES4-5C v66  StarMPS Totter - TDVPData											% LabBook 18/12/2015
% Purpose to Benchmark StarMPS Trotter
clear
defPlot(1,:) = {'20151218 - DPMES TDVP Benchmark - StarMPS Trotter',		[ 2: 5], {'xlim',[0,260],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v4-high-v2 non-disorder, No Trotter
foldPattern = '20151218-1624-57-DPMES4-5C-Star-v66-dk20D10dopt5L8';
filePattern = 'results-Till400Step0.1v66-OBBmax60-Dmax20-expvCustom700-1core.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('$\\Delta t = %g$, No Trotter',res(i).para.tdvp.deltaT));
	res(i) = res(i).setComment('Params-v4-high-v2, no Trotter');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-5:
foldPattern = '20151218-1626-59-DPMES4-5C-Star-v66-dk20D10dopt5L8';
filePattern = 'results-Till400Step.*v66-OBBmax60-Dmax20-expvCustom700-1core.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('$\\Delta t = %g$',res(i).para.tdvp.deltaT));
	res(i) = res(i).setComment('Params-v4-highv2');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
	set(ax,defPlot{fignum,3}{:});
	set(gca,'color','none');
	[~,ph] = res(1).plot('rhoii','-fsev','-unicol','k',f);	% overlay of old v3-high simulation
	legend(h(1).pl);
end
% export_fig(['img/20151218 - DPMES TDVP Benchmark - Dynamics'],'-transparent','-png','-m2',gca)

% further plots:
figure(3);clf;hold all; res.plot('calctime','-fsev',figure(3),'-xlim',[0,260]); ylabel('Total CPU time in h'); legend toggle;
% export_fig(['img/20151218 - DPMES TDVP Benchmark - Total t'],'-transparent','-png','-m2',gca)
f = figure(4);clf;f.Name = '20151218-DPMES-StarMPS-Trotter-Benchmark-tPerSweep'; resizePlot(f,2,2);hold all; 
res.plot('calctime-d-sec',figure(4),'-fsev','-xlim',[0,260],'-ylim',[0,119]); ylabel('CPU time/sweep/sec'); legend toggle; title('');resizePlot(f,2,2);
% export_fig(['img/',get(gcf,'name')],'-transparent','-painters','-png','-pdf','-m3')
%% Special convergence plot
% as taken from DPMES paper, reference in last position
% 				 png name														res#				ax prop															columnleg prop (#col, boxon,lineratio)
defPlot( 1,:) = {'20151218 - DPMES TDVP Benchmark - StarMPS Trotter',			[2:5,1],			{'xlim',[0,1.5],'ylim',[5e-3,1e-1],'yscale','log'},					{3,0,0.5}};
	% Define the Figure & Plot
h = struct();
h.f = figure(120); clf;
h.f.Name = defPlot{1,1};
h.ax = gca;
set(h.f,'DefaultAxesFontSize', 8,...
		'DefaultLineLineWidth',1);
resizePlot(h.f);
col = get(0,'defaultaxescolororder');		% not necessary?
hold all; box on;

	% Do Plots and copy
ii = 1; pick = defPlot{ii,2}; axes(h.ax(ii));
	% extract the TT populations
ttPop = cell2mat(arrayfun(@(x) TDVPData.subsample(reshape(real(x.getData('rhoii')),[],1),x.dt,1),res(pick)','UniformOutput',false));
	% estimate mean error
ttErr = abs(bsxfun(@minus,ttPop,ttPop(:,end)));
ttErr = bsxfun(@rdivide,ttErr,ttPop(:,end)); ttErr(isnan(ttErr)) = 0;		% makes error relative
ttMeanErr = (mean(ttErr,1));
	% get Runtimes
tSweep = arrayfun(@(x) 60*mean(diff(nonzeros(x.para.tdvp.calcTime))),res(pick));
h.pl = scatter(tSweep, ttMeanErr,'filled');
set(h.ax(ii), defPlot{ii,3}{:});
xlabel('$t_{CPU}$/Sweep/min');
ylabel('$\varepsilon(\hat\rho_{\mathrm{TT}})$');

	% shift text labels
offX = +ones(1,length(pick))*0.1; offY = ones(1,length(pick))*0;
% offX(5) = 0.2; offY(5) = 0;
% offX(7) = -1.0; offY(7) = -0.005;
% offX(8) = -0.1; offY(8) = 0.008;
% offX(9) = -0.6; offY(9) = 0.003;
% offY = offY;
for kk = 1:length(pick)
	t = text(tSweep(kk)+offX(kk),double(ttMeanErr(kk))+offY(kk),strrep(res(pick(kk)).LegLabel,', ','\\\\'),'FontSize',8);
end
set(gca,'color','none');
t = text(0.7,0.05,'$\varepsilon(\hat\rho_{\mathrm{TT}}) = \mathrm{mean}(\frac{\rho_{\mathrm{TT}}-\hat\rho_{\mathrm{TT}}}{\rho_{\mathrm{TT}}})$');
t = text(0.7,0.025,'$\rho_{\mathrm{TT}}: (120,30)$','FontSize',8);
t = text(0.7,0.013,'$t_{CPU}/\mathrm{Sweep} = 10.2$ min','FontSize',8);

%% DPMES4-5C v66-v72  StarMPS H correction - TDVPData								% LabBook 29/01/2016
% test v72 against v66 and effect of H_int correction
clear
defPlot(1,:) = {'20160129-DPMES4-5C-v66-v72-D10-4thParams-highv2-H',		[ 2: 3], {'xlim',[0,1e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v4-high-v2 non-disorder, wrong Hamiltonian factor
foldPattern = '20151123-1201-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('H bad v66'));
	res(i) = res(i).setComment('wrong H factor');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2: Params v4-high-v2 non-disorder, correct Hamiltonian factor
foldPattern = '20160129-1321-DPMES4-5C-Star-v66TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v66-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('H good v66'));
	res(i) = res(i).setComment('H corrected v66');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 3: Params v4-high-v2 v72 correct H factor
foldPattern = '20160129-1332-11-DPMES4-5C-Star-v72TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('H good v72'));
	res(i) = res(i).setComment('H corrected v72');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);
%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-resetColorOrder'), res(pick), 'UniformOutput', false);legend('TT','LE+','CT+','CT-')
% 	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = cellfun(@(x) x(1),ph,'UniformOutput',false);
% 	leg = legend([phArr{:}],res(pick).LegLabel,'location','Northwest');
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	res(1).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
end

%% DPMES4-5C v72  StarMPS CT shift - TDVPData										% LabBook 29/01/2016
% See effect of CT shift on dynamics
clear
defPlot(1,:) = {'20160129-DPMES4-5C-v72-D10-4thParams-highv2-CTshift',		[ 2: 5], {'xlim',[0,1e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v4-high-v2 v72 correct H factor
foldPattern = '20160129-1332-11-DPMES4-5C-Star-v72TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('CT unshifted'));
	res(i) = res(i).setComment('H corrected v72');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-5: Params v4-high-v2 v72 shifted CT
foldPattern = '20160129-1356-57-DPMES4-5C-Star-v72TCMde9-dk20D10dopt5L7Delta0..';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('CT +%s',res(i).para.folder(end-2:end)));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-resetColorOrder'), res(pick), 'UniformOutput', false);legend('TT','LE+','CT+','CT-')
% 	ph = arrayfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx)*0.658, abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
% 	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), abs(x.tresults.PPCWavefunction(1:x.tresults.lastIdx,[1,2,3])),...
% 				'DisplayName',sprintf('D%d',x.para.tdvp.maxBondDim(end))), res(pick,1), 'UniformOutput', false);
	axis tight;
	phArr = cellfun(@(x) x(1),ph,'UniformOutput',false);
% 	leg = legend([phArr{:}],res(pick).LegLabel,'location','Northwest');
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	res(1).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
end

%% DPMES5-7C v72  StarMPS LE- added - TDVPData										% LabBook 29/01/2016
% See effect of CT shift on dynamics
clear
defPlot(1,:) = {'20160129-DPMES5-7C-v72-D5-10-5param',				[ 2: 4], {'xlim',[0,1e3],'yscale','lin'}};
defPlot(2,:) = {'20160129-DPMES5-7C-v72-D5-10-5param-CTshift',		[ 5: 7], {'xlim',[0,1e3],'yscale','lin'}};
defPlot(3,:) = {'20160203-DPMES5-7C-v72-D5-10-5param-TT-CTshift',	[8:10], {'xlim',[0,1e3],'yscale','lin'}};
defPlot(4,:) = {'20160211-DPMES5-7C-v72-D5-10-5param-TT-CTshift',	[11:13], {'xlim',[0,1e3],'yscale','lin'}};

i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v4 v72, 4-5C
foldPattern = '20160129-1332-11-DPMES4-5C-Star-v72TCMde9-dk20D10dopt5L7';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('4-5C'));
	res(i) = res(i).setComment('H corrected v72');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-4: Params v5 v72 5-7C
foldPattern = '20160129-1417-29-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L7Delta0';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Ds.*-Dmax10-expvCustom700-1core-small';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	res(i) = res(i).setLegLabel(sprintf('D_s %g',res(i).para.tdvp.maxBondDim(1)));
	res(i) = res(i).setComment('5-7C');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 5-7: Params v5 v72 5-7C, CTshift
foldPattern = '20160201-1613-27-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L7Delta0.*';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Ds.*-Dmax10-expvCustom700-1core-small';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta(\w*.\w*)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s',temp{1}{1}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 8-10: Params v5 v72 5-7C, StateProj, TT init, CT shift
foldPattern = '20160203-1659-19-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L7Delta0.*State.';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta([0-9.]*).*State(\d)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s, State %s',temp{1}{:}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%11-13: Params v5 v72 5-7C, StateProj, TT init, CT shift - BETTER
foldPattern = '20160211-1819-38-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L8Delta0.*State1';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta([0-9.]*).*State(\d)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s, State %s',temp{1}{:}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum+10); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-resetColorOrder'), res(pick), 'UniformOutput', false);legend('TT','LE+','LE-','CT+','CT-')
	ph = res(pick).plot('rhoii','-fsev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
% 	ph = res(pick).plot('rhoii-osc-res','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	axis tight;
% 	leg = legend(ph(:,1),res(pick).LegLabel,'location','Northwest');		% leg for each series
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$t/fs$');
% 	xlabel('$t$');
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum+10,'twocolumn-single');
	set(gca,'color','none');
	grid on
	if fignum == 1
		res(1).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
	else
		res(3).plot('rhoii','-fsev','-unicol','k.');	% overlay of old v3-high simulation
	end
	drawnow
end
%%	Plot each dynamics, and FT of residual of following spectra
sets = [3,5,6,7,11,12,13];
% sets = 3
for fignum = 1:length(sets)
	f = figure(fignum*10+1); clf; hold all; ax = gca;
	pick = sets(fignum);
	ph = res(pick).plot('rhoii','-fsev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	formatPlot(f.Number,'twocolumn-single');
	set(gca,'color','none');
	grid on
	
	f = figure(fignum*10+2); clf; hold all; ax = gca;
	[ph,res(pick)] = res(pick).plot('rhoii-osc-res','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	formatPlot(f.Number,'twocolumn-single');
	set(gca,'color','none');
	ax.XLim = [0,1600];
	grid on

	f = figure(fignum*10+3); clf; hold all; ax = gca;
	ph = res(pick).plot('rhoii-osc-res','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	for kk = 1:numel(ph)
		ph(kk).YData = ph(kk).YData./max(ph(kk).YData);
	end
	formatPlot(f.Number,'twocolumn-single');
	set(gca,'color','none');
	ax.XLim = [0,1600];
	grid on

end


%% DPMES5-7C v72  StarMPS LinAbs - TDVPData											% LabBook 12/01/2016
% See effect of CT shift on dynamics
clear
%
defPlot(1,:) = {'20160212-DPMES5-7C-v72-D5-10-5param-LE-CTshift-linAbs',	[ 1:4], {'xlim',[400,800],'yscale','lin'}};
defPlot(2,:) = {'20160212-DPMES5-7C-v72-D5-10-5param-TT-CTshift-linAbs',	[ 5:7], {'xlim',[500,1e3],'yscale','lin'}};
%
i=0; cols = 5;
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
res = TDVPData();

% 1: Params v5 v72 5-7C, LE
foldPattern = '20160203-1659-19-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L7Delta0State2';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta([0-9.]*).*State(\d)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s, State %s',temp{1}{:}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 2-4: Params v5 v72 5-7C, LE, CT shift
foldPattern = '20160210-0109-3.-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L8Delta[-.0-9]*State2';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta([0-9.-]*).*State(\d)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s, State %s',temp{1}{:}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

% 5-7: Params v5 v72 5-7C, TT, CT shift - BETTER
foldPattern = '20160211-1819-38-DPMES5-7C-Star-v72TCMde9-dk20D5dopt5L8Delta0.*State1';
filePattern = 'results-Till1500Step0.1v72-OBBmax60-Dmax10-expvCustom700-1core-small.mat';
matches     = TDVPfolds(arrayfun(@(x) ~isempty(regexp(x.name,[foldPattern,'\\',filePattern],'once')),TDVPfolds));
res(i+size(matches,1),1) = TDVPData(); offset = i;
for file = {matches.name}
	file = file{1};
	i = i+1;
	res(i) = TDVPData(file);			% comment first!
	temp = regexp(res(i).para.folder,'Delta([0-9.]*).*State(\d)','tokens');	% contains strings extracted from folder name
	res(i) = res(i).setLegLabel(sprintf('CT +%s, State %s',temp{1}{:}));
	res(i) = res(i).setComment('CT shifted');
end
[y,I] = sort([res((offset+1):end).dt]);
res((offset+1):end,1) = res(offset+I,1);

%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	ph = res(pick).plot('linabs','-nmev');
	axis tight;
	leg = legend(ph(:,1),res(pick).LegLabel,'location','Northwest');		% leg for each series
 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
end
	
%%	find correlation between TT trace and "parameter distance"
chainPara = [];
for ii = 1:length(res)
	chainPara = [chainPara;res(ii).chainParaInfo()];

end
d = size(chainPara);

diffPara = chainPara;
diffParaNorm = chainPara;
for ii = 1:d(1)
	for jj = 1:d(2)
		for kk = {'xi','Gamma','epsilon','t'}
% 			diffPara(ii,jj).(kk{1}) = norm((chainPara(ii,jj).(kk{1})-chainPara(1,jj).(kk{1})));									% absolute diff
			diffPara(ii,jj).(kk{1}) = chainPara(ii,jj).(kk{1})./chainPara(1,jj).(kk{1})-1;	% normalised relative diff
			diffParaNorm(ii,jj).(kk{1}) = mean(abs(diffPara(ii,jj).(kk{1})));	% normalised relative diff
		end
	end
end

% norm of disorder for each chain / parameter set
for kk = {'xi','Gamma','epsilon','t'}
	diffParaMat.(kk{1}) = reshape([diffParaNorm.(kk{1})],d);
end

% find very good and very bad results:
for ii = 1:length(res)
	rhoii = res(ii).gettRhoiiSystem();
	resTT(ii).TT     = abs(rhoii(1:res(ii).tresults.lastIdx,1));
	resTT(ii).TTtail = resTT(ii).TT(res(ii).para.tdvp.t(1:res(ii).tresults.lastIdx)>255);
	resTT(ii).TTmean = mean(resTT(ii).TTtail);
	resTT(ii).TTstd  = std(resTT(ii).TTtail);
end

%% TDVP SBM2C Vtens-StarMPS Benchmark												% LabBook 12/10/2015
% Modified SBM2C to study pure dephasing
clear
defPlot(1,:) = {'20151012-SBM2C-v63-Vtens',					[ 1: 2], {'xlim',[0,20],'yscale','lin'}};
defPlot(2,:) = {'20151012-SBM2C-v64-StarMPS',				[ 3: 8], {'xlim',[0,20],'yscale','lin'}};
defPlot(3,:) = {'20151013-SBM2C-v64-StarMPS-s1-3',			[ 9:10], {'xlim',[0,500],'yscale','lin'}};
i=0; cols = 5;

%1-2: Old Vtens code
foldPattern = '20151012-1502-SpinBoson2C-VT-OrthPol-v64TCM73-alpha0.1delta0epsilon0.1dk20D5dopt5L50-art--sx';
filePattern = 'results-Till20Step*v64-OBBmax30-expvCustom0-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, dt%g',res{i,1}.para.tdvp.maxOBBDim,res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'Vtens, OBB30'};

%3-8: New StarMPS code
foldPattern = '20151012-1502-SpinBoson2C-Star-OrthPol-v64TCM73-alpha0.1delta0epsilon0.1dk20D5dopt5L50-art--sx';
filePattern = 'results-Till20Step*v64-OBBmax*0-expvCustom*-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, dt%g',res{i,1}.para.tdvp.maxOBBDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, OBB30'};

%9-10: New StarMPS code
foldPattern = '20151013-0114-SpinBoson2C-Star-OrthPol-v64TCM74-*alpha0.1delta0epsilon0.1dk100D5dopt5L250-art--sx';
filePattern = 'results-Till500Step0.2v64-OBBmax40-Dmax30-expvCustom700-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res{i+size(folds,1),cols} = []; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,2} = res{i,1}.para.tdvp.deltaT;
	res{i,3} = sprintf('OBB%d, dt%g',res{i,1}.para.tdvp.maxOBBDim, res{i,2});
	res{i,4} = res{i,1}.para.L;
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),[2 4]);
res(offset+1:i,5) = {'StarMPS, OBB30'};

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(x.tresults.lastIdx), res(pick,1)));
	plot([0,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.visibility(1:x.tresults.lastIdx)), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,3},'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum,'twocolumn-single');
% 	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
% 	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
	if fignum ~= 1
		ax.ColorOrderIndex = 1;
		ph = cellfun(@(x) plot(x.tresults.t(1:x.tresults.lastIdx), x.tresults.spin.visibility(1:x.tresults.lastIdx),'.black'), res(defPlot{1,2},1), 'UniformOutput', false);
	end
end

% plot time taken
figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% DPMES5-7C v73  TreeMPS from LE+ & TT, L18, L2 CTshift - TDVPData						% LabBook 15/03/2016
% See effect of CT shift on dynamics form LE+ and TT
clear
defPlot(1,:) = {'20160315-DPMES5-7C-v73-D5-20-5param-LE-CTshift',		[ 1: 6], {'xlim',[0,1e3],'yscale','lin'}};
% defPlot(2,:) = {'20160315-DPMES5-7C-v73-D5-20-5param-TT-CTshift',		[ 6: 9], {'xlim',[0,1e3],'yscale','lin'}};
% defPlot(3,:) = {'20160315-DPMES5-7C-v73-D5-L2-5param-LE-CTshift',		[10:15], {'xlim',[0,3.3e3],'yscale','lin'}};
% defPlot(4,:) = {'20160315-DPMES5-7C-v73-D5-L2-5param-TT-CTshift',		[16:19], {'xlim',[0,3.3e3],'yscale','lin'}};

%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1: from LE+, L=18 CT up
dirPat = '20160312-0432-5.-DPMES5-7C-Tree-v73TCMde9-L18CT.*LE.*';
filPat = 'results-Till1500Step0.1v73-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
matches = [matches; m(I)];

% 
% % 2: from TT, L=18
% dirPat = '20160312-0432-5.-DPMES5-7C-Tree-v73TCMde9-L18CT.*TT.*';
% filPat = 'results-Till1500Step0.1v73-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
% m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
% matches = [matches; m(I)];

% % 3: from LE, L=2
% dirPat = '20160315-2359-07-DPMES5-7C-Tree-v73TCMde10-L2CT.*LE.*';
% filPat = 'results-Till5000Step0.1v73-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
% m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
% matches = [matches; m(I)];
% 
% % 4: from TT, L=2
% dirPat = '20160315-2359-07-DPMES5-7C-Tree-v73TCMde10-L2CT.*TT.*';
% filPat = 'results-Till5000Step0.1v73-OBBmax60-Dmax20-expvCustom700-1core-small.mat';
% m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
% matches = [matches; m(I)];

% % 5: from LE, L=2, more observables dt=1
% dirPat = '20160407-2035-50-DPMES5-7C-Tree-v73TCMde9-L2CT.*LE\+';
% filPat = 'results-Till10000Step0.1v73-OBBmax60-Dmax\(5-20\)-expvCustom700-1core-small.mat';
% m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
% matches = [matches; m(I)];

% 6: from TT, L=2, more observables dt=1
% dirPat = '20160407-2035-4.-DPMES5-7C-Tree-v73TCMde10-L2CT.*TT';
% filPat = 'results-Till10000Step0.1v73-OBBmax60-Dmax\(5-20\)-expvCustom700-1core-small.mat';
% m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
% matches = [matches; m(I)];

% 7: 1-6: from LE+, L=18 CT down
dirPat = '20160516-2007-43-DPMES5-7C-Tree-v73TCMde11-L18CT.*LE.';
filPat = 'results-Till5000Step0.1v73-OBBmax60-Dmax\(5-20\)-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'CT([-.0-9]*)','tokens');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x{1}),tokens));
matches = [matches; m(I)];

res = TDVPData({matches.name});
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
% 	ph = arrayfun(@(x) x.plot('rhoii','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-unicol'), res(pick), 'UniformOutput', false);
% 	ph = arrayfun(@(x) x.plot('rhoii','-fsev','-resetColorOrder'), res(pick), 'UniformOutput', false);legend('TT','LE+','LE-','CT+','CT-')
	ph = res(pick).plot('rhoii','-fsev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
% 	ph = res(pick).plot('rhoii-osc-res','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	axis tight;
% 	leg = legend(ph(:,1),res(pick).LegLabel,'location','Northwest');		% leg for each series
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
% 	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	drawnow
	if fignum == 2
		ph = res(2).plot('rhoii','-fsev','-unicol','k.');
	end
end
	%% Correlate FT of pop and occ
	corrCT = zeros(5,7,5);				% states x chains x CT shifts
for ii = 1:5
	m        = res(ii).lastIdx;
	dt       = res(ii).t(2);
	t        = res(ii).t;
	rhoii    = abs(res(ii).getData('rhoii'));
	rhoiiAvg = TDVPData.movAvg(rhoii,170/dt);
	rhoiiRes = rhoii-rhoiiAvg;
	occRC    = squeeze(real(res(ii).occC(1:m,2,:))); % t x L x nChain
	occRCAvg = TDVPData.movAvg(occRC,760/dt);
	occRCRes = occRC-occRCAvg;
	
% 	f = figure(1);hold all; plot(t,rhoii); 	plot(t,rhoiiAvg);
% 	f = figure(2);hold all; plot(t,occRC);	plot(t,occRCAvg);
% 	
	% FFT of residuals
% 	window = hann(m,'periodic');
	window = kaiser(m,2.5);
	maxN = pow2(nextpow2(20*m));			% if > lastIdx -> zero padding
	freq = (0:maxN-1)/dt/maxN;
	freq = freq/0.658*4.135*8065.73;		% to cm
	rhoiiFT = fft(rhoiiRes.*(window*ones(1,size(rhoiiRes,2))),maxN,1);
	occRCFT = fft(occRCRes.*(window*ones(1,size(occRCRes,2))),maxN,1);
	
	f = figure(1);clf; hold all; plot(freq,abs(rhoiiFT)); set(gca,'xlim',[0,2000]);
	f = figure(2);clf; hold all; plot(freq,abs(occRCFT)); set(gca,'xlim',[0,2000]);
	
	% plot correlation, normalised spectra
	rhoiiFTnorm = bsxfun(@rdivide, abs(rhoiiFT), sqrt(sum(abs(rhoiiFT).^2)));
	occRCFTnorm = bsxfun(@rdivide, abs(occRCFT), sqrt(sum(abs(occRCFT).^2)));
	f = figure(3+ii); ax = gca;
	maxFreq = 2000;
	corr = abs(rhoiiFTnorm(freq<maxFreq,:))'*abs(occRCFTnorm(freq<maxFreq,:));
	corrCT(:,:,ii) = corr;
% 	pl = scatter3(reshape(repmat(1:5,1,7),[],1),reshape(repmat(1:7,5,1),[],1),reshape(corr,[],1),100,reshape(corr,[],1),'filled');
	pl = scatter3(reshape(repmat(1:5,1,7),[],1),reshape(repmat(1:7,5,1),[],1),reshape(corr+repmat(0:6,5,1),[],1),100,reshape(corr,[],1),'filled');
	view([0,0]);
	title(res(ii).folder);
	ax.XTick = 1:5; ax.XTickLabel = {'TT','LE^+','LE^-','CT^+','CT^-'};
	ax.YTick = 1:7; ax.YTickLabel = {'A_{1,1}','A_{1,2}','A_2','B_1','B_{2,1}','B_{2,2}','B_{2,3}'};
	ax.ZTick = (1:7); ax.ZTickLabel = {'A_{1,1}','A_{1,2}','A_2','B_1','B_{2,1}','B_{2,2}','B_{2,3}'};
end
%%
states = {'TT','$LE^+$','$LE^-$','$CT^+$','$CT^-$'};
for ii = 1:5
	f = figure(10+ii); clf; hold all;
	plot([-5,0,5,10,20],squeeze(corrCT(ii,:,:))')
% 	legend('$A_{1,1}$','$A_{1,2}$','$A_2$','$B_1$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$');
	xlabel('CT shift (\%)')
	ylabel(['Correlation of $\langle n\rangle$ with ',states{ii}]);
	formatPlot(f);
	export_fig(sprintf('img/%d',f.Number),'-transparent','-png','-m3', f);
end
%% DPMES5-7C v73  TreeMPS from LE+, L18 - D sys chain sweep - TDVPData					% LabBook 25/03/2016
% See effect of CT shift on dynamics form LE+ and TT
clear
defPlot(1,:) = {'20160325-DPMES5-7C-v73-D(5-10-5-60)-5param-LE',		[ 1: 9], {'xlim',[0,3.3e3],'yscale','lin'}};
defPlot(2,:) = {'20160325-DPMES5-7C-v73-D(5-(10-60))-5param-LE',		[ 1: 3], {'xlim',[0,3.3e3],'yscale','lin'}};
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1: from LE+, L=18
dirPat = '20160325-123.-..-DPMES5-7C-Tree-v73TCMde.*-L18CT0LE\+';
filPat = 'results-Till5000Step0.1v73-OBBmax60-Dmax\(.*-.*\)-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax.(?<Dsys>[0-9]*)-(?<Dchain>[0-9]*).-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.Dchain),tokens));
m = m(I)
tokens = regexp({m.name},'Dmax.(?<Dsys>[0-9]*)-(?<Dchain>[0-9]*).-','names');			% second sorting
[y,I] = sort(cellfun(@(x) str2double(x.Dsys),tokens));

matches = [matches; m(I)];

res = TDVPData({matches.name});
%% Legend Labels
legLab = cellfun(@(x) sprintf('D%s-%s',x.Dsys,x.Dchain),tokens(I),'UniformOutput',false)
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	
	[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
	set(ph(1:end,1),{'Displayname'},legLab(pick)');
	legend('TT','LE+','LE-','CT+','CT-')
% 	ph = res(pick).plot('rhoii-osc-res-med','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	axis tight;
% 	leg = legend(ph(:,1),res(pick).LegLabel,'location','Northwest');		% leg for each series
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	drawnow
end
%% Compare Dchain for Dsys=5
f=figure(11); clf; hold all; ax = gca; f.Name = 'DPMES Convergence for Dsys5';
pick = [1:3];
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% Compare Dchain for Dsys=7
f=figure(12); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dsys7';
pick = [4:6];
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% Compare Dchain for Dsys=10
f=figure(13); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dsys10';
pick = [7:9];
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% Compare Dsys for Dchain=5
f=figure(14); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dchain5';
pick = [1,4,7];
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% Compare Dsys for Dchain=20
f=figure(15); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dchain20';
pick = [1,4,7]+1;
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% Compare Dsys for Dchain=60
f=figure(16); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dchain60';
pick = [1,4,7]+2;
[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
%% CPU Time
f=figure(17); clf; hold all; ax = gca;f.Name = 'DPMES Convergence for Dchain Dsys - CPU time';
pick = 1:9;
[h,ph] = res(pick).plot('calctime-d','-fsev',f); ylim([0,5])
% set(ph(2,1:end),{'Color'},{[1 1 1]*0.4});
% set(ph(1,1:end),{'Color'},{[1 1 1]*0.7});
leg = legend(ph(:,1),legLab(pick));
title('')
set(gca,'color','none');
grid on
	

%% DPMES5-7C v73  TreeMPS from LE+ & TT, L2 CTshift 5ps - TDVPData					% LabBook 25/03/2016
% See effect of CT shift on dynamics form LE+ and TT
clear
defPlot(1,:) = {'20160407-DPMES5-7C-v73-L2-5param-LE-CTshift',		[ 1: 6], {'xlim',[0,6.6e3],'yscale','lin'}};
defPlot(2,:) = {'20160407-DPMES5-7C-v73-L2-5param-TT-CTshift',		[ 7:12], {'xlim',[0,6.6e3],'yscale','lin'}};
%%To update library:
%TDVPfolds = TDVPData.getTDVPLib();save('TDVPLib.mat','TDVPfolds');
%
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1-6: from LE+, L=2
dirPat = '20160407-2035-50-DPMES5-7C-Tree-v73TCMde9-L2CT.*LE\+';
filPat = 'results-Till10000Step0.1v73-OBBmax60-Dmax\(5-20\)-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'CT(?<CTshift>[0-9\-\.]*)LE','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.CTshift),tokens));
matches = [matches; m(I)];

% 7-12: from TT, L=2
dirPat = '20160407-2035-4.-DPMES5-7C-Tree-v73TCMde10-L2CT.*TT';
filPat = 'results-Till10000Step0.1v73-OBBmax60-Dmax\(5-20\)-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'CT(?<CTshift>[0-9\-\.]*)TT','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.CTshift),tokens));
matches = [matches; m(I)];

res = TDVPData({matches.name});
%% Legend Labels
legLab = cellfun(@(x) sprintf('CT %s',x.CTshift),tokens(I),'UniformOutput',false)
%%
for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all; ax = gca;
	f.Name = defPlot{fignum,1};
	pick = defPlot{fignum,2};			% plot all
	
	[h,ph] = res(pick).plot('rhoii','-fsev','-resetColorOrder',f);
	set(ph(1:end,1),{'Displayname'},legLab');
	legend('TT','LE+','LE-','CT+','CT-')
% 	ph = res(pick).plot('rhoii-osc-res-med','-cmev','-resetColorOrder');legend('TT','LE+','LE-','CT+','CT-')
	axis tight;
% 	leg = legend(ph(:,1),res(pick).LegLabel,'location','Northwest');		% leg for each series
% 	legend boxoff
% 	fs = 22;
% 	leg.FontSize = fs;
	set(ax,defPlot{fignum,3}{:});
	ylabel('$\rho_{ii} (t)$');
	fs = 22;
	formatPlot(fignum,'twocolumn-single');
	set(gca,'color','none');
	grid on
	drawnow
end

%% DPMES-Tree1 v74 TreeMPS Performance comparison - TDVPData						% LabBook 31/08/2016
clear
% x(1) = TDVPData('20160829-1710-46-DPMES5-7C-Tree-v73-L18CT0LE+\results-Till5000Step0.1v73-OBBmax20-Dmax(5-10)-expvCustom700-1core-small.mat');					% this has corrected H_I
% x(1) = TDVPData('Data\20160325-1238-25-DPMES5-7C-Tree-v73TCMde10-L18CT0LE+\results-Till5000Step0.1v73-OBBmax60-Dmax(10-20)-expvCustom700-1core-small.mat');		% this has still wrong H_I. Chain dimension agree with Tree1, should be closer than (5-10)
x(1) = TDVPData('Data\20160325-1238-25-DPMES5-7C-Tree-v73TCMde10-L18CT0LE+\results-Till5000Step0.1v73-OBBmax60-Dmax(5-20)-expvCustom700-1core-small.mat');		% this has still wrong H_I. sys-chain D agrees with above, chain D with below
% x(3) = TDVPData('20160831-1344-17-DPMES-Tree1-Tree-v73-L18CT0LE+\results-Till1000Step0.1v73-OBBmax20-Dmax(10-20)-expvCustom700-1core-small.mat');				% from OE-Rugor, too small node Dim
x(2) = TDVPData('20160831-1258-02-DPMES-Tree1-Tree-v73-L18CT0LE+\results-Till1000Step0.1v73-OBBmax20-Dmax(100-20)-expvCustom700-1core-small.mat');				% from OE-Rugor, big timestep, yet good!
% x(5) = TDVPData('20160929-1744-01-DPMES-Tree1-Tree-v74-L18CT0LE+\results-Till1000Step0.01v74-OBBmax20-Dmax(100-20)-expvCustom256-1core-small.mat');				% OE-Rugor, good timestep, too slow!
% x(6) = TDVPData('20160929-1745-52-DPMES-Tree1-Tree-v74-L18CT0LE+\results-Till1000Step0.01v74-OBBmax20-Dmax(100-50)-expvCustom256-1core-small.mat');				% OE-Rugor, Bigger chain Dmax, no difference
%x(5) = TDVPData('20161028-1350-50-DPMES-Tree2-Tree-v75-L18CT0LE+\results-Till40Step0.1v75-OBBmax20-Dmax(100-20)-expvCustom700-1core-small.mat');				% from Flop40, with Tree2! old order
% x(3) = TDVPData('20161028-1320-01-DPMES-Tree2-Tree-v75-L18CT0LE+\results-Till40Step0.1v75-OBBmax20-Dmax(100-100)-expvCustom700-1core-small.mat');				% from Flop40, with Tree2! new order: B22,B21
% x(4) = TDVPData('20161028-1404-41-DPMES-Tree2-Tree-v75-L18CT0LE+\results-Till100Step0.1v75-OBBmax20-Dmax(100-20)-expvCustom200-1core-small.mat');				% from Flop40, with Tree2! old order, Lanczos, SV[6,8]
% x(5) = TDVPData('20161028-1404-58-DPMES-Tree2-Tree-v75-L18CT0LE+\results-Till100Step0.1v75-OBBmax20-Dmax(50-20)-expvCustom200-1core-small.mat');				% from Flop40, with Tree2! old order, Lanczos, SV[6,8]
x(3) = TDVPData('20161031-1428-04-DPMES-Tree2-Tree-v76TCMde11-L18CT0LE+\results-Till5000Step0.5v76-OBBmax60-Dmax(70-20)-expvCustom270-1core-small.mat');		% from node11, with Tree2! SV[6,8], expvM = 5
close all
x.plot('rhoii','-fsev','-resetColorOrder',figure)
x.plot('calctime-d-sec',figure)
%%
x(1) = TDVPData('20161023-1556-43-DPMES-Tree1-Tree-v75-L18CT0LE+\results-Till40Step0.1v75-OBBmax20-Dmax(100-100)-expvCustom256-1core-small.mat');				% from Flop40, standard
x(2) = TDVPData('20161102-1814-31-DPMES-Tree1-Tree-v76-L18CT0LE+\results-Till40Step0.1v76-OBBmax20-Dmax(100-100)-expvCustom200-1core-small.mat');				% from Flop40, with expvM = 5
x(3) = TDVPData('20161028-1404-58-DPMES-Tree2-Tree-v75-L18CT0LE+\results-Till100Step0.1v75-OBBmax20-Dmax(50-20)-expvCustom200-1core-small.mat');				% from Flop40, with Tree2! old order, Lanczos, SV[6,8]
close all
x.plot('rhoii','-fsev','-resetColorOrder',figure)
x.plot('calctime-d-sec',figure)
%% DPMES-Tree1 v74 TreeMPS dt convergence - TDVPData								% LabBook 23/09/2016
clear
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1-6: from LE+, L=2
dirPat = '2016090.-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step0..*v74-OBBmax20-Dmax\(100-100\)-expvCustom700-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'CT(?<CTshift>[0-9\-\.]*)LE','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.CTshift),tokens));
matches = [matches; m(I)];

res = TDVPData({matches.name});

close all
%% plot
h = res.plot('rhoii','-fsev','-resetColorOrder',figure)
f = figure(); hold all; ax = gca; f.Name = 'Tree1-dt-convergence2-D(100-100)';
a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false)
ax.XLabel = h(1).ax(1).XLabel;
legend('0.1','0.05','0.02','0.01','0.005','0.002','0.001');
set([a{:}],{'Color'},mat2cell(ax.ColorOrder,ones(1,7),3));
ylabel('$\rho_{ii} (t)$');
set(gca,'color','none');
grid on

% x.plot('calctime-d-sec',figure)
%% DPMES-Tree1 v74 TreeMPS dt convergence finer 5fs - TDVPData								% LabBook 27/09/2016
clear
defPlot(1,:) = {'20160927-DPMES-Tree1-v73-dt-convergence1',						[ 1: 6], {'xlim',[0,2]}};
defPlot(2,:) = {'20160927-DPMES-Tree1-v74-dt-convergence2',						[ 7:12], {'xlim',[0,2]}};
defPlot(3,:) = {'20160927-DPMES-Tree1-v74-dt-convergence3-[120,257]',			[13:31], {'xlim',[0,2]}};
defPlot(4,:) = {'20160927-DPMES-Tree1-v74-dt-convergence4-[0,256]',				[32:50], {'xlim',[0,2]}};
defPlot(5,:) = {'20160927-DPMES-Tree1-v74-dt-convergence5-[0,0]',				[51:69], {'xlim',[0,2]}};
defPlot(6,:) = {'20160928-DPMES-Tree1-v74-dt-convergence6-[120,256]-[6,6.5]',	[70:88], {'xlim',[0,2]}};
defPlot(7,:) = {'20160928-DPMES-Tree1-v74-dt-convergence7-[120,257]-[7,9]',		[89:107], {'xlim',[0,2]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1-6: old
dirPat = '20160926-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step0.*v74-OBBmax20-Dmax\(100-100\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 7-12: old, should be v74 actually!
dirPat = '20160926-.*-DPMES-Tree1-Tree-v73-L18CT0LE\+';
filPat = 'results-Till3Step0.*v73-OBBmax20-Dmax\(50-100\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 13-31: done with [120,257]
dirPat = '20160927-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step.*v74-OBBmax20-Dmax\(50-100\)-expvCustom257-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 32-50: done with [0,256]
dirPat = '20160927-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step.*v74-OBBmax20-Dmax\(50-100\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 51-69: done with [0,0] ->  expvCustom only
dirPat = '20160927-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step.*v74-OBBmax20-Dmax\(50-100\)-expvCustom0-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 70-88: done with [120,256], SV: [6,6.5]
dirPat = '20160928-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step.*v74-OBBmax20-Dmax\(50-100\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

% 89-107: done with [120,257], SV: [7,9]
dirPat = '20160928-.*-DPMES-Tree1-Tree-v74-L18CT0LE\+';
filPat = 'results-Till3Step.*v74-OBBmax20-Dmax\(50-100\)-expvCustom257-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Step(?<dt>[0-9\-\.]*)v','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.dt),tokens));
matches = [matches; m(I)];

res = TDVPData({matches.name});
%% plot
for ii = 1:size(defPlot,1)
	f = figure(); hold all; ax = gca; f.Name = defPlot{ii,1};
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder',figure);
	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
	legs = strsplit(sprintf('%.3g ',[res(pick).dt]));
	legend(legs{1:end-1});
	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
	[a.LineStyle] = deal(linestyles{1:length(a)});
	ylabel('$\rho_{ii} (t)$');
	set(gca,'color','none');
	grid on
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree1 v76 TreeMPS dt convergence and Dnode - TDVPData								% LabBook 10/11/2016
clear
defPlot(1,:) = {'20161112-DPMES-Tree1-v76-dt0.5-BondConvergence-Dchain20',				[ 1: 6], {'xlim',[0,3300]}};
defPlot(2,:) = {'20161112-DPMES-Tree1-v76-dt0.5-BondConvergence-Dchain40',				[ 7: 9], {'xlim',[0,3300]}};
% defPlot(3,:) = {'20161112-DPMES-Tree1-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-6: Dchain 20
dirPat = '20161111-1701-4.-DPMES-Tree1-Tree-v76TCMde10-L18CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax40-Dmax\(.*-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 7-9: Dchain 40
dirPat = '20161111-1701-4.-DPMES-Tree1-Tree-v76TCMde10-L18CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax40-Dmax\(.*-40\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 13: Tree1 dt 0.1 for reference

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));

	%% plot full dynamics in grid
close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[2,ceil(length(pick)/2)]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		title(ax(kk),sprintf('$D_{node}=%s, \\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
% 	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
% 	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
% 	legend({res(pick).LegLabel});
% 	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
% 	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
% 	[a.LineStyle] = deal(linestyles{1:length(a)});
% 	ylabel('$\rho_{ii} (t)$');
% 	set(gca,'color','none');
% 	grid on
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree2 v76 TreeMPS dt convergence and Dnode - TDVPData								% LabBook 10/11/2016
clear
defPlot(1,:) = {'20161101-DPMES-Tree2-v76-dt0.1-BondConvergence1',						[ 1: 6], {'xlim',[0,3300]}};
defPlot(2,:) = {'20161101-DPMES-Tree2-v76-dt0.5-BondConvergence2',						[ 7:12], {'xlim',[0,3300]}};
% defPlot(3,:) = {'20161101-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-6: dt 0.1
dirPat = '20161031-1427-3.-DPMES-Tree2-Tree-v76TCMde11-L18CT0LE\+';
filPat = 'results-Till5000Step0.1v76-OBBmax60-Dmax\(.*-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 7-12: dt 0.5
dirPat = '20161031-1428-0.-DPMES-Tree2-Tree-v76TCMde11-L18CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax60-Dmax\(.*-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 13: Tree1 dt 0.1 for reference

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot CT convergence
for ii = 1:size(defPlot,1)
	f = figure(); hold all; ax = gca; f.Name = defPlot{ii,1};
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder',figure);
	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
	legend({res(pick).LegLabel});
	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
	[a.LineStyle] = deal(linestyles{1:length(a)});
	ylabel('$\rho_{ii} (t)$');
	set(gca,'color','none');
	grid on
end
% x.plot('calctime-d-sec',figure)

	%% plot full dynamics in grid
%close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[2,3]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		title(ax(kk),sprintf('$D_{node}=%s, \\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
% 	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
% 	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
% 	legend({res(pick).LegLabel});
% 	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
% 	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
% 	[a.LineStyle] = deal(linestyles{1:length(a)});
% 	ylabel('$\rho_{ii} (t)$');
% 	set(gca,'color','none');
% 	grid on
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree2 v76 TreeMPS Broad 20cm vs L - TDVPData								% LabBook 10/11/2016
clear
defPlot(1,:) = {'20161112-DPMES-Tree2-v76-dt0.5-A1-broad-vs-L',							[ 1: 4], {'xlim',[0,1500]}};
defPlot(2,:) = {'20161112-DPMES-Tree2-v76-dt0.5-All-broad-vs-L',						[ 5: 8], {'xlim',[0,1100]}};
% defPlot(3,:) = {'20161112-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-4: broad A1 only
dirPat = '20161111-1551-43-DPMES-Tree2-Tree-v76TCMde9-L.*CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'L(?<L>[0-9\.]*)CT','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.L),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 5-8: broad all
dirPat = '20161111-1551-52-DPMES-Tree2-Tree-v76TCMde9-L.*CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'L(?<L>[0-9\.]*)CT','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.L),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];


res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot full dynamics in grid
close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[2,2]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		title(ax(kk),sprintf('$L=%s, \\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree2 v76 TreeMPS Broad 20cm vs chains - TDVPData								% LabBook 10/04/2017
clear
defPlot(1,:) = {'20170410-DPMES-Tree2-v76-dt0.5-broad-vs-chains',							[ 1: 4], {'xlim',[0,1000]}};
% defPlot(3,:) = {'20161112-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-4: broad different chains
dirPat = '20170406-0944-..-DPMES-Tree2-Tree-v76TCMde9-L100CT0LE\+';
filPat = 'results-Till5000Step0.5v76-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'L(?<L>[0-9\.]*)CT','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.L),tokens));
matches = [matches; m(I)];
legLabels = {'B1','A1','B2','A2'};

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot full dynamics in grid
close all
for ii = 1:1%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[2,2]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		title(ax(kk),sprintf('Broadened %s',legLabels{kk}));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
end
% x.plot('calctime-d-sec',figure)

%% DPMESclust7-1 Star model															% LabBook 09/04/2017
clear
defPlot(1,:) = {'20170409-DPMESclust7-1-Star-v77-dt0.5-vs-D',					[ 1: 4], {'xlim',[0,670]}};
defPlot(2,:) = {'20170409-DPMESclust7-1-Star-v77-dt2-vs-D',						[ 6: 9], {'xlim',[0,670]}};
% defPlot(3,:) = {'20161112-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-5: dt = 0.5
dirPat = '20170408-0329-34-DPMESclust7-1-Tree-v77TCMde10-L100CT0LE\+';
filPat = 'results-Till1000Step0.5v77-OBBmax40-Dmax\(5-.*\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax.5-(?<D>[0-9]*)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.D),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 6-9: dt = 2
dirPat = '20170408-0329-34-DPMESclust7-1-Tree-v77TCMde10-L100CT0LE\+';
filPat = 'results-Till1000Step2v77-OBBmax40-Dmax\(5-.*\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax.5-(?<D>[0-9]*)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.D),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];


res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot full dynamics in grid
close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[2,2]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:min(length(ax),length(pick))
		title(ax(kk),sprintf('$D=%s, \\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
end

%% DPMES-Tree3/4 Comparison															% LabBook 09/04/2017
clear
defPlot(1,:) = {'20170409-DPMES-Tree3-v77-dt0.5',					[ 1: 2], {'xlim',[0,4500]}};
defPlot(2,:) = {'20170409-DPMES-Tree4-v77-dt0.5',					[ 3: 4], {'xlim',[0,4500]}};
% defPlot(3,:) = {'20161112-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1: Tree3 no broadening
dirPat = '20170409-1735-52-DPMES-Tree3-Tree-v77TCMde11-L100CT0LE\+';
filPat = 'results-Till5000Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
matches = [matches; m];
legLabels = [legLabels, {'Tree3'}];

% 2: Tree3 with 20cm broadening
dirPat = '20170409-1736-02-DPMES-Tree3-Tree-v77TCMde11-L100CT0LE\+';
filPat = 'results-Till5000Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
matches = [matches; m];
legLabels = [legLabels, {'Tree3 broad'}];

% 1: Tree4 no broadening
dirPat = '20170409-1735-52-DPMES-Tree4-Tree-v77TCMde11-L100CT0LE\+';
filPat = 'results-Till5000Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
matches = [matches; m];
legLabels = [legLabels, {'Tree4'}];

% 2: Tree4 with 20cm broadening
dirPat = '20170409-1736-02-DPMES-Tree4-Tree-v77TCMde11-L100CT0LE\+';
filPat = 'results-Till5000Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
matches = [matches; m];
legLabels = [legLabels, {'Tree4 broad'}];


res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot full dynamics in grid
close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid',[1,2]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:min(length(ax),length(pick))
		title(ax(kk),sprintf('%s, $\\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
end

%% DPMES-Tree2 Potential Dynamics													% LabBook 14/04/2017
clear
defPlot(1,:) = {'20170414-DPMES-Tree2-v77-dt0.5',					[ 1: 1], {'xlim',[0,450]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1: Tree2 no broadening
dirPat = '20170412-1742-25-DPMES-Tree2-Tree-v77-L18CT0LE\+';
filPat = 'results-Till1000Step0.1v77-OBBmax40-Dmax\(30-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
matches = [matches; m];
legLabels = [legLabels, {'noBroad'}];


res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));

	%% plot
close all
for ii = 1:1%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('heff-pop-diab','-fsev','-ylim',[1.5,4.5]);
	h(1).f.Name = defPlot{ii,1};
	ax = [h.f.Children];
	[ax.XLim] = deal(defPlot{ii,3}{2});
end

%% Heff-swapping: 20170412-1742 (not good)
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170412-1742-25-DPMES-Tree2-Tree-v77-L18CT0LE\+';
filPat = 'results-Till1000Step0.1v77-OBBmax40-Dmax\(30-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
	%% Define swaps
	% first find swaps 1 after the other (only do this once per dataset)
swaps = [  41 2 3;...
		  124 4 5;...
		  302 2 5;...
		  317 2 5;...
		  363 2 5;...
		  757 2 4;...	# soft change
		  924 1 3;...
		  958 4 5;...	# soft change
		 1120 2 5;...	# soft change
		 1473 1 4;...	# soft change
		 2248 2 5;...	# very soft		(could be removed)
		 2280 2 5;...	# very soft
		 2305 2 5;...	# very soft
		  ]		

% h = res.plot('heff-pop-diab-swap','-state',swaps,'-ylim',[1.8,4.2],'-rowheight',12,'-rowwidth',49,'-patchthickness',0.2);
%% Kin-Heff-swapping: 20170417-1449 Tree2
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170417-1449-08-DPMES-Tree2-Tree-v77-L18CT0LE\+';
filPat = 'results-Till1500Step0.1v77-OBBmax40-Dmax\(30-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   7 4 5;...
		   41 2 3;...
		  124 4 5;...
		  302 2	4;...
		  317 2 4;...
		  363 2 4;...
		  757 2 5;...
		  924 1 3;...
		  958 4 5;...
		 1120 2 4;...
		 1473 1 5;...
		 2247 2 4;...
		 2279 2 4;...
		 2307 2 4;...
		];
h=res(1).plot('heff-pop-diab-swap','-state',swaps,'-rowheight',12,'-rowwidth',48,'-ylim',[1.6,3.8],'fsev');
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170417-1449-heff-swap';
delete(h.ax(6));
xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% Kin-Heff-swapping: 20170417-1451 Tree2 (B1 broad)
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170417-1451-55-DPMES-Tree2-Tree-v77-L100CT0LE\+';
filPat = 'results-Till1500Step0.1v77-OBBmax40-Dmax\(30-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   7 4 5;...
		   38 2 3;...
		  130 4 5;... 
		  304 2 4;...
		  321 2 4;...
		  365 2 4;
		  810 2 5;
		 1057 1 3;
		 1296 4 5;
		 2090 1 5;
		];
h=res(1).plot('heff-pop-diab-swap','-state',swaps,'-rowheight',12,'-rowwidth',48,'-ylim',[1.6,3.8]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170417-1451-heff-swap';
delete(h.ax(6));
xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% Kin-Heff-swapping: 20170417-1536 Tree4
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170417-1536-29-DPMES-Tree4-Tree-v77-L100CT0LE\+';
filPat = 'results-Till1500Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   5 2 3;
		   58 4 5;
		  161 3 5;
		  326 1 5;
		  437 3 4;
		  457 3 4;
		  464 3 4;
		  608 3 4;
		  627 3 4;
		  676 3 4;
		  681 3 4;
		];
h=res(1).plot('heff-pop-diab-swap','-state',swaps,'-rowheight',12,'-rowwidth',48,'-ylim',[1.6,3.8]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170417-1536-heff-swap';
delete(h.ax(6));
xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% PES-Heff-swapping: 20170419-1218 Tree2
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170419-1218-25-DPMES-Tree2-Tree-v77-L18CT0LE\+';
filPat = 'results-Till1500Step0.2v77-OBBmax40-Dmax\(40-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [  22	2	3;
		   64	4	5;
		  152	2	5;
		  159	2	5;
		  183	2	5;
		  378	2	4;
		  463	1	3;
		  480	4	5;
		  563	2	5;
		  736	1	4;
		];
h=res(1).plot('heff-pop-diab-swap','-fsev','-state',swaps,'-rowheight',12,'-rowwidth',48,'-ylim',[1.5,3.4],'-xlim',[0,1000]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170419-1218-PES-swap';
delete(h.ax(6));
% xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% PES-Heff-swapping: 20170428-1244 Tree4 Sym Clean
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170428-1244-31-DPMES-Tree4-Tree-v77-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   5	2	3;
		   58	4	5;
		  162	2	5;
		  173	2	5;
		  324	1	3;
		  434	4	5;
		];
h=res(1).plot('heff-pop-diab-swap','-state',swaps,'-rowheight',12,'-rowwidth',48,'-ylim',[1.2,3.4],'-xlim',[0,900]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170428-1244-PES-swap';
delete(h.ax(6));
% xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% PES-Heff-swapping: 20170428-1348 Tree4 75 Modes
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170428-1348-57-DPMES-Tree4-Tree-v77-75Modes-L20CT0LE\+';
filPat = 'results-Till1500Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom256-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   5	2	3;
		  118	4	5;
		  170	2	5;
		  320	1	3;
		  331	2	5;
		  361	4	5;
		];
res(2) = res(1).getData('heff-swap','state',swaps);
h=res(2).plot('heff-pop-diab','-rowheight',12,'-rowwidth',48,'-ylim',[1.0,3.9],'-xlim',[0,1100]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170428-1348-PES-swap';
delete(h.ax(6));
% xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');
%% PES-Heff-swapping: 20170428-2117 Tree4 full
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];

% 1
dirPat = '20170428-2117-08-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(50-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
% tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
% [y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m];
% legLabels = [y];


res = TDVPData({matches.name});
%% swap
swaps = [   5	2	3;
		   58	4	5;
		  162	2	5;
		  173	2	5;
		  324	1	3;
		  434	4	5;
		];
res(2) = res(1).getData('heff-swap','state',swaps);
h=res(2).plot('heff-pop-diab','-rowheight',12,'-rowwidth',48,'-ylim',[1.2,3.4],'-xlim',[0,900]);
h.f.Position([1,2]) = [0,1];
h.f.Name = '20170428-1244-PES-swap';
delete(h.ax(6));
% xlabel(h.ax(2),'$t/fs$'); ylabel(h.ax(2),'$E/eV$');

%% DPMES-Tree4 v78 TreeMPS Dnode convergence  - TDVPData								% LabBook 30/05/2017
clear
defPlot(1,:) = {'20170428-DPMES-Tree4-v78-dt0.5-DNodeConvergence',						[ 1: 5], {'xlim',[0,1200]}};
% defPlot(3,:) = {'20161101-DPMES-Tree2-v76-dt0.5-BondConvergence3',						[ 1:12], {'xlim',[0,200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-5: dt 0.5
dirPat = '20170428-2116-06-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(.*-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
legLabels = [legLabels, y];

res = TDVPData({matches.name});
res = res.setLegLabel(mat2cell(legLabels,1,ones(1,length(legLabels))));
	%% plot single state convergence
state = 1;
for ii = 1:size(defPlot,1)
	f = figure(); hold all; ax = gca; f.Name = defPlot{ii,1};
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder',figure);
	a = cellfun(@(x) copyobj(x(state),ax),{h.pl},'UniformOutput',false); a=[a{:}];
	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
	legend({res(pick).LegLabel});
	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
	[a.LineStyle] = deal(linestyles{1:length(a)});
	ylabel('$\rho_{ii} (t)$');
	set(gca,'color','none');
	grid on
end
% x.plot('calctime-d-sec',figure)

	%% plot full dynamics in grid
%close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid');
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		title(ax(kk),sprintf('$D_{node}=%s, \\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
% 	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
% 	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
% 	legend({res(pick).LegLabel});
% 	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
% 	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
% 	[a.LineStyle] = deal(linestyles{1:length(a)});
% 	ylabel('$\rho_{ii} (t)$');
% 	set(gca,'color','none');
% 	grid on
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree4 v78 TreeMPS Dnode[80,100,120] convergence  + PES- TDVPData				% LabBook 30/05/2017
clear
defPlot(1,:) = {'20170504-DPMES-Tree4-v78-dt0.5-DConvergence',						[ 1: 6], {'xlim',[0,1200]}};
defPlot(2,:) = {'20170504-DPMES-Tree4-v78-dt0.2-DConvergence',						[ 7:13], {'xlim',[0,1200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-6: dt 0.5
dirPat = '20170504-1825-10-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(.*-.*\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);	
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 7:13 dt 0.2
dirPat = '20170504-1825-10-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.2v77-OBBmax40-Dmax\(.*-.*\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
matches = [matches; m(I)];
legLabels = [legLabels, y];

res = TDVPData({matches.name});
res = res.setLegLabel(legLabels);
	%% plot single state convergence
state = 1;
for ii = 1:size(defPlot,1)
	f = figure(); hold all; ax = gca; f.Name = defPlot{ii,1};
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder',figure);
	a = cellfun(@(x) copyobj(x(state),ax),{h.pl},'UniformOutput',false); a=[a{:}];
	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
	legend({res(pick).LegLabel});
	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
	[a.LineStyle] = deal(linestyles{1:length(a)});
	ylabel('$\rho_{ii} (t)$');
	set(gca,'color','none');
	grid on
end
% x.plot('calctime-d-sec',figure)

	%% plot full dynamics in grid
%close all
for ii = 1:2%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid');
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		t=title(ax(kk),sprintf('%s, $\\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		t.Interpreter = 'latex';
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
% 	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
% 	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
% 	legend({res(pick).LegLabel});
% 	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
% 	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
% 	[a.LineStyle] = deal(linestyles{1:length(a)});
% 	ylabel('$\rho_{ii} (t)$');
% 	set(gca,'color','none');
% 	grid on
end
% x.plot('calctime-d-sec',figure)

%% DPMES-Tree4 v78 TreeMPS Dnode[10:120] convergence - TDVPData							% LabBook 30/05/2017
clear
defPlot(1,:) = {'20170530-DPMES-Tree4-v78-DConvergence-overview',						[1:4,12,6:9,11], {'xlim',[0,1200]}};

load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-5: dt 0.5
dirPat = '20170428-2116-06-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(.*-20\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$(%s,%s)$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
legLabels = [legLabels, y];

% 6-11: dt 0.5
dirPat = '20170504-1825-10-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.5v77-OBBmax40-Dmax\(.*-.*\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);	
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$(%s,%s)$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
matches = [matches; m(I)];
legLabels = [legLabels, y];

% 12: dt 0.2 (50,30)
dirPat = '20170504-1825-10-DPMES-Tree4-Tree-v77TCMde9-L100CT0LE\+';
filPat = 'results-Till1700Step0.2v77-OBBmax40-Dmax\(50-30\)-expvCustom270-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
matches = [matches; m(I)];
legLabels = [legLabels, y];


res = TDVPData({matches.name});
res = res.setLegLabel(legLabels);

	%% plot full dynamics in grid (only testing)
%close all
for ii = 1:1%size(defPlot,1)
	
	pick = defPlot{ii,2};
	h = res(pick).plot('rhoii','-fsev','-resetColorOrder','-grid');
	h(1).f.Name = defPlot{ii,1};
	ax = [h.ax];
	[ax.XLim] = deal(defPlot{ii,3}{2});
	for kk = 1:length(ax)
		t=title(ax(kk),sprintf('%s, $\\frac{t_{CPU}}{sweep} = %.2g min$',res(pick(kk)).LegLabel, 60*mean(diff(nonzeros(res(pick(kk)).para.tdvp.calcTime)))));
		t.Interpreter = 'latex';
		if kk ~= length(ax)
			a = copyobj(h(end).pl,ax(kk));		% copy D=100 into other gridcells
			set(a,'Color',[1,1,1]*0.6);
		end
	end
% 	a = cellfun(@(x) copyobj(x(4),ax),{h.pl},'UniformOutput',false); a=[a{:}];
% 	ax.XLabel = h(1).ax(1).XLabel; close(h(1).f);
% 	legend({res(pick).LegLabel});
% 	set(a,{'Color'},mat2cell(ax.ColorOrder(mod(0:length(a)-1,7)+1,:),ones(1,length(a)),3));
% 	linestyles = reshape(repmat({'-','--',':'},7,1),[],1);
% 	[a.LineStyle] = deal(linestyles{1:length(a)});
% 	ylabel('$\rho_{ii} (t)$');
% 	set(gca,'color','none');
% 	grid on
end
% x.plot('calctime-d-sec',figure)

	%% plot TT error against time/sweep
%close all
clear h
for ii = 1:size(defPlot,1)
	pick = defPlot{ii,2};
	h.f = figure(ii); clf; h.ax = gca; h.f.Name = defPlot{ii,1}; hold all;
	% extract the TT populations
	ttPop = cell2mat(arrayfun(@(x) TDVPData.subsample(real(x.getData('rhoii','state',1)),x.dt,1),res(pick)','UniformOutput',false));
	% estimate mean error
	ttErr = bsxfun(@minus,-ttPop,-ttPop(:,end));
	ttErr = bsxfun(@rdivide,ttErr,ttPop(:,end)); ttErr(isnan(ttErr)) = 0;		% makes error relative
	ttMeanErr = (mean(ttErr,1))*100;
	
% 	ttErr = mean(abs(bsxfun(@minus,-ttPop,-ttPop(:,end))),1);
% 	ttMeanErr = ttErr./mean(abs(diff(ttPop(:,end))));
	% get Runtimes
	tSweep = arrayfun(@(x) 60*mean(diff(nonzeros(x.para.tdvp.calcTime))),res(pick));
	h.pl = scatter(tSweep, ttMeanErr,'filled');
	xlabel('$t/Sweep/min$');
	ylabel('Mean rel. error in $\rho_{\mathrm{TT}}$');
	box on;
	% shift text labels
	offX = -ones(1,length(pick))*0.9; offY = ones(1,length(pick))*0.04;
	offY(1) = -0.05;
	offX(4) = -2.1;offY(4) = 0;
	offX(6) = -1.5;offY(6) = 0.12;
	offX(7) = 0;offY(7) = 0.12;
	offX(8) = 0;
	offX(9) = -0.5;
	for kk = 1:length(pick)
		t = text(tSweep(kk)+offX(kk),double(ttMeanErr(kk))+offY(kk),strrep(res(pick(kk)).LegLabel,', ','\\\\'));
	end
	set(gca,'color','none');
% 	xlim([0,11.9]);
% 	ylim([-0.03,0.56]);
	plot(h.ax.XLim,[0,0],'k');
end
resizePlot(h.f)

%% DPMES-Tree4 v79 TreeMPS ProperPot - TDVPData											% LabBook 25/07/2017
clear
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-2: dt 0.5
dirPat = '20170717-.*-DPMES-Tree4-Tree-v79-ProperPot-L70CT0LE\+';
filPat = 'results-Till1700Step0.5v79-OBBmax40-Dmax.*-expvCustom240-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$(%s,%s)$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
legLabels = [legLabels, y];

res = TDVPData({matches.name});
res = res.setLegLabel(legLabels);
	%% Do Plots
	res.plot('rhoii','-fsev');
	%% PES, TES in grid
	defPlot(1,:) = {'20170725-DPMES-Tree4-v79-ProperPot-ACzoom',				[1],		{'xlim',[0, 200],'ylimTop',[2.04, 2.119],'ylimBot',[1.54, 2.08]}};
	defPlot(2,:) = {'20170725-DPMES-Tree4-v79-ProperPot',						[1],		{'xlim',[0,1130],'ylimTop',[1.6 , 3.8  ],'ylimBot',[1.5, 3.5]}};
	
	ii = 2;
	x = res(defPlot{ii,2});
	h = TDVPData.plotGrid(2,2,'rowheight',5,'rowwidth',16); resizePlot(h.f,1.6,2);	% side-by-side
h.f.Name = defPlot{ii,1};
xlimAll = defPlot{ii,3}{2}; ylimTop = defPlot{ii,3}{4}; ylimBot = defPlot{ii,3}{6}; thickness = 0.1; dotsize = 50;
htemp  = x.setHeffTo('tes').plot('heff-full-pop-diab-v2','-fsev'   ,'-xlim',xlimAll,'-ylim',ylimTop,'-patchthickness',thickness,h.ax(1)); htemp.pl(end).SizeData = dotsize;
htemp2 = x.plot('pes-from-tes-pop-diab-v2','-fsev','-xlim',xlimAll,'-ylim',ylimBot,'-patchthickness',thickness,h.ax(2)); htemp2.pl(end).SizeData = dotsize;
htemp2 = x.plot('tes-from-pes-pop-diab-v2','-fsev','-xlim',xlimAll,'-ylim',ylimTop,'-patchthickness',thickness,h.ax(3)); htemp2.pl(end).SizeData = dotsize;
htemp2 = x.setHeffTo('pes').plot('heff-full-pop-diab-v2','-fsev'   ,'-xlim',xlimAll,'-ylim',ylimBot,'-patchthickness',thickness,h.ax(4)); htemp2.pl(end).SizeData = dotsize;

ylabel(h.ax(3),'');
ylabel(h.ax(4),'');

legend(htemp.pl(1:5),{'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'}, 'location','North');
%
legText = {'a','b','c','d'};
legPosX = ones(1,5)*0.02;
legPosY = ones(1,4)*0.92;
for ii = 1:length(legText)
	t=text(h.ax(ii),legPosX(ii),legPosY(ii),sprintf('$\\mathbf{%s}$',legText{ii}),'sc','FontSize',8.8,'FontWeight','bold');
end

%% DPMES-Tree6 v79 TreeMPS divB1 - TDVPData												% LabBook 25/07/2017
clear
load('TDVPLib.mat');
%
TDVPfolds = TDVPfolds(arrayfun(@(x) ~isempty(strfind(x.name,'DPMES')),TDVPfolds));
matches = [];
legLabels = [];

% 1-2: 1: bad, 2: good calculation
dirPat = '20170720-1328-3.-DPMES-Tree6-Tree-v79-divB1-L70CT0LE\+';
filPat = 'results-Till.*Step.*v79-OBBmax40-Dmax.*-expvCustom240-1core-small.mat';
m = TDVPData.getMatches(TDVPfolds,dirPat,filPat);
tokens = regexp({m.name},'Dmax\((?<DNode>[0-9\.]*)-(?<DChain>[0-9\.]*)\)','names');			% start sorting
[y,I] = sort(cellfun(@(x) str2double(x.DNode),tokens));
matches = [matches; m(I)];
% y = arrayfun(@(ii) sprintf( '$D_{Node}=%s, D_{Chain}=%s$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
y = arrayfun(@(ii) sprintf( '$(%s,%s)$',tokens{ii}.DNode,tokens{ii}.DChain), I, 'UniformOutput', false);
legLabels = [legLabels, y];

res = TDVPData({matches.name});
res = res.setLegLabel(legLabels);
	%% Do Plots
	res.plot('rhoii','-fsev',gcf,'-resetColorOrder');
	%% PES, TES in grid
	defPlot(1,:) = {'20170725-DPMES-Tree6-v79-divB1-ACzoom',				[2],		{'xlim',[0, 200],'ylimTop',[2.04, 2.119],'ylimBot',[1.54, 2.08]}};
	defPlot(2,:) = {'20170725-DPMES-Tree6-v79-divB1',						[2],		{'xlim',[0, 600],'ylimTop',[1.6 , 3.8  ],'ylimBot',[1.5, 3.5]}};
	
	ii = 1;
	x = res(defPlot{ii,2});
	h = TDVPData.plotGrid(2,2,'rowheight',5,'rowwidth',16); resizePlot(h.f,1.6,2);	% side-by-side
h.f.Name = defPlot{ii,1};
xlimAll = defPlot{ii,3}{2}; ylimTop = defPlot{ii,3}{4}; ylimBot = defPlot{ii,3}{6}; thickness = 0.1; dotsize = 50;
htemp  = x.setHeffTo('tes').plot('heff-full-pop-diab-v2','-fsev'   ,'-xlim',xlimAll,'-ylim',ylimTop,'-patchthickness',thickness,h.ax(1)); htemp.pl(end).SizeData = dotsize;
htemp2 = x.plot('pes-from-tes-pop-diab-v2','-fsev','-xlim',xlimAll,'-ylim',ylimBot,'-patchthickness',thickness,h.ax(2)); htemp2.pl(end).SizeData = dotsize;
htemp2 = x.plot('tes-from-pes-pop-diab-v2','-fsev','-xlim',xlimAll,'-ylim',ylimTop,'-patchthickness',thickness,h.ax(3)); htemp2.pl(end).SizeData = dotsize;
htemp2 = x.setHeffTo('pes').plot('heff-full-pop-diab-v2','-fsev'   ,'-xlim',xlimAll,'-ylim',ylimBot,'-patchthickness',thickness,h.ax(4)); htemp2.pl(end).SizeData = dotsize;

ylabel(h.ax(3),'');
ylabel(h.ax(4),'');

legend(htemp.pl(1:5),{'TT','LE$^+$','LE$^-$','CT$^+$','CT$^-$'}, 'location','North');
%
legText = {'a','b','c','d'};
legPosX = ones(1,5)*0.02;
legPosY = ones(1,4)*0.92;
for ii = 1:length(legText)
	t=text(h.ax(ii),legPosX(ii),legPosY(ii),sprintf('$\\mathbf{%s}$',legText{ii}),'sc','FontSize',8.8,'FontWeight','bold');
end

%% TDVP SBM multi (1): Plot Visibility / Coherence
fignum = 3; figure(fignum); clf; hold all;
pick = [1:length(res)];			% plot all
% pick = [1:5];						% plot selective
xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
plot([1,xmax],[0,0],'black');
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
% ph{8}.LineStyle = ':';				%temp
set(gca,'ylim',[-1,1]);
% set(gca,'xlim',[0,325]);
% set(gca,'xlim',[1,xmax]);set(gca,'xscale','log');
xlabel('t');
ylabel('$\left<\sigma_z\right>$');
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum);

%% TDVP SBM multi (2): Plot Visibility / Coherence		Predefined
for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
% 	cellfun(@(x) set(x,'Color',col(1,:)), ph(1:5));
% 	cellfun(@(x) set(x,'Color',col(2,:)), ph(6:10));
% 	cellfun(@(x) set(x,'Color',col(3,:)), ph(11:15));
% 	ph{mod(1:15,3)}. = ':';				%temp
	set(gca,'ylim',[-1,1]);
	% set(gca,'xlim',[1,1e8]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
% 	leg = legend([ph{:}],strsplit(sprintf('%g,',res{pick,3}),','),'location','best');
% 	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM multi (3): Plot VMPS GS <n>
fignum = 2; figure(fignum); clf; hold all;
pick = [1:length(res)];			% plot all
% pick = [8:11];
% pick = [23,35];
% pick = [6,12,21,22,23,24,25,26];						% plot selective
% pick = [8,15,16,11,12,13,14];
ph = cellfun(@(x) plot(real(x.results.nx)), res(pick,1), 'UniformOutput', false);
% ph{2}.LineStyle = ':';				%temp
% ph{4}.LineStyle = ':';				%temp
set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<n_{k,VMPS}>$')
set(gca,'yscale','log')
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum)

%% TDVP SBM multi (4): Plot GS Energy convergence
fignum = 3; figure(fignum); clf; hold all;
% pick = [1:length(res)];			% plot all
% pick = [8,15,16,11,12,13,14];						% plot selective
% pick = [8:11];
% pick = [23,35];
ph = cellfun(@(x) plot(cell2mat(x.results.EvaluesLog)-min(cell2mat(x.results.EvaluesLog))), res(pick,1), 'UniformOutput', false);
% disp(sprintf('%.15e\n',cell2mat(cellfun(@(x) x.results.E, res(pick,1), 'UniformOutput', false))))
set(gca,'YScale','log');
% try
% title(sprintf('$E_0 = %.10g, \\Lambda =  %.2g, z =  %.2g$',results.E, para.Lambda, para.z));catch end
xlabel('Site$\cdot$Loop');
ylabel('$E-E_0$');
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum)
yLim = get(gca,'YLim');
for i = 1:para.loop
%     line([para.L*i para.L*i],yLim,'LineWidth',1,'Color','black');
end
if wantSave
    export_fig(sprintf('%s%s-MLSBM-Econvergence-Lambda%.2gz%.2gp16',saveto,para.filename(1:13),para.Lambda,para.z),'-transparent','-png','-painters')
end

%% TDVP SBM multi (5): plot slider <n> STAR
mode = 1;		% 0: lin, 1: log
fignum = 7;
f = figure(fignum); clf; x = res(1,:);
f.Position(3) = 840; pos = f.Position; ax = gca;
pl = surf(log10(x{1}.tresults.star.n(:,:,1)));
% setting auto-refresh
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(gcf));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$s = %g, \\alpha = %g$',pl.UserData{1}.para.chain{1}.s, pl.UserData{1}.para.chain{1}.alpha)));

pl.UserData=x;
cb = colorbar;cb.Title.Interpreter = 'latex';
prefix = 'pl.UserData{1}.tresults.star';
if mode
	zlabel('$log_{10}\left<n_k\right>$');
	pl.ZDataSource = sprintf('log10(%s.n(1:find(%s.t,1,last),:,1));',prefix,prefix);
% 	pl.ZDataSource = sprintf('log10(%s.n)-log10(ones(size(%s.n,1),1)*%s.n(1,:));',prefix,prefix,prefix);
else
	zlabel('$\left<n_k\right>$');
	pl.ZDataSource = sprintf('%s.n(1:find(%s.t,1,last),:,1);',prefix,prefix);
end
pl.XDataSource = [prefix,'.omega'];
pl.YDataSource = sprintf('%s.t(1:find(%s.t,1,last))',prefix,prefix);last = 'last';refreshdata;
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
cb.Title.String = ax.ZLabel.String;
% set(gca,'zscale','log');
set(gca,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;rotate3d on;axis tight;
if mode
% 	ax.ZLim = [0,10];
% 	ax.CLim = ax.ZLim;
end
% set(gca,'xlimmode','manual','xlim',[0,x{1}.tresults.star.omega(end)]);

% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(1);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'userdata',res(round(source.Value),:)));

%% TDVP SBM multi (6): plot slider <n> CHAIN - MOD MC
mode = 1;		% 0: lin, 1: log
fignum = 9; f = figure(fignum); f = figure(fignum);
del = findobj(f,'Type','hgjavacomponent'); del.delete;	% get rid of old sliders
f.Position(3) = 1.8*f.Position(4); pos = f.Position; hold all;
ax = gca; ax.UserData = 1; f.UserData = 1;
hPl = handle(ax); hProp = findprop(hPl,'UserData');
hPl2 = handle(f); hProp2 = findprop(hPl2,'UserData');
cb = colorbar;cb.Title.Interpreter = 'latex';
if mode
	zlabel('$log_{10}\left<n_k\right>$');
	pl = surf(log10(abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, 1))));
	hPl.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'zdata', log10(abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, f.UserData)))));
	hPl2.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'zdata', log10(abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, f.UserData)))));
else
	zlabel('$\left<n_k\right>$');
	pl = surf(abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, 1)));
	hPl.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'zdata', abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, f.UserData))));
	hPl2.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'zdata', abs(res{ax.UserData,1}.tresults.n(1:res{ax.UserData,1}.tresults.lastIdx, :, f.UserData))));
end
hPl.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'xdata', 1:res{ax.UserData,1}.para.L));
hPl.addlistener(hProp, 'PostSet', @(src, event) set(pl, 'ydata', res{ax.UserData,1}.tresults.t(1:res{ax.UserData,1}.tresults.lastIdx)*0.658));
cb.Title.String = ax.ZLabel.String;
xlabel('Site $k$');
% ylabel('Time $\omega_c t$');
ylabel('Time $t/fs$');
% set(gca,'zscale','log');
set(ax,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;
rotate3d on;axis tight;
if mode
	ax.ZLim = [-2,1];
	ax.CLim = ax.ZLim;
end
% slider definition for datasets:
sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(3);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax,'userdata',round(source.Value)));

% slider definition for chain number
sld2 = javax.swing.JScrollBar(0,1,1,1,res{1,1}.para.nChains+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld2, [pos(3)*0.1,5,70,15], gcf);
sld2.setUnitIncrement(1); sld2.setBlockIncrement(1);
hsld2 = handle(sld2,'CallbackProperties');
set(hsld2,'AdjustmentValueChangedCallback',@(source,callbackdata) set(f,'userdata',round(source.Value)));


% hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$s = %g, \\alpha = %g$',pl.UserData{1}.para.chain{1}.s, pl.UserData{1}.para.chain{1}.alpha)));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('%s',res{ax.UserData,1}.para.folder)));
%% TDVP SBM multi (7): plot slider POLARON
mode = 0;		% 0: lin, 1: log
fignum = 9;
f = figure(fignum); clf; x = res(1,:);
f.Position(3) = 840; pos = f.Position; ax = gca;
pl = surf(x{1}.tresults.star.x(:,:,1,1));
% setting auto-refresh
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(gcf));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$s = %g, \\alpha = %g$',pl.UserData{1}.para.chain{1}.s, pl.UserData{1}.para.chain{1}.alpha)));

pl.UserData=x;
cb = colorbar;cb.Title.Interpreter = 'latex';
if mode
	zlabel('$log_{10}\left<x_k\right>$');
	pl.ZDataSource = 'log10(pl.UserData{1}.tresults.star.x(:,:,1,1));';
else
	zlabel('$\left<x_k\right>$');
	pl.ZDataSource = 'pl.UserData{1}.tresults.star.x(:,:,1,1);';
end
pl.XDataSource = 'pl.UserData{1}.tresults.star.omega';
pl.YDataSource = 'pl.UserData{1}.tresults.star.t';refreshdata;
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
cb.Title.String = ax.ZLabel.String;
% set(gca,'zscale','log');
set(gca,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;rotate3d on;axis tight;
if mode
% 	ax.ZLim = [-30,0];
% 	ax.CLim = ax.ZLim;
end
% set(gca,'xlimmode','manual','xlim',[0,x{1}.tresults.star.omega(end)]);

% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(1);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'userdata',res(round(source.Value),:)));

%% TDVP SBM multi (8): plot slider CHAIN CURRENT		% most modern now.
mode = 0;		% 0: lin, 1: log
fignum = 10;
f = figure(fignum); clf; x = res(1,:);
f.Position(3) = 840; pos = f.Position; ax = gca;
pl = surf(x{1}.tresults.j);
% setting auto-refresh
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(pl));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$s = %g, \\alpha = %g$',pl.UserData{1}.para.chain{1}.s, pl.UserData{1}.para.chain{1}.alpha)));

pl.UserData=x;
cb = colorbar;cb.Title.Interpreter = 'latex';
p = 'pl.UserData{1}.tresults';							% prefix to shorten expressions
if mode
	zlabel('$log_{10}\left<j_k\right>$');
	pl.ZDataSource = sprintf('log10(abs(%s.j(1:%s.lastIdx,:)));',p,p);
else
	zlabel('$\left<j_k\right>$');
	pl.ZDataSource = sprintf('%s.j(1:%s.lastIdx,:);',p,p);
end
pl.XDataSource = '1:pl.UserData{1}.para.L-1';
pl.YDataSource = sprintf('%s.t(1:%s.lastIdx)',p,p); refreshdata;
xlabel('Site k');
ylabel('Time $\omega_c t$');
cb.Title.String = ax.ZLabel.String;
% set(gca,'zscale','log');
set(gca,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;rotate3d on;axis tight;
if mode
	ax.ZLim = [-30,0];
	ax.CLim = ax.ZLim;
end
% set(gca,'xlimmode','manual','xlim',[0,x{1}.tresults.star.omega(end)]);

% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(1);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'userdata',res(round(source.Value),:)));
%%
[f,ax,pl] = slider2DPlot(res, 10, 'tresults.j', 'para.L-1', 'tresults.t')
%% TEST SBM_genpara:
para1.chain.mapping = 'OrthogonalPolynomials';	para1.chain.spectralDensity = 'Leggett_Hard'; para1.chain.discrMethod = 'Analytic';
para1.chain.discretization = 'None';			para1.chain.method = 'Analytic';
para1.Lambda = 2; para1.z = 1; para1.L = 100; para1.s = 1; para1.alpha = 0.1;
para1 = SBM_genpara(para1);
para2.chain.mapping = 'OrthogonalPolynomials';	para2.chain.spectralDensity = 'Leggett_Hard'; para2.chain.discrMethod = 'Numeric';
para2.chain.discretization = 'LogZ';			para2.chain.method = 'Stieltjes';
para2.Lambda = 2; para2.z = 1; para2.L = 100; para2.s = 1; para2.alpha = 0.1;
para2 = SBM_genpara(para2);
%%
% clear;
para1.chain.mapping = 'OrthogonalPolynomials';	para1.chain.spectralDensity = 'Leggett_Hard'; para1.chain.discrMethod = 'None';
para1.chain.discretization = 'LogZ';			para1.chain.method = 'Analytic';
para1.Lambda = 1.1; para1.z = 1; para1.L = 50; para1.s = 1; para1.alpha = 0.1;
t1 = cputime; para1 = SBM_genpara(para1); para1.time = cputime-t1;
para2.chain.mapping = 'OrthogonalPolynomials';	para2.chain.spectralDensity = 'Leggett_Hard'; para2.chain.discrMethod = 'Analytic';
para2.chain.discretization = 'LogZ';			para2.chain.method = 'Stieltjes';
para2.Lambda = 1.001; para2.z = 1; para2.L = 50; para2.s = 1; para2.alpha = 0.1;
t1 = cputime; para2 = SBM_genpara(para2); para2.time = cputime-t1;
%%
hold all;
% plot(para1.t);plot(para1.epsilon);
% plot(para2.t);plot(para2.epsilon);
plot(para1.t-para2.t);plot(para1.epsilon-para2.epsilon);
% set(gca,'yscale','log');
legend('para1.t','para1.\epsilon','para2.t',num2str(para2.Lambda));
% legend('t Lanzcos','\epsilon','t Stieltjes','\epsilon');
formatPlot(1)

%% TDVP time extrapolation
n = 1:para.timeslice; clf;
% n=1:648;
n = n(para.tdvp.calcTime(n,1)>0);
plot(n,para.tdvp.calcTime(n,1)); hold on;
s = csaps(n,para.tdvp.calcTime(n,1));
sn = fnxtr(s);
fnplt(sn,[0,60],0.5,'r')
fnval(sn,60)

%% TDVP compare with Orth2010: 13a -- not really good. Extraction was not the best
% http://arohatgi.info/WebPlotDigitizer/app/?
orthData = cell(5,1);
orthData{1} = load('./../../Orth2010-13a-0.01.csv');
orthData{2} = load('./../../Orth2010-13a-0.05.csv');
orthData{3} = load('./../../Orth2010-13a-0.1.csv');
orthData{4} = load('./../../Orth2010-13a-0.15.csv');
orthData{5} = load('./../../Orth2010-13a-0.2.csv');

myData{1} = res{15,1};
myData{2} = res{16,1};
myData{3} = res{17,1};
myData{4} = res{18,1};
myData{5} = res{19,1};

fignum = 1; clf; hold all;
cellfun(@(x) plot(x(:,1),x(:,2)), orthData, 'UniformOutput', false)
cellfun(@(x) plot(x.para.tdvp.t(:),x.tresults.spin.sz(:)), myData, 'UniformOutput', false)
set(gca,'xlim',[0,325]);
formatPlot(1)

%% Format <n> plot for Poster:
ax = gca; f=gcf;
% set(gca, 'LooseInset', [0,0,0,0]); % not with colorbar!!
ax.FontSize = 14;
ax.TickDir = 'out';
% ax.Units = 'centimeters';
% for figure 840x420:
% oldPos = ax.Position;
% ax.Position = [oldPos(1),oldPos(2), 15, 9.8];
f.Units = 'centimeters';
oldPos = f.Position;
f.Position = [oldPos(1),oldPos(2), 19.35, 12.28];	% in cm also produces 15x9.8cm axes
% f.Units = 'pixels';
% f.Position = [oldPos(1),oldPos(2), 840, oldPos(4)];
% export_fig('test','-transparent','-png','-m4')
%% Save Surf plots for Poster:
figHandles = get(0,'Children');
for i = 1:size(figHandles)
	f=figure(figHandles(i).Number);ax = gca;
	% set(gca, 'LooseInset', [0,0,0,0]); % not with colorbar!!
	ax.FontSize = 14;
	ax.TickDir = 'out';
	% ax.Units = 'centimeters';
	% for figure 840x420:
	% oldPos = ax.Position;
	% ax.Position = [oldPos(1),oldPos(2), 15, 9.8];
	f.Units = 'centimeters';
	oldPos = f.Position;
% 	f.Position = [oldPos(1),oldPos(2), 19.35, 12.28];	% in cm also produces 15x9.8cm axes
	f.Position = [oldPos(1),oldPos(2), 19.35,  6.45];	% in cm half height!
	% f.Units = 'pixels';
	% f.Position = [oldPos(1),oldPos(2), 840, oldPos(4)];
	export_fig(sprintf('img/test%d',f.Number),'-transparent','-png','-m4');
% 	close(f);
end
% close(1);
%% inset for poster:
if exist('axIn','var')
	axIn.delete;
end
if exist('hInset','var')
	hInset.delete;
end
ax=gca; oldPos = ax.Position;
axIn = axes('Position',[oldPos(1)+0.32, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);		% use this if main with colorbar
% axIn = axes('Position',[oldPos(1)+0.36, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);	% use this if main without colorbar
axIn.TickDir = 'out';
axIn.FontSize=14; axIn.XColor=[0.7,0.7,0.7];axIn.YColor=axIn.XColor;

hold all;
% surface
% surf(tresults.star.omega,tresults.star.t, ax.Children.ZData);
% xlabel('Mode $\omega_k / \omega_c$');
% ylabel('Time $\omega_c t$');
% [M,I] = max(ax.Children.ZData(:,:),[],2);
% plot3(tresults.star.omega(I(2:end)),tresults.star.t(2:end),ones(1,length(tresults.star.t)-1).*200,':r')
% plot3([0.1,0.1],[0,500],[200,200],':black');
% axIn.YLim = [0,500];

% spectra
% tslices = [2,11,31,51, 164];
% plot(tresults.star.omega, ax.Children.ZData(tslices,:));
% legend(strsplit(sprintf('%g ',tresults.star.t(tslices))));
% xlabel('$\omega_k / \omega_c$');
% ylabel(ax.ZLabel.String);
% set(gca,'View',[0 90]);
% shading interp
% axis tight
% axIn.XLim = [0,0.2];

% FT inset via copy axes!
fMaster = 8;
fInset = 11;
axIn.delete;
figure(fInset); ax_Inset = gca;
figure(fMaster);
hInset = copyobj(ax_Inset,fMaster);
hInset.Position = [oldPos(1)+0.32, oldPos(2)+0.36, oldPos(3)*0.5, oldPos(4)*0.5];
hInset.XLim = [0,0.8]; hInset.XTick = sort([hInset.XTick,0.1]);
hInset.YLim = [0,1];
% hInset.CLim = ax.CLim;
%% Extract resonance of renormalized splitting
mode = 0;		% 0: lin, 1: log
figure(6); clf;
if mode
	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n));
% 	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n)-log10(ones(size(tresults.star.n,1),1)*tresults.star.n(1,:)));

	zlabel('$\log_{10}\left<n_k\right>$');
else
	surf(tresults.star.omega,tresults.star.t,tresults.star.n);
	zlabel('$\left<n_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = '$\left<n_k\right>$';
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1]);
% 	set(gca,'clim',[-4,1]);
	cb.Title.String = '$\log_{10}\left<n_k\right>$';
end
if exist('axIn','var')
	axIn.delete;
end
ax=gca; oldPos = ax.Position;
axIn = axes('Position',[oldPos(1)+0.32, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);		% use this if main with colorbar
% axIn = axes('Position',[oldPos(1)+0.36, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);	% use this if main without colorbar
axIn.TickDir = 'out';
axIn.FontSize=14; axIn.XColor=[0.7,0.7,0.7];axIn.YColor=axIn.XColor;

hold all;
surf(tresults.star.omega,tresults.star.t, ax.Children.ZData);
set(gca,'View',[0 90]);
shading interp
axis tight
axIn.XLim = [0,0.2]; axIn.YLim = [0,500];
[M,I] = max(ax.Children.ZData(:,1:200),[],2);
plot3(tresults.star.omega(I(2:end)),tresults.star.t(2:end),ones(1,length(tresults.star.t)-1).*1000,':r')
plot3([0.1,0.1],[0,350],[10,10],':black');

tresults.star.omega(I(2:end))

%% Analyze Chain Parameters SBM_genpara.m

para.s = 1; para.alpha=0.2; para.L=10; para.chain.spectralDensity = 'Leggett_Hard'; para.chain.mapping='OrthogonalPolynomials';para.chain.method='Analytic';
para1 = SBM_genpara(para);
% para.z=1; para.Lambda=1.01; para.chain.mapping = 'OrthogonalPolynomials';para.chain.discrMethod = 'Numerical';para.chain.discretization = 'Linear';para.chain.method = 'Stieltjes';
% para2 = SBM_genpara(para);
%%
figure(1);clf;hold on;
% scatter(para1.epsilon,para1.t); %set(gca,'yscale','log');
plot([para1.epsilon, para1.t]);

%% Map from chain to star:
obs = getObservable({'bath2correlators'},mps,Vmat,para);
fignum = 5;
figure(fignum);clf;
surf(real(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$Re(\left<a^\dagger_m a_n\right>)$')
fignum = fignum+1; figure(fignum);clf;
surf(imag(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$Im(\left<a^\dagger_m a_n\right>)$')
fignum = fignum+1; figure(fignum);clf;
surf(abs(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$\left|\left<a^\dagger_m a_n\right>\right|$')

%%
obs = getObservable({'staroccupation'},mps,Vmat,para);
figure(10);clf;
plot(obs(1,:),obs(2,:));

%% Slice current figure in x
f = get(gcf); ax = get(gca); hold on;
% plLine = plot3([0,0],ax.YLim,[1 1].*max(ax.ZLim),'r');
fnew = figure(100*f.Number);
pl = plot(ax.Children(end).YData, ax.Children(end).ZData(:,100));
pos = f.Position; f.Position = [pos(1),pos(2),840,pos(4)]; pos = f.Position;
pl.UserData=1;
pl.XDataSource = 'ax.Children(end).YData';
pl.YDataSource = 'ax.Children(end).ZData(:,pl.UserData)'; refreshdata;
xlabel(ax.YLabel.String);
ylabel(ax.ZLabel.String);
com = ', ';
% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(ax.Children(end).ZData,2)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(5);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'UserData',round(source.Value)));
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(fnew));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$ x = %g $',ax.Children(end).XData(pl.UserData))));
% hPl.addlistener(hProp,'PostSet',@(src,event) set(plLine,'XData',[1,1].*ax.Children(end).XData(pl.UserData)));

%% Slice current figure in y
f = get(gcf); ax = get(gca); hold on;
% plLine = plot3([0,0],ax.YLim,[1 1].*max(ax.ZLim),'r');
fnew = figure(101*f.Number);
pl = plot(ax.Children(end).XData, ax.Children(end).ZData(100,:));
pos = f.Position; f.Position = [pos(1),pos(2),840,pos(4)]; pos = f.Position;
pl.UserData=1;
pl.XDataSource = 'ax.Children(end).XData';
pl.YDataSource = 'ax.Children(end).ZData(pl.UserData,:)'; refreshdata;
xlabel(ax.XLabel.String);
ylabel(ax.ZLabel.String);
com = ', ';
% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(ax.Children(end).ZData,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(5);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'UserData',round(source.Value)));
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(fnew));
hPl.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$ y = %g $',ax.Children(end).YData(pl.UserData))));
% hPl.addlistener(hProp,'PostSet',@(src,event) set(plLine,'XData',[1,1].*ax.Children(end).XData(pl.UserData)));

%% Spectra of gcf
f = get(gcf); ax = get(gca); hold on;
fnew = figure(101*f.Number);

tslices = [2,11,31,51, 164];
pl = plot(tresults.star.omega, ax.Children.ZData(tslices,:));
legend(strread(num2str(tresults.star.t(tslices)),'%s'));
xlabel('$\omega_k / \omega_c$');
ylabel(ax.ZLabel.String);
xlim([0,0.2]);

%% move from mps{i} to mps{i+1}
para.sweepto = 'r'; j = para.sitej;
[mps{j}, Cn, para,results] = prepare_onesite_truncate(mps{j}, para,j,results);
mps{j+1} = contracttensors(Cn,2,2, mps{j+1},3,1);
para.sitej = j+1

%% Move from mps{i} to Vmat{i}
j = para.sitej;
[Amat,V] = prepare_onesiteAmat(mps{j},para,j);              % left-normalize A, SVD in n.
[BondDimLeft, BondDimRight, OBBDim]  = size(Amat);
Vmat{j} = Vmat{j} * transpose(V);							% set focus on Vmat: V_(n,n~)

%% prepare HOSVD test
% need to have focused Vmat{j}
j = para.sitej;
Vmat{j} = reshape(Vmat{j},[para.dk(:,j)',para.d_opt(j)]);
para.d_opt(3,:) = para.d_opt(1,:);
para.d_opt(1:2,:) = para.dk(:,:);
%% HOSVD test
[S U sv ,para1] = hosvd(Vmat{j},para,3);

%% calculate t-RDM
rdm = cell(size(tmps,1),1);
for i = 1:size(tmps,1)
	rdm{i} = getObservable({'rdm',[1 2]},tmps(i,:),tVmat(i,:),para);
	rdm{i} = reshape(rdm{i}, [4,4]);
end
%% Plot Slider-RDM
fignum = 10;
f = figure(fignum); clf; x = rdm;
f.Position(3) = 840; pos = f.Position; ax = gca;
pl=surf(real(x{1}));
% setting auto-refresh
hAx = handle(ax); hProp = findprop(hAx,'UserData');
hAx.addlistener(hProp,'PostSet',@(src,event) refreshdata(ax));
hAx.addlistener(hProp,'PostSet',@(src,event) title(sprintf('$t = %g$',para.tdvp.t(ax.UserData))));

ax.UserData = 1;
cb = colorbar; cb.Title.Interpreter = 'latex';
p = 'ax.UserData';							% prefix to shorten expressions
zlabel('$\rho$');
pl.ZDataSource = sprintf('real(x{%s})',p); refreshdata;
% refreshdata;
set(gca,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;rotate3d on;axis tight;
ax.ZLim = [-0.1,0.5]; ax.CLim = ax.ZLim;

% slider definition and UserData setting:
sldmax = length(x);
sld = javax.swing.JScrollBar(0,1,1,1,sldmax);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.05,5,200,15], f);
sld.setUnitIncrement(1); sld.setBlockIncrement(10);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(ax,'UserData',round(source.Value)));

%% RDM of site 2 from t-RDM of sites(1,2)
ncon({reshape(rdm{20},[2,2,2,2])},{[-1,1,-2,1]})
%% find process acting on site 2 in Pauli Basis
son = zeros(4,2,2);				% sigma ortho normal basis
son(4,:,:) = eye(2)/sqrt(2); son(1,:,:) = sx/sqrt(2);
son(2,:,:) = sy/sqrt(2);     son(3,:,:) = sz/sqrt(2);
EAm = zeros(2,2,4);
t = 3;		% which timeslice
rho = reshape(rdm{t},[2,2,2,2]);
for i = 1:4
	EAm(:,:,i) = ncon({rho,			squeeze(son(i,:,:))'},...
					  {[-1,2,-2,1], [1,2]})*2;
end
%% find process acting on site 2 in cell-basis
son = zeros(4,2,2);				% ortho normal basis
son(1,1,1) = 1; son(2,2,1) = 1;
son(3,1,2) = 1; son(4,2,2) = 1;
% EAm = zeros(2,2,4);
t = 50;		% which timeslice
rho = reshape(rdm{t},[2,2,2,2]);
for i = 1:4
	EAm(:,:,i) = ncon({rho,			squeeze(son(i,:,:))'},...
					  {[-1,2,-2,1], [1,2]})*2;
% 	FAm(:,:,i) = ncon({rho,			squeeze(son(i,:,:))'},...
% 					  {[2,-1,1,-2], [1,2]})*2;
end
lambda = reshape(EAm,[4,4]);		% this is the desired dynamical map epsilon!

%% Learn Epsilon for all Times
son = zeros(4,2,2);				% ortho normal basis
son(1,1,1) = 1; son(2,2,1) = 1;
son(3,1,2) = 1; son(4,2,2) = 1;
Epsilon = zeros(4,4,length(rdm));
for t = 1:length(rdm)
	rho = reshape(rdm{t},[2,2,2,2]);
	for i = 1:4
		EAm(:,:,i) = ncon({rho,			squeeze(son(i,:,:))'},...
						  {[-1,2,-2,1], [1,2]})*2;
	end
	Epsilon(:,:,t) = reshape(EAm,[4,4]);		% this is the desired dynamical map epsilon!
end
%% Derive the Tranfer Tensors
n = size(Epsilon,3);
T = zeros(4,4,n);		% T_n,m
T(:,:,1) = Epsilon(:,:,2);
for i = 2:n-1
	T(:,:,i) = Epsilon(:,:,i+1);
	for j = 2:i
		T(:,:,i) = T(:,:,i) - T(:,:,i+1-j)*Epsilon(:,:,j);
	end
	anorm(i) = norm(T(:,:,i));
end

%% reconstruct Dynamics using the TTM: (DEPRECATED. Better method above!)
finalT = 300;
% x = res{1,1}; tresults = x.tresults; para = x.para;
[sx,sy,sz] = spinop('Z');
tic;
n = round(finalT/para.tdvp.deltaT)+1;
rhoT = zeros(length(tresults.TTM.T)*4,1);
Esigma = zeros(n,3);
T = reshape(tresults.TTM.T, 4,[]);			% creates [T(1) T(2) T(3) ...]
for i = 1:n
	if i == 1
		rho = [1,0,0,0]';
	else
		rho = T*rhoT;
	end
	rhoT = [rho; rhoT(1:end-4)];					% prepend new vector rho(i)
	rho = reshape(rho,[2,2]);						% reshape rho(i) for observables
	Esigma(i,1) = trace(sx*rho);
	Esigma(i,2) = trace(sy*rho);
	Esigma(i,3) = trace(sz*rho);
end
tresults.spin.sx = real(Esigma(:,1));
tresults.spin.sy = real(Esigma(:,2));
tresults.spin.sz = real(Esigma(:,3));
tresults.t       = 0:para.tdvp.deltaT:finalT;
figure(1);clf;hold all;
plot(tresults.t,tresults.spin.sx);
plot(tresults.t,tresults.spin.sy);
plot(tresults.t,tresults.spin.sz);
set(gca,'ylim',[-1,1]);
toc

%% Plot Memory Kernel ( elements of T)
fignum = 10; f = figure(fignum); clf; hold all; ax = gca;
% x = res{1,1}; tresults = x.tresults; para = x.para;
t = (1:length(tresults.TTM.T)).*para.tdvp.deltaT;
leg = cell(0,0);
fun = @real;
plot(t(2:end),fun(squeeze(tresults.TTM.T(1,1,2:end)./(para.tdvp.deltaT.^2)))); leg = [leg,{'11 \rightarrow 11'}]; % 11->11 = - 11->22
plot(t(2:end),fun(squeeze(tresults.TTM.T(2,1,2:end)./(para.tdvp.deltaT.^2)))); leg = [leg,{'11 \rightarrow 21'}];% 11->21 = 11->12
plot(t(2:end),fun(squeeze(tresults.TTM.T(2,3,2:end)./(para.tdvp.deltaT.^2)))); leg = [leg,{'12 \rightarrow 21'}];% 12->21
plot(t(2:end),fun(squeeze(tresults.TTM.T(2,2,2:end)./(para.tdvp.deltaT.^2)))); leg = [leg,{'21 \rightarrow 21'}];% 12->21
% plot(t(2:end),fun(squeeze(x.tresults.TTM.T(1,2,2:end)./(x.para.tdvp.deltaT.^2)))); leg = [leg,{'21 \rightarrow 11'}];% 12->21 basically 0
legend(leg);

% ax.YScale = 'log';

%% Convert tresults to single precision
% tresults.star.x = single(tresults.star.x);
tresults.star.n = single(tresults.star.n);
tresults.star.t = single(tresults.star.t);
tresults.star.omega = single(tresults.star.omega);
tresults.spin.sx = single(tresults.spin.sx);
tresults.spin.sy = single(tresults.spin.sy);
tresults.spin.sz = single(tresults.spin.sz);
% tresults.j = single(tresults.j);
% tresults.t = single(tresults.t);
tresults.nx = single(real(tresults.nx));
save(para.tdvp.filenameSmall, 'para','tresults');

%% Diagonalise Chain
n = para.L-1; beta = 40;
A = full(sparse(1:n,1:n,para.chain{1}.epsilon,n,n)+sparse(2:n,1:n-1,para.chain{1}.t(2:end),n,n)+sparse(1:n-1,2:n,para.chain{1}.t(2:end),n,n));
[V,D] = eig(A);
occDiag = 1./(exp(beta.*diag(D))-1);
% plot(diag(D),occDiag);set(gca,'xscale','log');		% display values
Aocc = V*diag(occDiag)*V.';
% plot(diag(A),diag(Aocc));set(gca,'xscale','log');
plot(2:para.L,diag(Aocc))
%% Diagonalise Chain with spin
n = para.L; beta = 40;
sigmaZ= [1,-1]; ii = 1; Aocc = zeros(n,n,2);
for sz = sigmaZ
	A = full(sparse(1:n,1:n,[-para.hz*sz/2; para.chain{1}.epsilon],n,n)+...
		sparse(2:n,1:n-1,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n)+...
		sparse(1:n-1,2:n,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n));
	[V,D] = eig(A);
	occDiag = 1./(exp(beta.*diag(D))-1); occDiag(1) = 0;
	% plot(diag(D),occDiag);set(gca,'xscale','log');		% display values
	Aocc(:,:,ii) = V*diag(occDiag)*V.'; Aocc(1,1,ii) = 0;				% 1-site is meaningless!
	ii = ii+1;
	% plot(diag(A),diag(Aocc));set(gca,'xscale','log');
% plot(1:para.L,diag(Aocc))
end
pdown = (1-tresults.spin.sz(end))/2; pup = 1-pdown;
% pdown = 0.5; pup = 1-pdown;
occN = [diag(Aocc(:,:,1)),diag(Aocc(:,:,2))];
% plot(occN); l=legend('$\left< n \right> @ \left| \uparrow \right>$','$\left< n \right> @ \left| \downarrow \right>$'); l.Interpreter = 'latex';
hold on;
% plot(occN(:,1)*pup + occN(:,2)*pdown);

%% Plot analytic Ohmic renormalization
figure(1); clf; hold all;
a = 0:0.01:1; d = 0.1;		% alpha delta
p1 = plot(a,d.^(1./(1-a)));
p2 = plot(a,((gamma(1-2.*a).*cos(pi.*a)).^(1./(2-2.*a))).*d.^(1./(1-a))); % this is from NIBA and much worse!

%% Plot NIBA ohmic SBM evolution
figure(1); clf; hold all;
% x = res{11}; para = x.para; tresults = x.tresults;
% a = para.chain{1}.alpha; d = 0.1; t = 0:0.1:180;
a = 0.01; d = 0.1; t = 0:0.1:2300;
% calculate the sum to eps precision
summ  = 0;
summP = inf;
n = 0;
while norm(abs(summP-summ)) > 1e-16
% while n<70
    summP = summ;
    summ = summ + ((-1).^n./gamma(1+2.*(1-a).*n).*(d.^(1./(1-a)).*t ).^(2.*n.*(1-a)));
    n = n + 1;
end
% plot the result
p1 = plot(t,summ);
xlim([0,2300]);
ylim([-1,1]);
% plot(tresults.t, tresults.spin.sz);

%% Plot iSBM multi-Chain Analytic
t = 0:0.2:200;%para.tdvp.tmax;
para.nChains = 1; para.chain{1}.alpha = 0.1; para.chain{1}.s = 1;
exponent = 0;
w_cutoff = 1;
for ii = 1:para.nChains
	a = para.chain{ii}.alpha; s = para.chain{ii}.s; wc = w_cutoff(ii);
	syms w
	f = int(2.*a.* w.^(s-2).* wc.^(1-s).*(1-cos(w.*t)),w,0,wc);
	exponent = exponent + vpa(f,3);
end
plot(t,exp(-exponent));
%%
plot(t,exponent./t)

%% Plot iSBM pure dephasing decoherence
t = 0:0.2:200;
a = 0.1; s=1; wc=1; wcut=inf;
Jsmooth = @(w) 2.*a.* w.^(s-2).* wc.^(1-s).*exp(-w/wc);
Jhard = @(w) 2.*a.* w.^(s-2).* wc.^(1-s);
exponent = integral(@(w) Jhard(w).*(1-cos(w.*t)),0,wc,'ArrayValued',true);

figure(2); clf; hold all;
plot(t,exp(-exponent));
plot(t(2:end),diff(exponent)./t(2))
%% Plot iSBM2C pure dephasing decoherence
t = 0:0.2:20;
a = 0.1; s=1; wc=[1,5]; wcut=inf;
exponent = 0;
for ii = 1:2
% 	Jsmooth = @(w) 2.*a.* w.^(s-2).* wc(ii).^(1-s).*exp(-w/wc(ii));
	Jhard = @(w) 2.*a.* w.^(s-2).* wc(ii).^(1-s);
	exponent = exponent+ integral(@(w) Jhard(w).*(1-cos(w.*t)),0,wc(ii),'ArrayValued',true);
end
figure(2); clf; hold all;
plot(t,exp(-exponent)*0.5);
plot(t(2:end),diff(exponent)./t(2))
%% Create ideal 2D <n> for Thermal cooling
n = para.L; beta = 40; dt = 0.1;
sigmaZ= [1,-1]; ii = 1; Aocc = zeros(n,n,2);
occN = zeros(round(beta/dt)-3,n,2);
jj = 1;
for t = 4*dt:dt:beta
	ii = 1;
	for sz = sigmaZ
		A = full(sparse(1:n,1:n,[-para.hz*sz/2; para.chain{1}.epsilon],n,n)+...
				sparse(2:n,1:n-1,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n)+...
				sparse(1:n-1,2:n,[para.chain{1}.t(1)*sz; para.chain{1}.t(2:end)],n,n));
		[V,D] = eig(A);
		occDiag = 1./(exp(t.*diag(D))-1); occDiag(1) = 0;
		Aocc(:,:,ii) = V*diag(occDiag)*V.'; Aocc(1,1,ii) = 0;				% 1-site is meaningless!
		ii = ii+1;
	end
	occN(jj,:,:) = [diag(Aocc(:,:,1)),diag(Aocc(:,:,2))];
	jj = jj+1;
end
clf;hold on;
surf(1:n,4*dt:dt:beta,log10(occN(:,:,1)));
axis tight; rotate3d on; shading interp;

%% Export to CSV
x = tresults.t; x = reshape(x,numel(x),1);
y = tresults.spin.visibility; y = reshape(y,numel(y),1);
csvwrite('visibility02.dat',[x,y]);

%% Plot series of figures and save
%% rhoii-osc-res-med
% f = figure;
for ii = [1]
% 	res(ii).lastIdx = 5001;
% h=res(ii).plot('rhoii-osc-res-med','-cmev','-dist','-norm','-resetcolororder',f);
h=res(ii).plot('rhoii-osc-res-med','-cmev','-dist','-norm');
% h.f.Visible = 'off';
% [h.ax.Units] = deal('pixels');
% [h.ax.ActivePositionProperty] = deal('outerposition');
% h.ax(end).ActivePositionProperty = 'outerposition';
% h.f.Units = 'centimeters'; h.f.Position(3) = h.f.Position(3)*0.66; h.f.Units = 'pixels';
if length(h.ax) > 1
	arrayfun(@(x) set(x,'YTick',[0,x.YTick(end)/2,x.YTick(end)]),h.ax);
else
	h.ax.YTick = 0:5;
end
[h.ax.XLim] = deal([0,1999]);
[h.ax.Color] = deal('none');
set(h.pl,{'Displayname'},{'TT','LE+','LE-','CT+','CT-'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) legend(x,'boxoff'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);

formatPlot(h.f,'twocolumn-single');
% h.f.Units = 'centimeters'; h.f.Position(3) = 2.5*h.f.Position(4); h.f.Units = 'pixels';
% export_fig(sprintf('img/%d',h.f.Number),'-transparent','-png','-m2', h.f);
end
%% chain-n-site2-ft
f = figure;
for ii = [4:5]
h = res(ii).plot('chain-n-site2-ft','-cmev','-dist','-norm','-resetcolororder',f);
% h = res(ii).plot('chain-n-site2-ft','-cmev','-dist','-norm');
% h.f.Units = 'centimeters'; h.f.Position(3) = h.f.Position(3)*0.66; h.f.Units = 'pixels';
set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) legend(x,'boxoff'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
if length(h.ax) > 1
	arrayfun(@(x) set(x,'YTick',[0,x.YTick(end)/2,x.YTick(end)]),h.ax);
else
	h.ax.YTick = 0:7;
end
ylabel('$|FT (n_2)|$');
[h.ax.XLim] = deal([0,1999]);
[h.ax.Color] = deal('none');
% formatPlot(h.f);
drawnow;
formatPlot(h.f,'twocolumn-single');
% h.f.Units = 'centimeters'; h.f.Position(3) = 2.5*h.f.Position(4); h.f.Units = 'pixels';
% export_fig(sprintf('img/%d',h.f.Number),'-transparent','-png','-m2', h.f);
end
%% chain-x-site2-ft
for ii = [1:1]
h = res(ii).plot('chain-x-site2-ft','-cmev','-dist','-norm');
% h.f.Units = 'centimeters'; h.f.Position(3) = h.f.Position(3)*0.66; h.f.Units = 'pixels';
set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) legend(x,'boxoff'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
if length(h.ax) > 1
	arrayfun(@(x) set(x,'YTick',[0,x.YTick(end)/2,x.YTick(end)]),h.ax);
else
	h.ax.YTick = 0:7;
end
ylabel('$|FT (x_2)|$');
[h.ax.XLim] = deal([0,1999]);
[h.ax.Color] = deal('none');
% formatPlot(h.f);
drawnow;
formatPlot(h.f,'twocolumn-single');
% h.f.Units = 'centimeters'; h.f.Position(3) = 2.5*h.f.Position(4); h.f.Units = 'pixels';
% export_fig(sprintf('img/%d',h.f.Number),'-transparent','-png','-m2', h.f);
end
%% chain-x2-site2-ft
for ii = [1:1]
h = res(ii).plot('chain-x2-site2-ft','-cmev','-dist','-norm');
% h.f.Units = 'centimeters'; h.f.Position(3) = h.f.Position(3)*0.66; h.f.Units = 'pixels';
set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) legend(x,'boxoff'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
if length(h.ax) > 1
	arrayfun(@(x) set(x,'YTick',[0,x.YTick(end)/2,x.YTick(end)]),h.ax);
else
	h.ax.YTick = 0:7;
end
ylabel('$|FT (x^2_2)|$');
[h.ax.XLim] = deal([0,1999]);
[h.ax.Color] = deal('none');
% formatPlot(h.f);
drawnow;
formatPlot(h.f,'twocolumn-single');
% h.f.Units = 'centimeters'; h.f.Position(3) = 2.5*h.f.Position(4); h.f.Units = 'pixels';
% export_fig(sprintf('img/%d',h.f.Number),'-transparent','-png','-m2', h.f);
end
%% rhoij-imag
for ii = [1]
h = res(ii).plotSld1DFT('rhoij-imag','-cmev','-dist');
% h.f.Units = 'centimeters'; h.f.Position(3) = h.f.Position(3)*0.66; h.f.Units = 'pixels';
% set(h.pl,{'Displayname'},{'$A_{1,1}$','$A_{1,2}$','$A_{2}$','$B_{1}$','$B_{2,1}$','$B_{2,2}$','$B_{2,3}$'}')
arrayfun(@(x) legend(x,'show'),h.ax);
arrayfun(@(x) legend(x,'boxoff'),h.ax);
arrayfun(@(x) grid(x,'on'),h.ax);
ylabel('$|FT Im(\rho_{ij})|$');
[h.ax.XLim] = deal([0,15000]);
[h.ax.Color] = deal('none');
% formatPlot(h.f);
drawnow;
formatPlot(h.f,'twocolumn-single');
% h.f.Units = 'centimeters'; h.f.Position(3) = 2.5*h.f.Position(4); h.f.Units = 'pixels';
% export_fig(sprintf('img/%d',h.f.Number),'-transparent','-png','-m2', h.f);
end

%% add stems into norm dist reaction coord plot - L18
x = res(1);
% figure
hold all
for nc = 1:7
% 	stem(x.para.chain{nc}.dataPoints(:,1)*8065.6,nc-1+abs(x.para.chain{nc}.dataPoints(:,2))/max(abs(x.para.chain{nc}.dataPoints(:,2))),'BaseValue',nc-1);
% or:
	dat = abs(x.para.chain{nc}.dataPoints);
	dat(:,1) = dat(:,1)*8065.6;							% ev to cm
	dat(:,2) = dat(:,2)/max(dat(:,2));					% normalise
	for ii = 1:size(dat,1)
		line([1 1]*dat(ii,1),[0,1]*dat(ii,2)-nc+7,'color',[0 0 0],'linewidth',2);
	end
end
%% add stems into norm dist reaction coord plot - L2
x = res(1);
% figure
hold all
for nc = 1:7
% 	stem(x.para.chain{nc}.dataPoints(:,1)*8065.6,nc-1+abs(x.para.chain{nc}.dataPoints(:,2))/max(abs(x.para.chain{nc}.dataPoints(:,2))),'BaseValue',nc-1);
% or:
	dat = [x.para.chain{nc}.epsilon(1)*8065.6, 1];
	line([1 1]*dat(1),[0,1]-nc+7,'color',[0 0 0],'linewidth',2);
end

%%
h=res(1).plotSld1D('chain-n-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end
%%
h=res(1).plotSld1D('chain-n-d-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
legend('$TT$','$LE^+$','$LE^-$','$CT^+$','$CT^-$')
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end
%%
h=res(1).plotSld1D('chain-x-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end
%%
h=res(1).plotSld1D('chain-x-d-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
legend('$TT$','$LE^+$','$LE^-$','$CT^+$','$CT^-$')
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end
%%
h=res(1).plotSld1D('chain-x2-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end
%%
h=res(1).plotSld1D('chain-x2-d-t','-fsev');
formatPlot(h.f);
h.sld{1}.setValue(2);
legend('$TT$','$LE^+$','$LE^-$','$CT^+$','$CT^-$')
for ii = 1:7
	h.sld{2}.setValue(ii);
	export_fig(sprintf('img/%s-%d',get(gcf,'Name'),ii),'-transparent','-png',h.ax);
end

%% Analyse TreeMPS node for inter-chain entanglement
% which could be reduced by chain combination
%x = TDVPData('H:\Documents\Theory\schroederflorian-vmps-tdvp\TDVP-Git\Data\20160325-1237-26-DPMES5-7C-Tree-v73TCMde11-L18CT0LE+\results-Till5000Step0.1v73-OBBmax60-Dmax(7-60)-expvCustom700-1core.mat');
% load('Data\20160325-1237-26-DPMES5-7C-Tree-v73TCMde11-L18CT0LE+\results-Till5000Step0.1v73-OBBmax60-Dmax(7-60)-expvCustom700-1core.mat','treeMPS','para'); % used for ER-Tree2
load('..\cacheComputations\20170408-0329-34-DPMESclust7-1-Tree-v77TCMde10-L100CT0LE+\results-Till1000Step0.5v77-OBBmax40-Dmax(5-7)-expvCustom270-1core.mat','treeMPS','para');
%% Get vNE for each partition
% lower vNE: lower entanglement!
res = getTensorVNE(squeeze(treeMPS.mps{1}));
%% Get vNE Adjacency of TreeMPS Root node
res = TDVPData.getTensorVNEAdjacency(squeeze(treeMPS.mps{1}));
% Get graph for adjacency
g = graph(res,'upper');
%% Plot the Graph
figure(30);
h = plot(g,  'LineWidth', 5*g.Edges.Weight/max(g.Edges.Weight), 'Layout', 'circle'); %		'EdgeLabel',g.Edges.Weight,
%% [5,9] and [6,7] SV:
A = treeMPS.mps{1}; dA = size(A);
[B, dB] = tensShape(A,'unfoldiso', [5,9],dA);
[~,S,~] = svd2(B);


%% Export adiabatic RC displacement
% only takes the strongest adiabatic state!
% obj = TDVPData('20161031-1428-05-DPMES-Tree2-Tree-v76TCMde11-L18CT0LE+\results-Till5000Step0.5v76-OBBmax60-Dmax(100-20)-expvCustom270-1core-small.mat');
obj = TDVPData('20170409-1735-52-DPMES-Tree4-Tree-v77TCMde11-L100CT0LE+\results-Till5000Step0.5v77-OBBmax40-Dmax(50-20)-expvCustom270-1core-small.mat');
maxN = 440;
x = zeros(maxN,7);		% t x nChain
for ii = 1:7
	rcCol = find(~all(squeeze(obj.xCa(1:maxN,:,1,ii))==0),1);
	x(:,ii) = squeeze(real(obj.xCa(1:maxN,rcCol,1,ii)));		% t x L x state x nChain
end
csvwrite('adiabX.dat',x);

%% Analyse adiabatic potential
adiabEnv = 1;
out = real(cell2mat(arrayfun(@(i) sort(eig(squeeze(tresults.Heff(i,adiabEnv,:,adiabEnv,:))))',[1:tresults.lastIdx-1]','UniformOutput',false)));
figure;
plot(tresults.t(2:tresults.lastIdx),out(1:end,:))

%% try to simultaneously extract distribution of wavefunction
%
for ii = 1:tresults.lastIdx
	tresults.system.state
end
%% Save all currently opend figures
f_handles = get(0,'children');
for ii = 1:length(f_handles)
% 	export_fig(['img/',f_handles(ii).Name],'-transparent','-png','-pdf','-m2', f_handles(ii));
% 	export_fig(['img/',f_handles(ii).Name],'-transparent','-png','-m2', f_handles(ii));
	export_fig(sprintf('img/%d',f_handles(ii).Number),'-transparent','-painters','-png','-m2', f_handles(ii));
% 	figure(f_handles(ii));
% 	export_fig(sprintf('img/%d',f_handles(ii).Number),'-transparent','-png','-m2', gca);
end

%% Deserialise all Variables
Vars = whos;
for ii = 1:size(Vars)
if strcmp(Vars(ii).class,'uint8')
eval(sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
end
end

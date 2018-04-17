function PlotData = ExtractAllData()
% Used to get useful data from results.mat of each subfolder.
% saves extractes data into file:
% Makes Plot to visualize data
%
% Revision: 23/06/2014 15:41
%
%%
% Overwrite extracted.mat ?
overwrite = 0;
global loadDataFromFile
loadDataFromFile = '/results.mat';
global saveDataToFile
saveDataToFile = '/extracted.mat';
global saveDataToFileError
saveDataToFileError = '/extractedError.mat';

regexpLINUX = ['[^:;]*'];
regexpWIN = ['[^;]*'];

doPlot = 0;
doAverage = 0;          % set = 1 for average over static disorder! 0 else
hasDisorder = 0;
loadResults = 1;

if strfind(system_dependent('getos'),'Windows')
    regexpUSE = regexpWIN;
elseif strfind(system_dependent('getos'),'Linux')
    regexpUSE = regexpLINUX;
else
    disp('OS not recognised');
end
dirs = regexp(genpath(pwd),regexpUSE,'match');      % get all subfolders of current directory
dirs = dirs(2:end);                                 % remove current folder

%%
if ~loadResults
    for i = 1:1:length(dirs)
        if (overwrite == 0 && exist([dirs{i},saveDataToFile],'file')) || ~exist([dirs{i},loadDataFromFile],'file')
            continue
        end
        % Extract Data and save to File, per folder
        ExtractOneData(dirs{i});
    end
end
%% Load already existing PlotData??

% generate Parameter and Essentials list
PlotData = zeros(length(dirs),16);   % (s, alpha, log(shifterror))
PPCWavefunction = zeros(length(dirs),16);   % only applies to Rs.Molischianum
BosonChainOccupation = zeros(length(dirs),47);
for i = 1:1:length(dirs)
    try
        if loadResults
            load([dirs{i},'/results.mat']);         % for MLSBM
        else
            load([dirs{i},'/extracted.mat']);
        end

        %disp(dirs{i})                          % only disp if data is corrupt or not readably.
        if ~strcmp(para.model,'MLSpinBoson')
            PlotData(i,1) = para.s;
            PlotData(i,2) = para.alpha;
            PlotData(i,3) = -log10(extracted.ShiftMaxError);
            PlotData(i,4) = -log10(extracted.NMaxError);
            PlotData(i,5) = extracted.timeMins;
            PlotData(i,6) = extracted.Sx;       % <sx>
            PlotData(i,7) = para.hx;            % = delta
            PlotData(i,8) = results.spin.sz;    % <sz>
            if max(para.dk) < 1000
                PlotData(i,9) = norm(abs(extracted.aaCreationCorrelator(1:end-1)-extracted.aaCreationCorrelator2));
            else
                PlotData(i,9) = -1;     % dimension too high--> no calculation
            end
            PlotData(i,10) = max(abs(results.nx)); % Max(<n>) to see changes in chain population
            PlotData(i,11) = para.hz;           % = -epsilon
            PlotData(i,12) = para.loop;
            PlotData(i,13) = extracted.entropy; % Entanglement entropy (cf. Le Hur 2008)
            PlotData(i,14) = results.Amat_vNE(1);
            PlotData(i,15) = results.Vmat_vNE(2);
            PlotData(i,16) = results.Vmat_sv{1, 2}(1);  % find relationship
            BosonChainOccupation(i,:) = results.nx;      % save the entire chain occupation
        else
            if ~isfield(results,'time');
                disp(['Incomplete: ',dirs{i}]);
                continue;
            end
            PlotData(i,1) = para.Lambda;
            PlotData(i,2) = para.z;
            PlotData(i,3) = results.participation;  % extracted.participation; = inverse participation of Cogdell 2006
            PlotData(i,4) = max(abs(results.nx));   % Max(<n>) to see changes in chain population
            PlotData(i,5) = results.E;
            PlotData(i,6) = para.MLSB_p;            % extracted.period;
            if isfield(para,'MLSB_etaFactor')
                PlotData(i,7) = para.MLSB_etaFactor;
            else
                PlotData(i,7) = 1;      % case before MLSB_etaFactor was introduced.
            end
            PlotData(i,8) = para.loop;
            if isfield(results,'tunnelEnergy')
                PlotData(i,9) = results.tunnelEnergy;
            else
                PlotData(i,9) = getObservable({'tunnelenergy',op},mps,Vmat,para);
            end
            PlotData(i,10) = -log10(results.E- min(cell2mat(results.EvaluesLog)));   % Deviation from lowest energy
            if isfield(para,'MLSB_disSDV')
                PlotData(i,11) = para.MLSB_disSDV;      % strength of disorder
            else
                PlotData(i,11) = 0;
            end
            PlotData(i,12) = results.Amat_vNE(1);       % prop to entanglement entropy
            % copy the wavefunction of the ring
            PPCWavefunction(i,:) = diag(getObservable({'rdm',1},mps,Vmat,para));
            BosonChainOccupation(i,:) = results.nx;      % save the entire chain occupation
        end
    catch

    end

    clear('para','extracted','results');
end
%dlmwrite('20140310-Delta-Alpha-Scan.txt',sortrows(PlotData,7));
save([datestr(now,'yyyymmdd'),'-ExtractedResults.mat'],'PlotData','PPCWavefunction','BosonChainOccupation');  % for MLSBM
%save([datestr(now,'yyyymmdd'),'-ExtractedResults.mat'],'PlotData','BosonChainOccupation');

%% Execute this to AVERAGE over DISORDER
%   Use Variables later by assigning PlotData = PlotDataAvg;
if doAverage && hasDisorder
    uniqueCombinations = unique([PlotData(:,6), PlotData(:,7),PlotData(:,11)],'rows');
    PlotDataAvg = zeros(size(uniqueCombinations,1),size(PlotData,2));
    PPCWavefunctionAvg = zeros(size(uniqueCombinations,1),size(PPCWavefunction,2));
    % uncomment if not working:
    BosonChainOccupationAvg = zeros(size(uniqueCombinations,1),size(BosonChainOccupation,2));

    for i = 1:size(uniqueCombinations,1)
        id = PlotData(:,6)==uniqueCombinations(i,1) & PlotData(:,7)==uniqueCombinations(i,2) ...
            &PlotData(:,11)==uniqueCombinations(i,3);
        PlotDataAvg(i,:) = mean(PlotData(id,:));
        PPCWavefunctionAvg(i,:) = mean(PPCWavefunction(id,:));
        BosonChainOccupationAvg(i,:) = mean(BosonChainOccupation(id,:));
    end
    PlotDataAll = PlotData;
    PPCWavefunctionAll = PPCWavefunction;
    BosonChainOccupationAll = BosonChainOccupation;
end

if ~doPlot
	return;
end

%% Make Plots if doPlot = 1;

colormap jet
figure(1)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,3),40,PlotData(:,3),'fill')
title('Deviation of Shifts');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$log_{10}(\Delta \delta)$');
figure(2)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,4),40,PlotData(:,4),'fill')
title('Deviation of $<n>$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$log_{10}(\Delta <n>)$');
figure(3)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,5),40,PlotData(:,5),'fill')
title('Time needed');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('t in min');
%% Select with respect to certain conditions:
% Delta = 7; alpha = 2
figure(50)
hold all
select = PlotData(:,2)==0.5;
pd = sortrows(PlotData(select,:),-11);
delta = pd(:,7); alpha = pd(:,2);
reference = - 4./pi.*sqrt(delta.^(1./(1-alpha)).*(4/pi).^(1./alpha -1)).*log(delta.^(1./(1-alpha)).*(4/pi).^(1./alpha -1));
%plot(-pd(:,11),(reference-pd(:,6))./reference)
plot(-pd(:,11), pd(:,6));
plot(-pd(:,11),reference)
xlabel('$\varepsilon$')
ylabel('$<\sigma_x>$')
set(gca,'XScale','log');
formatPlot(50)
%%
figure(51)
hold all
select = PlotData(:,2)==0.5;
pd = sortrows(PlotData(select,:),-11);
reference = -1+pd(:,7).^2.*log(-pd(:,11));
plot(-pd(:,11),(reference-pd(:,6))./reference)
set(gca,'XScale','log');
%% Plot <sz> vs delta, alpha
figure(4)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,6),30,PlotData(:,7),'fill')
title('$<\sigma_x>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_x>$');
%% Plot <sz> vs delta, alpha
figure(5)
subplot(1,2,1)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,8),10,PlotData(:,8),'fill')
title('$<\sigma_z>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
subplot(1,2,2)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,8),10,PlotData(:,8),'fill')
set(gca,'View',[0 90]);
title('$<\sigma_z>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
%% Plot <sx> vs delta
select =  PlotData(:,2)==0.5;
IndexedMatrix = [PlotData(select,7) PlotData(select,6)];
IndexedMatrix = sortrows(IndexedMatrix,1);                      % sort with respect to y axis
figure(14);
plot(IndexedMatrix(:,1),IndexedMatrix(:,2));
xlabel('$\Delta$');
ylabel('$<\sigma_x>$')
formatPlot(14)
%% Plot <sz> vs delta
select =  PlotData(:,2)==0.5;
IndexedMatrix = [PlotData(select,7) PlotData(select,8)];
IndexedMatrix = sortrows(IndexedMatrix,1);                      % sort with respect to y axis
figure(15);
plot(IndexedMatrix(:,1),IndexedMatrix(:,2));
xlabel('$\Delta$');
ylabel('$<\sigma_z>$')
formatPlot(15)
%% Plot <sz> vs alpha, epsilon   USED
%   Similar to Fig2 in Le Hur 2013
defEpsilon = [1e-10,1e-8,1e-6,1e-4,1e-3,5e-3];
colors=hsv(length(defEpsilon));
figure(16);
hold all
i=1;                % count colors
for epsilon = defEpsilon
    select = PlotData(:,11) == -epsilon;          % select as binary array. can add other conditions
    Indexed = [PlotData(select,2) -abs(PlotData(select,8))];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    i=i+1;
end
xlabel('$\alpha$');
ylabel('$<\sigma_z>$')
formatPlot(16)
%legend({'$10^{-10}$','$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-3}$','$0.005$'})
%export_fig('20140601-Ohmic-sz-epsilon-alpha','-transparent','-pdf','-png','-painters')
%% Plot <sz> vs alpha, delta
%   Similar to ??
alphaC = [0.74,0.81,0.84,0.88,0.91,0.93,0.94,0.95,0.96,0.97,0.98,0.99];
defDelta = [1e-3,5e-3,10e-3,20e-3,30e-3,40e-3,50e-3,55e-3,65e-3,75e-3,85e-3,95e-3];
figure(17);
hold all
plot(log(defDelta),log(alphaC),'Marker','*','LineStyle','none');
%set(gca,'YScale','log','XScale','log');
xlabel('$\Delta$');
ylabel('$\alpha_c$')
%text(0.03,0.9,'$\Delta/\omega_c$');
formatPlot(17)
%legend(leg,'Location','West');
%export_fig('20140602-Ohmic-sz-delta-alpha','-transparent','-pdf','-png','-painters');
%% Plot <sx> vs alpha, epsilon
%   Similar to Fig2 in Le Hur 2013
defEpsilon = [1e-10,1e-8,1e-6,1e-4,1e-3,5e-3];
colors=hsv(length(defEpsilon));
figure(17);
hold all
i=1;                % count colors
for epsilon = defEpsilon
    select = PlotData(:,11) == -epsilon;          % select as binary array. can add other conditions
    Indexed = [PlotData(select,2) PlotData(select,6)];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    i=i+1;
end
xlabel('$\alpha$');
ylabel('$<\sigma_x>$')
text(0.5,0.9,'$\Delta/\omega_c=0.01$');
text(0.85,0.8,'$\varepsilon/\omega_c$');
formatPlot(17)
legend({'$10^{-10}$','$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-3}$','$0.005$'},'Location','East')
%export_fig('20140601-Ohmic-sx-epsilon-alpha','-transparent','-pdf','-png','-painters')

%% Plot <sx> vs epsilon at alpha == 1
%   As mentioned in LeHur2008, Eq.37
figure(17);
select = PlotData(:,2)==1;     % select as binary array. can add other conditions
Indexed = [-PlotData(select,11) PlotData(select,6)];          % dataset [x y]
Indexed = sortrows(Indexed,1);
plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none');
set(gca,'XScale','log');
xlabel('$\varepsilon$');
ylabel('$<\sigma_x>$')
formatPlot(17)

%% Plot \delta <sx> vs alpha, epsilon  --diff to reference
%  Plot error against formula.
%  As mentioned in LeHur2008, Eq.37
defEpsilon = [1e-10,1e-8,1e-6,1e-4,1e-3,5e-3];
colors=hsv(length(defEpsilon));
figure(17);
hold all
i=1;                % count colors
for epsilon = defEpsilon
    select = PlotData(:,11) == -epsilon & PlotData(:,2)>0.5;          % select as binary array. can add other conditions
    reference = PlotData(select,7)./(2.*PlotData(select,2)-1);      % delta/(2*alpha-1)
    Indexed = [PlotData(select,2) (reference-PlotData(select,6))./reference];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    i=i+1;
end
%set(gca,'YScale','log');
xlabel('$\alpha$');
ylabel('$\frac{<\sigma_{x,lit}>-<\sigma_{x,VMPS}>}{<\sigma_{x,lit}>}$')
text(0.7,0.9,'$\Delta/\omega_c=0.01$');
text(0.92,0.84,'$\varepsilon/\omega_c$');
formatPlot(17)
legend({'$10^{-10}$','$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-3}$','$0.005$'},'Location','East')
%export_fig('20140601-Ohmic-sx-epsilon-alpha','-transparent','-pdf','-png','-painters')

%% Plot <sx> vs alpha, delta 2D
%   Similar to ??
defDelta = [1e-3,5e-3,10e-3,20e-3,30e-3,40e-3,55e-3,65e-3,75e-3,85e-3,95e-3]; %50e-3,
colors=hsv(length(defDelta));
figure(17);
hold all
i=1;                % count colors
leg = cell(1);
for delta = defDelta
    select = PlotData(:,7) == delta;          % select as binary array. can add other conditions
    sum(select)
    Indexed = [PlotData(select,2) PlotData(select,6)];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    leg{i} = sprintf('%0.3f',delta);
    i=i+1;
end
xlabel('$\alpha$');
ylabel('$<\sigma_x>$')
text(0.85,0.9,'$\Delta/\omega_c$');
formatPlot(17)
legend(leg,'Location','East');
%export_fig('20140602-Ohmic-sx-delta-alpha','-transparent','-pdf','-png','-painters')
%% Plot \delta <sx> vs alpha, delta    --diff to reference
%   Similar to ??
defDelta = [1e-3,5e-3,10e-3,20e-3,30e-3,40e-3,55e-3,65e-3,75e-3,85e-3,95e-3]; %50e-3,
colors=hsv(length(defDelta));
figure(19);
hold all
i=1;                % count colors
leg = cell(1);
for delta = defDelta
    select = PlotData(:,7) == delta & PlotData(:,2)>0.7;          % select as binary array. can add other conditions
    reference = PlotData(select,7)./(2.*PlotData(select,2)-1);      % delta/(2*alpha-1)
    Indexed = [PlotData(select,2) (reference-PlotData(select,6))./reference];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    leg{i} = sprintf('%0.3f',delta);
    i=i+1;
end
%set(gca,'YScale','log');
xlabel('$\alpha$');
ylabel('$\frac{<\sigma_{x,lit}>-<\sigma_{x,VMPS}>}{<\sigma_{x,lit}>}$')
text(0.87,0.42,'$\Delta/\omega_c$');
formatPlot(19)
legend(leg,'Location','NorthEast');
%export_fig('20140602-Ohmic-sx-delta-alpha','-transparent','-pdf','-png','-painters')
%% Plot loop, time, Energy, Entropy, SV vs alpha, epsilon
%   Similar to Fig2 in Le Hur 2013
defEpsilon = [1e-10,1e-8,1e-6,1e-4,1e-3,5e-3];
colors=hsv(length(defEpsilon));
figure(18);
% 12: loop, 5: time, 10: max(<n>), 13: Entropy, 14: Amat_vNE(1), 15: Vmat_vNE(2)
%  6: sx, 8: sz, 16: max SV
PlotDat=PlotData(:,13);
hold all
i=1;                % count colors
for epsilon = defEpsilon
    select = PlotData(:,11) == -epsilon;          % select as binary array. can add other conditions
    Indexed = [PlotData(select,2) PlotDat(select)];          % dataset [x y]
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    i=i+1;
end
xlabel('$\alpha$');
ylabel('$S$')
formatPlot(18)
%legend({'$10^{-10}$','$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-3}$','$0.005$'})
export_fig('20140601-Ohmic-S-epsilon-alpha','-transparent','-pdf','-png','-painters')
%% Plot the chain occupation for increasing alpha 3D
select = PlotData(:,11)==-1e-10;  % epsilon
IndexedMatrix = [PlotData(select,2) BosonChainOccupation(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);
figure(5);
surf(1:49,IndexedMatrix(:,1),IndexedMatrix(:,2:end));
xlabel('Site $k$');
ylabel('$\alpha$')
zlabel('$<n>(k)$')
%% Plot the chain occupation for increasing alpha 2D
defEps = 1e-10;
defAlpha = [4e-1, 5e-1, 89e-2, 99e-2];
colors=lines(length(defAlpha));
figure(19);
hold all
select = PlotData(:,11)==-defEps;  % epsilon
IndexedMatrix = [PlotData(select,2) BosonChainOccupation(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);
for a = defAlpha
    plot(IndexedMatrix(IndexedMatrix(:,1)==a,2:end)','Marker','*','LineStyle','none');
end
xlabel('$\alpha$')
ylabel('$<n>(k)$')

%% Plot error in <a_i a_i+1> vs delta, alpha
figure(6)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,9),10,PlotData(:,9),'fill')
title('$<a_i^\dagger a_{i+1}^\dagger> - (f_i f_{i+1}-<a_i^\dagger><a_{i+1}^\dagger>)$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
%% Plot Max(<n>) vs delta, alpha
figure(7)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,10),10,PlotData(:,10),'fill')
title('$Max(<n>)$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$Max(<n>)$');

%% MLSBM from here:
%  PlotData: 1 = Lambda; 2 = z; 3 = participation; 4 = Max(<n>); 5 = E0; 6 = period;
%            7 = etaFactor;     8 = loop;          9 = tunnelE ; 10 = log(E0-Emin);
%           11 = Disorder SDV; 12 = Entropy;

%% Plot Participation, Max <n>, E0 vs \Lambda, z    (scatter)
figure(1)
scatter3(PlotData(:,1),PlotData(:,2),PlotData(:,3),10,PlotData(:,3),'fill')
title('$Max(<n>)$ vs $\Lambda$ and z');
xlabel('$\Lambda$');
ylabel('$z$');
zlabel('$Max(<n>)$');
%% Plot Max <n>, E0, vs. p, eta     (scatter)
% intuitive scatter plot
% no save
figure(2); rotate3d off
j = 9;
selectObs = {3,'$P$'; 4,'$Max(<n>)$'; 5,'$E_0/eV$'; 8,'Loop'; 9,'$<H_{s,int}>/eV$';...
    10,'$log_{10}(E_0-E_{min})$'; 11, '$\sigma_{\Delta E}$'; 12, 'S'};
PlotDat=PlotData(:,j);
scatter3(PlotData(:,6),PlotData(:,7),PlotDat,50,PlotDat,'fill')
xlabel('$p$');
ylabel('$\eta$')
zlabel(selectObs{cell2mat(selectObs(:,1)) == j,2});
formatPlot(2);
rotate3d on

%% Plot Max <n>, E0, vs. eta, p     (2D)
j = 12;
selectObs = {3,'$P$'; 4,'$Max(<n>)$'; 5,'$E_0/eV$'; 8,'Loop'; 9,'$<H_{s,int}>/eV$';...
    10,'$log_{10}(E_0-E_{min})$'; 11, '$\sigma_{\Delta E}$'; 12, 'S'};
criteria = ones(size(PlotData,1),1);
if hasDisorder
    if doAverage
        PlotData = PlotDataAvg;
    else
        PlotData = PlotDataAll;
    end
    defSDV = unique(PlotData(:,11))';
    criteria = criteria & round(real(PlotData(:,11))) == round(real(defSDV(2)));
end
defPeriods = unique(PlotData(:,6))';     % obtains all existing values
colors=hsv(length(defPeriods));
figure(3)
hold on
i=1;
legLabels = cell(1);
for periods = defPeriods
    select = PlotData(:,6) == periods & PlotData(:,7)>=0 & criteria;     % period == 4 or 16
    if sum(select) == 0
        continue
    end
    Indexed = [PlotData(select,7) real(PlotData(select,j))];
    fprintf('Sum of Imaginary values: %.5e\n',sum(imag(PlotData(select,j))));
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','Color',colors(i,:));
    legLabels{i} = sprintf('%d',periods);
    i=i+1;
end

xlabel('$\eta$');
ylabel(selectObs{cell2mat(selectObs(:,1)) == j,2});
formatPlot(3);

legend(legLabels,'Location','SouthWest')


%% Plot Participation vs. Period
figure(4)
%select = PlotData(:,6)==4;          % period == 4
select = PlotData(:,7)==1;          % etaFactor == 4
plot(PlotData(select,6),PlotData(select,3))
xlabel('$p$');
ylabel('$P$');
%% Plot Participation vs. eta, p    2D  -- pdf, png export!
defPeriods = unique(PlotData(PlotData(:,6)>0,6))';     % obtains all existing values
colors=hsv(length(defPeriods));
figure(4)
hold all
i=1;
legLabels = cell(1);
for periods = defPeriods
    select = PlotData(:,6)==periods & PlotData(:,7)<4.4;          % period == 4 or 16
    Indexed = [PlotData(select,7) PlotData(select,3)];
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','Color',colors(i,:));
    legLabels{i} = sprintf('%d',periods);
    i=i+1;
end
xlabel('$\eta$');
ylabel('$P$');
formatPlot(4);
legend(legLabels,'Location','SouthWest')
line([0 5],[16 16],'LineWidth',1,'Color','black');
%% Plot Tunneling Energy vs. eta, p     2D
defPeriods = unique(PlotData(:,6))';     % obtains all existing values
colors=hsv(length(defPeriods));
figure(5)
hold all
i=1;
legLabels = cell(1);
for periods = defPeriods
    select = PlotData(:,6)==periods & PlotData(:,7)>=0;          % period == 4 or 16
    Indexed = [PlotData(select,7) PlotData(select,9)];
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','Color',colors(i,:));
    legLabels{i} = sprintf('%d',periods);
    i=i+1;
end
xlabel('$\eta$');
ylabel('$E_{tunnel}/eV$');
formatPlot(5);
legend(legLabels,'Location','NorthWest')
%% Plot the wavefunction for different strengths
% modify the required period

period = 9;        % change to select wavefunction
% variable to add criteria for select
criteria = ones(size(PlotData,1),1);
if hasDisorder
    if doAverage
        PlotData = PlotDataAvg; PPCWavefunction = PPCWavefunctionAvg;
    else
        PlotData = PlotDataAll; PPCWavefunction = PPCWavefunctionAll;
    end
    defSDV = unique(PlotData(:,11))';
    criteria = criteria & round(real(PlotData(:,11))) == round(real(defSDV(2)));
end

select =  PlotData(:,6)==period & criteria;
IndexedMatrix = [PlotData(select,7) PPCWavefunction(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);                      % sort with respect to y axis
if hasDisorder && ~doAverage
    IndexedMatrix(:,1) = 1:sum(select);                         % plot against increment
end
figure(6); rotate3d off; clf;
surf(1:size(IndexedMatrix,2)-1,IndexedMatrix(:,1),IndexedMatrix(:,2:end));
%shading interp
xlabel('Site $k$');
ylabel('$\eta$')
zlabel('$|\Psi(k)|^2$')
formatPlot(6)
axis tight
rotate3d on

%% Plot the wavefunction for different periods --- can be removed?
select = PlotData(:,7)==3; %& PlotData(:,7) >=0;
IndexedMatrix = [PlotData(select,6) PPCWavefunction(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);
figure(4);
rotate3d off
surf(1:16,IndexedMatrix(:,1),IndexedMatrix(:,2:end));
xlabel('Site $k$');
ylabel('$p$')
zlabel('$|\Psi(k)|^2$')
formatPlot(4)
rotate3d on

%% Plot the Boson Occupation for different strengths
% modify the required period

period = 8;        % change to select wavefunction
% variable to add criteria for select
criteria = ones(size(PlotData,1),1);
if hasDisorder
    if doAverage
        PlotData = PlotDataAvg; BosonChainOccupation = BosonChainOccupationAvg;
    else
        PlotData = PlotDataAll; BosonChainOccupation = BosonChainOccupationAll;
    end
    defSDV = unique(PlotData(:,11))';
    criteria = criteria & round(real(PlotData(:,11))) == round(real(defSDV(2)));
end

select =  PlotData(:,6)==period & criteria;
IndexedMatrix = [PlotData(select,7) BosonChainOccupation(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);                      % sort with respect to y axis
if hasDisorder && ~doAverage
    IndexedMatrix(:,1) = 1:sum(select);                         % plot against increment
end
figure(6); rotate3d off; clf;
surf(1:size(IndexedMatrix,2)-1,IndexedMatrix(:,1),real(IndexedMatrix(:,2:end)));
%shading interp
set(gca,'View',[47 47]);
xlabel('Site $k$');
ylabel('$\eta$')
zlabel('$|\Psi(k)|^2$')
formatPlot(6)
axis tight;
rotate3d on

%% Scatter interpolation
%% using griddata
figure
[xi, yi] = meshgrid(1:0.5:29, -1:0.1:5);
x = PlotData(:,6);
y = PlotData(:,7);
z = PlotData(:,3);
zi = griddata(x,y,z, xi,yi);
surf(xi,yi,zi);
shading interp
xlabel('$p$'), ylabel('$\eta$'), zlabel('$P$')
rotate3d on
%% using Interpolant class
figure
F = scatteredInterpolant([x y],z)
[Xq,Yq] = meshgrid(2:0.5:16, 0:0.1:5);
F.Method = 'natural';
Vq = F(Xq,Yq);
surf(Xq,Yq,Vq);
xlabel('p'), ylabel('$\eta$');
zlabel('P');
shading interp
formatPlot(1)
rotate3d on
end

function ExtractOneData(path)
% Extracts Data and saves it to file
% Use as:
% ExtractOneData(dirs{i});
global loadDataFromFile
global saveDataToFile
global saveDataToFileError
try
	load([path,loadDataFromFile]);
catch
	fprintf('%s : file does not exist\n',path);
	return;
end
E = struct();
if isfield(results,'time');
	E.timeMins = results.time / 60;
else
	E.aborted = 1;
    disp(['aborted: ',path]);
    return;
end
E.loops = para.loop;
%% Calculate Shift analytically for iSBM:
    t = para.t; e = para.epsilon;
    A = gallery('tridiag',t(2:end),e(1:end),t(2:end));      % creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
    B = zeros(para.L-1,1);
    if isfield(results,'spin')                              % if spin was calculated in the end
        B(1) = -t(1)*(results.spin.sz);
    else
        B(1) = -t(1)*1;
    end
    % E.ShiftCalc = results.shift{1};                        % this is analytical solution for spin.
    E.ShiftCalc = [0; A\B.*sqrt(2)]';
%% Compare analytical and simulated iSBM shift
    E.ShiftSim = para.shift;
    E.ShiftDeviation = E.ShiftCalc./E.ShiftSim-1;
    E.ShiftMaxError = max(abs(E.ShiftDeviation));           % Max absolute Error
    E.ShiftMeanError = mean(abs(E.ShiftDeviation(2:end)));  % Mean absolute Error
    E.ShiftMax = max(abs(E.ShiftCalc));                     % Max Shift
%% Compare Occupation numbers
    E.NCalc = (E.ShiftCalc.*E.ShiftCalc./2);
    try
        E.NSim = results.nx;
        E.converged = 1;
    catch
        E.converged = 0;
        fprintf('%s is not converged\n',para.filename(1:58));
        return;
        if max(para.dk) < 1000
            E.NSim = calbosonocc_SBM1(mps,Vmat,para,results);   % still calculate from approximate state, takes quite long in bad cases!
        else
            fprintf('%s has dk = %u\n',para.filename(1:58),max(para.dk));
            return;
        end
    end
    E.NDeviation = E.NCalc./E.NSim-1;
    E.NMaxError = max(abs(E.NDeviation));
    E.NMeanError = mean(abs(E.NDeviation(2:end)));
    E.NMaxCalc = max(abs(E.NCalc));
    E.NMaxSim = max(abs(E.NSim));
%% Save spin <sx> and calculate entanglement entropy (see Le Hur 2008) for Plot against alpha
% try because not applicable to MLSBM
    if E.converged
        try
            E.Sx  = results.spin.sx;
            E.Sz  = results.spin.sz;
            pUp   = (1+sqrt(E.Sx^2+E.Sz^2))/2;  % sy = 0
            pDown = (1-sqrt(E.Sx^2+E.Sz^2))/2;
            E.entropy = -pUp*log2(pUp)-pDown*log2(pDown);
        catch
        end
    else
        try
            E.spin = calspin(mps,Vmat,para,results);
            E.Sx = E.spin.sx;
        catch
        end
    end
%% Calculate <a_i a_i+1>, call aaCreationCorrelator
    if max(para.dk) < 1000
        E.aaCreationCorrelator = calboson2siteCreationCorrelator_SBM1(mps,Vmat,para);   % calculated directly
        E.aCreationCorrelator = calboson1siteCreationCorrelator_SBM1(mps,Vmat,para);    % <a^+_i>
        E.aaCreationCorrelator2 = para.shift(1:end-1).*para.shift(2:end) - E.aCreationCorrelator(1:end-1).*E.aCreationCorrelator(2:end);            % calculated as f_i f_i+1 - <a^+_i><a^+_i+1>
    else
        fprintf('%s has dk = %u\n',para.filename(1:58),max(para.dk));
        return;
    end

%% Results specifically for MLSBM
    if strcmp(para.model,'MLSpinBoson')
        try
            E.participation = results.participation;
        catch
            E.participation = calParticipation(calReducedDensity(mps,Vmat,para,1));
        end
		E.EvaluesLog = cell2mat(results.EvaluesLog);
		E.EnergyConvergence = E.EvaluesLog-min(E.EvaluesLog);					% Plot against log-axis therefore -min
        E.E0 = results.E;
        E.period = para.MLSB_p;
    end

%% Save to file
extracted = E;
if E.converged
	save([path,saveDataToFile],'extracted','para','results');
else
	save([path,saveDataToFileError],'extracted','para','results');
end
end

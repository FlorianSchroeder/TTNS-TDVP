%test decompose function:

para = struct();
para.SVDmethod = 'qr';
para.parity = 'n';
%% setup time arrays
tAQR = zeros(20,20);
tAQRave = tAQR;
tASVD = zeros(20,20);
tASVDave = tASVD;
%%
runs = 40;
for k = 1:1:runs
for j = 2:2:40       % d_opt
for i = 2:2:40        %  D or d_k
    A = randn(i,i,j);
    tStart = tic;
    [Q,V] = decompose(A,2,para);
    tAQR(i/2,j/2) = toc(tStart);
    clear('Q','V')
    tStart = tic;
    [Q,V,vNE,sv] = decompose(A,2,para);
    tASVD(i/2,j/2) = toc(tStart);
    clear('Q','V')
end
end
%
tAQRave = tAQRave+tAQR;
tASVDave = tASVDave+tASVD;
end
tAQRave = tAQRave./runs;
tASVDave = tASVDave./runs;
%%
tAQRave = medfilt2(tAQRave);
tASVDave = mediflt2(tASVDave);
%% Plot Data for V
saveto = '..\Presentations\20140206 - Summary\';
surf(tAQRave)
hold on
surf(tASVDave)
%shading interp
%set(gca,'View',[0 90]);
title('QR vs SVD of V comparison including reshape')
ylabel('$\frac{d_{k}}{2}$')
xlabel('$\frac{d_{opt}}{2}$')
%%
export_fig(sprintf('%sQR vs SVD for V',saveto),'-transparent','-png','-eps')

%% Plot Data for A
saveto = '..\Presentations\20140206 - Summary\';
surf(tAQRave)
hold on
surf(tASVDave)
set(gca,'View',[50 30]);
title('QR vs SVD of A comparison including reshape')
ylabel('$\frac{D}{2}$')
xlabel('$\frac{d_{opt}}{2}$')
%% print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
export_fig(sprintf('%sQR vs SVD for A',saveto),'-transparent','-png','-eps')
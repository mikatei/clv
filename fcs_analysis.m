%% Attune FCS file reader
close all
clear all
clc

Ha = @(x,n,K) x.^n./(x.^n+K^n);

ec = [1,1,1,1,20,20,20,20,1,1,0,0]*0.01;
bs = [1,1,20,20,1,1,20,20,0,0,1,1]*0.01;
cumate = [1 2/3 (2/3)^2 (2/3)^3 (2/3)^4 (2/3)^5 (2/3)^6 0]*100; % cumate (uM)
Mv = Ha(cumate,4,35);

%% Column numbers for corresponding traits
%FCS 3.1 analysis
iFSCA = 3;
iFSCH = 13;
iSSCA = 4;
iSSCH = 14;
iGFPA = 5;
iGFPH = 15;
iRFPA = 6;
iRFPH = 16;
fbins = [0 5];
sbins = [0 5];
gbins = [0 5];
rbins = [0 5];
gfbins = [-5 5];
rfbins = [-5 5];

%% Load
FSC_data_all = nan(96,3.3e5);
SSC_data_all = nan(96,3.3e5);
GFP_data_all = nan(96,3.3e5);
RFP_data_all = nan(96,3.3e5);
VL1_data_all = nan(96,3.3e5);
load plate1
for i = 1:4
    for j = 1:8
        ind = 12*(j-1)+2*(i-1)+1;
        FSC_data_all(ind,1:length(FSC_data(4*(j-1)+i,:))) = FSC_data(4*(j-1)+i,:);
        SSC_data_all(ind,1:length(SSC_data(4*(j-1)+i,:))) = SSC_data(4*(j-1)+i,:);
        GFP_data_all(ind,1:length(GFP_data(4*(j-1)+i,:))) = GFP_data(4*(j-1)+i,:);
        RFP_data_all(ind,1:length(RFP_data(4*(j-1)+i,:))) = RFP_data(4*(j-1)+i,:);
        VL1_data_all(ind,1:length(VL1_data(4*(j-1)+i,:))) = VL1_data(4*(j-1)+i,:);
    end
end
load plate2
for i = 1:4
    for j = 1:8
        ind = 12*(j-1)+2*i;
        FSC_data_all(ind,1:length(FSC_data(4*(j-1)+i,:))) = FSC_data(4*(j-1)+i,:);
        SSC_data_all(ind,1:length(SSC_data(4*(j-1)+i,:))) = SSC_data(4*(j-1)+i,:);
        GFP_data_all(ind,1:length(GFP_data(4*(j-1)+i,:))) = GFP_data(4*(j-1)+i,:);
        RFP_data_all(ind,1:length(RFP_data(4*(j-1)+i,:))) = RFP_data(4*(j-1)+i,:);
        VL1_data_all(ind,1:length(VL1_data(4*(j-1)+i,:))) = VL1_data(4*(j-1)+i,:);
    end
end
load plate3
for i = 1:4
    for j = 1:8
        ind = 12*(j-1)+8+i;
        FSC_data_all(ind,1:length(FSC_data(4*(j-1)+i,:))) = FSC_data(4*(j-1)+i,:);
        SSC_data_all(ind,1:length(SSC_data(4*(j-1)+i,:))) = SSC_data(4*(j-1)+i,:);
        GFP_data_all(ind,1:length(GFP_data(4*(j-1)+i,:))) = GFP_data(4*(j-1)+i,:);
        RFP_data_all(ind,1:length(RFP_data(4*(j-1)+i,:))) = RFP_data(4*(j-1)+i,:);
        VL1_data_all(ind,1:length(VL1_data(4*(j-1)+i,:))) = VL1_data(4*(j-1)+i,:);
    end
end
FSC_data = FSC_data_all;
SSC_data = SSC_data_all;
GFP_data = GFP_data_all;
RFP_data = RFP_data_all;
VL1_data = VL1_data_all;
L = 96;

%% data processing
T = nan(8,12); % total cells
G = T; % GFP+ cells
for i = 1:96
    SSC = SSC_data(i,:);
    VL1 = VL1_data(i,:);
    good = and(SSC>10^3.3,VL1>10^2);
    GFP = GFP_data(i,:);
    if mod(i,12) == 0
        T(ceil(i/12),12) = sum(good);
        G(ceil(i/12),12) = sum(and(GFP>10^2.5,good))./sum(and(GFP>0,good));
    else
        T(ceil(i/12),mod(i,12))= sum(good);
        G(ceil(i/12),mod(i,12)) = sum(and(GFP>10^2.5,good))./sum(and(GFP>0,good));
    end
end

%% plot
figure('Position',[0 0 1200 800],'Name','fcs result')
for j = 1:12
    subplot(3,4,j)
    hold all
    for i = 1:8
        plot(8-i,G(i,j),'bo','linewidth',3)
    end
    axis([1 8 0 1])
    axis square
    ylabel('vioABE fraction')
    title(sprintf('Ec %.2f, Bs %.2f',ec(j),bs(j)))
    set(gca,'fontsize',18)
    set(gcf,'PaperPositionMode','auto')
end

%% Mean
AM = nan(8,6);
mgfp = nan(8,1);
sgfp = mgfp;
figure('Position',[0 0 1200 800],'Name','mean fcs result')
for j = 1:6
    subplot(2,3,j)
    hold all
    for i = 1:8
        mgfp(i) = nanmean([G(i,2*j-1),G(i,2*j)]);
        sgfp(i) = nanstd([G(i,2*j-1),G(i,2*j)]);
        errorbar(Mv(i),mgfp(i),sgfp(i),'bo','linewidth',2)
    end
    axis([-0.1 1.1 0 1.1])
    axis square
    title(sprintf('Ec %.2f, Bs %.2f',ec(2*j),bs(2*j)))
    ylabel('vioABE fraction')
    set(gca,'fontsize',18)
    set(gcf,'PaperPositionMode','auto')
end
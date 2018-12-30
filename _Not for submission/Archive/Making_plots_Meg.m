% Plots for paper with Sarah
% last edited by Megan, March 21, 2018
clear all
close all
load('Mombasa_TandP.mat','P0','Pij','Tij','T0')

% yearly values at Mombasa
for y = 1:length(P0)/12
    CRU.T(y) = mean(T0(12*(y-1)+1:12*y));
    CRU.P(y) = mean(P0(12*(y-1)+1:12*y));   
end

for y = 1:size(Pij,1)/12
    GCM.T(y,:) = mean(Tij(12*(y-1)+1:12*y,:),1);
    GCM.P(y,:) = mean(Pij(12*(y-1)+1:12*y,:),1);
end
clearvars -except CRU GCM
% filtering out high freq noise (10 year moving window)
windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
GCM.Tsmooth = filter(b,a, GCM.T); %xfilt is the low frequency signal
GCM.Psmooth = filter(b,a, GCM.P); %xfilt is the low frequency signal
CRU.Tsmooth = filter(b,a, CRU.T); %xfilt is the low frequency signal
CRU.Psmooth = filter(b,a, CRU.P); %xfilt is the low frequency signal

for y = 1:length(CRU.T)-39
    CRU.deltaT(y) = mean(CRU.T(y+20:y+39))-mean(CRU.T(y:y+19));
    CRU.deltaP(y) = (mean(CRU.P(y+20:y+39))-mean(CRU.P(y:y+19)))/mean(CRU.P(y:y+19)); 
end

for y = 1:size(GCM.T,1)-39
    GCM.deltaT(y,:) = mean(GCM.T(y+20:y+39,:),1)-mean(GCM.T(y:y+19,:),1);
    GCM.deltaP(y,1:21) = (mean(GCM.P(y+20:y+39,:),1)-mean(GCM.P(y:y+19,:),1))./mean(GCM.P(y:y+19,:),1); 
end

%% Plots of basic climate projections
fig1 = figure; 
subplot(1,2,1) 
plot([1951:1:2100],GCM.Tsmooth(51:end,:),'LineWidth',1);
hold on;
plot([1951:1:2012],CRU.Tsmooth(51:end),'k','LineWidth',3);title('Temperature')
set(gca,'FontSize', 6); xlabel('Year'); ylabel('T')

subplot(1,2,2) 
plot([1951:1:2100],GCM.Psmooth(51:end,:),'LineWidth',1);title('Precipitation')
hold on;
plot([1951:1:2012],CRU.Psmooth(51:end),'k','LineWidth',3);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('P')
figure_width = 7; % in inches
figure_height = 2.50; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'TandP';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)
%%
fig1 = figure; 
subplot(1,2,1) 
plot([1921:1:2081],GCM.deltaT(1:end,:),'LineWidth',1);
hold on;
plot([1921:1:1993],CRU.deltaT(1:end),'k','LineWidth',3);
title('\Delta T'); xlim([1921,2081]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta T over 20-years');


subplot(1,2,2) 
plot([1921:1:2081],GCM.deltaP(1:end,:),'LineWidth',1);title('% \Delta P')
hold on;
plot([1921:1:1993],CRU.deltaP(1:end),'k','LineWidth',3);xlim([1921,2081]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta P over 20-years');

figure_width = 7; % in inches
figure_height = 2.50; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'DeltaTandP';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)

%% IPCC method
GCM.Tchange = GCM.T-repmat(mean(GCM.T(86:105,:),1),200,1);
GCM.Pchange = (GCM.P-repmat(mean(GCM.P(86:105,:),1),200,1));
for y = 1:size(GCM.Tchange,1)
    IPCC.meanT(y) = double(mean(GCM.Tchange(y,:)));
    IPCC.stdT(y) = double(std(GCM.Tchange(y,:)));
    IPCC.meanP(y) = double(mean(GCM.Pchange(y,:)));
    IPCC.stdP(y) = double(std(GCM.Pchange(y,:)));
end

IPCC.meanT = filter(b,a, IPCC.meanT);
IPCC.meanP = filter(b,a, IPCC.meanP);

IPCC.stdT = filter(b,a, IPCC.stdT);
IPCC.stdP = filter(b,a, IPCC.stdP);
%%
fig1 = figure; 
subplot(1,2,1) 
boundedline([2005:1:2100]',IPCC.meanT(105:200)',[1.64*IPCC.stdT(105:200)' 1.64*IPCC.stdT(105:200)'],'alpha','cmap',[0.7,0,0]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanT(11:105)',[1.64*IPCC.stdT(11:105)' 1.64*IPCC.stdT(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on; 
plot([1911:1:2005]',CRU.Tsmooth(11:105)-mean(CRU.T(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); 
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta T (relative to 1986-2005)');

subplot(1,2,2) 
boundedline([2005:1:2100]',IPCC.meanP(105:200)',[1.64*IPCC.stdP(105:200)' 1.64*IPCC.stdP(105:200)'],'alpha','cmap',[0,0,0.7]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanP(11:105)',[1.64*IPCC.stdP(11:105)' 1.64*IPCC.stdP(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on; 
plot([1911:1:2005]',CRU.Psmooth(11:105)-mean(CRU.P(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); ylabel('\Delta P (relative to 1986-2005)');xlabel('Year')
set(gca,'FontSize', 6)
figure_width = 7; % in inches
figure_height = 3; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'IPCCbounded';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)

%%
fig1 = figure; 
subplot(1,2,1) 
boundedline([2005:1:2100]',IPCC.meanT(105:200)',[1.64*IPCC.stdT(105:200)' 1.64*IPCC.stdT(105:200)'],'alpha','cmap',[0.7,0,0]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanT(11:105)',[1.64*IPCC.stdT(11:105)' 1.64*IPCC.stdT(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on;
plot([1911:1:2100],GCM.Tsmooth(11:end,:)-repmat(mean(GCM.T(86:105,:),1),190,1),'LineWidth',1);
hold on; 
plot([1911:1:2005]',CRU.Tsmooth(11:105)-mean(CRU.T(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); 
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta T')

subplot(1,2,2) 
boundedline([2005:1:2100]',IPCC.meanP(105:200)',[1.64*IPCC.stdP(105:200)' 1.64*IPCC.stdP(105:200)'],'alpha','cmap',[0,0,0.7]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanP(11:105)',[1.64*IPCC.stdP(11:105)' 1.64*IPCC.stdP(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on; 
plot([1911:1:2100],GCM.Psmooth(11:end,:)-repmat(mean(GCM.P(86:105,:),1),190,1),'LineWidth',1);
hold on; 
plot([1911:1:2005]',CRU.Psmooth(11:105)-mean(CRU.P(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); 
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta P')

figure_width = 7; % in inches
figure_height = 3; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'IPCCbounded_withGCMlines';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)

%% Plots of learning BMA approach.
load('BMA_results_T0istoday21-Mar-2018 12:10:37.mat')
Nsamples = 10^4;
Ntsteps = 5;
GN = randi(1000, Nsamples, Ntsteps,2);
DT = [squeeze(MUT(GN(:,1,1))), squeeze(NUT(GN(:,:,1)))];
DT = cumsum(DT,2); % this is the cumulative 

GN = randi(1000, Nsamples, Ntsteps,3);
DP = [squeeze(MUP(GN(:,1,1))), squeeze(NUP(GN(:,:,1)))];
DP = ones(size(DP))+DP;
percDP = DP(:,1);
for year = 2:6
    percDP(:,year) = DP(:,year).*percDP(:,year-1);
end

DT = sort(DT,1);
DP = sort(percDP,1);

medInd = Nsamples/2;
perc5 = floor(5*Nsamples/100);
perc95 = floor(95*Nsamples/100);

fig1 = figure; 
subplot(1,2,1) 
boundedline([1980:20:2080]',CRU.Tsmooth(80)-DT(medInd,1)+DT(medInd,:)',[(DT(medInd,:)-DT(perc5,:))' (DT(perc95,:)-DT(medInd,:))'],'alpha','cmap',[0.7,0,0]); 
plot([1911:1:2012]',CRU.Tsmooth(11:112),'k','LineWidth',2); xlim([1920,2080]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('T')

subplot(1,2,2) 
boundedline([1980:20:2080]',CRU.Psmooth(80)*DP(medInd,:)',CRU.Psmooth(90)*[(DP(medInd,:)-DP(perc5,:))' (DP(perc95,:)-DP(medInd,:))'],'alpha','cmap',[0,0,0.7]); 
plot([1911:1:2012]',CRU.Psmooth(11:112),'k','LineWidth',2); xlim([1920,2080]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('P')

figure_width = 7; % in inches
figure_height = 2.5; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'BMA_nolearning';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)

%% BMA with GCMs overtop

GCM.Tchange2 = GCM.T-repmat(mean(GCM.T(61:80,:),1),200,1);
GCM.Pchange2 = (GCM.P-repmat(mean(GCM.P(61:80,:),1),200,1));
GCM.Tsmooth2 = filter(b,a, GCM.Tchange2); %xfilt is the low frequency signal
GCM.Psmooth2 = filter(b,a, GCM.Pchange2); %xfilt is the low frequency signal


fig1 = figure; 
subplot(1,2,1) 
plot([1911:1:2100],CRU.Tsmooth(80)+GCM.Tsmooth2(11:end,:),'LineWidth',1);
hold on
boundedline([1980:20:2080]',CRU.Tsmooth(80)-DT(medInd,1)+DT(medInd,:)',[(DT(medInd,:)-DT(perc5,:))' (DT(perc95,:)-DT(medInd,:))'],'cmap',[0.7,0,0]); 
plot([1911:1:2012]',CRU.Tsmooth(11:112),'k','LineWidth',2); xlim([1920,2080]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('T')

subplot(1,2,2) 
plot([1911:1:2100],CRU.Psmooth(80)+GCM.Psmooth2(11:end,:),'LineWidth',1);
hold on
boundedline([1980:20:2080]',CRU.Psmooth(80)*DP(medInd,:)',CRU.Psmooth(90)*[(DP(medInd,:)-DP(perc5,:))' (DP(perc95,:)-DP(medInd,:))'],'cmap',[0,0,0.7]); 
plot([1911:1:2012]',CRU.Psmooth(11:112),'k','LineWidth',2); xlim([1920,2080]);
set(gca,'FontSize', 6); xlabel('Year'); ylabel('P')

figure_width = 7; % in inches
figure_height = 2.5; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'BMA_nolearningwithGCMstrings';
figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)


%% Just looking at Deltas and not absolute values.


%%
load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')


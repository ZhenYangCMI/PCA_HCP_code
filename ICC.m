clear
clc
close all

numSeed=4;
numSub=22;
numCommonState=5;
numStateSes1=5;
numStateSes2=6;

dataLength='all_10min';
load('MyColormaps2','mycmap2')
load('MyColormaps','mycmap')


resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session1/'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/',dataLength,filesep, 'session2/'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/',dataLength, filesep, 'correlTwoSessions/'];

% compute ICC and perform t-test on totoal number of transition
tmp1=load([resultDir1, 'transitions_session1.mat']);
transit1=tmp1.transitions;
tmp2=load([resultDir2, 'transitions_session2.mat']);
transit2=tmp2.transitions;

figure(1)
imagesc(transit1)
colorbar
title('Number of transitions ')
xlabel('Seeds')
ylabel('Subjects')
caxis([0 15])
set(gca,'xTick',1:4);
xlim([0.5 4.5]);
set(figure(1),'Colormap',mycmap2)
saveas(figure(1),[figDir,'TotNumTransitSession1.png'])

figure(2)
imagesc(transit2)
colorbar
title('Number of transitions ')
xlabel('Seeds')
ylabel('Subjects')
caxis([0 15])
set(gca,'xTick',1:4);
xlim([0.5 4.5]);
set(figure(2),'Colormap',mycmap2)
saveas(figure(2),[figDir,'TotNumTransitSession2.png'])

numPC=10;


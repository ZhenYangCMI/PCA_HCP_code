% phase randomization

clear
clc
close all

sessionList={'rfMRI_REST1_RL','rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};
subList=load('/home/data/Projects/HCP/subListDynamic.txt');
numROI=88;
numPC=10;
numSession=length(sessionList);

numVol=1200;
winSize=61;
step=3;
numWinPerSub=floor((numVol-winSize)/step)+1;
numSub=length(subList);
numWinAllSub=numWinPerSub*numSub;

resultDir=['/home/data/Projects/HCP/results/phaseRand/'];
figDir = ['/home/data/Projects/HCP/figs/phaseRand/'];

allSesFeature=[];
%for j=1:length(sessionList)
for j=1
    session=char(sessionList{j})
    
    for i=1:length(subList)
        sub=num2str(subList(i));
        dataDir=['/home/data/Projects/HCP/data/', session, '/', sub, '/'];
        loadFile=load([dataDir, '/ROISignals_HCP.mat']);
        
        ROIsignal=loadFile.ROISignals;
        for k=1:numROI
            [ out_ts_arr ] = autocorr_resample_ycg( squeeze(ROIsignal(:,k)), 1 );
            ROISignalRand(:,k)=out_ts_arr;
        end
        save([dataDir, 'ROISignals_HCP_rand.mat'], 'ROISignalRand')
    end
end

for j=1:length(sessionList)
    
    session=char(sessionList{j})
    
    ROIoutputDir=['/home/data/Projects/HCP/data/', session, '/'];
    [ fullCorrWin zFullCorrWin ] = FCwinCreation( subList, ROIoutputDir, winSize, numWinPerSub, step);
    [featureWin, zFeatureWin] = featureExtract(zFullCorrWin);
    [ norm1FeatureWin norm2FeatureWin ] = standardizeFeature( subList, featureWin, numWinPerSub );
    
end
allSesNormFeature(1+numWinPerSes*(j-1):numWinPerSes*j, :);

[COEFF,SCORE,latent, tsquare] = princomp(allSesNormFeature);

a=COEFF(:,1:10);

allSesW=zeros(numSub*4,numWinPerSub, numPC);
for i=1:numSub*4
    
    eachSub= allSesNormFeature(1+(i-1)*numWinPerSub: i*numWinPerSub, :)';
    W=a'*eachSub;
    for j=1:numPC
        allSesW(i,:,j)=W(j,:);
    end
    
end

allRandW(:,:,:,t)=allSesW;


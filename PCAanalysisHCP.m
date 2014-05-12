% PCA analysis
clear
clc

% Global variable
sessionList={'rfMRI_REST1_RL','rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};
%sessionList={'rfMRI_REST2_RL'};
subList=load('/home/data/Projects/HCP/sublist95sub.txt');

maskDir='/home/data/Projects/HCP/mask/';
standardBrainMask=[maskDir,'BrainMask_05_61x73x61.img'];
ROImask=[maskDir,'AAL_61x73x61_YCG_reduced.nii'];

numVol=1200;
winSize=61;
step=3;
numWinPerSub=floor((numVol-winSize)/step)+1;
numSub=length(subList);
numWinAllSub=numWinPerSub*numSub;

norm='norm2';

%% 0. mask creation AAL mask with ROI 1 to 88

% mask=[mask, '/AAL_61x73x61_YCG.nii']
% [MaskData,MaskVox,MaskHead]=rest_readfile(mask);
%
% MaskData(find(MaskData>88))=0;
% reducedMask=MaskData;
%
% fileName='/home/data/Projects/microstate/DPARSF_preprocessed/fullBrainFC/mask/AAL_61x73x61_YCG_reduced.nii';
% rest_WriteNiftiImage(reducedMask,MaskHead,fileName);


%% 1. Extract the TS from ROIs

% for j=1:length(sessionList)
%     session=char(sessionList{j})
%
%     ROIoutputDir=['/home/data/Projects/HCP/data/', session, '/'];
%
%     [ counter ] = extractTC(subList, ROIoutputDir, standardBrainMask, ROImask);
% end
%
% disp ('ROI time series extraction done!')

%% 2. FC window creation, extract the feature from the lowertriangle

% for j=1:length(sessionList)
%     
%     session=char(sessionList{j})
%     
%     ROIoutputDir=['/home/data/Projects/HCP/data/', session, '/'];
%     
%     if ~exist(['/home/data/Projects/HCP/results/', session], 'dir')
%         mkdir(['/home/data/Projects/HCP/results/', session])
%     end
%     
%     if ~exist(['/home/data/Projects/HCP/figs/', session], 'dir')
%         mkdir(['/home/data/Projects/HCP/figs/', session])
%     end
%     
%     resultDir=['/home/data/Projects/HCP/results/', session, '/'];
%     figDir = ['/home/data/Projects/HCP/figs/', session, '/'];
%     
%     % FC window creation
%     disp('Create FC windows')
%     [ fullCorrWin zFullCorrWin ] = FCwinCreation( subList, ROIoutputDir, winSize, numWinPerSub, step);
%     save([resultDir, 'fullCorrWin_', session, '.mat'], 'fullCorrWin')
%     save([resultDir, 'zFullCorrWin_', session, '.mat'], 'zFullCorrWin')
%     
%     % extract the lower triangle and z standardize each window
%     disp('Extract lower triangle feature')
%     [featureWin, zFeatureWin] = featureExtract(zFullCorrWin);
%     save([resultDir, 'feature_', session, '.mat'], 'featureWin')
%     save([resultDir, 'zNormWinFeature_', session, '.mat'], 'zFeatureWin')
%     
%     % standardize each subject's data by the global variable and row-wise demean (norm2) or not row-wise demean (norm1)
%     disp('Standardize the dyanmic FC matric')
%     [ norm1FeatureWin norm2FeatureWin ] = standardizeFeature( subList, featureWin, numWinPerSub );
%     save([resultDir,'/norm1FeatureWin_',session,'.mat'],'norm1FeatureWin')
%     save([resultDir,'/norm2FeatureWin_',session,'.mat'],'norm2FeatureWin')
%     
%     %% 3 run PCA
%     disp ('Run PCA')
%     normFeature=eval([norm, 'FeatureWin']);
%     [COEFF,SCORE,latent, tsquare] = princomp(normFeature);
%     cumVar=cumsum(latent)./sum(latent);
%     save([resultDir,'/eigenvector_',norm, '_',session,'.mat'],'COEFF')
%     save([resultDir,'/eigenvalue_',norm, '_',session,'.mat'],'latent')
%     save([resultDir,'/cumVar_', norm,  '_',session,'.mat'],'cumVar')
%     save([resultDir,'/tsquare_', norm,  '_',session,'.mat'],'tsquare')
%     
%     x=1:length(latent);
%     numEigenValue=[length(latent), 100, 20];
%     close all
%     figure(1)
%     for m=1:length(numEigenValue)
%         num=numEigenValue(m);
%         subplot(3,1,m)
%         [AX, H1, H2]=plotyy(x(1:num),latent(1:num), x(1:num), cumVar(1:num));
%         set(get(AX(1),'Ylabel'),'String','Eigen Value')
%         set(get(AX(2),'Ylabel'),'String','Cummulative Variance')
%         xlabel('Number of Eigen Values')
%         title([num2str(num),  'Eigen Values'])
%         set(H1,'LineStyle','-', 'LineWidth', 2)
%         set(H2,'LineStyle','-', 'LineWidth', 2)
%     end
%     saveas(figure(1), [figDir, 'eigenValue_', norm, '_', session, '.jpg'])
%     
% end

%% 4 run PCA on concatenated data
disp('Run PCA on concatenated data')
allSesNormFeature=[];
for j=1:length(sessionList)

    session=char(sessionList{j})

    dataDir=['/home/data/Projects/HCP/results/', session, '/'];
    loadFile=load([dataDir,'/', norm, 'FeatureWin_',session,'.mat'])
    normFeature=loadFile.([norm, 'FeatureWin']);
    numWinPerSession=size(normFeature,1);
    allSesNormFeature(1+numWinPerSession*(j-1):numWinPerSession*j, :)=normFeature;
end

session='all';
if ~exist(['/home/data/Projects/HCP/results/', session], 'dir')
        mkdir(['/home/data/Projects/HCP/results/', session])
    end

    if ~exist(['/home/data/Projects/HCP/figs/', session], 'dir')
        mkdir(['/home/data/Projects/HCP/figs/', session])
    end

resultDir=['/home/data/Projects/HCP/results/all/'];
figDir = ['/home/data/Projects/HCP/figs/all/'];


[COEFF,SCORE,latent, tsquare] = princomp(allSesNormFeature);
cumVar=cumsum(latent)./sum(latent);
save([resultDir,'/eigenvector_',norm, '_',session,'.mat'],'COEFF')
save([resultDir,'/eigenvalue_',norm, '_',session,'.mat'],'latent')
save([resultDir,'/cumVar_', norm,  '_',session,'.mat'],'cumVar')
save([resultDir,'/tsquare_', norm,  '_',session,'.mat'],'tsquare')

x=1:length(latent);
numEigenValue=[length(latent), 100, 20];
close all
figure(1)
for m=1:length(numEigenValue)
    num=numEigenValue(m);
subplot(3,1,m)
[AX, H1, H2]=plotyy(x(1:num),latent(1:num), x(1:num), cumVar(1:num));
set(get(AX(1),'Ylabel'),'String','Eigen Value')
set(get(AX(2),'Ylabel'),'String','Cummulative Variance')
xlabel('Number of Eigen Values')
title([num2str(num),  'Eigen Values'])
set(H1,'LineStyle','-', 'LineWidth', 2)
set(H2,'LineStyle','-', 'LineWidth', 2)
end
saveas(figure(1), [figDir, 'eigenValue_', norm, '_', session, '.jpg'])



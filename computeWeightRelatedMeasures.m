close all
clear
clc

sessionList={'rfMRI_REST1_RL','rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};
subList=load('/home/data/Projects/Zhen/PCA_HCP/sublist95sub.txt');

norm='norm2';
numPC=10;
numSession=length(sessionList);

allSesNormFeature=[];
for j=1:length(sessionList)
    
    session=char(sessionList{j})
    
    dataDir=['/home/data/Projects/Zhen/PCA_HCP/results/', session, '/'];
    loadFile=load([dataDir,'/', norm, 'FeatureWin_',session,'.mat'])
    normFeature=loadFile.([norm, 'FeatureWin']);
    numWinPerSession=size(normFeature,1);
    allSesNormFeature(1+numWinPerSession*(j-1):numWinPerSession*j, :)=normFeature;
end


resultDir=['/home/data/Projects/Zhen/PCA_HCP/results/all/'];
figDir = ['/home/data/Projects/Zhen/PCA_HCP/figs/poster_figs/'];

session='all';

tmp=load([resultDir,'/eigenvector_',norm, '_',session,'.mat']);
PCs=tmp.COEFF;

a=PCs(:,1:10);

% projects the eigenVector to subject dyanmic FC matrix

numVol=1200;
winSize=61;
step=3;
numWinPerSub=floor((numVol-winSize)/step)+1;
numSub=length(subList);
numWinAllSub=numWinPerSub*numSub;

allW=zeros(numSub*4,numWinPerSub, numPC);
for i=1:numSub*4
    
    % compute the weight by regression
    eachSub= allSesNormFeature(1+(i-1)*numWinPerSub: i*numWinPerSub, :)';
    W=a'*eachSub;
    
    for j=1:numPC
        
        allW(i,:,j)=W(j,:);
        % calc the % of pos weights
        prctPosW(i,j)=length(find(W(j,:)>0))/numWinPerSub*100;
        
        % calc the average magnitude of pos weights
        avgPosW(i,j)=mean(W(j,find(W(j,:)>0)));
        
        % calc the skewness of the W TS
        skewW(i,j)=skewness(W(j,:));
        
        % calc the std of thw W TS
        stdW(i,j)=std(W(j,:));
        
        % calc the transitions from - to + or + to -
        count=0;
        for k=1:size(W,2)-1
            if W(j,k)*W(j,k+1)<0
                count=count+1;
            end
        end
        transition(i,j)=count;
    end
    
end

durationPosW=prctPosW./transition;

% plot the weights for an exemplar subject
%color='r','b','g','y','k','c','m', 'r--', 'b--','g--']
color=[1 215/255 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0;205/255 92/255 92/255;85/255 107/255 47/255; 0 128/255 128/255];
%for k=1:numSub
%sub=k;  % sub <=95;
sub=47
%figure(k)

for j=1:numSession % numSession
    session=char(sessionList{j})
    figure(j)
    for i=1:10
        
        plot(allW(sub+numSub*(j-1),:,i), 'Color', color(i,:), 'LineWidth', 2)
        hold on
    end
    if mod(j,2)==1
        ses='1';
    else
        ses='2';
    end
    %title(['Session ', num2str(ceil(j/2)), ' Run', ses])
    ylim([-50 70])
    %xlabel('Time windows')
    %ylabel('Time-dependent Weights')
    set(gca, 'Linewidth',2, 'FontName','Arial', 'fontsize', 32)
    saveas(figure(j), [figDir, 'sub', num2str(sub), '_', session, '_weightsTC.jpg'])
end
close all
% end

Y=zeros(numPC,1);

for j=1:numPC
    meanW(j,1)=mean(mean(allW(:,:,j)));
    globalStdW(j,1)=std(std(allW(:,:,j)));
    pos=reshape(llW(:,:,j), [], 1);
    Y(j,1)= prctile(pos(find(pos>0)), 20);
end

for i=1:numSub*4
    for j=1:numPC
        tmpW=allW(i,:,j);
        prctSigPosW(i,j)=length(find(tmpW>Y(j)))/numWinPerSub*100;
        avgSigPosW(i,j)=mean(tmpW(find(tmpW>Y(j))));
        medianSigPosW(i,j)=median(tmpW(find(tmpW>Y(j))));
        medianPosW(i,j)=median(tmpW(find(tmpW>0)));
        flag=0;
        count=0;
        for m=1:size(tmpW,2)
            if flag==0 && tmpW(m)>=Y(j)
                count=count+1;
                flag=1;
            elseif tmpW(m)<Y(j)
                flag=0;
            end
        end
        
        timesSigPosW(i,j)=count;
        
        durationSigPosW(i,j)=prctSigPosW(i,j)/timesSigPosW(i,j);
        
    end
end

metricList={'prctPosW', 'avgPosW', 'stdW',  'skewW', 'transition', 'durationPosW'}
%metricList={'prctPosW', 'avgPosW', 'medianPosW',  'prctSigPosW', 'avgSigPosW', 'medianSigPosW', 'stdW', 'timesSigPosW', 'durationSigPosW', 'transition', 'skewW'}



% plot and evaluate the ICC of these measures
ICCREMLAllMeasure=zeros(length(metricList), numPC);
for m=1:length(metricList)
    metricName=char(metricList{m})
    metric=eval(metricList{m});
    
    figure(1)
    imagesc(metric)
    colorbar
    ylabel('sub 4 sessions')
    xlabel('PCs')
    metricName
    if strcmp(metricName, 'transition')
        caxis([0 50])
    elseif strcmp(metricName, 'avgPosW')|| strcmp(metricName, 'medianPosW')|| strcmp(metricName, 'medianSigPosW') || strcmp(metricName, 'medianPosW')
        caxis([0 10])
    elseif strcmp(metricName, 'skewW')
        caxis([-1.5 1.5])
    elseif strcmp(metricName, 'durationSigPosW')
        caxis([0 8])
    elseif strcmp(metricName, 'prctPosW') || strcmp(metricName, 'prctSigPosW')
        caxis([30 70])
    elseif strcmp(metricName, 'timesSigPosW')
        caxis([0 20])
    else
        caxis([0 10])
    end
    saveas(figure(1), [figDir, metricName, '.jpg'])
    
    
    %compute the ICC within an between session
    numDay=2;
    %ICCTypeList={'withinDay1', 'withinDay2', 'between'};
    ICCTypeList={'between'}
    numICCType=length(ICCTypeList);
    ICCREML=zeros(numPC, numICCType);
    ICCIPN=zeros(numPC,numICCType);
    pValue=zeros(numPC,numICCType);
    
    day1ses1=metric(1:numSub, :);
    day1ses2=metric(1+numSub:2*numSub, :);
    day2ses1=metric(1+2*numSub:3*numSub, :);
    day2ses2=metric(1+3*numSub:4*numSub, :);
    day1avg=(day1ses1+day1ses2)/2;
    day2avg=(day2ses1+day2ses2)/2;
    
    for j=1:numICCType
        ICCType=char(ICCTypeList{j});
        if strcmp(ICCType,'withinDay1')
            a1=day1ses1;
            a2=day1ses2;
        elseif strcmp(ICCType,'withinDay2')
            a1=day2ses1;
            a2=day2ses2;
        else
            a1=day1avg;
            a2=day2avg;
        end
        
        for i=1:numPC
            disp(['Work on PC ', num2str(i)])
            session1Data=squeeze(a1(:,i));
            session2Data=squeeze(a2(:,i));
            
            dataConcate=[session1Data;session2Data];
            time = [ones(numSub,1);2*ones(numSub,1)];
            sID=[[1:numSub]';[1:numSub]'];
            [ ICC1, idx_fail] = do_ICC(dataConcate, time, [], [], sID);
            ICCREML(i,j)=ICC1;
            
            %             data=[session1Data,session2Data];
            %             ICC2 = IPN_icc(data,1,'single');
            %             ICCIPN(i,j)=ICC2;
            %             [h,p]=ttest(session1Data, session2Data);
            %             pValue(i,j)=p;
        end
    end
    ICCREMLAllMeasure(m ,:)=ICCREML';
    %ICCIPN
    %pValue
    save([resultDir, 'ICCREMLwithinAndbetwDays_', metricName, '.mat'], 'ICCREMLAllMeasure')
    
    figure(2)
    imagesc(ICCREMLAllMeasure)
    colorbar
    %if strcmp(metricName, 'avgPosW')
    %   caxis([-1 1])
    %else
    caxis([-0.8 0.8])
    %end
    set(gca, 'Linewidth',2, 'FontName','Arial', 'fontsize', 32, 'XTick', 1:10, 'YTick', 1:5)
    saveas(figure(2), [figDir, 'ICCREMLbetwSession_allMeasure.jpg'])
end

% relate weight time course based measure to behavior

% for i=1:length(metricList)
%     metric=char(metricList{i})
%     metricData=eval(metric);
%     day1=(metricData(1:numSub,:)+metricData(1+numSub:2*numSub, :))/2;
%     day2=(metricData(1+numSub*2:3*numSub,:)+metricData(1+3*numSub:4*numSub, :))/2;
%
%     allMeasureDay1(:, 1+(i-1)*numPC:i*numPC)=day1;
%     allMeasureDay2(:, 1+(i-1)*numPC:i*numPC)=day2;
%
% end
% % save([resultDir, 'allMeasureDay1.txt'],  '-ascii', '-tabs', 'allMeasureDay1')
% % save([resultDir, 'allMeasureDay2.txt'],  '-ascii', '-tabs', 'allMeasureDay2')
% %
% % % load motion
% % tmp=load('/home/data/Projects/Zhen/PCA_HCP/data/motionRelMeanRMS.mat')
% % motion=tmp.motion;
% % motionday1=(motion(:,1)+motion(:,2))/2;
% % motionday2=(motion(:,3)+motion(:,4))/2;
% %
% % % load behavior measure
% [number,txt,raw]= xlsread('/home/data/Projects/Zhen/PCA_HCP/data/NPHCP_1_15_14.xlsx');
% NP=[number(:, 7:8)]; % flanker is dropped due to extreme high correlation with Processing speed
%
% [r1 p1]=partialcorr([allMeasureDay1(:, 1:2),allMeasureDay1(:, 5:9)], NP(1:95, :), [number(1:95, 14)]);
% [r2 p2]=partialcorr([allMeasureDay2(:, 1:2),allMeasureDay2(:, 5:9)] , NP(1:95, :), [number(1:95, 15)]);
% %
% [h p]=ttest2(allMeasureDay1(1:23,:), allMeasureDay1(24:49, :));
% [h2 p2]=ttest2(allMeasureDay2(1:23,:), allMeasureDay2(24:49, :));

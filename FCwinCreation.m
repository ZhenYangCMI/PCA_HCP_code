function [ fullCorrWin zFullCorrWin ] = FCwinCreation( subList, ROIoutputDir, winSize, numWinPerSub, step)
%This function create full connectivity windows

% create the convolved window
[finalWin]=winCreation(winSize,0);
finalWin=finalWin';
numWinAllSub=numWinPerSub*length(subList);

for i=1:length(subList)
    sub=num2str(subList(i));
    disp (['Working on sub ', sub,' ......'])
    subDir=[ROIoutputDir,sub];
% for real data    
%ROISignals = load([subDir,'/ROISignals_HCP.mat']);
    %TC=ROISignals.ROISignals;
% for surrogate data
ROISignals = load([subDir,'/ROISignals_HCP.mat']);
    TC=ROISignals.ROISignals;
    
    % apply the sliding window to the time series and segment the TS into TS windows
    asize = size(TC); % numROIs
    
    %win=zeros(winSize,asize(2),numWinAllSub);
    for q=1+numWinPerSub*(i-1):numWinPerSub*i;
        for n=((q-1)-(i-1)*numWinPerSub)*step+1:winSize+((q-1)-(i-1)*numWinPerSub)*step
            for m=1:asize(2)
                win(n-((q-1)-(i-1)*numWinPerSub)*step,m,q)=TC(n,m)*finalWin((n-((q-1)-(i-1)*numWinPerSub)*step),1);
            end
        end
    end
end
disp (['Window applying done for all subjects. All windows of time series are saved in one matrix!'])

% generate full correlation for each window
disp (['Compute the full correlation for each window'])
for q=1:numWinAllSub
    
    % compute the full correlation: Pearson's r
    fullCorrWin(:,:,q)=corrcoef(win(:,:,q));
end

% Fisher z tranform the correlations
zFullCorrWin=0.5*log((1+fullCorrWin)./(1-fullCorrWin));

disp('Full correlation are computed for each window.')

end


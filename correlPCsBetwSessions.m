
clear
clc
close all

sessionList={'rfMRI_REST1_RL','rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};
figDir = ['/home/data/Projects/HCP/figs/'];
norm='norm2';
numPC=10;
numSession=length(sessionList);

for j=1:length(sessionList)
    
    session=char(sessionList{j})
    
    resultDir=['/home/data/Projects/HCP/results/', session, '/'];
    
    tmp=load([resultDir,'/eigenvector_', norm, '_', session,'.mat']);
    eigenVector=tmp.COEFF;
    allSession(:,(1+10*(j-1)):j*10)=eigenVector(:,1:10);
end

correl=zeros(size(allSession, 2), size(allSession, 2));
for i=1:size(allSession,2)
    for j=1:size(allSession,2)
        [d z]=procrustes(allSession(:,i), allSession(:,j));
        correl(i,j)=corr(allSession(:,i),z);
    end
end

% plot correlation matrix for transformed and untransformed PCs
figure(1)
imagesc(correl)
colorbar
caxis([-1 1])

[r, p]=corr(allSession);

figure(2)
imagesc(r)
colorbar
caxis([-1 1])

 saveas(figure(1), [figDir, 'correlPCsOf4SessionsProcrustes.jpg'])
 saveas(figure(2), [figDir, 'correlPCsOf4Sessions.jpg'])

% plot the transformed and matched PCs
allMatch=zeros(numPC, numSession);
Z=[];
for j=1:numSession
    ses1=allSession(:, 1:10);
    a=correl(1:10, 1+10*(j-1):10*j);
    b=allSession(:,1+10*(j-1):10*j);
    for i=1:numPC
        match=find(a(i,:)==max(a(i,:)));
        allMatch(i,j)=match;
        [d1 z1]=procrustes(ses1(:,i),b(:,match));
        Z1(:, i, j)=z1;
    end
end

% due to the same PC in session 2 and 3 are aligned with one PC in session
% 1, these two PCs were manually matched based on the correlation values.
allMatch(8,2)=8;
allMatch(6,3)=5;
allMatch(9,3)=10;


% after rematch, transform the matched PCs again
for j=1:numSession
    ses1=allSession(:, 1:10);
        b=allSession(:,1+10*(j-1):10*j);
    for i=1:numPC
        [d2 z2]=procrustes(ses1(:,i),b(:,allMatch(i,j)));
        Z2(:, i, j)=z2; 
        correl2(i,j)=corr(ses1(:,i), z2);
    end
end


close all
figure(1)
imagesc(correl2(:,2:end))
colorbar
caxis([0.3 1])
set(gca, 'XTick', 1:3)
saveas(figure(1), [figDir, 'correlbetwEachSesProcrustestransformed.png'])
% plot the tranformed PCs
% PCs in Z2 are transformed and matched already
close all
sessionList={'rfMRI_REST1_RL','rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};
for j=1:length(sessionList)
    sessionPC=squeeze(Z2(:,:,j));
        
    %figure(j)
    for i=1:numPC
        PC=squareform(sessionPC(:,i));
        %subplot(3,4,i)
        figure(i)
        imagesc(PC)
        colorbar
        if i==1
            caxis([-0.025 0.025])
        else
            caxis([-0.03 0.03])
        end
        title(['PC', num2str(i)])
        axis square
         saveas(figure(i), [figDir, 'PC_', num2str(i), '_', norm, '_', char(sessionList{j}) '_matched.jpg'])
    end
   
end

% compute the correlations between the matched PCs
ses1=allSession(:, 1:10);
for j=2:numSession
   otherSes=Z2(:,:,j);
            for i=1:numPC
          [r, p]=corr(ses1(:,i), otherSes(:, i));
          rMatchedPC(i,j-1)=r;
      end
end
   

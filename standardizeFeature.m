function [ norm1FeatureWin norm2FeatureWin ] = standardizeFeature( subList, featureWin, numWinPerSub )
%This function standardize each subject's iFC windows by global variables
%and do demean for each feature across time windows
norm1FeatureWin=zeros(size(featureWin,1), size(featureWin,2));
norm2FeatureWin=zeros(size(featureWin,1), size(featureWin,2));
for i=1:length(subList)
    sub=num2str(subList(i))
    featureSub=featureWin(1+numWinPerSub*(i-1):numWinPerSub*i, :);
    featureSub1D=reshape(featureSub, [],1);
    featureSubNorm1=(featureSub-repmat(mean(featureSub1D), size(featureSub,1), size(featureSub,2)))./repmat(std(featureSub1D), size(featureSub,1), size(featureSub,2));
    featureSubNorm2=featureSubNorm1-repmat(mean(featureSubNorm1),size(featureSub, 1), 1);
    norm1FeatureWin(1+numWinPerSub*(i-1):numWinPerSub*i, :)=featureSubNorm1;
    norm2FeatureWin(1+numWinPerSub*(i-1):numWinPerSub*i, :)=featureSubNorm2;
end

end


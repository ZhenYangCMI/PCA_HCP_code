clear
clc
close all

addpath /home/data/Projects/Zhen/commonCode/BrainConnectivityToolbox20121204
addpath /home/data/Projects/Zhen/commonCode/BrainNetViewer1.43

session='all';

resultDir=['/home/data/Projects/Zhen/PCA_HCP/results/', session, '/'];

tmp=load([resultDir,'eigenvector_norm2_', session,'.mat']);
eigenVector=tmp.COEFF;

PCs=eigenVector(:,1:10);

% compute the threshold for edges
Alledges=reshape(PCs, [],1);
pos=Alledges(find(Alledges>0));
neg=Alledges(find(Alledges<0));

%for i=1:size(PCs, 2)
PC=PCs(:, 1);
PCThresh=PC;
for j=1:size(PC,1)
    if PC(j)>0 && PC(j)<prctile(pos,50) 
        PCThresh(j)=0;
    elseif PC(j)<0 && PC(j)>prctile(neg, 50)
        PCThresh(j)=0;
    else
        PCThresh(j)=PC(j);
    end
end 

NodeCoordinates=load('/home/data/Projects/Zhen/PCA_HCP/mask/Node_AAL88.txt');
EdgeMatrix=squareform(PC);
EdgeMatrixThreshed=squareform(PCThresh);
EdgeThreshold=0;
NodeWeight=sum(EdgeMatrix)';

[Ci Q] = modularity_louvain_und_sign(EdgeMatrix);
NodeColor=Ci';
NodeColorMap=[0 1 0;0 0 1];

for i=1:88
    NodeLabel{i,1}={{'-'}};
end
NodeLabel{1,1}={{'L PFC'}};
NodeLabel{2,1}={{'R PFC'}};
ModularIndex=unique(Ci)

viewtype='MediumView';
SurfFileName='/home/data/Projects/Zhen/commonCode/BrainNetViewer1.43/Data/SurfTemplate/BrainMesh_ICBM152.nv';

%H_BrainNet = y_CallBrainNetViewer_NodeEdge(NodeCoordinates,EdgeMatrix,EdgeThreshold,NodeWeight,1,NodeColor,NodeColorMap,NodeLabel,ModularIndex,viewtype,SurfFileName)
H_BrainNet = y_CallBrainNetViewer_NodeEdge(NodeCoordinates,EdgeMatrixThreshed)

clear
clc
close all

session='all';

resultDir=['/home/data/Projects/HCP/results/', session, '/'];

tmp=load([resultDir,'/eigenvector_norm2_', session,'.mat']);
eigenVector=tmp.COEFF;

PCs=eigenVector(:,1:10);

%for i=1:size(PCs, 2)
PC=PCs(:, 1);
%         a=prctile(PC,2)
%         b=prctile(PC, 98)
%         PC(find(PC>a & PC<b))=0;

NodeCoordinates=load('/home/data/Projects/HCP/mask/Node_AAL88.txt');
EdgeMatrix=squareform(PCs(:,1));
EdgeThreshold=0.02;
NodeWeight=sum(EdgeMatrix)';
NodeColor=repmat([1;2;3;4],22,1);
NodeColor(1:8,1)=[3;4;1;2;3;4;1;2];
NodeColorMap=[0 1 0;0 0 1; 1 0 0;1 0 1];
for i=1:88
    NodeLabel{i,1}={{'-'}};
end
NodeLabel{1,1}={{'L PFC'}};
NodeLabel{2,1}={{'R PFC'}};
ModularIndex=[1,2,3,4];
viewtype='MediumView'
SurfFileName='BrainMesh_ICBM152.nv';
%H_BrainNet = y_CallBrainNetViewer_NodeEdge(NodeCoordinates,EdgeMatrix,EdgeThreshold,NodeWeight,1,NodeColor,NodeColorMap,NodeLabel,ModularIndex,viewtype,SurfFileName)
H_BrainNet = y_CallBrainNetViewer_NodeEdge(NodeCoordinates,EdgeMatrix,EdgeThreshold,NodeWeight,1,NodeColor,NodeColorMap,NodeLabel,ModularIndex,viewtype)

function [ tag ] = extractTC(subList, ROIoutputDir, standardBrainMask, ROImask)
%This function standardize the volume and extract the TS of ROIs of a given
%mask

tag=0;
for i=1:length(subList)
    sub=num2str(subList(i))
    tag=tag+1;
    if ~exist([ROIoutputDir, sub], 'dir')
        mkdir(ROIoutputDir,sub)
    end
    subDir=[ROIoutputDir,sub];
    
    data=[subDir, '/funImg3mm.nii'];
    
    [AllVolume, VoxelSize, ImgFileList, Header1, nVolumn] =rest_to4d(data);
    [nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
    brainSize = [nDim1 nDim2 nDim3];
    
    % remove the regions outside of the brain and convert data into 2D
    MaskData=rest_loadmask(nDim1, nDim2, nDim3, standardBrainMask);
    MaskData =logical(MaskData);%Revise the mask to ensure that it contain only 0 and 1
    AllVolume=reshape(AllVolume,[],nDimTimePoints)';
    MaskDataOneDim=reshape(MaskData,1,[]);
    MaskIndex = find(MaskDataOneDim);
    AllVolume=AllVolume(:,MaskIndex);
    
    
    % Z_norm the time series for each voxel
    AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);
    AllVolume(isnan(AllVolume))=0;
    
    
    % Convert 2D file back into 4D
    AllVolumeBrain = single(zeros(nDimTimePoints, nDim1*nDim2*nDim3));
    AllVolumeBrain(:,MaskIndex) = AllVolume;
    AllVolumeBrain=reshape(AllVolumeBrain',[nDim1, nDim2, nDim3, nDimTimePoints]);
    
    
    % write 4D file as a nift file
    NormAllVolumeBrain=[subDir,'/','zFunImg.nii'];
    rest_Write4DNIfTI(AllVolumeBrain,Header1,NormAllVolumeBrain)
    
    disp ('Time series of each voxel was Z-score normalized.')
    
    
    % extract time series for seeds and ROIs
    
    y_ExtractROISignal(NormAllVolumeBrain, ...
        {ROImask},[subDir,'/HCP'],MaskData,1);
end

end


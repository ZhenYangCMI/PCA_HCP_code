% HCP data analysis

% organize the data

subList=load('/home/data/Projects/HCP/subjectList_97sub.txt');
scanList={'rfMRI_REST1_RL', 'rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};

motion=zeros(length(subList), length(scanList));

for i=1:length(subList)
    sub=num2str(subList(i))
    
    
    if exist(['/home/data/Originals/HCP/hcp1/', sub], 'dir')
        sourceDir=['/home/data/Originals/HCP/hcp1/', sub];
    else
        sourceDir=['/home/data/Originals/HCP/hcp2/', sub];
    end
    
    for j=1:length(scanList)
        scan=char(scanList{j})
        
        %copy over the image file
        sourceData=[sourceDir, '/MNINonLinear/Results/', scan, '/', scan, '_hp2000_clean.nii.gz'];
        
                if ~exist(['/home/data/Projects/HCP/data/', scan, '/', sub], 'dir')
                    mkdir(['/home/data/Projects/HCP/data/', scan, '/', sub])
                end
        
                destinationDir=['/home/data/Projects/HCP/data/', scan, '/', sub, '/'];
        
                copyfile(sourceData, [destinationDir, 'funImg.nii.gz'])
        
%         % read in and save the motion file
%         if strcmp(sub, '116120') | strcmp(sub, '207628')
%             meanRMS=0;
%         else
%             meanRMS=load([sourceDir, '/MNINonLinear/Results/', scan, '/Movement_RelativeRMS_mean.txt']);
%         end
%         
%         motion(i,j)=meanRMS;
    end
end

save(['/home/data/Projects/HCP/data/motion_full139sub.txt'],  '-ascii', '-tabs', 'motion')
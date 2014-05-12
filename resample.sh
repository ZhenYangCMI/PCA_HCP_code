
## 1. resample functional data 
# subList="cat /home/data/Projects/HCP/subjectList_97sub.txt" 

for scan in rfMRI_REST1_RL rfMRI_REST1_LR rfMRI_REST2_RL rfMRI_REST2_LR; do
echo $scan

	for sub in `cat /home/data/Projects/HCP/sublist97_test.txt`; do
	echo $sub
		if [ ! -e /home/data/Projects/HCP/data/${scan}/${sub}/funImg3mm.nii ]; then
	 	cmd="3dresample -dxyz 3.0 3.0 3.0 -prefix /home/data/Projects/HCP/data/${scan}/${sub}/funImg3mm.nii -inset /home/data/Projects/HCP/data/${scan}/${sub}/funImg.nii.gz"
		echo $cmd
		eval $cmd
	fi
	done
done




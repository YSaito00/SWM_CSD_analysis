# These scripts are adapted from the original code used in the study and have been modified for public release.
# They may depend on specific aspects of the developer’s computing environment; users should verify compatibility before running.
# For inquiries regarding the code or its usage, please contact Yoshito Saito (yoshitos@student.unimelb.edu.au).

# Directory settings
home= # home directory
csv_file="${home}/output.csv"
case_file="${home}/caselist.txt"
T1wdir= # directory where T1w images are stored
DWIdir= # directory where preprocessed FOD images are stored

# Create header row of the output CSV file
echo -n "Subject," > "${csv_file}"
skip_lines=0
tail -n +$((skip_lines + 1)) ${home}/brainregions_des.txt | while IFS=',' read -r brainregion col2; do # brainregions_des.txt contains ROI numbers and ROI names separated by commas on each line
    echo -n "${col2}_GMWM_fd,${col2}_L1_fd,${col2}_GMWM_FA,${col2}_L1_FA," >> "${csv_file}"
done
echo "" >> "${csv_file}"

while read case; do
	# Run 5ttgen and extract segmentation maps for GM, WM, and CSF
	echo -n "${case}," >> "${csv_file}"
	mkdir ${home}/${case}
	mrconvert ${T1wdir}/${case}_01_MR/T1w/T1w_acpc_dc_restore.nii.gz ${home}/${case}/T1w.mif -force
	5ttgen fsl ${home}/${case}/T1w.mif ${home}/${case}/5tt_nocoreg.mif -force
	mrtransform ${home}/${case}/5tt_nocoreg.mif -linear ${DWIdir}/${case}/diff2struct_mrtrix.txt -inverse ${home}/${case}/5tt_coreg.mif -force
	mrconvert ${home}/${case}/5tt_coreg.mif -coord 3 0 -axes 0,1,2 ${home}/${case}/5tt_coreg_GM.mif  -force
	mrconvert ${home}/${case}/5tt_coreg.mif -coord 3 1 -axes 0,1,2 ${home}/${case}/5tt_coreg_sub.mif  -force
	mrconvert ${home}/${case}/5tt_coreg.mif -coord 3 2 -axes 0,1,2 ${home}/${case}/5tt_coreg_WM.mif  -force
	mrconvert ${home}/${case}/5tt_coreg.mif -coord 3 3 -axes 0,1,2 ${home}/${case}/5tt_coreg_CSF.mif  -force
	mrconvert ${home}/${case}/5tt_coreg_GM.mif ${home}/${case}/5tt_coreg_GM.nii.gz -force
	mrconvert ${home}/${case}/5tt_coreg_WM.mif ${home}/${case}/5tt_coreg_WM.nii.gz -force
	mrconvert ${home}/${case}/5tt_coreg_sub.mif ${home}/${case}/5tt_coreg_sub.nii.gz -force
	mrconvert ${home}/${case}/5tt_coreg_CSF.mif ${home}/${case}/5tt_coreg_CSF.nii.gz -force
	fslmaths ${home}/${case}/5tt_coreg_WM.nii.gz -sub ${home}/${case}/5tt_coreg_GM.nii.gz -bin ${home}/${case}/5tt_coreg_WM-GM_bin.nii.gz
	fslmaths ${home}/${case}/5tt_coreg_WM.nii.gz -sub ${home}/${case}/5tt_coreg_CSF.nii.gz -bin ${home}/${case}/5tt_coreg_WM-CSF_bin.nii.gz
	fslmaths ${home}/${case}/5tt_coreg_WM.nii.gz -sub ${home}/${case}/5tt_coreg_sub.nii.gz -bin ${home}/${case}/5tt_coreg_WM-sub_bin.nii.gz
	# Create a mask containing only the WM component
	fslmaths ${home}/${case}/5tt_coreg_WM-sub_bin.nii.gz -mul ${home}/${case}/5tt_coreg_WM-GM_bin.nii.gz -mul ${home}/${case}/5tt_coreg_WM-CSF_bin.nii.gz ${home}/${case}/5tt_coreg_WMmask_bin.nii.gz
	mrtransform ${home}/${case}/5tt_coreg_WMmask_bin.nii.gz -template ${DWIdir}/${case}/${case}_wmfod_norm.mif ${home}/${case}/5tt_coreg_WMmask_inFODspace.nii.gz
	# Extract the GMWM boundary
	5tt2gmwmi ${home}/${case}/5tt_coreg.mif ${home}/${case}/gmwmSeed_coreg.mif -force
	mrconvert ${home}/${case}/gmwmSeed_coreg.mif  ${home}/${case}/gmwmSeed_coreg.nii.gz -force
	fslmaths ${home}/${case}/gmwmSeed_coreg.nii.gz -bin ${home}/${case}/gmwmSeed_coreg_bin.nii.gz
	mrtransform ${home}/${case}/gmwmSeed_coreg_bin.nii.gz -template ${DWIdir}/${case}/${case}_wmfod_norm.mif ${home}/${case}/gmwmSeed_coreg_inFODspace.nii.gz
	# Register FA image to FOD space
	convert_xfm -omat ${DWIdir}/${case}/struct2FOD.mat -inverse ${DWIdir}/${case}/diff2struct_fsl.mat
	flirt -in ${T1wdir}/${case}_01_MR/T1w/Diffusion/dti_FA.nii.gz -ref ${DWIdir}/${case}/${case}_wmfod_norm_b0.nii.gz -applyxfm -init ${DWIdir}/${case}/struct2FOD.mat -out ${home}/${case}/FA_inFODspace.nii.gz

	# Extract GMWM (SWM) and sub-SWM layers adjacent to each region
	tail -n +$((skip_lines + 1)) ${home}/brainregions_des.txt | while IFS=',' read -r brainregion col2; do
		fslmaths ${T1wdir}/${case}_01_MR/T1w/aparc.a2009s+aseg.nii.gz -thr ${brainregion} -uthr ${brainregion} -bin ${home}/${case}/aparc+aseg_${brainregion}.nii.gz
	 	mrtransform ${home}/${case}/aparc+aseg_${brainregion}.nii.gz -template ${DWIdir}/${case}/${case}_wmfod_norm.mif -linear ${DWIdir}/${case}/diff2struct_mrtrix.txt -inverse ${home}/${case}/aparc+aseg_${brainregion}_coreg.nii.gz -force
	 	fslmaths ${home}/${case}/aparc+aseg_${brainregion}_coreg.nii.gz -thr 0.5 -bin ${home}/${case}/aparc+aseg_${brainregion}_coreg.nii.gz
		# Extract GMWM layer adjacent to each region (SWM)
	 	fslmaths ${home}/${case}/aparc+aseg_${brainregion}_coreg.nii.gz -kernel sphere 2 -fmean -thr 0.0001 -bin ${home}/${case}/aparc+aseg_${brainregion}_coreg_mask.nii.gz
	 	fslmaths ${home}/${case}/gmwmSeed_coreg_inFODspace.nii.gz -mas ${home}/${case}/aparc+aseg_${brainregion}_coreg_mask.nii.gz ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz
	 	# Extract sub-SWM layer containing only WM elements adjacent to GMWM of each region
	 	fslmaths ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz -kernel sphere 2 -fmean -thr 0.0001 -bin ${home}/${case}/gmwmSeed_coreg_${brainregion}_L1.nii.gz
	 	fslmaths ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz -binv ${home}/${case}/inverse_mask_${brainregion}.nii.gz
	 	fslmaths ${home}/${case}/gmwmSeed_coreg_${brainregion}_L1.nii.gz -mas ${home}/${case}/inverse_mask_${brainregion}.nii.gz -mas ${home}/${case}/5tt_coreg_WMmask_inFODspace.nii.gz ${home}/${case}/gmwmSeed_coreg_${brainregion}_L1.nii.gz
	done

	# Compute the number of fixels per voxel
	mrcalc ${DWIdir}/${case}/fixel_in_subject_space/fd.mif 0 -ge ${home}/${case}/fixel_ROI/fd_bin.mif -force
	fixel2voxel ${home}/${case}/fixel_ROI/fd_bin.mif sum ${home}/${case}/fixel_ROI/fixel_vox.nii.gz -force
	# Compute the total FD value per voxel
	fixel2voxel ${DWIdir}/${case}/fixel_in_subject_space/fd.mif sum ${home}/${case}/fixel_ROI/fd_sum.nii.gz -force

	tail -n +$((skip_lines + 1)) ${home}/brainregions_des.txt | while IFS=',' read -r brainregion col2; do
		# Calculate FD in SWM of each region
		fslmaths ${home}/${case}/fixel_ROI/fixel_vox.nii.gz -mul ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz ${home}/${case}/fixel_ROI/fixel_${brainregion}_GMWM_vox.nii.gz
		# t = total number of fixels in that region’s SWM, y = total FD
		t=$(fslstats ${home}/${case}/fixel_ROI/fixel_${brainregion}_GMWM_vox.nii.gz -M -V | awk '{print $1*$2}')
		fslmaths ${home}/${case}/fixel_ROI/fd_sum.nii.gz -mul ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz ${home}/${case}/fixel_ROI/fd_gmwm_${brainregion}.nii.gz
		y=$(fslstats ${home}/${case}/fixel_ROI/fd_gmwm_${brainregion}.nii.gz -M -V | awk '{print $1*$2}')
		y_decimal=$(printf "%.4f" $y)
		t_decimal=$(printf "%.4f" $t)
		y=$(echo "scale=4; $y_decimal / $t_decimal" | bc)
		echo -n "$y," >> "${csv_file}"
		# Calculate FD in sub-SWM of each region
		fslmaths ${home}/${case}/fixel_ROI/fixel_vox.nii.gz -mul ${home}/${case}/gmwmSeed_coreg_${brainregion}_L1.nii.gz ${home}/${case}/fixel_ROI/fixel_${brainregion}_L1_vox.nii.gz
		t=$(fslstats ${home}/${case}/fixel_ROI/fixel_${brainregion}_L1_vox.nii.gz -M -V | awk '{print $1*$2}')
		fslmaths ${home}/${case}/fixel_ROI/fd_sum.nii.gz -mul ${home}/${case}/gmwmSeed_coreg_${brainregion}_L1.nii.gz ${home}/${case}/fixel_ROI/fd_L1_${brainregion}.nii.gz
		y=$(fslstats ${home}/${case}/fixel_ROI/fd_L1_${brainregion}.nii.gz -M -V | awk '{print $1*$2}')
		y_decimal=$(printf "%.4f" $y)
		t_decimal=$(printf "%.4f" $t)
		y=$(echo "scale=4; $y_decimal / $t_decimal" | bc)
		echo -n "$y," >> "${csv_file}"

		# Calculate FA in SWM of each region
		t=$(fslstats ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz -M -V | awk '{print $1*$2}')
		fslmaths ${home}/${case}/FA_inFODspace.nii.gz -mul ${home}/${case}/gmwmSeed_coreg_${brainregion}.nii.gz ${home}/${case}/fixel_ROI/FA_inFODspace_in${brainregion}.nii.gz
		y=$(fslstats ${home}/${case}/fixel_ROI/FA_inFODspace_in${brainregion}.nii.gz -M -V | awk '{print $1*$2}')
		y_decimal=$(printf "%.4f" $y)
		t_decimal=$(printf "%.4f" $t)
		y=$(echo "scale=4; $y_decimal / $t_decimal" | bc)
		echo -n "$y," >> "${csv_file}"
		# Calculate FA in sub-SWM of each region
		fslmaths ${home}/${case}/fixel_ROI/fixel_${brainregion}_L1_vox.nii.gz -bin ${home}/${case}/fixel_ROI/gmwmSeed_coreg_${brainregion}_L1.nii.gz
		t=$(fslstats ${home}/${case}/fixel_ROI/gmwmSeed_coreg_${brainregion}_L1.nii.gz -M -V | awk '{print $1*$2}')
		fslmaths ${home}/${case}/FA_inFODspace.nii.gz -mul ${home}/${case}/fixel_ROI/gmwmSeed_coreg_${brainregion}_L1.nii.gz ${home}/${case}/fixel_ROI/FA_inFODspace_in${brainregion}_L1.nii.gz
		y=$(fslstats ${home}/${case}/fixel_ROI/FA_inFODspace_in${brainregion}_L1.nii.gz -M -V | awk '{print $1*$2}')
		y_decimal=$(printf "%.4f" $y)
		t_decimal=$(printf "%.4f" $t)
		y=$(echo "scale=4; $y_decimal / $t_decimal" | bc)
		echo -n "$y," >> "${csv_file}"
	done
	echo "" >> "${csv_file}"
done < "${case_file}"

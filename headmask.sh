# script to get the headmask using the UTE template and the raw IMAGE

flirt -in res_ute.nii.gz -ref /media/sf_shared/regute2/ref_ute.nii.gz -dof 12 -omat ute_2_template.m -out ute_affine.nii.gz
fnirt --ref=ref_ute.nii.gz --in=ute_affine.nii.gz --iout=warped.nii.gz --fout=field.nii.gz --cout=coeffs.nii.gz --regmod=membrane_energy 
invwarp -w field.nii.gz -o inv_field.nii.gz -r ute_affine.nii.gz
applywarp -i /media/sf_shared/regute2/finalmask.nii.gz -r ./ute_affine.nii.gz -o mask_affine.nii.gz -w inv_field.nii.gz --interp=nn
convert_xfm -omat inv_affine.m -inverse ute_2_template.m
flirt -in mask_affine.nii.gz -ref res_ute.nii.gz -applyxfm -init inv_affine.m -interp nearestneighbour -out mask.nii.gz

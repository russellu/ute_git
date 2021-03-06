README for electrode hand labeling and segmentation based off of UTE image intensity.
author: Russell Butler (russell.buttler@usherbrooke.ca) 


basic instructions: 

step 1) run get_headmask.sh in a directory containing the subject's raw UTE (1 subject per directory)
	after setting the appropriate paths
	this will output a bunch of files, including a head mask in that subject's native space, which
	is to be used for the subsequent hand labeling and segmentation

step 2) run hand_labeling.m, and label the electrodes by placing the crosshair over each electrode, and
	clicking once on each electrode. a label will show up, marking the electrode's position, and name.
	once you have finished all the hand labeling, the program will generate the 30 random points on the
	scalp, and close.

step 3) run segment_electrodes.m, which will output the final coordinates in "xyz.nii.gz" 


utedata:

mask.nii.gz 
	=> manually segmented head mask
	=> voxels are 0s,1s, and 2s
		=> 0s=air, 1s=head, 2s=neck (used to trim pancake in hand_labeling.m)
		=> when applying the warp to this mask, make sure to use nearest neighbour interpolation

ref_ute.nii.gz 
	=> raw UTE template (non gradient), created by averaging 14 affinely registered single subjects
	=> for the affine transformation

ref_ute_gradient.nii.gz 
	=> UTE gradient template, for the nonlinear transformation




utescripts:

get_headmask.sh
	=> script to get the head mask, using the template and the original UTE

hand_labeling.m
	=> script to hand label electrodes, outputs hand labeled points in 3d MRI space

segment_electrodes.m
	=> script to segment electrodes, outputs center of mass of segmented electrodes

                                                                                                              
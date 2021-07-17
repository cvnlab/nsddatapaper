% we used this to determine where invalid data might be in the NSD dataset (so that we can report it).

% edit as necessary
cd ~/Dropbox/nsddata/ppdata/subj08/func1mm/

% look at all
itksnap -g T1_to_func1mm.nii.gz -o valid.nii.gz valid_nsdimagery.nii.gz valid_nsdsynthetic.nii.gz valid_prffloc.nii.gz mean.nii.gz brainmask.nii.gz

% look at which session in particular?
itksnap -g T1_to_func1mm.nii.gz -o valid.nii.gz valid_session*.nii.gz

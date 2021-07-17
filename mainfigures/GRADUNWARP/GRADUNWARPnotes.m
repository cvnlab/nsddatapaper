%% really quick and dirty to get screenshots of the gradunwarp

bet DATECENSORED_125001T1wMPRs022a1001.nii.gz test.nii.gz -m

fslmaths test_mask.nii.gz -kernel box 9 -dilM test_mask_dil.nii.gz

res is from itksnap quick and dirty

fslmaths DATECENSORED_125001T1wMPRs022a1001.nii.gz -mas res.nii.gz out1_masked.nii.gz -odt input
fslmaths DATECENSORED_125001T1wMPRs022a1001_gradunwarped.nii.gz -mas res.nii.gz out2_masked.nii.gz -odt input

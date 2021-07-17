% Create .nii versions of FreeSurfer's hippocampal segmentation
for p=1:8
  dir0 = sprintf('~/nsddata/freesurfer/subj%02d/mri',p);
  unix_wrapper(sprintf('mri_convert %s/lh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.mgz %s/lh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.nii.gz',dir0,dir0));
  unix_wrapper(sprintf('mri_convert %s/rh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.mgz %s/rh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.nii.gz',dir0,dir0));
  a1 = load_untouch_nii(sprintf('%s/lh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.nii.gz',dir0));
  a2 = load_untouch_nii(sprintf('%s/rh.hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.nii.gz',dir0));
  a1.img = repmat(sum(a1.img+a2.img,1),[size(a1.img,1) 1 1]);
  save_untouch_nii(a1,sprintf('%s/hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace.MIP.nii.gz',dir0));
  unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s',dir0));
end

% View in ITK-SNAP:
% - maximize / layout / move window to upper left
% - zoom x 2
% - select along midline
% - colormap jet
% - 0.8 for middle control point
% - no crosshair
% - save snapshots
% - create file "hippocampus slice planning.ai"
% - this is used to prescribe the high-res T2 slices

% see analysis_transforms.m regarding transformation of "hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace"

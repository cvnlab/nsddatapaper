% define
dirs = {'func1mm' 'func1pt8mm'};
res = {[1 1 1] [1.8 1.8 1.8]};
trs = [1 4/3];

% do it
for dd=1:length(dirs)
  for subjix=1:8

    % define
    file0  = sprintf('%s/ppdata/subj%02d/%s/mean.nii.gz',nsddatadir,subjix,dirs{dd});
    file0b = sprintf('%s/ppdata/subj%02d/%s/aseg.nii.gz',nsddatadir,subjix,dirs{dd});

    % load data
    a1  = load_untouch_nii(file0);
    a1b = load_untouch_nii(file0b);

    % homogenize mean volume
    cd ~/Dropbox/KKTEMP/
    for deg=5  % 1:9
      [newvol,brainmask,polymodel] = homogenizevolumes(double(a1.img),[99 1/10 deg 10],[],ismember(a1b.img,[3 42 8 47]));  % cortical gray matter and cerebellum gray matter
      imwrite(uint8(255*makeimagestack(newvol(:,:,1:10:end),[0 3])),gray(256),sprintf('bc_subj%02d.png',subjix));
      imwrite(uint8(255*makeimagestack(brainmask(:,:,1:10:end),1)),gray(256),sprintf('mask_subj%02d.png',subjix));
      imwrite(uint8(255*makeimagestack(polymodel(:,:,1:10:end),[0 1000])),hot(256),sprintf('poly_subj%02d.png',subjix));
    end
  
    % save coilbias and bc
    nsd_savenifti(single(polymodel),res{dd},sprintf('%s/ppdata/subj%02d/%s/mean_coilbias.nii.gz',nsddatadir,subjix,dirs{dd}),trs(dd));
    nsd_savenifti(single(newvol),   res{dd},sprintf('%s/ppdata/subj%02d/%s/mean_bc.nii.gz',      nsddatadir,subjix,dirs{dd}),trs(dd));
  
  end
end

%% IMPORT AND VET

% DO THIS
id = 'thalamus';
sourcename = 'combinedthalamus';

% ALSO DO THIS
id = 'MTL';
sourcename = 'MTLinT2';

% do it on stone
allcounts = [];
for subjix=1:8

  % define
  lhfile = sprintf('~/nsd/ppdata/%s/subj%02d_%s_lh.nii.gz',id,subjix,sourcename);
  rhfile = sprintf('~/nsd/ppdata/%s/subj%02d_%s_rh.nii.gz',id,subjix,sourcename);

  % load (512 x 512 x 512, 0.5 mm, uint8 format; LPI; exactly the same as our official T1 0.5mm format)
  lhroi = load_untouch_nii(lhfile);
  rhroi = load_untouch_nii(rhfile);
  
  % do some checks (are all values non-negative integers from 0 to N?)
  v1 = double(flatten(union(lhroi.img(:),[])));
  assert(isequal(v1,0:length(v1)-1));
  v2 = double(flatten(union(rhroi.img(:),[])));
  assert(isequal(v2,0:length(v2)-1));
  assert(isequal(v1,v2));
  assert(all(flatten(lhroi.img>0 & rhroi.img>0)==0));  % check lh and rh are mutually exclusive

  % do some counts for checking purposes
  allcounts(subjix,:,1) = countinstances(lhroi.img);
  allcounts(subjix,:,2) = countinstances(rhroi.img);
  
  % inspect
  img = max(lhroi.img,[],3);
  img = img + 10*max(rhroi.img,[],3);
  figureprep;
  imagesc(img); colormap(jet); axis image tight;
  figurewrite(sprintf('subj%02d',subjix),[],[],'~/Dropbox/KKTEMP');
  
  % write master version (0.5-mm resolution)
  outputfileLH =   sprintf('%s/ppdata/subj%02d/anat/roi/other/%s.%s_0pt5.nii.gz',nsd_datalocation,subjix,'lh',id);
  nsd_savenifti(int16(lhroi.img),            [0.5 0.5 0.5],outputfileLH);
  outputfileRH =   sprintf('%s/ppdata/subj%02d/anat/roi/other/%s.%s_0pt5.nii.gz',nsd_datalocation,subjix,'rh',id);
  nsd_savenifti(int16(rhroi.img),            [0.5 0.5 0.5],outputfileRH);
  outputfileBOTH = sprintf('%s/ppdata/subj%02d/anat/roi/other/%s_0pt5.nii.gz',nsd_datalocation,subjix,id);
  nsd_savenifti(int16(lhroi.img + rhroi.img),[0.5 0.5 0.5],outputfileBOTH);

  % make other versions (0.8-mm and 1-mm resolutions for anat)
  filestodo = {outputfileLH outputfileRH outputfileBOTH};
  strstodo1 = {['lh.' id]         ['rh.' id]          id};
  strstodo2 = {['lh.' id '_1pt0'] ['rh.' id '_1pt0'] [id '_1pt0']};
  for zz=1:length(filestodo)
    out0 = sprintf('%s/ppdata/subj%02d/anat/roi/%s.nii.gz',      nsd_datalocation,subjix,strstodo1{zz});
    nsd_mapdata(subjix,'anat0pt5','anat0pt8',filestodo{zz},'wta',NaN,out0,'int16');
    out0 = sprintf('%s/ppdata/subj%02d/anat/roi/other/%s.nii.gz',nsd_datalocation,subjix,strstodo2{zz});
    nsd_mapdata(subjix,'anat0pt5','anat1pt0',filestodo{zz},'wta',NaN,out0,'int16');
  end

  % make other versions (1.0-mm resolution for func, computed based on 0.8 anat version)
  filestodo = {sprintf('%s/ppdata/subj%02d/anat/roi/%s.%s.nii.gz',nsd_datalocation,subjix,'lh',id) ...
               sprintf('%s/ppdata/subj%02d/anat/roi/%s.%s.nii.gz',nsd_datalocation,subjix,'rh',id) ...
               sprintf('%s/ppdata/subj%02d/anat/roi/%s.nii.gz',   nsd_datalocation,subjix,id)};
  strstodo = {['lh.' id]      ['rh.' id]      id};
  for zz=1:length(filestodo)
    out0 = sprintf('%s/ppdata/subj%02d/func1mm/roi/%s.nii.gz',nsd_datalocation,subjix,strstodo{zz});
    nsd_mapdata(subjix,'anat0pt8','func1pt0',filestodo{zz},'wta',NaN,out0,'int16');
  end

  % make other versions (1.8-mm resolution for func, computed based on 1.0 anat version)
  filestodo = {sprintf('%s/ppdata/subj%02d/anat/roi/other/%s.%s_1pt0.nii.gz',nsd_datalocation,subjix,'lh',id) ...
               sprintf('%s/ppdata/subj%02d/anat/roi/other/%s.%s_1pt0.nii.gz',nsd_datalocation,subjix,'rh',id) ...
               sprintf('%s/ppdata/subj%02d/anat/roi/other/%s_1pt0.nii.gz',   nsd_datalocation,subjix,id)};
  strstodo = {['lh.' id]      ['rh.' id]      id};
  for zz=1:length(filestodo)
    out0 = sprintf('%s/ppdata/subj%02d/func1pt8mm/roi/%s.nii.gz',nsd_datalocation,subjix,strstodo{zz});
    nsd_mapdata(subjix,'anat1pt0','func1pt8',filestodo{zz},'wta',NaN,out0,'int16');
  end

end

% check
allcounts

% do visual sanity checks for all subjects!

% in the worst possible case (1.8mm) check that lh and rh are still mutually exclusive
for subjix=1:8
  lhfile = sprintf('%s/ppdata/subj%02d/func1pt8mm/roi/lh.%s.nii.gz',nsd_datalocation,subjix,id);
  rhfile = sprintf('%s/ppdata/subj%02d/func1pt8mm/roi/rh.%s.nii.gz',nsd_datalocation,subjix,id);
  lh = load_untouch_nii(lhfile);
  rh = load_untouch_nii(rhfile);
  assert(all(flatten(lh.img>0 & rh.img>0)==0));
end

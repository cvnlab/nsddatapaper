%% %%%%%%%% Manually define parameters

nsdsetup;

%% %%%%%%%% Use nsd_mapdata to bring brain mask to EPI space

% run on stone [anat->func]
for subjix=1:8

  % load the source 1pt0 brain mask
  sourcedata = sprintf('%s/ppdata/%s-structurals/brainmask_1pt0.nii.gz',nsddir,allnsdids{subjix});  

  % map it to the functional 1pt0 space
  outputfile = sprintf('%s/ppdata/subj%02d/func1mm/brainmask.nii.gz',nsddatadir,subjix);
  nsd_mapdata(subjix,'anat1pt0','func1pt0',sourcedata,'linear',0,outputfile,'int16');

  % map it to the functional 1pt8 space
  outputfile = sprintf('%s/ppdata/subj%02d/func1pt8mm/brainmask.nii.gz',nsddatadir,subjix);
  nsd_mapdata(subjix,'anat1pt0','func1pt8',sourcedata,'linear',0,outputfile,'int16');

end

%% %%%%%%%% Use nsd_mapdata to bring T1/T2/SWI/TOF to func spaces

% run on stone [anat->func]
strs = {'T1' 'T2' 'SWI' 'TOF'};
for subjix=1:8

  for p=1:length(strs)

    % map it to the functional 1pt0 space
    sourcedata = sprintf('%s/ppdata/subj%02d/anat/%s_0pt5_masked.nii.gz',nsddatadir,subjix,strs{p});  
    outputfile = sprintf('%s/ppdata/subj%02d/func1mm/%s_to_func1mm.nii.gz',nsddatadir,subjix,strs{p});
    nsd_mapdata(subjix,'anat0pt5','func1pt0',sourcedata,[],0,outputfile,'int16');

    % map it to the functional 1pt8 space
    sourcedata = sprintf('%s/ppdata/subj%02d/anat/%s_1pt0_masked.nii.gz',nsddatadir,subjix,strs{p});  
    outputfile = sprintf('%s/ppdata/subj%02d/func1pt8mm/%s_to_func1pt8mm.nii.gz',nsddatadir,subjix,strs{p});
    nsd_mapdata(subjix,'anat1pt0','func1pt8',sourcedata,[],0,outputfile,'int16');

  end

end

%% %%%%%%%% Use nsd_mapdata to bring 1mm T1/T2/SWI/TOF (masked) to MNI space

% run on stone [anat->MNI]
strs = {'T1' 'T2' 'SWI' 'TOF'};
for subjix=1:8
  for p=1:length(strs)
    sourcedata = sprintf('%s/ppdata/subj%02d/anat/%s_1pt0_masked.nii.gz',nsddatadir,subjix,strs{p});
    outputfile = sprintf('%s/ppdata/subj%02d/anat/%s_to_MNI.nii.gz',nsddatadir,subjix,strs{p});
    nsd_mapdata(subjix,'anat1pt0','MNI',sourcedata,'cubic',0,outputfile);
  end
end

%% %%%%%%%% Use nsd_mapdata to bring EPI 1mm to anat1pt0 and MNI

% run on stone
strs = {'MNI' 'anat1pt0'};
for subjix=1:8
  for p=1:length(strs)
    sourcedata = sprintf('%s/ppdata/subj%02d/func1mm/meanFIRST5.nii.gz',nsddatadir,subjix);
    outputfile = sprintf('%s/ppdata/subj%02d/anat/EPI_to_%s.nii.gz',nsddatadir,subjix,strs{p});
    nsd_mapdata(subjix,'func1pt0',strs{p},sourcedata,'cubic',0,outputfile,'int16');
  end
end

%% %%%%%%%% Use nsd_mapdata to prepare FS-related outputs

% Confirm the FS transformations:
%
% for subjix=1:8
% 
%   sourcedata = sprintf('%s/freesurfer/subj%02d/mri/T1.mgz',nsddatadir,subjix);
%   vol = cvnloadmgz(sourcedata);
% 
%   sourcedata = sprintf('%s/ppdata/subj%02d/anat/T1_0pt8_masked.nii.gz',nsddatadir,subjix);
%   vol2 = load_untouch_nii(sourcedata);
%   
%   imwrite(uint8(makeimagestack(double(vol2.img(:,:,30:20:end)),1)*255),gray(256),sprintf('~/Dropbox/KKTEMP/%d_1.png',subjix));
% 
%   volA = flipdim(flipdim(permute(vol,[1 3 2]),3),1);
%   volB = zeros(size(volA));
%   volB(2:end,:,2:end) = volA(1:end-1,:,1:end-1);
%   imwrite(uint8(makeimagestack(double(volB(:,:,30:20:end)),1)*255),gray(256),sprintf('~/Dropbox/KKTEMP/%d_2.png',subjix));

todo = {'aseg' 'surfaceimperfections' 'hippoSfLabels-T1-HST1T2.v10.FSvoxelSpace'};
todoout = {'aseg' 'surfaceimperfections' 'hippoSfLabels'};
for pp=1:length(todo)
  for subjix=1:8

    % load aseg
    sourcedata = matchfiles(sprintf('%s/freesurfer/subj%02d/mri/?h.%s.mgz',nsddatadir,subjix,todo{pp}));
    if isempty(sourcedata)
      sourcedata = matchfiles(sprintf('%s/freesurfer/subj%02d/mri/%s.mgz',nsddatadir,subjix,todo{pp}));
    end
    vol = 0;
    for p=1:length(sourcedata)
      vol = vol + cvnloadmgz(sourcedata{p});
    end
    figure; hist(double(vol(:)),100)

    % bring it to our anat0pt8 space
    vol = flipdim(flipdim(permute(vol,[1 3 2]),3),1);
    volB = zeros(size(vol));
    volB(2:end,:,2:end) = vol(1:end-1,:,1:end-1);
  
    % map it to the functional 1pt0 space
    outputfile = sprintf('%s/ppdata/subj%02d/func1mm/%s.nii.gz',nsddatadir,subjix,todoout{pp});
    nsd_mapdata(subjix,'anat0pt8','func1pt0',volB,'wta',0,outputfile,'int16');

    % map it to the functional 1pt8 space
    outputfile = sprintf('%s/ppdata/subj%02d/func1pt8mm/%s.nii.gz',nsddatadir,subjix,todoout{pp});
    nsd_mapdata(subjix,'anat0pt8','func1pt8',volB,'wta',0,outputfile,'int16');

    % save it in our anat0pt8 space
    outputfile = sprintf('%s/ppdata/subj%02d/anat/%s_0pt8.nii.gz',nsddatadir,subjix,todoout{pp});
    nsd_savenifti(int16(volB),[.8 .8 .8],outputfile);

    % map it to the anat1pt0 space
    outputfile = sprintf('%s/ppdata/subj%02d/anat/%s_1pt0.nii.gz',nsddatadir,subjix,todoout{pp});
    nsd_mapdata(subjix,'anat0pt8','anat1pt0',volB,'wta',0,outputfile,'int16');

    % map it to the anat0pt5 space
    outputfile = sprintf('%s/ppdata/subj%02d/anat/%s_0pt5.nii.gz',nsddatadir,subjix,todoout{pp});
    nsd_mapdata(subjix,'anat0pt8','anat0pt5',volB,'wta',0,outputfile,'int16');

  end
end

%% %%%%%%%% Use nsd_mapdata to transform surface imperfections to surface format

% define
hemis = {'lh' 'rh'};
layers = {'pial' 'layerB1' 'layerB2' 'layerB3' 'white'};

% loop over subject
for subjix=1:8

  % load 0.8mm surface imperfections volume
  file0 = sprintf('%s/ppdata/subj%02d/anat/surfaceimperfections_0pt8.nii.gz',nsddatadir,subjix);
  vol = load_untouch_nii(file0);
  
  % smooth with FWHM 2mm
  volB = smoothvolumes(vol.img==2,[0.8 0.8 0.8],[2 2 2]);

  % for each hemisphere, interpolate onto surfaces using cubic and take the max across depth
  for hh=1:2
    data = [];
    for zz=1:length(layers)
      data(:,zz) = nsd_mapdata(subjix,'anat0pt8',sprintf('%s.%s',hemis{hh},layers{zz}),volB,'cubic',[],[],'single');
    end
    data = max(data,[],2);
    nsd_savemgz(data,sprintf('%s/freesurfer/subj%02d/label/%s.surfaceimperfections.mgz',nsddatadir,subjix,hemis{hh}), ...
      sprintf('%s/freesurfer/subj%02d',nsddatadir,subjix));
  end

end

%% %%%%%%%% Use nsd_mapdata to transform some func-related measures to surface format [inherited from examples_funcviz.m]

% define
dirs = {sprintf('%s/ppdata/subj%%02d/func1mm/',nsddatadir) ...
        sprintf('%s/ppdata/subj%%02d/func1mm/betas_fithrf_GLMdenoise_RR/',nsd_datalocation('betas'))};
files =     {{'mean' 'valid' 'R2' 'signaldropout'} {'R2'}};
filenames = {{'mean' 'valid' 'R2' 'signaldropout'} {'b3R2'}};
hemis = {'lh' 'rh'};
fsdir = [nsddatadir '/freesurfer/subj%02d'];
nsubj = 8;
ndepth = 3;

% do it
for dd=1:length(dirs)
  for ff=1:length(files{dd})
    for hh=1:length(hemis)
      for ss=1:nsubj

        % load data
        sourcedata = sprintf('%s/%s.nii.gz',sprintf(dirs{dd},ss),files{dd}{ff});
        a1 = getfield(load_untouch_nii(sourcedata),'img');

        % map data to subject-native surfaces (cubic, conversion to single) and average across depth.
        data1 = [];
        for ii=1:ndepth
          data1(:,:,ii) = nsd_mapdata(ss,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),a1,'cubic',[],[],'single');
        end
        data1 = mean(data1,3);

        % save
        nsd_savemgz(data1,sprintf('%s/label/%s.%s.mgz',sprintf(fsdir,ss),hemis{hh},filenames{dd}{ff}),sprintf(fsdir,ss));

      end
    end
  end
end

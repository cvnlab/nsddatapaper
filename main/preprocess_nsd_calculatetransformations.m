function preprocess_nsd_calculatetransformations(subjix,mode)

% function preprocess_nsd_calculatetransformations(subjix,mode)
% 
% <subjix> is 1-8
% <mode> is
%   1 means EPI to ANAT
%   2 means ANAT to EPI
%   3 means ANAT to SURF
%   4 means FSAVERAGE TO SURF
%   5 means SURF TO FSAVERAGE
%   6 means ANAT to MNI
%   7 means MNI to ANAT
%   8 means EPI to SURF, MNI to SURF
%   9 means EPI to MNI, MNI to EPI

%% Theory

% "anat1pt0": common anatomical space for T1 and T2 (256 x 256 x 256 at 1 mm)
% "anat0pt8": common anatomical space for T1 and T2 (320 x 320 x 320 at 0.8 mm)
% "anat0pt5": common anatomical space for T1 and T2 (512 x 512 x 512 at 0.5 mm)
% "func1pt0": upsampled EPI space (1 mm)
% "func1pt8": regular EPI space (1.8 mm)
% "white": FreeSurfer white-matter surface
% "pial":  FreeSurfer pial surface
% "layerB1": FreeSurfer gray-matter surface at 25% depth (outer, close to pial)
% "layerB2": FreeSurfer gray-matter surface at 50% depth (middle)
% "layerB3": FreeSurfer gray-matter surface at 75% depth (inner, close to white)
% "fsaverage": FreeSurfer's surface-based common anatomical space
% "MNI": MNI space, determined by nonlinearly warping a subject's T1 to an MNI template

% All the "anat*" spaces are in alignment, and they have the same origin and
%   field-of-view. The origin lies at the exact center of each image slab.
%   Note that because the volumes have even numbers of voxels along each
%   dimension, the origin lies exactly in-between voxels.
% The "func*" spaces are not aligned with the "anat*" spaces. However, we have
%   determined the necessary warping between the "func*" and "anat*" spaces.
% The "func1pt0" space has an origin at the center of its slab.
% The "func1pt8" space also has an origin at the center of its slab. However, note
%   that the func1pt0 and func1pt8 spaces (due to imperfections in how the integers
%   work out) do not have exactly the same fields of view. Instead, the two spaces
%   are aligned such that the center of the "first corner voxel" of each of these two
%   spaces is identical. Thus, the two spaces are in register, but the overall field of 
%   view and location of the origin are slightly different.

% We aim to calculate and save coordinates (or nearest-neighbor indices) for 
% all possible space-to-space transformations.
%
% There should be a convenient routine that takes these coordinates and some data,
% does the required interpolation (or indexing), and then returns the result (possibly
% saving to disk). This will be nsd_mapdata.m.
%
% Note that for transformations involving fsaverage, for simplicity, we only save out the nearest-
% neighbor indexing version (not continuous coordinates, since that would require interpolation on
% spherical surfaces).

% Here are some examples of space-to-space transformations:
% Example: func1pt0 -> anat1pt0 means give me func1pt0 coordinates in an anat1pt0-like volume
% Example: func1pt0 -> layerB1 means give me func1pt0 coordinates associated with each layerB1 vertex
% Example: anat1pt0 -> layerB2 means give me anat1pt0 coordinates associated with each layerB2 vertex
% Example: func1pt0 -> MNI means give me func1pt0 coordinates in an MNI-like volume
% Example: fsaverage -> white means give me the index of the fsaverage vertex associated with each white vertex
% Example: white -> fsaverage means give me the index of the white vertex associated with each fsaverage vertex
% Example: func1pt0 -> fsaverage means give me func1pt0 coordinates associated with each fsaverage vertex

% Note that this function does not attempt to go from surfaces to volumes.

% Here are some example files that we will save:
% Files: func1pt0-to-anat1pt0.nii.gz
% Files: anat1pt0-to-func1pt0.nii.gz
% Files: func1pt0-to-layerB1.mgz
% Files: white-to-fsaverage.mgz

% Note that 9999 is a special value that means no valid data are available at that location.

% In general, we always write LPI NIFTIs and our coordinates are always matrix coordinates
% (1,1,1) is the center of the first voxel. Note that the MNI Template is in RPI order.
% Transformations involving MNI still respect the RPI load order (thus, 1 in the first 
% dimension means the center of the first voxel on the right side of the volume).
% When we write MNI NIFTI files, to be consistent we write out LPI NIFTIs. This means
% that the stored ordering is not the same as the ordering of the MNI Template NIFTIs, 
% but as long as the loading software understands and interprets correctly the ordering, 
% everything should be exactly consistent.

%% Setup

% setup
nsdsetup;
structdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});
tdir = sprintf('%s/ppdata/subj%02d/transforms',nsddatadir,subjix);
mkdirquiet(tdir);

% define
epifile =     sprintf('%s/ppdata/subj%02d/func1mm/meanFIRST5.nii.gz',nsddatadir,subjix);  % thing we warped in ANTS
t2file =      sprintf('%s/T2_1pt0_masked.nii.gz',          structdir);  % ref used in ANTS
trans1file =  sprintf('%s/EPIsyn_1Warp.nii.gz',            structdir);  % warp file
trans2file =  sprintf('%s/EPIsyn_0GenericAffine.mat',      structdir);  % affine file
itrans1file = sprintf('%s/EPIsyn_1InverseWarp.nii.gz',     structdir);  % inverse of warp
itrans2file = sprintf('"[%s/EPIsyn_0GenericAffine.mat,1]"',structdir);  % inverse of affine

%% Calculate transformations

switch mode





%% EPI TO ANAT

case 1

% Heavy lifting (map a grid).
% determine transform from EPI (1pt0) to T2 (1pt0).
% we do this by mapping a grid using linear interpolation.
% we obtain decimal EPI coordinates in the space of the T2.
% invalid T2 voxels (no EPI data) are assigned 9999.
a1 = load_untouch_nii(epifile);  % load EPI volume
coords = single([]);  % 256 x 256 x 256 x 3
for p=1:3
  a2 = a1;
  a2.img = fillmatrix(1:size(a1.img,p),size(a1.img),p);  % modulate along one dimension of EPI volume
  infile = [tempname '.nii'];
  save_untouch_nii(a2,infile);
  outfile = [tempname '.nii'];
  unix_wrapper(sprintf('antsApplyTransforms --dimensionality 3 --input %s --reference-image %s --output %s --interpolation Linear --transform %s --transform %s --default-value 9999',infile,t2file,outfile,trans1file,trans2file));
  coords(:,:,:,p) = getfield(load_untouch_nii(outfile),'img');
  delete(infile); delete(outfile);
end
nsd_savenifti(single(coords),[1 1 1],sprintf('%s/func1pt0-to-anat1pt0.nii.gz',tdir));

% Prepare finer resolution through interpolation.
% based on anat 1pt0 results (coords), deal with anat 0pt5 and 0pt8 results.
% our strategy is simply to perform linear interpolation at appropriate points in space.
% there may be some "bleed out" at the edges of the volume, but that can be dealt with later.
% ensure 9999s in the 1pt0 results propagate as NaNs.
voxelsizes =   [0.5 0.8];
voxelnumbers = [512 320];
voxellabels =  {'anat0pt5' 'anat0pt8'};
for p=1:length(voxelsizes)
  temp = resamplingindices(1,256,voxelnumbers(p));
  [xx,yy,zz] = ndgrid(temp,temp,temp);
  newcoords = [flatten(xx); flatten(yy); flatten(zz);];
  vals = single([]);  % N x N x N x 3
  for q=1:3
    vol = coords(:,:,:,q);
    vol(vol==9999) = NaN;  % don't let 9999 create bogus values
    vals(:,:,:,q) = nanreplace(reshape(ba_interp3_wrapper(vol,newcoords,'linear'),size(xx)),9999);
  end
  nsd_savenifti(single(vals),repmat(voxelsizes(p),[1 3]),sprintf('%s/func1pt0-to-%s.nii.gz',tdir,voxellabels{p}));
end

% Change the source file through simple deterministic coordinate transforms.
% deal with func 1pt8 resolution.
% simply take the 1pt0 files and do some massaging.
voxelsizes =  [1.0 0.5 0.8];
voxellabels = {'anat1pt0' 'anat0pt5' 'anat0pt8'};
for p=1:length(voxelsizes)

  % load the results from func1pt0
  a1 = load_untouch_nii(sprintf('%s/func1pt0-to-%s.nii.gz',tdir,voxellabels{p}));
  
  % the functional stuff is prepared in LPI (first voxel is left, posterior, inferior).
  % in our pre-processing, the fixed point for the 1.0 and 1.8 data is at Anterior, Right, Inferior.
  vals = a1.img;
  vals(:,:,:,1) = (vals(:,:,:,1) - nsdmatrixsize{subjix}(2)) / 1.8 + nsdmatrixsizeLOW{subjix}(2);
  vals(:,:,:,2) = (vals(:,:,:,2) - nsdmatrixsize{subjix}(1)) / 1.8 + nsdmatrixsizeLOW{subjix}(1);
  vals(:,:,:,3) = (vals(:,:,:,3) -                        1) / 1.8 + 1;
  
  % save
  nsd_savenifti(single(vals),repmat(voxelsizes(p),[1 3]),sprintf('%s/func1pt8-to-%s.nii.gz',tdir,voxellabels{p}));

end







%% ANAT TO EPI

case 2

% Heavy lifting (map a grid).
% determine transform from T2 (1pt0) to EPI (1pt0).
% we do this by mapping a grid using linear interpolation.
% we obtain decimal T2 coordinates in the space of the EPI.
% invalid EPI voxels (no T2 data) are assigned 9999.
a1 = load_untouch_nii(t2file);  % load T2 volume
coords = single([]);  % X x Y x Z x 3
for p=1:3
  a2 = a1;
  a2.img = fillmatrix(1:size(a1.img,p),size(a1.img),p);  % modulate along one dimension of T2 volume
  infile = [tempname '.nii'];
  save_untouch_nii(a2,infile);
  outfile = [tempname '.nii'];
  unix_wrapper(sprintf('antsApplyTransforms --dimensionality 3 --input %s --reference-image %s --output %s --interpolation Linear --transform %s --transform %s --default-value 9999',infile,epifile,outfile,itrans1file,itrans2file));
  coords(:,:,:,p) = getfield(load_untouch_nii(outfile),'img');
  delete(infile); delete(outfile);
end
nsd_savenifti(single(coords),[1 1 1],sprintf('%s/anat1pt0-to-func1pt0.nii.gz',tdir));

% Prepare finer resolution through interpolation.
% based on func 1pt0 results (coords), deal with func 1pt8 results.
% our strategy is simply to perform linear interpolation at appropriate points in space.
% there may be some "bleed out" at the edges of the volume, but that can be dealt with later.
% ensure 9999s in the 1pt0 results propagate as NaNs.
voxelsizes =   [1.8];
voxellabels =  {'func1pt8'};
for p=1:length(voxelsizes)
  [xx,yy,zz] = ndgrid(fliplr(linspacefixeddiff(size(coords,1),-1.8,nsdmatrixsizeLOW{subjix}(2))), ...  % voodoo
                      fliplr(linspacefixeddiff(size(coords,2),-1.8,nsdmatrixsizeLOW{subjix}(1))), ...
                      linspacefixeddiff(1,1.8,nsdmatrixsizeLOW{subjix}(3)));
  newcoords = [flatten(xx); flatten(yy); flatten(zz);];
  vals = single([]);  % N x N x N x 3
  for q=1:3
    vol = coords(:,:,:,q);
    vol(vol==9999) = NaN;  % don't let 9999 create bogus values
    vals(:,:,:,q) = nanreplace(reshape(ba_interp3_wrapper(vol,newcoords,'linear'),size(xx)),9999);
  end
  nsd_savenifti(single(vals),repmat(voxelsizes(p),[1 3]),sprintf('%s/anat1pt0-to-%s.nii.gz',tdir,voxellabels{p}));
end

% Change the source file through simple deterministic coordinate transforms.
% deal with anat 0pt5 and 0pt8 resolutions.
% simply take the 1pt0 files and do some massaging.
anatvoxelsizes = [0.5 0.8];
anatvoxellabels = {'anat0pt5' 'anat0pt8'};
voxelsizes =  [1.0 1.8];
voxellabels = {'func1pt0' 'func1pt8'};
for q=1:length(anatvoxelsizes)
  for p=1:length(voxelsizes)

    % load the results from anat1pt0
    a1 = load_untouch_nii(sprintf('%s/anat1pt0-to-%s.nii.gz',tdir,voxellabels{p}));
  
    % simple rescaling...
    vals = a1.img;
    vals(:,:,:,1) = (vals(:,:,:,1) - 0.5) / anatvoxelsizes(q) + 0.5;
    vals(:,:,:,2) = (vals(:,:,:,2) - 0.5) / anatvoxelsizes(q) + 0.5;
    vals(:,:,:,3) = (vals(:,:,:,3) - 0.5) / anatvoxelsizes(q) + 0.5;
  
    % save
    nsd_savenifti(single(vals),repmat(voxelsizes(p),[1 3]),sprintf('%s/%s-to-%s.nii.gz',tdir,anatvoxellabels{q},voxellabels{p}));

  end
end







%% ANAT (T1 or T2) TO SURF

case 3

% At the coordinates of each non-distorted surface, give me the coordinates into a
% given anat space.

% define voxel-to-world matrices for the 1.0, 0.8, and 0.5 mm anatomical data
VTW = { ...
    [ 1  0  0 -127.5;
      0  1  0 -127.5;
      0  0  1 -127.5;
      0  0  0      1;] ...
    [.8  0  0 -127.6;
      0 .8  0 -127.6;
      0  0 .8 -127.6;
      0  0  0      1] ...
    [.5  0  0 -127.75;
      0 .5  0 -127.75;
      0  0 .5 -127.75;
      0  0  0       1]};
anatlabels = {'anat1pt0' 'anat0pt8' 'anat0pt5'};

% define FS stuff 
%  unix(sprintf('mri_info --vox2ras-tkr ~/nsd/nsddata/freesurfer/subj%02d/mri/T1.mgz',p))
%  unix(sprintf('mri_info --vox2ras ~/nsd/nsddata/freesurfer/subj%02d/mri/T1.mgz',p))
Torig = [-.8 0 0 128; 0 0 .8 -128; 0 -.8 0 128; 0 0 0 1];         % vox2ras-tkr
Norig = [-.8 0 0 128.4; 0 0 .8 -127.6; 0 -.8 0 128.4; 0 0 0 1];   % vox2ras

% load each FS surface
surfs = {'white' 'pial' 'layerB1' 'layerB2' 'layerB3'};
prefixes = {'lh' 'rh'};
for p=1:length(prefixes)
  fsmgh = [];  % this is cached for each hemisphere
  
  for q=1:length(surfs)

% ORIGINAL APPROACH. BUGGY.  
%     % load in the surface coordinates
%     vertices = freesurfer_read_surf_kj(sprintf('%s/subj%02d/surf/%s.%s',nsdfsdir,subjix,prefixes{p},surfs{q}));
%     vertices = bsxfun(@plus,vertices',[128; 129; 128]);  % 3 x V [LPI, 1-256 mm] [consistent with fstoint.m of FS T1.nii.gz]
%     vertices(1,:) = vertices(1,:) + 0.8;  % our input anatomical volume is mismatched with FS's volume. to compensate, we need to shift 0.8-mm (1 voxel) towards Right.
%     vertices(3,:) = vertices(3,:) + 0.8;  % our input anatomical volume is mismatched with FS's volume. to compensate, we need to shift 0.8-mm (1 voxel) towards Superior.
%     vertices(4,:) = 1;           % now: 4 x V
%     vertices = VTW{1}*vertices;  % map surface vertices to world coordinates. We are now in millimeter coordinates and same space as our official anatomical volumes!!
%     
%     % process each anat space (we have to save coordinates with respect to the 3 anatomical resolutions)
%     for res=1:length(VTW)
%     
%       % map world coordinates to voxel coordinates (4 x V)
%       temp = inv(VTW{res})*vertices;
% 
%       % save a file like lh.anat0pt5-to-white.mgz
%       fsmgh = cvnwritemgz(sprintf('%s/subj%02d',nsdfsdir,subjix),sprintf('%s-to-%s',anatlabels{res},surfs{q}),temp(1:3,:),prefixes{p},tdir,[],fsmgh);
%     
%     end

    % load in the surface coordinates
    vertices = freesurfer_read_surf_kj(sprintf('%s/subj%02d/surf/%s.%s',nsdfsdir,subjix,prefixes{p},surfs{q}));
    vertices = vertices';  % 3 x V
    vertices(4,:) = 1;     % 4 x V
    vertices = Norig*inv(Torig)*vertices;  % map from rastkr to vox; map from vox to ras
    
    % process each anat space (we have to save coordinates with respect to the 3 anatomical resolutions)
    for res=1:length(VTW)
    
      % map world coordinates to voxel coordinates (4 x V)
      temp = inv(VTW{res})*vertices;  % 4 x V
      temp = temp(1:3,:)+1;           % 3 x V. Notice the +1. This is because in our convention, 1 is the center of the first voxel.

      % save a file like lh.anat0pt5-to-white.mgz
      fsmgh = cvnwritemgz(sprintf('%s/subj%02d',nsdfsdir,subjix),sprintf('%s-to-%s',anatlabels{res},surfs{q}),temp,prefixes{p},tdir,[],fsmgh);
    
    end

  end
  
end






%% FSAVERAGE TO SURF
% For each vertex of a subject's native surface, give me the index of the
% nearest fsaverage vertex.

%% SURF TO FSAVERAGE
% For each vertex of fsaverage, give me the index of the vertex of the
% nearest subject's native surface vertex.

case {4 5}

% loop
prefixes = {'lh' 'rh'};
for p=1:length(prefixes)
  fsmgh = [];  % this is cached for each hemisphere

  % where are the files?
  fsavgsurffile =  sprintf('%s/fsaverage/surf/%s.sphere',nsdfsdir,prefixes{p});
  nativesurffile = sprintf('%s/subj%02d/surf/%s.sphere.reg',nsdfsdir,subjix,prefixes{p});
  
  % load surfaces (note that we skip the post-processing of vertices and faces since unnecessary for what we are doing)
  clear fsavgsurf;
  [fsavgsurf.vertices, fsavgsurf.faces] =  freesurfer_read_surf_kj(fsavgsurffile);
  clear nativesurf;
  [nativesurf.vertices,nativesurf.faces] = freesurfer_read_surf_kj(nativesurffile);
  
  % do the calculations
  if mode==4

    % for nativesurf, find fsaverage indices
    vals = 1:size(fsavgsurf.vertices,1);
    interptype = 'nearest';
    f = griddata(fsavgsurf.vertices(:,1),fsavgsurf.vertices(:,2),fsavgsurf.vertices(:,3),vflatten(vals), ...
                 nativesurf.vertices(:,1),nativesurf.vertices(:,2),nativesurf.vertices(:,3),interptype);  % N x 1
    assert(all(isint(f) & f>=1));

    % save a file like lh.fsaverage-to-white.mgz
    fsmgh = cvnwritemgz(sprintf('%s/subj%02d',nsdfsdir,subjix),sprintf('%s-to-%s','fsaverage','white'),f',prefixes{p},tdir,[],fsmgh);

  else

    % for fsaverage, find nativesurf indices
    vals = 1:size(nativesurf.vertices,1);
    interptype = 'nearest';
    f = griddata(nativesurf.vertices(:,1),nativesurf.vertices(:,2),nativesurf.vertices(:,3),vflatten(vals), ...
                 fsavgsurf.vertices(:,1),fsavgsurf.vertices(:,2),fsavgsurf.vertices(:,3),interptype);  % N x 1
    assert(all(isint(f) & f>=1));

    % save a file like lh.white-to-fsaverage.mgz
    fsmgh = cvnwritemgz(sprintf('%s/subj%02d',nsdfsdir,subjix),sprintf('%s-to-%s','white','fsaverage'),f',prefixes{p},tdir,[],fsmgh);

  end

end






%% ANAT TO MNI
% For each voxel in MNI space, give me coordinates in the native subject anatomical space.

case 6

% define
Input = sprintf('%s/templates/MNI152_T1_1mm.nii.gz',nsddatadir);           % MNI template
WD = sprintf('%s/MNI',structdir);

% load existing file. we will clobber.
RefMask2 = sprintf('%s/T1_1pt0.nii.gz',structdir);                         % master T1 (1mm)
a1 = load_untouch_nii(RefMask2);

% Heavy lifting (map a grid).
% determine transform from ANAT (1mm) to MNI (1mm).
% we do this by mapping a grid using linear interpolation.
% we obtain decimal 1mm ANAT coordinates in MNI space.
coords = single([]);  % X x Y x Z x 3
for p=1:3
  a2 = a1;
  a2.img = fillmatrix(1:size(a1.img,p),size(a1.img),p);  % modulate along one dimension of volume
  infile = [tempname '.nii'];
  save_untouch_nii(a2,infile);
  outfile = [tempname '.nii'];
  unix_wrapper(sprintf('applywarp --datatype=float --rel --interp=trilinear --in="%s" --ref="%s" -w "%s"/str2standard.nii.gz -o "%s"',infile,Input,WD,outfile));
  coords(:,:,:,p) = getfield(load_untouch_nii([outfile '.gz']),'img');
  delete(infile); delete([outfile '.gz']);
end
nsd_savenifti(single(coords),[1 1 1],sprintf('%s/anat1pt0-to-MNI.nii.gz',tdir));

% Change the source file through simple deterministic coordinate transforms.
% deal with anat 0pt5 and 0pt8 resolutions.
% simply take the 1pt0 files and do some massaging.
anatvoxelsizes = [0.5 0.8];
anatvoxellabels = {'anat0pt5' 'anat0pt8'};
for q=1:length(anatvoxelsizes)

  % load the results from anat1pt0
  a1 = load_untouch_nii(sprintf('%s/anat1pt0-to-MNI.nii.gz',tdir));

  % simple rescaling...
  vals = a1.img;
  vals(:,:,:,1) = (vals(:,:,:,1) - 0.5) / anatvoxelsizes(q) + 0.5;
  vals(:,:,:,2) = (vals(:,:,:,2) - 0.5) / anatvoxelsizes(q) + 0.5;
  vals(:,:,:,3) = (vals(:,:,:,3) - 0.5) / anatvoxelsizes(q) + 0.5;

  % save
  nsd_savenifti(single(vals),[1 1 1],sprintf('%s/%s-to-MNI.nii.gz',tdir,anatvoxellabels{q}));

end









%% MNI TO ANAT
% For each voxel in a native subject anatomical space, give me coordinates into MNI (1mm) space.

case 7

% define
Input = sprintf('%s/T1_1pt0.nii.gz',structdir);                         % master T1 (1mm)
WD = sprintf('%s/MNI',structdir);

% load existing file. we will clobber.
RefMask2 = sprintf('%s/templates/MNI152_T1_1mm_brain_mask_dil_dilM.nii.gz',nsddatadir);
a1 = load_untouch_nii(RefMask2);

% Heavy lifting (map a grid).
% determine transform from MNI (1mm) to ANAT (1mm).
% we do this by mapping a grid using linear interpolation.
% we obtain decimal MNI coordinates in the 1mm subject anatomical space.
coords = single([]);  % X x Y x Z x 3
for p=1:3
  a2 = a1;
  a2.img = fillmatrix(1:size(a1.img,p),size(a1.img),p);  % modulate along one dimension of volume
  infile = [tempname '.nii'];
  save_untouch_nii(a2,infile);
  outfile = [tempname '.nii'];
  unix_wrapper(sprintf('applywarp --datatype=float --rel --interp=trilinear --in="%s" --ref="%s" -w "%s"/standard2str.nii.gz -o "%s"',infile,Input,WD,outfile));
  coords(:,:,:,p) = getfield(load_untouch_nii([outfile '.gz']),'img');
  delete(infile); delete([outfile '.gz']);
end
nsd_savenifti(single(coords),[1 1 1],sprintf('%s/MNI-to-anat1pt0.nii.gz',tdir));

% Prepare finer resolution through interpolation.
% based on anat 1pt0 results (coords), deal with anat 0pt5 and 0pt8 results.
% our strategy is simply to perform linear interpolation at appropriate points in space.
% there may be some "bleed out" at the edges of the volume, but that can be dealt with later.
% ensure 9999s in the 1pt0 results propagate as NaNs.
voxelsizes =   [0.5 0.8];
voxelnumbers = [512 320];
voxellabels =  {'anat0pt5' 'anat0pt8'};
for p=1:length(voxelsizes)
  temp = resamplingindices(1,256,voxelnumbers(p));
  [xx,yy,zz] = ndgrid(temp,temp,temp);
  newcoords = [flatten(xx); flatten(yy); flatten(zz);];
  vals = single([]);  % N x N x N x 3
  for q=1:3
    vol = coords(:,:,:,q);
    vol(vol==9999) = NaN;  % don't let 9999 create bogus values
    vals(:,:,:,q) = nanreplace(reshape(ba_interp3_wrapper(vol,newcoords,'linear'),size(xx)),9999);
  end
  nsd_savenifti(single(vals),repmat(voxelsizes(p),[1 3]),sprintf('%s/MNI-to-%s.nii.gz',tdir,voxellabels{p}));
end





%% EPI/MNI TO SURF
% At the coordinates of each non-distorted surface, give me the coordinates into the
% EPI or MNI space.  This is accomplished by EPI/MNI -> T2 -> SURF.

case 8

voxellabels = {'func1pt0' 'func1pt8' 'MNI'};
surfs = {'white' 'pial' 'layerB1' 'layerB2' 'layerB3'};
prefixes = {'lh' 'rh'};
for q=1:length(prefixes)
  fsmgh = [];  % this is cached for each hemisphere
  for p=1:length(voxellabels)

    % load XYZ-to-anat (where is each T2 anat voxel in XYZ space?) [some are 9999]
    a1 = load_untouch_nii(sprintf('%s/%s-to-anat1pt0.nii.gz',tdir,voxellabels{p}));

    for r=1:length(surfs)
      
      % load anat-to-surface (where is each surface vertex in anat space?)
      coord = cvnloadmgz(sprintf('%s/%s.anat1pt0-to-%s.mgz',tdir,prefixes{q},surfs{r}));  % V x 1 x 1 x 3

      % interpolate to determine: where is each surface vertex in XYZ space?
      temp = [];  % 3 x V
      for s=1:3
        vol = a1.img(:,:,:,s);
        vol(vol==9999) = NaN;  % don't let 9999 create bogus values
        temp(s,:) = nanreplace(ba_interp3_wrapper(vol,permute(coord,[4 1 2 3]),'cubic'),9999);
      end
      
      % save a file like lh.func1pt0-to-white.mgz
      fsmgh = cvnwritemgz(sprintf('%s/subj%02d',nsdfsdir,subjix),sprintf('%s-to-%s',voxellabels{p},surfs{r}),temp,prefixes{q},tdir,[],fsmgh);

    end
  end
end




%% EPI/MNI TO MNI/EPI
% At the coordinates of each MNI voxel, give me coordinates into EPI.
% This is accomplished by EPI -> T2 -> MNI.
% At the coordinates of each EPI voxel, give me coordinates into MNI.
% This is accomplished by MNI -> T2 -> EPI.

case 9

EPIv = {'func1pt0' 'func1pt8'};
MNIv = {'MNI'};
alltodo = {{EPIv MNIv} {MNIv EPIv}};
for p=1:length(alltodo)
  for q=1:length(alltodo{p}{1});
    for r=1:length(alltodo{p}{2});
      source0 = alltodo{p}{1}{q};
      target0 = alltodo{p}{2}{r};
      
      % load SOURCE-to-anat (where is each anat voxel in SOURCE space?)
      a1 = load_untouch_nii(sprintf('%s/%s-to-anat1pt0.nii.gz',tdir,source0));
      
      % load anat-to-TARGET (where is each TARGET voxel in anat space?)
      a2 = load_untouch_nii(sprintf('%s/anat1pt0-to-%s.nii.gz',tdir,target0));
      
      % interpolate
      newcoords = squish(a2.img,3)';  % 3 x ABC
      vals = single([]);  % X x Y x Z x 3
      for s=1:3
        vol = a1.img(:,:,:,s);
        vol(vol==9999) = NaN;  % don't let 9999 create bogus values
        vals(:,:,:,s) = nanreplace(reshape(ba_interp3_wrapper(vol,newcoords,'linear'),sizefull(a2.img,3)),9999);
      end
      nsd_savenifti(single(vals),a2.hdr.dime.pixdim(2:4),sprintf('%s/%s-to-%s.nii.gz',tdir,source0,target0));

    end
  end
end






end

% allow Dropbox update
unix_wrapper(sprintf('touch %s',tdir));

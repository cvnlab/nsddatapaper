%% %%%%%%%%%%%%%%%%%% EXPORT FILES

% notes:
% - you must run on stone to get the Dropbox to update! and might even need to touch the directory
% - for betas, we should run on surly (for speed).
%   in fact, i need surly in order to write nocompression betas .mat
% - we are zeroing-out betas! and zeroing-out timeseries volumes.
% - and we haven't done the screening session yet! [prf floc]
% - apparently, our strategy saves out LPI nifti files. this is straightforward Identity matrix.
%
% - we set origin at the center of the image slab.
% - telling make_nii the origin makes it do some magic so that 0,0,0 in world coordinates is that spot in the image matrix.
% - the matrix transforms the voxel coords such that center of image slab is 0,0,0.
% - note that for even numbers of voxels, there is no central voxel!
%
% - note that our integer saving scheme supports up to approximately -100/100% PSC maximum
%
% - in our local pp (and raw DICOM loading), we are:
%   matlab matrix: dim1 is A->P, dim2 is R->L, dim3 is I->S
% - thus to become our exported LPI, we have to swap the first two and reverse the first two.
% - our local pp nifti files are just 1.8 (or 1) and 1.33317 (or 1).  but the origin is not set and the flipping is weird. (and no masking yet)
%
% - mri_convert --out_orientation RAS actually creates a file that looks like LPI in ITK-SNAP.
%   this must mean that RAS corresponds to which direction is in POSITIVE direction.
% - compare FS input and FS output is the same raw data, just padded, flipped/oriented, and some minor intensity stuff
% - so, if we pad and flip/orient up front, we should be able to just rely on our input nifti.
% - so, LPI all the way; with identity-like matrix in the header (as from make_nii)
% - it does turn out that if you mri_convert T1.mgz T1.nii.gz, that it is RSP in ITK-SNAP but the underlying image content does not move.
%
% - dcm2nii from raw dicoms is coming with simple RPI (in itksnap) and identity headers. so we need to flip the first dimension.
% - interesting, MNI templates from FSL are also in RPI.
%
% - IMPORTANT: we don't want to use mri_convert to change resolution since it has origin issues.
%
% - in ITK-SNAP, in yoking multiple applications together, sometimes the same spot is not actually the same world coordinates.
%
% - testing make_nii for different resolutions:
%   0.5 is master, itk says -127.75
%   going to 1mm, itk says -127.5
%   going to .8mm, itk says -127.6 
%   thus, the origin numbers in the NIFTI *must* change in order for the center of the image slab to be the same.
%
% - note that for diffusion data, there is AP vs PA discrepancy in gradunwarp outputs (16bit vs 32bit). this is weird.
%
% - FS's internal T1.nii.gz is in RSP, and not only that, there is a weird 1-voxel offset relative to the
%   input volume. So, the underlying data is unchanged, but there are flipping, permuting, and voxel shifting happening.
%   See interpretation below. You can shift, flip, permute to match our canonical 0pt8 LPI space.
%
% - MNI template is [182 218 182]. World units (ITK) 0,0,0 corresponds to voxel coordinates (i.e. "Image units") (91,127,73) in RPI.
%
% remember difference wrt exported betas:
% - brain mask, int16, volfun orientation

%% %%%%%%% FREESURFER

% after we are done processing, let's move the subject directories to the official location.

% run on stone: [this is old/obsolete]
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj01/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj01/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj02/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj02/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj03/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj03/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj04/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj04/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj05/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj05/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj06/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj06/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj07/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj07/
rsync -av /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj08/* /home/stone/generic/Dropbox/nsddata/freesurfer/subj08/
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage          /home/stone/generic/Dropbox/nsddata/freesurfer/

% The final official edited FS outputs are at /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/
% The original FS outputs are at /home/surly-raid4/kendrick-data/nsd/nsddata_other/freesurferoriginals/

%% %%%%%%% EXPERIMENTS

% prf experiment
copy official RETBAR* and RETWEDGERINGMASH* files from /research/papers/2018\ hcp7tret/OSF 

%% %%%%%%% EXPERIMENTS (NSDSYNTHETIC)

% NOTE!!! this has the same filename as an internal version of this file.

cd /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD400-nsdsynthetic/mat_files_from_scan

masterordering = [];  % 1 x 744 with sequence of trials, indices relative to the 284
stimpattern = zeros(1,8,107);  % when did stimulus trials occur

a = matchfiles('*nsdsyn*.mat');
for p=1:length(a)
  a1 = load(a{p});
  cntRUN = 1;
  for q=1:size(a1.trialpattern,1)
    ix = find(a1.trialpattern(q,:));
    if ~isempty(ix)
      masterordering = [masterordering a1.curordering(cntRUN)];
      stimpattern(1,p,q) = 1;
      cntRUN = cntRUN + 1;
    end
  end
end

save('~/nsddata/experiments/nsdsynthetic/nsdsynthetic_expdesign.mat','masterordering','stimpattern');

%% %%%%%%% STRUCTURAL / ANATOMICAL DATA

% setup
nsdsetup;

% constants
wantsubj = 1:8;
typ = {'T1_????_masked' 'T2_????_masked' 'brainmask_????' 'T1*rep*masked' 'T2*rep*masked' 'SWI_????_masked' 'TOF_????_masked'};   %'EPI*Warped'

% loop
for subjix=wantsubj

  for p=[1:3 6 7]  %length(typ)
    unix_wrapper(sprintf('rsync -av %s/ppdata/%s-structurals/%s.nii.gz %s/ppdata/subj%02d/anat/',nsddir,allnsdids{subjix},typ{p},nsddatadir,subjix));
  end
  unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s/ppdata/subj%02d/anat/',nsddatadir,subjix));

  for p=4:5
    unix_wrapper(sprintf('rsync -av %s/ppdata/%s-structurals/%s.nii.gz /home/surly-raid4/kendrick-data/nsd/nsddata_other/subj%02d/',nsddir,allnsdids{subjix},typ{p},subjix));
  end

end

%% %%%%%%% BEHAVIORAL DATA (PRF + FLOC)

% setup
nsdsetup;

% loop over subjects
prffinal = [];
flocfinal = [];
for subjix=1:8
  datadir0 = nsdscreeningfun2(subjix);
  ppdir = regexprep(datadir0,'rawdata','ppdata');

  a1 = load([ppdir '/behavioralresults_prf.mat']);
  results = a1.results(1:3,:);
  save(sprintf('~/nsddata/bdata/prf/prf_subj%02d.mat',subjix),'results');
  prffinal(subjix) = mean(  (results(2,:)-results(3,:)) ./ results(1,:)  ) * 100;

  a1 = load([ppdir '/behavioralresults_floc.mat']);
  results = a1.results(1:3,:);
  save(sprintf('~/nsddata/bdata/floc/floc_subj%02d.mat',subjix),'results');
  flocfinal(subjix) = mean(results(2,:)./results(1,:))*100;

end

min(prffinal)
max(prffinal)
% ans =
% 
%           93.5017032074026
% 
% 
% ans =
% 
%           98.8611727886358
% 

min(flocfinal)
max(flocfinal)
%           90.8333333333333
% 
% 
% ans =
% 
%                       97.5
% 

%% %%%%%%% BEHAVIORAL DATA

% setup
nsdsetup;

% define
masterdir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
headers = {'SUBJECT' 'SESSION' 'RUN' 'TRIAL' '73KID' '10KID' 'TIME' ...
           'ISOLD' 'ISCORRECT' 'RT' 'CHANGEMIND' 'MEMORYRECENT' 'MEMORYFIRST' ...
           'ISOLDCURRENT' 'ISCORRECTCURRENT' 'TOTAL1' 'TOTAL2' 'BUTTON' 'MISSINGDATA'};

% loop over subjects
for subjix=setdiff(1:8,[1 2 3 4 5 6 8])

  % get all nsd session directories
  sessions = nsdallsessfun2(subjix);

  % loop over sessions, loading in the prepared behavioral data
  trialinfo = [];
  for ss=1:length(sessions)
    ppdir = regexprep(sessions{ss},'rawdata','ppdata');
    file0 = [ppdir '/behavioralresults_nsd.mat'];
    if ~exist(file0,'file')
      fprintf('MISSING: %s\n',file0);
      continue;
    end
    a1 = load(file0);
% NO LONGER NEEDED:
%     if size(a1.trialinfo,2)==15
%       a1.trialinfo(:,16) = 0;  % default is NOT missing data
%     elseif any(a1.trialinfo(:,16))
%       fprintf('DETECTED MISSING DATA for %s\n',sessions{ss});
%     end
    trialinfo = cat(1,trialinfo,a1.trialinfo);
  end
  
  % ensure that the times are relative (de-identification)
  trialinfo(:,7) = trialinfo(:,7) - floor(trialinfo(1,7));  % midnight of first day is time 0

  % write it out
  dir0 = sprintf('%s/ppdata/subj%02d/behav',masterdir,subjix);
  mkdirquiet(dir0);
  file0 = sprintf('%s/responses.tsv',dir0);
  savetext(file0,sprintf([repmat('%s\t',[1 18]) '%s'],headers{:}));  % 19 total
  dlmwrite(file0,trialinfo,'delimiter','\t','precision',20,'-append');
  unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s',dir0));

end

%% %%%%%%% BEHAVIORAL DATA (NSDSYNTHETIC)

% setup
nsdsetup;

% define
masterdir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
headers = {'SUBJECT' 'RUN' 'TRIAL' '284ID' 'TIME'};

% loop over subjects
for subjix=1:8

  % load in the prepared behavioral data
  ppdir = regexprep(datadirs{nsdsyntheticsession(subjix)},'rawdata','ppdata');
  file0 = [ppdir '/behavioralresults_nsdsynthetic.mat'];
  if ~exist(file0,'file')
    fprintf('MISSING: %s\n',file0);
    continue;
  end
  a1 = load(file0);
  trialinfo = a1.trialinfo;
  
  % ensure that the times are relative (de-identification)
  trialinfo(:,5) = trialinfo(:,5) - floor(trialinfo(1,5));  % midnight of first day is time 0

  % write it out
  file0 = sprintf('~/Dropbox/KKTEMP/subj%02d_nsdsynthetic.tsv',subjix);
  savetext(file0,sprintf([repmat('%s\t',[1 5-1]) '%s'],headers{:}));
  dlmwrite(file0,trialinfo,'delimiter','\t','precision',20,'-append');

end

%% CONTINUE

% loop over subjects
for subjix=1:8
  ppdir = regexprep(datadirs{nsdsyntheticsession(subjix)},'rawdata','ppdata');

  a1 = load([ppdir '/behavioralresults_nsdsynthetic.mat']);
  a1 = rmfield(a1,{'trialinfo' 'userkeys' 'userkeycounts' 'totaldur' 'badtimes'});
  a1
  save(sprintf('~/Dropbox/nsddata/bdata/nsdsynthetic/nsdsynthetic_subj%02d.mat',subjix),'-struct','a1');

end

%% %%%%%%% META [just for Word document]

% Ran on July 6 2019.

nsdsetup;
for p=1:8
  [leftN,rightN] = cvnreadsurface(sprintf('subj%02d',p),{'lh' 'rh'},'white','orig','justcount',true);
  fprintf('Subject %d\t%s\t%s\t%d\t%d\n',p,mat2str(nsdmatrixsize{p}([2 1 3])),mat2str(nsdmatrixsizeLOW{p}([2 1 3])),leftN,rightN);
end

%% %%%%%%% MOTION

% setup
nsdsetup;

% define
masterdir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
ppdirnames = {'preprocessVER1SECONDSTAGEfigures' 'preprocessVER1SECONDSTAGELOWfigures'};
funcnames = {'func1mm' 'func1pt8mm'};

% do it
for subjix=1:8
  for ver=1:length(funcnames)

    % get all nsd session directories (rawdata)
    sessions = nsdallsessfun2(subjix);    % paths
    sessionixs = nsdallsessfun3(subjix);  % positive integers
    
    % loop
    for ss=1:length(sessions)+3
    
      % don't re-do files!
      if 0
        if ~(sessionixs(ss)>=225 && sessionixs(ss)<=239)
          continue;
        end
      end

      % calc stuff
      if ss==length(sessions)+1
        session0 = datadirs{screeningmatrix(subjix)};
        suffix0 = '_prffloc';
      elseif ss==length(sessions)+2
        session0 = datadirs{nsdsyntheticsession(subjix)};
        suffix0 = '_nsdsynthetic';
      elseif ss==length(sessions)+3
        session0 = datadirs{nsdimagerysession(subjix)};
        suffix0 = '_nsdimagery';
      else
        session0 = sessions{ss};
        suffix0 = sprintf('_session%02d',ss);
      end
      ppdir = regexprep(session0,'rawdata','ppdata');

      % load
      a1 = load(sprintf([ppdir '/%s/record.mat'],ppdirnames{ver}),'mparams');
      
      for rr=1:length(a1.mparams)
        file0 = sprintf('%s/ppdata/subj%02d/%s/motion/motion%s_run%02d.tsv', ...
                        nsddatatimeseriesdir,subjix,funcnames{ver},suffix0,rr);
        mkdirquiet(stripfile(file0));
          % convert to "nicer" export format
        tr1 = spm_matrix(a1.mparams{rr}(1,:));  % ref voxels to world space
        exportparams = [];
        for p=1:size(a1.mparams{rr},1)-1
          tr2 = spm_matrix(a1.mparams{rr}(1+p,:));  % target voxels to world space
          temp0 = spm_imatrix(tr2/tr1);
          exportparams(p,:) = temp0(1:6);
        end
%        dlmwrite(file0,a1.mparams{rr}(2:end,1:6),'delimiter','\t','precision',20);
          % write them out
        rmdirquiet(file0);
        dlmwrite(file0,exportparams,'delimiter','\t','precision',20);
      end
    end
  end
end

% POSTHOC for nsd core sessions:
% fix subj08-nsd02.
% this was a frankenstein session that needed adjustments.
% see datafix_script3b.m for details.

%% %%%%%%% FUNCTIONAL DATA

% setup
nsdsetup;

% define
wantsubj = 1:8;
wantvars = [1 1 1 1 1];  % mean volume, R2 (and onoffbeta), single-trial betas, valid, timeseries
wantvers = [1 2];    % 1 is high res, 2 is low res

% OBSOLETE:
%   % both versions for some
%   wantvars = [1 1 0 1];
%   wantvers = [1 2];
% 
%   % only low res for betas
%   wantvars = [0 0 1 0];
%   wantvers = [2];
% 
%   % only high res for betas
%   wantvars = [0 0 1 0];
%   wantvers = [1];
% 
%   % both versions for all [FINAL]
%   wantvars = [1 1 1 1];
%   wantvers = [1 2];
  
% which ones to do
whichbetas = [1 2 3 4 5];

% NOTE: use this for nsdimagery:
wantvars = [1 0 1 1 1];
wantvers = [1 2];    % 1 is high res, 2 is low res
whichbetas = [5];

% NOTE: use this for nsdsynthetic:
wantvars = [1 0 1 1 1];
wantvers = [1 2];    % 1 is high res, 2 is low res
whichbetas = [6 7];

% constants
volfun = @(x) flipdim(flipdim(permute(x,[2 1 3 4 5 6 7 8 9 10]),2),1);
volfun2 = @(x) flipdim(rotatematrix(x,1,2,-1),1);
ivolfun = @(x) permute(flipdim(flipdim(x,1),2),[2 1 3 4 5 6 7 8 9 10]);
masterdir =       '/home/surly-raid4/kendrick-data/nsd/nsddata/';
masterdir_betas = '/home/surly-raid4/kendrick-data/nsd/nsddata_betas/';
betanames = {'betas_assumehrf' 'betas_fithrf' 'betas_fithrf_GLMdenoise_RR' 'restingbetas_fithrf' 'nsdimagerybetas_fithrf' 'nsdsyntheticbetas_fithrf' 'nsdsyntheticbetas_fithrf_GLMdenoise_RR'};

% loop
for subjix=wantsubj





%   % THIS IS OLD, OBSOLETE:
%   % this is for prior to nsdimagery (the main sessions)
%   if 0
%   if 0
%     % get all nsd session directories
%     sessions = nsdallsessfun2(subjix);    %alldirs0 = [nsdscreeningfun2(subjix) ];
%   end
% 
%     sessmatrixTEMP = sessmatrix;
% 
%   %%%% HERE IS WHERE WE EDIT!!!!!!   [THIS IS ALL ONE-OFFS. NSDIMAGERY WILL REVISIT]
%   % %1  sessmatrixTEMP(sessmatrixTEMP>93) = 0;
%   %   sessmatrixTEMP(sessmatrixTEMP<=93) = 0;
% 
%   %  sessmatrixTEMP(sessmatrixTEMP>152) = 0;
%   
%   %  sessmatrixTEMP(sessmatrixTEMP>197) = 0;  % beyond 197 DO NOT DO
%   %  sessmatrixTEMP(sessmatrixTEMP>224) = 0;  % beyond 224 DO NOT DO
%     sessmatrixTEMP(sessmatrixTEMP>239) = 0;  % beyond 239 DO NOT DO
%     sesstodo = filterout(sessmatrixTEMP(subjix,:),0)  % a vector of positive integers
%     sessions = datadirs(sesstodo)                     % cell vector of paths
%   %  sesstoexport = sesstodo >= 153                    % binary vector (same size as sesstodo) indicating these are new things to export/save
%   %  sesstoexport = sesstodo > 0;
%   %  sesstoexport = sesstodo >= 198                    % binary vector (same size as sesstodo) indicating these are new things to export/save
%     sesstoexport = sesstodo >= 225                    % binary vector (same size as sesstodo) indicating these are new things to export/save
% 
%     isstandalone = 0;
%   end
  
  
 
  % this is for nsdimagery
  sesstodo = nsdimagerysession(subjix);
  sessions = datadirs(sesstodo);
  sesstoexport = ones(1,length(sesstodo));
  isstandalone = 1;


  % this is for nsdsynthetic
  sesstodo = nsdsyntheticsession(subjix);
  sessions = datadirs(sesstodo);
  sesstoexport = ones(1,length(sesstodo));
  isstandalone = 1;

  % nsdcore
  sesstodo = nsdallsessfun3(subjix);
  sessions = datadirs(sesstodo);
  sesstoexport = ones(1,length(sesstodo));
  isstandalone = 0;


  % NOTE: DO NOT RUN IF NO NEW SESSION TO EXPORT
  if ~any(sesstoexport)
    continue;
  end


  % loop
  ppdirnames = {'preprocessVER1SECONDSTAGE' 'preprocessVER1SECONDSTAGELOW'};
  glmdirnames = {'glmdata' 'glmdataLOW'};
  funcnames = {'func1mm' 'func1pt8mm'};
  voxsize = {[1 1 1] [1.8 1.8 1.8]};
  for ver=wantvers  % 1 is high res, 2 is low res
  
    % init and load
    meanvol =  single([]);
    meanR2 =   single([]);
    onoffbeta = single([]);
    meanvalid = single([]);
    for ss=1:length(sessions)
      if ismember(sesstodo(ss),nsdimagerysession)
        tr = [1 1];
        suffix0 = '_nsdimagery';
      elseif ismember(sesstodo(ss),nsdsyntheticsession)
        tr = [1 4/3];
        suffix0 = '_nsdsynthetic';
      else
        tr = [1 4/3];
        suffix0 = sprintf('_session%02d',ss);
      end
        
      % figure out where various directories are
      glmdir = regexprep(sessions{ss},'rawdata',glmdirnames{ver});
      ppdir =  regexprep(sessions{ss},'rawdata','ppdata');
    
      % load mean volume
      if wantvars(1)
        a1 = load_untouch_nii(sprintf('%s/%s/mean.nii',ppdir,ppdirnames{ver}));  % stored as int16
        meanvol0 = int16(a1.img);
        sz = sizefull(meanvol0,3);
            %a1 = load_untouch_nii(sprintf([ppdir '/%s/valid.nii'],ppdirnames{ver}))  % int16
        meanvol(:,:,:,ss) = single(meanvol0);
      end
      
      % load R2
      if wantvars(2)
        a1 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'R2','modelmd');  % stored as single
        varexp0 = single(a1.R2);
        onoffbeta0 = single(a1.modelmd);
        sz = sizefull(varexp0,3);
        meanR2(:,:,:,ss) = varexp0;
        onoffbeta(:,:,:,ss) = onoffbeta0;
      end

      % load valid
      if wantvars(4)
        a1 = load_untouch_nii(sprintf('%s/%s/valid.nii',ppdir,ppdirnames{ver}));  % stored as int16 ?
        valid0 = int16(a1.img);
        sz = sizefull(valid0,3);
        meanvalid(:,:,:,ss) = single(valid0);
      end

      % save mean volume
      if wantvars(1) && sesstoexport(ss)
        file0 = sprintf('%s/ppdata/subj%02d/%s/mean%s.nii', ...
                        masterdir,subjix,funcnames{ver},suffix0);
        mkdirquiet(stripfile(file0));
        save_nii(settr_nii(make_nii(volfun(meanvol0),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
        gzip(file0); delete(file0);
      end

      % save R2
      if wantvars(2) && sesstoexport(ss)
        file0 = sprintf('%s/ppdata/subj%02d/%s/R2%s.nii', ...
                        masterdir,subjix,funcnames{ver},suffix0);
        mkdirquiet(stripfile(file0));
        save_nii(settr_nii(make_nii(volfun(varexp0),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
        gzip(file0); delete(file0);

        file0 = sprintf('%s/ppdata/subj%02d/%s/onoffbeta%s.nii', ...
                        masterdir,subjix,funcnames{ver},suffix0);
        mkdirquiet(stripfile(file0));
        save_nii(settr_nii(make_nii(volfun(onoffbeta0),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
        gzip(file0); delete(file0);
      end

      % save valid
      if wantvars(4) && sesstoexport(ss)
        file0 = sprintf('%s/ppdata/subj%02d/%s/valid%s.nii.gz', ...
                        masterdir,subjix,funcnames{ver},suffix0);
        nsd_savenifti(volfun(valid0),voxsize{ver},file0,tr(ver));
      end

      % deal with time-series
      if wantvars(5) && sesstoexport(ss)
        for runii=1:1000
          tic;

          % load run
          runfile0 = sprintf('%s/%s/run%02d.nii',ppdir,ppdirnames{ver},runii);
          if ~exist(runfile0,'file')
            break;
          end
          a1 = load_untouch_nii(runfile0);  % stored as int16

          % get mask
          mask0 = logical(ivolfun(getfield(load_untouch_nii(gunziptemp(sprintf('%s/ppdata/subj%02d/%s/brainmask.nii.gz', ...
                    nsddatadir,subjix,funcnames{ver}))),'img')));
                              
          % zero-out
          a1.img(repmat(~mask0,[1 1 1 size(a1.img,4)])) = 0;

          % save
          file0 = sprintf('%s/ppdata/subj%02d/%s/timeseries/timeseries%s_run%02d.nii.gz', ...
                          nsddatatimeseriesdir,subjix,funcnames{ver},suffix0,runii);
          nsd_savenifti(volfun(a1.img),voxsize{ver},file0,tr(ver));
          
          toc;
        end
      end

    end

    %% these are all grand averages
    
    % save
    if wantvars(1) && any(sesstoexport) && ~isstandalone
      file0 = sprintf('%s/ppdata/subj%02d/%s/mean.nii', ...
                      masterdir,subjix,funcnames{ver});
      mkdirquiet(stripfile(file0));
      save_nii(settr_nii(make_nii(volfun(single(mean(meanvol,4))),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
      gzip(file0); delete(file0);
    end
    
    % save
    if wantvars(1) && 0 && ~isstandalone   % NO NEED TO DO THIS AGAIN
      file0 = sprintf('%s/ppdata/subj%02d/%s/meanFIRST5.nii', ...
                      masterdir,subjix,funcnames{ver});
      mkdirquiet(stripfile(file0));
      save_nii(settr_nii(make_nii(volfun(single(mean(meanvol(:,:,:,1:5),4))),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
      gzip(file0); delete(file0);
    end

    % save
    if wantvars(2) && any(sesstoexport) && ~isstandalone
      file0 = sprintf('%s/ppdata/subj%02d/%s/R2.nii', ...
                      masterdir,subjix,funcnames{ver});
      mkdirquiet(stripfile(file0));
      save_nii(settr_nii(make_nii(volfun(single(nanmean(meanR2,4))),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
      gzip(file0); delete(file0);

      file0 = sprintf('%s/ppdata/subj%02d/%s/onoffbeta.nii', ...
                      masterdir,subjix,funcnames{ver});
      mkdirquiet(stripfile(file0));
      save_nii(settr_nii(make_nii(volfun(single(nanmean(onoffbeta,4))),voxsize{ver},([1 1 1]+sz([2 1 3]))/2),tr(ver)),file0);
      gzip(file0); delete(file0);
    end

    % save
    if wantvars(4) && any(sesstoexport) && ~isstandalone
      file0 = sprintf('%s/ppdata/subj%02d/%s/valid.nii.gz', ...
                      masterdir,subjix,funcnames{ver});
      nsd_savenifti(volfun(mean(meanvalid,4)),voxsize{ver},file0,tr(ver));
    end
    
    %% end.

    % touch just to ensure that synchronization utilities (e.g. Dropbox) notice
    unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s/ppdata/subj%02d/%s/',masterdir,subjix,funcnames{ver}));

    %%%%%% GO TO BETAS

    if wantvars(3)

    for betaii=whichbetas
      betaname = betanames{betaii};
    
      % init and load
      meanbeta = single([]);
      R2 =       single([]);
      r2ix = [];
      for ss=1:length(sessions)
      
        isrestingstate = betaii==4;
        if isrestingstate
          if ~(sessiontypes(sesstodo(ss))==8)   % if not a resting-state session, continue
            continue;
          end
        end

        isimagery = betaii==5;

        issynthetic = ismember(betaii,[6 7]);

        if ismember(sesstodo(ss),nsdimagerysession)
          tr = [1 1];
          suffix0 = '_nsdimagery';
        elseif ismember(sesstodo(ss),nsdsyntheticsession)
          tr = [1 4/3];
          suffix0 = '_nsdsynthetic';
        else
          tr = [1 4/3];
          suffix0 = sprintf('_session%02d',ss);
        end
  
        % figure out where various directories are
        glmdir = regexprep(sessions{ss},'rawdata',glmdirnames{ver});
        ppdir =  regexprep(sessions{ss},'rawdata','ppdata');

        % only do beta stuff if want to export
        if sesstoexport(ss)

          % load single-trial betas and other
          switch betaii
          case 1
            a2 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrial.mat',glmdir), ...
                      'modelmd','R2','R2run');
          case {2 5 6}
            a2 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir), ...
                      'modelmd','HRFindex','HRFindexrun','R2','R2run');
          case {3 7}
            a2a = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir), ...
                      'HRFindex','HRFindexrun');
  %            a8 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'R2');
            a2 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrfGLMdenoiseFracridge.mat',glmdir), ...    %OLD GLMdenoise_nsdBASICsingletrialfithrfGLMdenoiseRR.mat
                      'modelmd','R2','R2run','FRACvalue','scaleoffset');    % 'totalbadness'  % OLD: LAMBDAindex
            a2.HRFindex =    a2a.HRFindex;
            a2.HRFindexrun = a2a.HRFindexrun;
  %            xvalcurve = median(a2.totalbadness(a8.R2(:)>5,:),1);
            clear a2a;
          case 4
            a2 = load(sprintf('%s/restingGLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir), ...
                      'modelmd','R2','R2run');
          end
          mask0 = logical(ivolfun(getfield(load_untouch_nii(gunziptemp(sprintf('%s/ppdata/subj%02d/%s/brainmask.nii.gz', ...
                    nsddatadir,subjix,funcnames{ver}))),'img')));
          beta0 = int16(a2.modelmd * 300);                  % NOTE use integer
          beta0(repmat(~mask0,[1 1 1 size(beta0,4)])) = 0;  % NOTE zero-out
          meanbeta(:,:,:,ss) = mean(single(beta0),4);
          R2(:,:,:,ss) = single(a2.R2);
          r2ix = [r2ix ss];
        
        else
        
          if ~isrestingstate
            file0 = sprintf('%s/ppdata/subj%02d/%s/%s/meanbeta%s.nii.gz', ...
                            masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
            meanbeta(:,:,:,ss) = ivolfun(getfield(load_untouch_nii(file0),'img'));
          end

          file0 = sprintf('%s/ppdata/subj%02d/%s/%s/R2%s.nii.gz', ...
                          masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
          R2(:,:,:,ss) = ivolfun(getfield(load_untouch_nii(file0),'img'));
          r2ix = [r2ix ss];

        end

        % do actual costly saves if needed
        if sesstoexport(ss)

          % save single-trial betas nifti
          file0 = sprintf('%s/ppdata/subj%02d/%s/%s/betas%s.nii.gz', ...
                          masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
          nsd_savenifti(volfun(int16(beta0)),voxsize{ver},file0,tr(ver));

          % save single-trial betas hdf5
          file0 = sprintf('%s/ppdata/subj%02d/%s/%s/betas%s.hdf5', ...
                          masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
          betas = volfun2(beta0);
          delete(file0);
          h5create(file0,'/betas',size(betas),'Datatype','int16','ChunkSize',[1 1 1 size(betas,4)]);
          h5write(file0,'/betas',betas);
          %%save(file0,'betas','-nocompression','-v7.3');

          % clean up
          clear betas beta0;
        
          % R2
          file0 = sprintf('%s/ppdata/subj%02d/%s/%s/R2%s.nii.gz', ...
                          masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
          nsd_savenifti(volfun(single(a2.R2)),voxsize{ver},file0,tr(ver));               % FOR INCREMENTAL UPDATES, NEED TO LOAD SAVED FILE

          % R2run
          file0 = sprintf('%s/ppdata/subj%02d/%s/%s/R2run%s.nii.gz', ...
                          masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
          nsd_savenifti(volfun(single(a2.R2run)),voxsize{ver},file0,tr(ver));

          % resting-state does not get these
          if ~isrestingstate

            % meanbeta_session
            file0 = sprintf('%s/ppdata/subj%02d/%s/%s/meanbeta%s.nii.gz', ...
                            masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
            nsd_savenifti(volfun(single(meanbeta(:,:,:,ss))),voxsize{ver},file0,tr(ver));  % FOR INCREMENTAL UPDATES, NEED TO LOAD SAVED FILE

            % HRFindex
            if ismember(betaii,[2 3 5 6 7])
              file0 = sprintf('%s/ppdata/subj%02d/%s/%s/HRFindex%s.nii.gz', ...
                              masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
              nsd_savenifti(volfun(int16(a2.HRFindex)),voxsize{ver},file0,tr(ver));
            end

            % HRFindexrun
            if ismember(betaii,[2 3 5 6 7])
              file0 = sprintf('%s/ppdata/subj%02d/%s/%s/HRFindexrun%s.nii.gz', ...
                              masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
              nsd_savenifti(volfun(int16(a2.HRFindexrun)),voxsize{ver},file0,tr(ver));
            end

%             % LAMBDAindex
%             if ismember(betaii,[3 7])
%               file0 = sprintf('%s/ppdata/subj%02d/%s/%s/LAMBDAindex%s.nii.gz', ...
%                               masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
%               nsd_savenifti(volfun(int16(a2.LAMBDAindex)),voxsize{ver},file0,tr(ver));
%             end
            % FRACvalue
            if ismember(betaii,[3 7])
              file0 = sprintf('%s/ppdata/subj%02d/%s/%s/FRACvalue%s.nii.gz', ...
                              masterdir_betas,subjix,funcnames{ver},betaname,suffix0);
              nsd_savenifti(volfun(single(a2.FRACvalue)),voxsize{ver},file0,tr(ver));
            end

          end
        
        end

      end

      %% These are grand averages:
      
      if ~isrestingstate && ~isimagery && ~issynthetic && any(sesstoexport)

        % save mean across sessions of meanbeta
        file0 = sprintf('%s/ppdata/subj%02d/%s/%s/meanbeta.nii.gz', ...
                        masterdir_betas,subjix,funcnames{ver},betaname);
        nsd_savenifti(volfun(single(mean(meanbeta,4))),voxsize{ver},file0,tr(ver));

      end

      % save mean across sessions of R2
      if ~isempty(r2ix) && ~isimagery && ~issynthetic && any(sesstoexport)  % some subjects might not have any yet
        file0 = sprintf('%s/ppdata/subj%02d/%s/%s/R2.nii.gz', ...
                        masterdir_betas,subjix,funcnames{ver},betaname);
        nsd_savenifti(volfun(single(nanmean(R2(:,:,:,r2ix),4))),voxsize{ver},file0,tr(ver));
      end
      
      %% End.

    end

    % touch just to ensure that synchronization utilities (e.g. Dropbox) notice
    unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s/ppdata/subj%02d/%s/',masterdir_betas,subjix,funcnames{ver}));
    
    end

  end

end

%% %%%%%%% FUNCTIONAL DATA (prffloc timeseries)

% setup
nsdsetup;

% define
wantsubj = 1:8;
wantvers = [1 2];    % 1 is high res, 2 is low res

% constants
volfun = @(x) flipdim(flipdim(permute(x,[2 1 3 4 5 6 7 8 9 10]),2),1);
volfun2 = @(x) flipdim(rotatematrix(x,1,2,-1),1);
ivolfun = @(x) permute(flipdim(flipdim(x,1),2),[2 1 3 4 5 6 7 8 9 10]);

% more
ppdirnames = {'preprocessVER1SECONDSTAGE' 'preprocessVER1SECONDSTAGELOW'};
funcnames = {'func1mm' 'func1pt8mm'};
voxsize = {[1 1 1] [1.8 1.8 1.8]};
tr = [1 4/3];
suffix0 = '_prffloc';

% loop
for subjix=wantsubj
  datadir0 = nsdscreeningfun2(subjix);
  ppdir = regexprep(datadir0,'rawdata','ppdata');

  for ver=wantvers  % 1 is high res, 2 is low res

    %% deal with mean and valid

    a2 = load_untouch_nii(sprintf('%s/%s/mean.nii',ppdir,ppdirnames{ver}));  % stored as int16
    a3 = load_untouch_nii(sprintf('%s/%s/valid.nii',ppdir,ppdirnames{ver}));  % stored as int16 ?

    file0 = sprintf('%s/ppdata/subj%02d/%s/mean%s.nii.gz', ...
                    nsddatadir,subjix,funcnames{ver},suffix0);
    nsd_savenifti(volfun(int16(a2.img)),voxsize{ver},file0,tr(ver));
    file0 = sprintf('%s/ppdata/subj%02d/%s/valid%s.nii.gz', ...
                    nsddatadir,subjix,funcnames{ver},suffix0);
    nsd_savenifti(volfun(int16(a3.img)),voxsize{ver},file0,tr(ver));

    % touch just to ensure that synchronization utilities (e.g. Dropbox) notice
    unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s/ppdata/subj%02d/%s/',nsddatadir,subjix,funcnames{ver}));

    %% deal with timeseries

    for runii=1:1000

      tic;

      % load run
      runfile0 = sprintf('%s/%s/run%02d.nii',ppdir,ppdirnames{ver},runii);
      if ~exist(runfile0,'file')
        break;
      end
      a1 = load_untouch_nii(runfile0);  % stored as int16

      % get mask
      mask0 = logical(ivolfun(getfield(load_untouch_nii(gunziptemp(sprintf('%s/ppdata/subj%02d/%s/brainmask.nii.gz', ...
                nsddatadir,subjix,funcnames{ver}))),'img')));
                          
      % zero-out
      a1.img(repmat(~mask0,[1 1 1 size(a1.img,4)])) = 0;

      % save
      file0 = sprintf('%s/ppdata/subj%02d/%s/timeseries/timeseries%s_run%02d.nii.gz', ...
                      nsddatatimeseriesdir,subjix,funcnames{ver},suffix0,runii);
      nsd_savenifti(volfun(a1.img),voxsize{ver},file0,tr(ver));
      
      toc;

    end
    
  end
  
end

%% %%%%%%% FUNCTIONAL DATA TO FSAVERAGE

% see [analysis_transformfsaverage.m]

%% %%%%%%% FUNCTIONAL DATA TO MNI

% see [analysis_transformMNI.m]

%% %%%%%%% PERMISSIONS ISSUES

% what are our final opinions on permissions for raw data?

% Sep 4 2019:
find /home/surly-raid4/kendrick-data/nsd -type d -exec chmod g-w {} ";"
find /home/surly-raid3/kendrick-data/nsd -type d -exec chmod g-w {} ";"
find /home/surly-raid4/kendrick-data/nsd/nsddata_betas -type d -exec chmod g+w {} ";"

% Sep 17 2019:  (as root@surly)   [Hmm, this still leaves nsddata exposed...]
find /home/surly-raid1/kendrick-data/nsd -type f -exec chmod g-w {} ";"
find /home/surly-raid2/kendrick-data/nsd -type f -exec chmod g-w {} ";"
find /home/surly-raid3/kendrick-data/nsd -type f -exec chmod g-w {} ";"
find /home/surly-raid4/kendrick-data/nsd -type f -exec chmod g-w {} ";"
find /home/surly-raid4/kendrick-data/nsd/nsddata -type f -exec chmod g+w {} ";"
find /home/surly-raid4/kendrick-data/nsd/nsddata -type d -exec chmod g+w {} ";"

% Had to do this manually
find . -type f -user logan -delete

%% %%%%%%% CLEAN UP ISSUES

% Sep 5 2019 [weird ._ files]
[kendrick@surly ~/nsd/rawdata]$
rm */dicom/extra_scans/*localizer/._*

%% %%%%%%% DISK SPACE ISSUES

% RAN APPROXIMATELY AUG 20 2019 + OCT 17 2019.
% we want to remove these .nii files (since we export).
% but one difference is that the exported files have the brain mask applied! (oh well)

% get some candidates
ls -1 /home/surly-raid1/kendrick-data/nsd/ppdata/*-nsd??/preprocessVER1SECONDSTAGE/run*nii > ~/Dropbox/KKTEMP/markfordeletion.txt
ls -1 /home/surly-raid1/kendrick-data/nsd/ppdata/*-nsd??/preprocessVER1SECONDSTAGE/run*nii > ~/Dropbox/KKTEMP/markfordeletion2.txt
ls -1 /home/surly-raid1/kendrick-data/nsd/ppdata/*-nsd??/preprocessVER1SECONDSTAGE/run*nii > ~/Dropbox/KKTEMP/markfordeletion3.txt

% check the list, remove too recent ones, add rm, run as sh.

% archive the record at markfordeletion?.txt

%% %%%%%%% HRFs

% get HRFs for 1mm
a1 = load('~/nsd/ppdata/hrfmanifold.mat');
hrfs = a1.hrfs;  % 31 time points x 20 different timecourses
file0 = sprintf('%s/templates/hrfs_func1mm.mat',nsddatadir);
save(file0,'hrfs','-v7.3');

% get HRFs for 1pt8mm
a1 = load('~/nsd/ppdata/hrfmanifold_3pertrial.mat');
hrfs = a1.hrfs;  % 23 time points x 20 different timecourses
file0 = sprintf('%s/templates/hrfs_func1pt8mm.mat',nsddatadir);
save(file0,'hrfs','-v7.3');

%% %%%%%%% DESIGN MATRICES (NSD)

% setup
nsdsetup;

% define
masterdir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
ppdirnames = {'preprocessVER1SECONDSTAGEfigures' 'preprocessVER1SECONDSTAGELOWfigures'};
funcnames = {'func1mm' 'func1pt8mm'};  % 1-s, 4/3-s
suffixes = {'' '_3pertrial'};

% do it
for subjix=1:8, subjix
  for ver=1:length(funcnames), ver

    % get all nsd session directories (rawdata)
    sessions = nsdallsessfun2(subjix);    % paths
    sessionixs = nsdallsessfun3(subjix);  % positive integers
    length(sessions)
    
    % loop
    for ss=1:length(sessions)
      ppdir = regexprep(sessions{ss},'rawdata','ppdata');
      
      % load
      a1 = load(sprintf('%s/designmatrixFULL%s_nsd.mat',ppdir,suffixes{ver}),'stimulus');
      
      % deal
      assert(length(a1.stimulus)==12);
      allruns = {};
      for rr=1:length(a1.stimulus)           % for each run
        for qq=1:size(a1.stimulus{rr},1)     % for each TR
          ix = find(a1.stimulus{rr}(qq,:));
          if isempty(ix)
            ix = 0;
          end
          allruns{rr}(qq) = ix;
        end
      end
      
      % deal with resting-state
      if sessiontypes(sessionixs(ss))==8
        allruns{13} = zeros(1,length(allruns{1}));
        allruns{14} = zeros(1,length(allruns{1}));
        allruns = allruns([13 1:12 14]);
      end
    
      % write
      for rr=1:length(allruns)
        file0 = sprintf('%s/ppdata/subj%02d/%s/design/design_session%02d_run%02d.tsv', ...
                        nsddatatimeseriesdir,subjix,funcnames{ver},ss,rr);
        mkdirquiet(stripfile(file0));
        dlmwrite(file0,allruns{rr}(:),'delimiter','\t','precision',20);
      end
      
    end
  end
end

% - Note that trailing volumes of fMRI data must be dropped in order to match!

% SANITY CHECK (ONLY):
% subj 5, sess 32, first NSD run
a1 = load('/Users/kendrick/Dropbox/nsddata/experiments/nsd/nsd_expdesign.mat');
test = a1.subjectim(5,a1.masterordering(750*31+(1:63)));
a2=load('~/Dropbox/KKTEMP/1.tsv');
test2 = flatten(a2(a2~=0));
isequal(test,test2)   % check expdesign.mat 73k IDs is same as in .tsv
test3 = flatten(find(a2~=0));
test4 = find(upsamplematrix(flatten(a1.stimpattern(32,1,:)),4,2,[],0));
isequal(test3,test4)  % check trial onsets in tsv is same as in expdesign.mat

%% %%%%%%% DESIGN MATRICES (FLOC, NSDSYNTHETIC, NSDIMAGERY)

% history:
% - 2021/07/23 - rerun the nsdsynthetic to fix a bug. regenerated the .tsv files.

% ran 2, but revisit??? [separate by task??]
% 3 is a work in progress!!

% setup
nsdsetup;

% define
masterdir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
ppdirnames = {'preprocessVER1SECONDSTAGEfigures' 'preprocessVER1SECONDSTAGELOWfigures'};
funcnames = {'func1mm' 'func1pt8mm'};  % 1-s, 4/3-s
suffixes = {'' '_3pertrial'};
simpleids = {'floc' 'nsdsynthetic' 'nsdimagery'};

% do it
for ttt=1:2  %3%1:3  % floc, nsdsynthetic

for subjix=1:8, subjix
  for ver=1:length(funcnames), ver

    switch ttt
    case 1
      ppdir = regexprep(nsdscreeningfun2(subjix),'rawdata','ppdata');
      a1 = load(sprintf('%s/designmatrix%s_floc.mat',ppdir,suffixes{ver}),'stimulus');
      assert(length(a1.stimulus)==6);
    case 2
      ppdir = regexprep(datadirs{nsdsyntheticsession(subjix)},'rawdata','ppdata');
      a1 = load(sprintf('%s/designmatrixSINGLETRIAL%s_nsdsynthetic.mat',ppdir,suffixes{ver}));  % stimulus, stimorder, stimorderALT
      assert(length(a1.stimulus)==8);
    case 3
% 1-s preprocessing.
%       ppdir = regexprep(datadirs{nsdimagerysession(subjix)},'rawdata','ppdata');
%       a1 = load(sprintf('%s/designmatrixSINGLETRIAL%s_nsdimagery.mat',ppdir,suffixes{1}));  % NOTE: 1-s prep!! ?????
%       assert(length(a1.stimulus)==12);
    end
    
    % deal
    switch ttt
    case 1
      allruns = {};
      for rr=1:length(a1.stimulus)           % for each run
        for qq=1:size(a1.stimulus{rr},1)     % for each TR
          ix = find(a1.stimulus{rr}(qq,:));
          if isempty(ix)
            ix = 0;
          end
          allruns{rr}(qq) = ix;
        end
      end
    case 2
      allruns = {};
      cnt = 1;
      for rr=1:length(a1.stimulus)         % for each run
        for qq=1:size(a1.stimulus{rr},1)   % for each TR
          ix = find(a1.stimulus{rr}(qq,:));
          if isempty(ix)
            ix = 0;
          else
            ix = a1.stimorder(cnt);  % 1-284   (unique distinct images!)
            cnt = cnt + 1;
          end
          allruns{rr}(qq) = ix;
        end
      end
    end
    
    % write
    for rr=1:length(allruns)
      file0 = sprintf('%s/ppdata/subj%02d/%s/design/design_%s_run%02d.tsv', ...
                      nsddatatimeseriesdir,subjix,funcnames{ver},simpleids{ttt},rr);
      mkdirquiet(stripfile(file0));
      dlmwrite(file0,allruns{rr}(:),'delimiter','\t','precision',20);
    end
    
  end
end

end

%% %%%%%%% EXPORT PCNUM

% setup
nsdsetup;
glmdirs = {'glmdata' 'glmdataLOW'};

% loop
alldata = NaN*zeros(8,40,2);
alldataB = NaN*zeros(8,40,2,11);
for gg=1:length(glmdirs)
  for subjix=1:8
    todo = sessmatrix(subjix,:);
    for ii=1:length(todo)
      zz = todo(ii);
      if zz==0
        continue;
      end
      glmdir = regexprep(datadirs{zz},'rawdata',glmdirs{gg});
      a2 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrfGLMdenoiseFracridge.mat',glmdir),'pcnum','totalbadness');
      a8 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'R2');
      totalbadness = squish(a2.totalbadness,3);
      xvaltrend = -median(totalbadness(a8.R2(:)>5,:),1);  % NOTE: sign flip so that high is good
      alldata(subjix,ii,gg) = a2.pcnum;
      alldataB(subjix,ii,gg,:) = xvaltrend;
    end
  end
end

% write it out
funcdirs = {'func1mm' 'func1pt8mm'};
for p=1:2
  file0 = sprintf('~/nsd/nsddata/information/b3pcnum_%s.tsv',funcdirs{p});
  dlmwrite(file0,alldata(:,:,p)','delimiter','\t','precision',5,'-append');
end

%% %%%%%%% DEAL WITH PHYSIO DATA

% setup
nsdsetup;

% do it
beginslop = [];
for subjix=1:8, subjix

  % get all nsd session directories (rawdata)
  sessions = nsdallsessfun2(subjix);    % paths
  sessionixs = nsdallsessfun3(subjix);  % positive integers
  
  % loop
  for ss=21:30, ss   % NOTE: physio was recorded in nsd21-nsd30
  
    % check total number of files
    files0 = matchfiles(sprintf('%s/physio/*',sessions{ss}));
    assert(length(files0)==60,'A');
    
    % get the .puls files
    files0 = matchfiles(sprintf('%s/physio/*.puls',sessions{ss}));
    assert(length(files0)==15,'B');
    allphysiofiles = files0(2:end);

    % get the .resp files
    files0 = matchfiles(sprintf('%s/physio/*.resp',sessions{ss}));
    assert(length(files0)==15,'C');
    allphysiofiles = [allphysiofiles; files0(2:end)];
    
    % get the behavioral .mat files
    files1a = matchfiles(sprintf('%s/mat_files_from_scan/*_resting_*',sessions{ss}));
    files1b = matchfiles(sprintf('%s/mat_files_from_scan/*_nsd_*',sessions{ss}));
    assert(length(files1a)==2,'D');
    assert(length(files1b)==12,'E');
    allbehavfiles = [files1a(1); files1b; files1a(2)];
    
    % check dicom
    files2 = matchfiles(sprintf('%s/dicom/*cmrr_mbep2d_bold*',sessions{ss}));
    assert(length(files2)==14,'EE');
        
%     % read behavioral files
%     absfirstfmri = [];  % 1 x 14 with absolute times
%     for rr=1:length(allbehavfiles), fprintf('rr=%d,',rr);
%       a0 = load(allbehavfiles{rr},'timekeys');
%       assert(isequal(a0.timekeys{1,2},'absolutetimefor0'),'F');
%       assert(isequal(a0.timekeys{2,2},'trigger'),'G');
%       
%       % compute absolute time for first fMRI volume as detected by stimulus computer (seconds since midnight)
%       absfirstfmri(rr) = rem(a0.timekeys{1,1},1)*24*60*60 + a0.timekeys{2,1};
%     end
%     assert(all(diff(absfirstfmri) > 300),'H');
%     assert(all(diff(absfirstfmri) < 2000),'I');
    
    % read dicom files
    absfirstfmriALT = [];  % 1 x 14 with absolute times (THIS IS THE NEW VERSION)
    for rr=1:length(files2)
      files2a = matchfiles(sprintf('%s/*.dcm',files2{rr}));
      a33 = dicominfo(files2a{1});
      temp = a33.AcquisitionTime;
      temp2 = str2double(temp(1:2))*60*60 + ...
              str2double(temp(3:4))*60 + ...
              str2double(temp(5:end));
      absfirstfmriALT(rr) = temp2;  % seconds since midnight
%       fprintf('discrepancy = %.4f\n',absfirstfmriALT(rr) - absfirstfmri(rr));  % STRANGE: detection of 5 and AcquisitionTime is OFF.
    end

    % load and process each physio file
    for ff=1:length(allphysiofiles), fprintf('ff=%d,',ff);
      data = loadtext(allphysiofiles{ff});  % cell of strings
      fileext = allphysiofiles{ff}(end-3:end);  % 'puls' or 'resp'
      
      % sanity check
      assert(isequal(data{2},'ACQ FINISHED'),'J');
      
      % look for 5002 and go on
      temp = data{1};
      ix1 = regexp(temp,'5002'); assert(length(ix1)==2,'K');
      ix2 = regexp(temp,'6002'); assert(length(ix2)==1,'L');
      temp(ix1:(ix2+3))=[];
      ix1 = regexp(temp,'5002'); assert(length(ix1)==1,'M');
      temp = temp(1:ix1-1);
      
      % add more data
      temp2 = data{3};
      ix1 = regexp(temp2,'6002'); assert(length(ix1)==1,'N');
      ix2 = regexp(temp2,'5003'); assert(length(ix2)==1,'O');
      temp2 = temp2(ix1+4:ix2-1);
      
      % now we have all the numbers!
      valsA = eval(['[' temp ']']);
      valsB = eval(['[' temp2 ']']);
      valsA = valsA(~ismember(valsA,[5000 5002 5003 6002]));
      valsB = valsB(~ismember(valsB,[5000 5002 5003 6002]));

      % get rid of first four numbers
      if isequal(fileext,'puls')
        assert(isequal(valsA(1:4),[1 2 40 280]),'P');
      else
        assert(isequal(fileext,'resp'));
        assert(isequal(valsA(1:4),[1 2 20 2]),'Q');
      end
      valsA = valsA(5:end);
      totalsamplebeforeacq = length(valsA);
      vals = [valsA(:); valsB(:)];
      
      % calc start and stop times (seconds since midnight)
      temp = regexp(data{13},'LogStartMDHTime.+?(\d+)','tokens');
      start0 = str2double(temp{1}{1})/1000;
      temp = regexp(data{14},'LogStopMDHTime.+?(\d+)','tokens');
      stop0 = str2double(temp{1}{1})/1000;
      
      % sanity checks
      assert((stop0-start0) > 340 & (stop0-start0) < 345,'R');
      srate = (length(vals)-1) / (stop0-start0);  % samples per second
      assert(srate > 49.5 & srate < 49.8,'S');
      
      % figure
      if 0  % experimental...
        figure; hold on;
        time0 = linspace(start0,stop0,length(vals));  % assume equally spaced
        plot(time0,vals);
        ideal0 = absfirstfmriALT(ff) + linspacefixeddiff(0,1.606425,188);
        straightline(ideal0,'v','r-');
      end
      
      % more sanity checks
      beginslop(subjix,ss-20,ff,1) = absfirstfmriALT(mod2(ff,14)) - start0;                  % fMRI-start minus physio-start
      beginslop(subjix,ss-20,ff,2) = stop0 - (absfirstfmriALT(mod2(ff,14)) + 188*1.606425);  % physio-end minus idealized-fMRI-end
      beginslop(subjix,ss-20,ff,3) = stop0 - start0;                                         % physio-end minus physio-start

      % final check
      check1 = absfirstfmriALT(mod2(ff,14)) + 188*1.606425;   % acqtime + idealized time
      check2 = start0 + (1/srate)*totalsamplebeforeacq;       % sample-based stop
      beginslop(subjix,ss-20,ff,4) = check2 - check1;         % [physio's end time] - [dicom's guessed end time] in seconds
                                                              %       this checks: physio's self timing and ACQ detection
                                                              %                    against the DICOM header information

      % report
      fprintf('\nbegin = %.10f, end = %.10f, dur = %.10f, chk = %.4f\n', ...
                                                            beginslop(subjix,ss-20,ff,1), ...
                                                            beginslop(subjix,ss-20,ff,2), ...
                                                            beginslop(subjix,ss-20,ff,3), ...
                                                            beginslop(subjix,ss-20,ff,4));
      % the MDH duration is a little variable already, which is weird slop.
      % however, the begin and end times are quite good, suggesting that we are interpreting things correctly!!

      % pick the samples to include
      time0 = linspace(start0,stop0,length(vals));  % assume equally spaced!!
      ix1 = firstel(find(absfirstfmriALT(mod2(ff,14)) - time0 < 0));                    % first sample AFTER fmri starts
      ix2 = firstel(find(absfirstfmriALT(mod2(ff,14)) + 188*1.606425 - time0 < 0));     % first sample AFTER fmri ends
      vals = vals(ix1:ix2);    % crop the physio data
      time0 = time0(ix1:ix2);  % crop the time0 times
      
      % save files
      file0 = sprintf('%s/ppdata/subj%02d/physio/physio_session%02d_run%02d_%s.tsv', ...
                      nsddatatimeseriesdir,subjix,ss,mod2(ff,14),fileext);
      mkdirquiet(stripfile(file0));
         % REMOVE!!! vflatten(time0([1 end]) - absfirstfmriALT(mod2(ff,14))); 
         % note that after subtraction, the units are seconds and t=0 corresponds to the AcquisitionTime of the first DICOM
      dlmwrite(file0,[vals(:)],'delimiter','\t','precision',20);
      % new interpretation: the start and stop in the file basically corresponds to the actual total fmri acquisition

    end
    
  end
end

% final sanity plots (see nsd/physio)
figureprep;hist(flatten(beginslop(:,:,1:14,1)),100);figurewrite('puls1',[],[],'~/Dropbox/KKTEMP');    % puls, begin period: 40.15 s - 40.28 s
figureprep;hist(flatten(beginslop(:,:,15:end,1)),100);figurewrite('resp1',[],[],'~/Dropbox/KKTEMP');  % resp, begin period: 40.12 s - 40.26 s
figureprep;hist(flatten(beginslop(:,:,1:14,2)),100);figurewrite('puls2',[],[],'~/Dropbox/KKTEMP');    % puls, end extra:     1.28 s - 1.315 s
figureprep;hist(flatten(beginslop(:,:,15:end,2)),100);figurewrite('resp2',[],[],'~/Dropbox/KKTEMP');  % resp, end extra:     1.78 s - 1.82  s
figureprep;hist(flatten(beginslop(:,:,1:14,3)),100);figurewrite('puls3',[],[],'~/Dropbox/KKTEMP');    % puls, duration:      342.24 - 342.38 s
figureprep;hist(flatten(beginslop(:,:,15:end,3)),100);figurewrite('resp3',[],[],'~/Dropbox/KKTEMP');  % resp, duration:      342.7  - 342.86 s
figureprep;hist(flatten(beginslop(:,:,1:14,4)),100);figurewrite('puls4',[],[],'~/Dropbox/KKTEMP');    % puls, check:      
figureprep;hist(flatten(beginslop(:,:,15:end,4)),100);figurewrite('resp4',[],[],'~/Dropbox/KKTEMP');  % resp, check:
save('~/Dropbox/KKTEMP/physiochecking.mat','beginslop');

% number of samples in the files I write are:  (obsolete)
% min 14941
% max 14982
% 41 is about 0.8 s.
% not sure why.

% some rethinking from a sample file:
%
% ext is 343.170 s long.
% dicom acquisitiontime is 63281.5975 seconds since midnight
% log started at 63241.532
% time = linspace(63241.532,63584.702,68359);  % whole-file interpretation
% 68155 samples until the ACQ.
% samplerate = 343.170 / 68359;    % seconds per sample (whole-file interpretation)
% 63281.5975 + 188*1.606425 = 63583.6054 s   is the supposed acquisition stopping.
% 63241.532 + samplerate*68155 = 63583.6778966632 s is the sample-based stop.
% 
% puls is 342.262 s long.
% dicom acquisitiontime is 63281.5975 seconds since midnight
% log started at 63241.435
% time = linspace(63241.435,63583.697,17490);  % whole-file interpretation
% 17484 samples until the ACQ.
% samplerate = 342.262 / 17490;    % seconds per sample (whole-file interpretation)
% 63281.5975 + 188*1.606425 = 63583.6054 s   is the supposed acquisition stopping.
% 63241.435 + samplerate*17484 = 63583.5795859348 s is the sample-based stop.

%%% empirically measure dicoms

a = matchfiles('*cmrr*bold*');
rec = [];
for p=1:length(a)
  b = matchfiles([a{p} '/*dcm']);
  c = dicominfo(b{1});
  d = dicominfo(b{end});
  e1 = str2double(c.AcquisitionTime(1:2))*60*60 + str2double(c.AcquisitionTime(3:4))*60 + str2double(c.AcquisitionTime(5:end));
  e2 = str2double(d.AcquisitionTime(1:2))*60*60 + str2double(d.AcquisitionTime(3:4))*60 + str2double(d.AcquisitionTime(5:end));
  rec(p) = (e2-e1)/187;
end
mean(rec)

1.606425  * 188  (mean)
1.6064267        (another mean)
1.60642475       (another mean)
1.60642491       (another mean)

final estimate:
mean([1.606425 1.6064267 1.60642475 1.60642491])
***1.606425****   this is the official number (according to scanner how long TRs are)

%% %%%%%%%% Do setup

nsdsetup;

%% %%%%%%%% Use nsd_mapdata to get betas into MNI space

% Ran on Aug 18 2019 (1-224) [includes restingbetas and going forward]
% Update on Sep 3 2019 (225-239).
% DEC 29 2019: update and run for nsdimagery√
% DEC 29 2019: update and run for nsdsynthetic

% define
nsubj = 8;
betadirs = {'betas_fithrf' 'restingbetas_fithrf' 'nsdimagerybetas_fithrf' 'nsdsyntheticbetas_fithrf'};  % only do these versions; don't do 'betas_assumehrf' 'betas_fithrf_GLMdenoise_RR' 'nsdsyntheticbetas_fithrf_GLMdenoise_RR'

% do it
for subjix=1:8, subjix
  for bb=1:length(betadirs), bb
    for sess=1:40, sess

      if bb==3
        if sess > 1
          continue;
        end
        suffix0 = '_nsdimagery';
      elseif bb==4
        if sess > 1
          continue;
        end
        suffix0 = '_nsdsynthetic';
      else
        % if sess is more than we've done OR it has not been collected OR it's not in the range we want to do
        if sess > size(sessmatrix,2) || sessmatrix(subjix,sess)==0   %%%|| sessmatrix(subjix,sess) < 225 || sessmatrix(subjix,sess) > 239
          continue;
        end
        suffix0 = sprintf('_session%02d',sess);
      end

      % calc
      dir0 = sprintf('%s/ppdata/subj%02d/func1mm/%s',nsd_datalocation('betas'),subjix,betadirs{bb});        % 1mm beta volume directory
      file0 = sprintf('%s/betas%s.mat',dir0,suffix0);                                                       % 1mm beta file
      file1 = sprintf('%s/ppdata/subj%02d/func1mm/valid%s.nii.gz',nsd_datalocation,subjix,suffix0);         % 1mm valid file
      file2 = sprintf('%s/ppdata/subj%02d/func1mm/brainmask.nii.gz',nsd_datalocation,subjix);               % 1mm brain mask
                
      % if session doesn't exist, abort
      if ~exist(file0,'file')
        fprintf('warning %d %d\n',subjix,bb);
        continue;
      end
      
      % load data
      a1 = load(file0);  % betas
      a2 = load_untouch_nii(file1);  % valid
      mask0 = logical(getfield(load_untouch_nii(file2),'img'));  % brainmask
      
      % prep data
      a1.betas = single(a1.betas);                                         % convert to single so that we can use NaN. Note: still in *300 format
      a1.betas(repmat(~a2.img | ~mask0,[1 1 1 size(a1.betas,4)])) = NaN;   % set any voxel that is invalid OR outside-brainmask to NaN
      
      % transform to MNI space
      outputfile = sprintf('%s/ppdata/subj%02d/MNI/%s/betas%s.nii.gz', ...
        nsd_datalocation('betas'),subjix,betadirs{bb},suffix0);
      mkdirquiet(stripfile(outputfile));
      data1 = nsd_mapdata(subjix,'func1pt0','MNI',a1.betas,'cubic',NaN,outputfile,'int16');  % NOTE: in LPI, and NaNs are saved as 0s, we save as *300 integers

      % here we "reverse compute" valid
      valid0 = ~all(data1==0,4);  % valid means not all values are 0
      outputfile = sprintf('%s/ppdata/subj%02d/MNI/%s/valid%s.nii.gz', ...
        nsd_datalocation('betas'),subjix,betadirs{bb},suffix0);
      nsd_savenifti(int16(valid0),[1 1 1],outputfile,[],[183-91 127 73]);  % NOTE THIS HACK AND HARD-CODED CONSTANTS
      
    end
  end
end

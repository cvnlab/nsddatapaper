%% %%%%%%%% Do setup

nsdsetup;

%% %%%%%%%% Use nsd_mapdata to get betas onto fsaverage

% Ran on July 12 2019.
% Finished on July 14 2019.
%   And caught up resting-state to that point on Aug 5 2019 (that's through 197).
% Freshened on Aug 19 2019 [processed 198-224].
% Freshened on Sep 4 2019 [processed 225-239].
% DEC 29 2019: update and run for nsdimagery√
% DEC 29 2019: update and run for nsdsynthetic
% FEB 10 2020: add in saving of nativesurface and re-run for that component (hdf5 format).
% APR 19 2020: RE-RUN ALL nativesurface and fsaverage due to FS update
% May 16 2020: RE-RUN [3 7] beta versions because of the b3 update.
% Mar  2 2024: RUN 8 for nsdimagery b3.
%
% Notes:
% - It is possible that some vertices have NaNs for all of their 750 betas.

% define
hemis = {'lh' 'rh'};
nsubj = 8;
ndepth = 3;
betadirs = {'betas_assumehrf' 'betas_fithrf' 'betas_fithrf_GLMdenoise_RR' 'restingbetas_fithrf' 'nsdimagerybetas_fithrf' 'nsdsyntheticbetas_fithrf' 'nsdsyntheticbetas_fithrf_GLMdenoise_RR' 'nsdimagerybetas_fithrf_GLMdenoise_RR'};
fsavgdir = [nsd_datalocation '/freesurfer/fsaverage'];
fsdir = [nsd_datalocation '/freesurfer/subj%02d'];

% do it
for subjix=1:8, subjix
  for bb=1:length(betadirs), bb
    for sess=1:40, sess
    
      if ismember(bb,[5 8])
        if sess > 1
          continue;
        end
        suffix0 = '_nsdimagery';
      elseif ismember(bb,[6 7])
        if sess > 1
          continue;
        end
        suffix0 = '_nsdsynthetic';
      else
        % if sess is more than we've done OR it has not been collected OR it's not in the range we want to do
        if sess > size(sessmatrix,2) || sessmatrix(subjix,sess)==0     % || sessmatrix(subjix,sess) < 225 || sessmatrix(subjix,sess) > 239
          continue;
        end
        suffix0 = sprintf('_session%02d',sess);
      end

      % calc
      dir0 = sprintf('%s/ppdata/subj%02d/func1mm/%s',nsd_datalocation('betas'),subjix,betadirs{bb});  % 1mm beta volume directory
      file0 = sprintf('%s/betas%s.hdf5',dir0,suffix0);                                                % 1mm beta file
      file1 = sprintf('%s/ppdata/subj%02d/func1mm/valid%s.nii.gz',nsd_datalocation,subjix,suffix0);   % 1mm valid file

      % if session doesn't exist, abort
      if ~exist(file0,'file')
        fprintf('warning %d %d\n',subjix,bb);
        continue;
      end
      
      % load data
% OLD FORMAT:
%       a1 = load(file0);
%       a1.betas = single(a1.betas)/300;                            % convert to PSC
%       a1.betas(repmat(~a2.img,[1 1 1 size(a1.betas,4)])) = NaN;   % ensure invalid betas are NaN
      a2 = load_untouch_nii(file1);
      betas = h5read(file0,'/betas',[1 1 1 1],[Inf Inf Inf Inf]);
      betas = single(betas)/300;                            % convert to PSC
      betas(repmat(~a2.img,[1 1 1 size(betas,4)])) = NaN;   % ensure invalid betas are NaN
      
      % for each hemisphere, write out an fsaverage .mgh file with the betas
      for hh=1:length(hemis)

        % map via cubic interpolation to each depth and then average across depth
        data1 = single([]);  % V x 750 x 3depth
        for ii=1:ndepth
          data1(:,:,ii) = nsd_mapdata(subjix,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),betas,'cubic',NaN,[],'single');
        end
        data1 = mean(data1,3);  % NaNs propagate!
        
        % write a nativesurface version
        outputfile = sprintf('%s/ppdata/subj%02d/nativesurface/%s/%s.betas%s.hdf5', ...
          nsd_datalocation('betas'),subjix,betadirs{bb},hemis{hh},suffix0);
        mkdirquiet(stripfile(outputfile));
        rmdirquiet(outputfile);
        h5create(outputfile,'/betas',size(data1),'Datatype','int16','ChunkSize',[1 size(data1,2)]);
        h5write(outputfile,'/betas',int16(data1*300));  % note that NaNs become 0

        % map via nearest neighbor to fsaverage and then write file
        outputfile = sprintf('%s/ppdata/subj%02d/fsaverage/%s/%s.betas%s.mgh', ...
          nsd_datalocation('betas'),subjix,betadirs{bb},hemis{hh},suffix0);
        mkdirquiet(stripfile(outputfile));
        nsd_mapdata(subjix,sprintf('%s.white',hemis{hh}),'fsaverage',data1,[],NaN,outputfile,'single',fsavgdir);

      end
      
    end
  end
end

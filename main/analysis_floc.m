%% %%%%%%%% Perform GLM analysis of floc data

% do the analysis
for zz=1:8
  glm_floc(zz,'preprocessVER1SECONDSTAGE',   1,  [3 4 7 8 11 12],'',          'glmdata');
  glm_floc(zz,'preprocessVER1SECONDSTAGELOW',4/3,[3 4 7 8 11 12],'_3pertrial','glmdataLOW');
end

%% %%%%%%%% Load results and write out official .nii.gz results [this replaces export.m!]

% run on surly is okay?
nsdsetup;

% define some constants
volfun = @(x) flipdim(flipdim(permute(x,[2 1 3 4 5 6 7 8 9 10]),2),1);
voxsize = {[1 1 1] [1.8 1.8 1.8]};
tr = [1 4/3];
glmdirnames = {'glmdata' 'glmdataLOW'};
varstoload = {'R2' 'meanvol' 'pcR2final' 'modelmd' 'modelmd'};  % 5 is special t-values and angle-values
filenames = {'R2' 'meanvol' 'R2xval' 'betas' []};
funcnames = {'func1mm' 'func1pt8mm'};

% floc
stimix = {};
for p=1:5
  stimix{p} = {};
  stimix{p}{1} = (p-1)*2+[1 2];
  stimix{p}{2} = setdiff(1:10,stimix{p}{1});
end
for p=1:10
  stimix{5+p} = {p setdiff(1:10,(ceil(p/2)-1)*2 + (1:2))};
end
stimlabels = {'characters' 'bodies' 'faces' 'places' 'objects' ...
              'word' 'number' 'body' 'limb' 'adult' 'child' 'corridor' 'house' 'car' 'instrument'};

% load results and write niftis
for ver=1:length(voxsize)
  for subjix=1:8

    % calc
    ppdir = regexprep(nsdscreeningfun2(subjix),'rawdata','ppdata');
    glmdir = regexprep(ppdir,'ppdata',glmdirnames{ver});

    % load brainmask (it is in LPI; we are in ARI)
    a1 = load_untouch_nii(gunziptemp(sprintf('%s/ppdata/subj%02d/%s/brainmask.nii.gz',nsddatadir,subjix,funcnames{ver})));
    brainmask = volfun(double(a1.img));
    clear a1;

    % load all results
    a1 = load(sprintf('%s/GLMdenoise_flocFULL.mat',glmdir));
    
    % process each quantity
    for qq=1:length(varstoload)
      data0 = a1.(varstoload{qq});

      % notes:
      % - R2.  some are nan (missing data).
      % - meanvol is well behaved, no nans.
      % - pcR2final range is well behaved, but some are nans.
      % - modelmd has extreme! values. and also NaNs.
      % - the t-values will have NaNs propagating

      if qq==1
        ensurebad = isnan(data0);
      end

      if qq==5
        filenamestowrite = {};
        datatowrite = {};        
        for ss=1:length(stimix)

          A = stimix{ss}{1};
          B = stimix{ss}{2};
  
          beta = reshape(data0,[sizefull(data0,3) 6 10]);
          mainbeta =  mean(mean(beta(:,:,:,:,A),4),5);  % mean of all trials in A
          otherbeta = mean(mean(beta(:,:,:,:,B),4),5);  % mean of all trials in B
  
          temp = squish(permute(beta(:,:,:,:,A),[4 5 1 2 3]),2);
          mainbetase = permute(std(temp,[],1)/sqrt(size(temp,1)),[2 3 4 1]);

          temp = squish(permute(beta(:,:,:,:,B),[4 5 1 2 3]),2);
          otherbetase = permute(std(temp,[],1)/sqrt(size(temp,1)),[2 3 4 1]);

          vals = (mainbeta-otherbeta) ./ sqrt(mainbetase.^2 + otherbetase.^2);

          filenamestowrite{end+1} = [stimlabels{ss} 'tval'];
          datatowrite{end+1} = vals;

          % A (mainbeta) and B (otherbeta) are the x- and y-axes.
          % atan2 makes it such that 0 is A preference and pi/2 is B preference.
          % after circulardiff, pi/4 is A preference and -pi/4 is B preference.
          % after conversion, A preference is from 0 to 180 and B preference is from 0 to -180.
          vals = circulardiff(pi/4,atan2(otherbeta,mainbeta),2*pi)/pi*180;

          filenamestowrite{end+1} = [stimlabels{ss} 'anglemetric'];
          datatowrite{end+1} = vals;

        end
      else
        filenamestowrite = {filenames{qq}};
        datatowrite = {data0};
      end

      % write it out
      for dd=1:length(filenamestowrite)
        % enforce nans [NaN means either we are hiding (outside of brainmask) OR missing data due to coverage]
        datatowrite{dd}(~logical(brainmask) | ensurebad) = NaN;
        file0 = sprintf('%s/ppdata/subj%02d/%s/floc_%s.nii.gz', ...
                        nsddatadir,subjix,funcnames{ver},filenamestowrite{dd});
        nsd_savenifti(volfun(single(datatowrite{dd})),voxsize{ver},file0,tr(ver));
        unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s',stripfile(file0)));
      end

    end
    
  end
end

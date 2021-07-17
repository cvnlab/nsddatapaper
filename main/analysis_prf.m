%% %%%%%%%% Perform prf analysis

% NOTE: we want to use the vanilla analyzePRF, so we do this and check our version
rmpath(genpath('/home/stone/kendrick/code/analyzePRF/'));
which analyzePRF -all

% using parpool?
parpool(16);
parpool(12);
parpool(8);

% do the FULL analysis [possibly farm out to different machines!]
for zz=1:8
  glm_prf(zz,'preprocessVER1SECONDSTAGE',   1,  [1 2 5 6 9 10],'func1mm',   'small',      'glmdata');
  glm_prf(zz,'preprocessVER1SECONDSTAGELOW',4/3,[1 2 5 6 9 10],'func1pt8mm','small4div3s','glmdataLOW');
end

%%% interlude

% see how consistent findtailthreshold is
for p=1:8
  file0 = sprintf('/home/stone-ext4/generic/Dropbox/nsddata/ppdata/subj%02d/func1pt8mm/prf_R2.nii.gz',p);
  a1 = load_untouch_nii(file0);
  vals(p) = findtailthreshold(a1.img(:));
  sum(a1.img(:)>10.1)/numel(a1.img)*100   %vals: NaN -1000
end
% 10.1%

%%% proceed [EXPERIMENTAL, NOT FOR RELEASE]

% do the SEP analysis [possibly farm out to different machines!]
  % surly and stone
for zz=2:2:8
  glm_prfSEP(zz,'preprocessVER1SECONDSTAGELOW',4/3,[1 2 5 6 9 10],'func1pt8mm','small4div3s','glmdataLOW', 10.1);
end
for zz=2:2:8
  glm_prfSEP(zz,'preprocessVER1SECONDSTAGE',   1,  [1 2 5 6 9 10],'func1mm',   'small',      'glmdata',    10.1);
end

%% %%%%%%%% Load in the files and write out official .nii.gz results [this replaces export.m!]

% run on stone!
nsdsetup;

% define some constants
volfun = @(x) flipdim(flipdim(permute(x,[2 1 3 4 5 6 7 8 9 10]),2),1);
voxsize = {[1 1 1] [1.8 1.8 1.8]};
tr = [1 4/3];
glmdirnames = {'glmdata' 'glmdataLOW'};
varstoload = {'R2' 'ecc' 'ang' 'expt' 'gain' 'meanvol' 'rfsize'};
filenames = {'R2' 'eccentricity' 'angle' 'exponent' 'gain' 'meanvol' 'size'};
funcnames = {'func1mm' 'func1pt8mm'};
pxtodeg = 8.4/200;

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
    
    % load each quantity
    alldata = {};
    fprintf('\n\nsubjix = %d, %s\n\n',subjix,funcnames{ver});
    for qq=1:length(varstoload)
      data0 = loadmulti(sprintf('%s/analyzePRF_prfFULL_chunk????.mat',glmdir),varstoload{qq},1);
      alldata{qq} = data0;

      % check
      fprintf('\n**** var is %s\n',varstoload{qq});
      fprintf('min is %.4f\n',min(data0));
      fprintf('max is %.4f\n',max(data0));
      fprintf('number that are nan? %d\n',sum(isnan(data0)));
      fprintf('number that are not finite? %d\n',sum(~isfinite(data0)));
      fprintf('number that are less than 0? %d\n',sum(data0 < 0));
      fprintf('number that are equal to 0? %d\n',sum(data0 == 0));
      fprintf('prctile from 5 to 95? %s\n',mat2str(prctile(data0,[5 95])));
    end

    % process each quantity
    for qq=1:length(varstoload)
    
      % massage
      switch qq
      case 1      % R2 is 0-100%, roughly. However, there are NaNs due to missing data!
        ensurebad = isnan(alldata{qq});   % keep track of missing data voxels!
        alldata{qq}(alldata{qq} < -1000) = -1000;
      case 2      % ecc is >=0 to a really large number. no NaNs!
        ecczero = alldata{qq}==0;
        alldata{qq} = alldata{qq} * pxtodeg;  % convert to dva
        alldata{qq}(alldata{qq} > 1000) = 1000;
      case 3      % ang is 0-360. no NaNs!
        alldata{qq}(ecczero) = NaN;  % force illdefined angles to NaN
      case 4      % expt is >=0 to really large. no NaNs!
        exptzero = alldata{qq}==0;   % when is there infinite compression?
        alldata{qq}(alldata{qq} > 1000) = 1000;
      case 5      % gain is >=0 to really large. no NaNs!
        assert(sum(alldata{6}==0 & ~ensurebad)==0);          % confirm that no meanvol is exactly 0 (for non-bad)
        alldata{qq} = alldata{qq} ./ abs(alldata{6}) * 100;  % convert to PSC
        alldata{qq}(alldata{qq} > 1000) = 1000;
      case 6      % meanvol is -45 to 3193... no NaNs.
      case 7      % size is >0 to Inf.  no NaNs (but Inf).
        assert(sum(alldata{qq} == 0) == 0);   % confirm none are exactly zero
        alldata{qq} = alldata{qq} * pxtodeg;  % convert to dva
        alldata{qq}(exptzero) = 1000;   % don't allow Inf, but just set to 1000
        alldata{qq}(alldata{qq} > 1000) = 1000;
        assert(all(isfinite(alldata{qq})));
      end
      
      % reconstruct volume [NaN means either not analyzed (outside of brainmask) OR missing data due to coverage]
      data0 = copymatrix(NaN*ones(size(brainmask)),logical(brainmask),copymatrix(alldata{qq},ensurebad,NaN));

      % write it out
      file0 = sprintf('%s/ppdata/subj%02d/%s/prf_%s.nii.gz', ...
                      nsddatadir,subjix,funcnames{ver},filenames{qq});
      nsd_savenifti(volfun(single(data0)),voxsize{ver},file0,tr(ver));

    end
    
  end
end

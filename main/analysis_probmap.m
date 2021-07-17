% Here, we create some simple "probabilistic ROIs"
% whereby we transfer manually defined ROIs in each of the 8 NSD 
% subjects to fsaverage and then compute the fraction of subjects
% that are labeled at fsaverage vertex.
%
% We save both fsaverage .mgz files and also create some visualizations.

%% %%%%% Compute and visualize probmaps

% setup
setenv('SUBJECTS_DIR','/home/surly-raid3/kendrick-data/nsd/nsddata/freesurfer');
allrois = {'prf-visualrois' 'prf-eccrois' 'floc-faces' 'floc-words' 'floc-places' 'floc-bodies'};
hemis = {'lh' 'rh'};

% get ROI names
roinames = {};
for p=1:length(allrois)
  file0 = sprintf('~/Dropbox/nsddata/freesurfer/subj01/label/%s.mgz.ctab',allrois{p});
  fid = fopen(file0);
  C = textscan(fid,'%d%s');%,'Delimiter',' ');
  fclose(fid);
  assert(isequal(double(C{1}),(1:length(C{1}))'-1));
  roinames{p} = C{2}(2:end);
end
allroinames = cat(1,roinames{:});

% transfer ROIs to fsaverage
vals = zeros(163842,2,0);  % 163842 vertices x 2 hemis x ROIs with the final probmap values
cnt = 1;
Lookup = [];
for ff=1:length(allrois)   % for each ROI group

  % use nearest neighbor to map to fsaverage
  temp = [];    % 163842 vertices x 8 subjects x 2 hemis with the individual subject ROI integers
  for hh=1:length(hemis)       % for each hemisphere
    for subjix=1:8
      a5 = cvnloadmgz(sprintf('~/Dropbox/nsddata/freesurfer/subj%02d/label/%s.%s.mgz',subjix,hemis{hh},allrois{ff}));
      assert(~isempty(a5));
      temp(:,subjix,hh) = nsd_mapdata(subjix,sprintf('%s.white',hemis{hh}),'fsaverage',a5(:));
    end
  end
    
  % for each ROI, compute fractions of subjects
  for rr=1:length(roinames{ff})

    for hh=1:length(hemis)
  
      % compute fraction and record
      temp2 = mean(temp(:,:,hh)==rr,2);
      vals(:,hh,cnt) = temp2;
    
      % save to disk
      nsd_savemgz(temp2,sprintf('~/Dropbox/nsddata/freesurfer/fsaverage/label/%s.probmap_%s.mgz',hemis{hh},roinames{ff}{rr}), ...
                  '~/Dropbox/nsddata/freesurfer/fsaverage');

    end

    % visualize surface map [0 to 1]
    extraopts = {'rgbnan',1,'hemibordercolor',[1 1 1],'scalebarcolor','k'};  %% NOTE
    [rawimg,Lookup,rgbimg,himg] = cvnlookup('fsaverage',13,vflatten(vals(:,:,cnt)),[0 1],copper(256),1/8/2,[],0,extraopts);
    imwrite(rgbimg,sprintf('~/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_probmap_%s.png',roinames{ff}{rr}));
    
    % increment
    cnt = cnt + 1;
  
  end

end

% notes: inherit some from examples_fullanalysis2.m

%%%%%%%%%%%%%%%%%%%%%%%%% Analyze, looping over subjects.

% define
totalnsdsessions = [40 40 32 30 40 32 40 30];
hemis = {'lh' 'rh'};

% do it
alltvals = {};
for subjix=1:8
  nn = totalnsdsessions(subjix);

  % load betas (b2, in PSC units)
  alldata = zeros(327684,750*nn,'single');
  for sess=1:nn, sess
    data = cvnloadmgz(sprintf('~/nsddata_betas/ppdata/subj%02d/fsaverage/betas_fithrf_GLMdenoise_RR/*.betas_session%02d.mgh',subjix,sess));  % 327684 x 1 x 1 x 750
    alldata(:,(sess-1)*750+(1:750)) = permute(data,[1 4 2 3]);
  end
  clear data;

  % load behavioral data
  fid = fopen(sprintf('~/nsddata/ppdata/subj%02d/behav/responses.tsv',subjix));
  C = textscan(fid,repmat('%f',[1 19]),'Delimiter','\t','HeaderLines',1);
  fclose(fid);
  
  % check
  assert(size(C{1},1)==size(alldata,2));

  % compute t-test for hits > correct rejections
  tvals = [];  % 327k x (N+1pool+1permute) sessions
  todo = [num2cell(1:nn) {1:nn} {[]}];  % individual sessions as well as pooled-across-sessions
  for ss=1:length(todo), ss
  
    % handle special permute case
    isspecial = 0;
    if isempty(todo{ss})
      isspecial = 1;
      todo{ss} = 1:nn;
    end
  
    % match subject, match session, key on behavior
    A = find(C{1}==subjix & ismember(C{2},todo{ss}) & C{8}==1 & C{9}==1);  % old + correct = hit
    B = find(C{1}==subjix & ismember(C{2},todo{ss}) & C{8}==0 & C{9}==1);  % new + correct = correct rejection
    
    % handle permutation control
    if isspecial
      temp = permutedim([A; B]);
      A = temp(1:length(A));
      B = temp(length(A)+1:end);
    end

    % compute t-value
    mn1 = mean(alldata(:,A),2);
    se1 = std(alldata(:,A),[],2)/sqrt(length(A));
    mn2 = mean(alldata(:,B),2);
    se2 = std(alldata(:,B),[],2)/sqrt(length(B));
    tvals(:,ss) = (mn1 - mn2) ./ sqrt(se1.^2 + se2.^2);

  end
  
  % record
  alltvals{subjix} = tvals;

end

% in rare cases, no data exist in a given session. set these t-values to 0.
alltvals = cellfun(@(x) nanreplace(x,0,3),alltvals,'UniformOutput',0);

%%%%%%%%%%%%%%%%%%%%%%%%% Make maps

% define and prep
Lookup = [];
eo = {'hemibordercolor' 'w' 'rgbnan' 1 'drawscalebar' false};
mkdirquiet('~/Dropbox/KKTEMP/maps');

% loop
grandrec = [];
for subjix=1:8+1
  if subjix<=8
    whsessions = [5:5:20 totalnsdsessions(subjix)+(1:2)];  % individual sessions, pooled, and control
  end
  for ss=1:length(whsessions)
    if subjix==9
      temp = mean(grandrec(:,ss,:),3);  % simple mean across subjects
    else
      temp = alltvals{subjix}(:,whsessions(ss));
      grandrec(:,ss,subjix) = temp;
    end
    if ss<=20/5
      rng0 = [-8 8];
    else
      rng0 = [-40 40];
    end
    cvnlookup('fsaverage',13,temp,rng0,cmapsign4(256),3*j,Lookup,0,eo);
    imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/maps/subj%02d_ss%02d.png',subjix,ss));

    if subjix==9
      cvnlookup('fsaverage',{'medial' 'inflated' 0 1000 0 [1 1]},temp,rng0,cmapsign4(256),3*j,[],0,eo);
      imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/maps/INFLATED_subj%02d_medial_ss%02d.png',subjix,ss));
      
      cvnlookup('fsaverage',{'lateral' 'inflated' 0 1000 0 [1 1]},temp,rng0,cmapsign4(256),3*j,[],0,eo);
      imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/maps/INFLATED_subj%02d_lateral_ss%02d.png',subjix,ss));
    end

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Crop maps (LH only)

% duplicate maps to mapscrop!

% cd to mapscrop directory.

% crop
thecrop = imcropfiles('*.png',[1 1500 1 1662]);

% consolidate images into a large row (for each subject)
for p=1:8+1
  imwrite(concatimages(sprintf('subj%02d*.png',p),{[1 2 3 4 5 6]},0),sprintf('../allsubj%02d.png',p));
end

% one big one
imwrite(concatimages(sprintf('../allsubj*.png'),{1:9},1),sprintf('../consolidated.png'));

===============================================

notes:
- could be done as a movie [eh]

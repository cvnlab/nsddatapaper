%% %%%%%%%%%%%%%%%%%%%%%%% RUN ANALYZEPRF (VERY EXPENSIVE)

% Matlab (matlab2b) run some parpools on:
% $stone: 2 x 12
% atlas8: 16
% $surly: 16
% atlas9: 12
% $atlas10: 16
% $atlas11: 12
% surly: 12
%   $atlas13: 32
%   -atlas12: 2 x 32
%   -atlas2: 16
%   -atlas4: 16
%   -atlas13: 32
%   -atlas3: 12
% ~256 cores ... but now less

% created renice.sh script that renices all of kendrick's MATLABs.
%
% to change chunk size, if all nsdfaceprf.m jobs are COMPLETE, then
% edit LOADDATA1/2 files and we are good to go.

% start each machine with a queue daemon:
kkqueuestart(19);

% each queue script should start with:
pareasy(gethostparpoolsize); 

% HINTS:
% [pos/neg] x [contrast,foreground,background,face,body] x [1,2,3,4,5]

% face semi-automated model using the full data
for p=1:1, nsdfaceprf(p,1,2,1,1:3);, end

% contrastgrid model using the full data
for p=1:8, nsdfaceprf(p,3,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,4,1,1,1:3);, end  % NEW

% face automated model using the full data
for p=1:1, nsdfaceprf(p,1,1,1,1:3);, end

% word semi-automated model using the full data
for p=1:1, nsdfaceprf(p,2,2,1,1:3);, end

% word automated model using the full data
for p=1:1, nsdfaceprf(p,2,1,1,1:3);, end

% at this point ran all 8 subjects maybe?

% proceed to 1.8mm
for p=1:8, nsdfaceprf(p,3,1,1,4);, end
for p=1:8, nsdfaceprf(p,4,1,1,4);, end

% split half!!! face auto
for p=1:8, nsdfaceprf(p,1,1,2,1:3);, end
for p=1:8, nsdfaceprf(p,1,1,3,1:3);, end

% bodies model (automated)
for p=1:8, nsdfaceprf(p,5,1,1,1:3);, end

% foreground model (automated)
for p=1:8, nsdfaceprf(p,6,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,7,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,8,1,1,1:3);, end

% foreground model (automated) 1.8mm volume version both positive and negative
for p=1:8, nsdfaceprf(p,6,1,1,4);, end
for p=1:8, nsdfaceprf(p,7,1,1,4);, end

% contrastgrid model fixed exponent and baseline
for p=1:8, nsdfaceprf(p,9,1,1,1:3);, end

% contrastgrid model using contrastgridNEW preparation
for p=1:8, nsdfaceprf(p,10,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,11,1,1,1:3);, end

% face automated NEGATIVE  [but these results are weird because of offset issues?]
for p=1:8, nsdfaceprf(p,12,1,1,1:3);, end

% cerebellum  [did tk find these interesting?]
for p=1:8
  nsdfaceprf(p,1,1,1,5);   % face auto
  nsdfaceprf(p,2,1,1,5);   % word auto
  nsdfaceprf(p,5,1,1,5);   % bodies
  nsdfaceprf(p,6,1,1,5);   % foreground
  nsdfaceprf(p,8,1,1,5);   % background
  nsdfaceprf(p,10,1,1,5);  % contrastgridNEW
  nsdfaceprf(p,13,1,1,5);  % salience $$$
end

% salience
for p=1:8, nsdfaceprf(p,13,1,1,1:3);, end

% salience, 1.8mm volume version both positive and negative  $$$
for p=1:8, nsdfaceprf(p,13,1,1,4);, end
for p=1:8, nsdfaceprf(p,14,1,1,4);, end

%% %%%%%%%%%%%%%%%%%%%%%%% POST-PROCESS THE RESULTS (this includes saving various files)

% category model (note that rr is [] since we save all rrs together. also, this is a hack to help maintain old results)
for p=1:8
  nsdfaceprffinish(p,1,1,[],1:3,1);
  nsdfaceprffinish(p,1,2,[],1:3,1);
  nsdfaceprffinish(p,2,1,[],1:3,1);
  nsdfaceprffinish(p,2,2,[],1:3,1);
end

% fullprf model
for p=1:1
%  nsdfaceprffinish(1,1,2,1,1:3);
%  nsdfaceprffinish(1,3,1,1,1:3);
%  nsdfaceprffinish(1,1,1,1,1:3);
  nsdfaceprffinish(1,2,1,1,1:3);  
  nsdfaceprffinish(1,2,2,1,1:3);
end
for p=2:8
%  nsdfaceprffinish(p,1,1,1,1:3);
%  nsdfaceprffinish(p,2,1,1,1:3);
%  nsdfaceprffinish(p,1,2,1,1:3);
%  nsdfaceprffinish(p,2,2,1,1:3);
  nsdfaceprffinish(p,3,1,1,1:3);
end
for p=1:8
  nsdfaceprffinish(p,4,1,1,1:3);
end
for p=1:8
  nsdfaceprffinish(p,3,1,1,4);
  nsdfaceprffinish(p,4,1,1,4);
end
for p=1:8
  nsdfaceprffinish(p,1,1,2,1:3,[],0);
  nsdfaceprffinish(p,1,1,3,1:3);
end

% NEXT SET: ran on July 1 2020.
for p=1:8
  nsdfaceprffinish(p,5,1,1,1:3);
end
for p=1:8   % reran this to reinstute standard method
  nsdfaceprffinish(p,3,1,1,1:3);
end

% NEXT SET: ran on July 5+12 2020.
for p=1:8
  nsdfaceprffinish(p,6,1,1,1:3);
  nsdfaceprffinish(p,8,1,1,1:3);
end

% NEXT SET: ran on July 9 2020.
for p=1:8
  nsdfaceprffinish(p,6,1,1,4);
  nsdfaceprffinish(p,7,1,1,4);
end

% july 21 2020:
for p=1:8
  nsdfaceprffinish(p,7,1,1,1:3);
end

% sep 19 2020:
for p=1:8
  nsdfaceprffinish(p,9,1,1,1:3);
end

% sep 23-26 2020:
for p=1:8
  nsdfaceprffinish(p,10,1,1,1:3);
  nsdfaceprffinish(p,11,1,1,1:3);
end

% oct 4 2020:
for p=1:8
  nsdfaceprffinish(p,12,1,1,1:3);
end

% oct 13 2020
for p=1:8
  for wh=[1 2 5 6 8 10]
    nsdfaceprffinish(p,wh,1,1,5);
  end
end

% nov 1 2020:
for p=1:8
  nsdfaceprffinish(p,13,1,1,1:3);
end

% record of what we have:
% categoryresults_subj01-8.mat   - simple face vs non-face
% XXX prfresults_subj01-8.mat    - quick mode fitting for faces and words
% fullprfresults_subj01-8.mat    - full pRF fit for contrast grid model, face/word models
%
% also, we have in predictions/   [actually, we stopped writing these files out]
%   categorypred_subj01-8_lh-_p1-2_q1-2_rr1-4.hdf5
%   prfpred10x10_subj01-8_lh-_p1-2_q1-2_rr1-3.hdf5
%   prfpred_subj01-8_lh-_p1-2_q1-2_rr1-3.hdf5
%   and so on...

%% %%%%%%%%%%%%%%%%%%%%%%% VISUALIZE

% NOTE: we don't do the test-retest here.

for subjix=1:8

% define
modelname = 'fullprf';

% load
finalfile0 = sprintf('~/ext/figurefiles/nsd/%sresults_subj%02d.mat',modelname,subjix);
load(finalfile0);

% init
Lookup = [];

%%%%% EASY

if isequal(modelname,'category')

% define
common = sprintf('~/Dropbox/KKTEMP/figures.science/FIRSTLEVEL/%smodel/subj%02d',modelname,subjix);

% floc-faces
data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/?h.floc-faces.mgz',subjix));
cvnlookup(sprintf('subj%02d',subjix),13,data0,[0 5],jet(256),0.5,Lookup,0); 
cvnlookupwrite(rgbimg,rawimg,common,'floc-faces');

% floc-words
data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/?h.floc-words.mgz',subjix));
cvnlookup(sprintf('subj%02d',subjix),13,data0,[0 5],jet(256),0.5,Lookup,0); 
cvnlookupwrite(rgbimg,rawimg,common,'floc-words');

% curvature
data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/surf/?h.curvature.mgz',subjix));
cvnlookup(sprintf('subj%02d',subjix),13,data0<0,[-1 2],gray(256),[],Lookup,0); 
cvnlookupwrite(rgbimg,rawimg,common,'curvature');

end

%%%%% CATEGORY MODEL

if isequal(modelname,'category')

todo = [1 4];
for zz=1:length(todo)
  for p=1:length(types)
    for q=1:length(labeltypes)
    
      if todo(zz)==4
        if ~(p==1 && q==1)
          continue;
        end
      end

      temp = cat(1,allresults{p,q,1:2,todo(zz)});
      common = sprintf('~/Dropbox/KKTEMP/figures.science/FIRSTLEVEL/%smodel/type%d_subj%02d_p%d_q%d',modelnames{mtodo},todo(zz),subjix,p,q);

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,1),[-100 100],cmapsign4(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'tval');
      
      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,1),[-30 30],cmapsign4(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'tvalALT');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,2),[-90 90],cmapsign4(256),[],Lookup,0);   % TO REGEN
      cvnlookupwrite(rgbimg,rawimg,common,'ang');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,3),[0 5],hot(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'mag');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,4),[0 50],hot(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'magt');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,5),[-100 100],cmapsign4(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'allt');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,6),[0 100],jet(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'predR2');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,7),[-1 1],jet(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'predr');

      cvnlookup(sprintf('subj%02d',subjix),13,temp(:,7),[-.5 .5],jet(256),[],Lookup,0); 
      cvnlookupwrite(rgbimg,rawimg,common,'predrALT');

    end
  end
end

end

%%%%% PRF MODEL (QUICK AND FULL)

if ismember(modelname,{'prf' 'fullprf'})

for p=13%  12  %11   %10  %9  %[7]%4]%[8] %6 %[1:3 5]
  for q=1 %:2
  
    temp = cat(1,allresults{p,q,1:2,1});  % only resampling #1
    if isempty(temp)
      continue;
    end
    mkdirquiet(sprintf('~/Dropbox/KKTEMP/%smodel',modelname));
    common = sprintf('~/Dropbox/KKTEMP/%smodel/subj%02d_p%d_q%d',modelname,subjix,p,q);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,1),[0 360],hsv(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'ang',1);
    
    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,2).^.5,[0 6].^.5,jet(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'ecc',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,3),[0 1],copper(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'expt',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,4).^.5,[0 4].^.5,jet(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'rfsize',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,6),[0 3],hot(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'gain',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,8),[-90 90],cmapsign4(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'prfanglemetric',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,9),[0 100],jet(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'predR2',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,10),[-1 1],jet(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'predr',1);

    cvnlookup(sprintf('subj%02d',subjix),13,temp(:,10),[-.5 .5],jet(256),[],Lookup,0); 
    cvnlookupwrite(rgbimg,rawimg,common,'predrALT',1);


  end
end

end

end

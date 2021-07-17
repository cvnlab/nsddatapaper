% Expects subjix to be defined.
% Creates a number of variables in the workspace (does not include the data).

% define
hemis = {'lh' 'rh'};
datanames = {'lh' 'rh' 'subcortex' 'func1pt8mm' 'cerebellum'};
imageres = 50;  % resolution for the prf model (8.4/50 = 0.168 deg)
pxtodeg = 8.4/imageres;

% calc
fsdir = sprintf('~/nsd/nsddata/freesurfer/subj%02d',subjix);
fsdirlabel = sprintf('/home/stone/generic/Dropbox/nsdfaceprf/freesurfer/subj%02d/label/',subjix);

% load data
load(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d.mat',subjix));

% define
matrixsize1mm = {
 [145 186 148]
 [146 190 150]
 [145 190 146]
 [152 177 143]
 [141 173 139]
 [152 202 148]
 [139 170 145]
 [143 184 139]
};
matrixsize1pt8mm = {
[81 104 83]
[82 106 84]
[81 106 82]
[85 99 80]
[79 97 78]
[85 113 83]
[78 95 81]
[80 103 78]
};

% calc
subcortexorigin = (1+matrixsize1mm{subjix})/2 - ([d1(1) d2(1) d3(1)]-1);  % special origin for small subcortex volume files

% define
nsdstimfile = '~/nsd/nsddata_stimuli/stimuli/nsd/nsd_stimuli.hdf5';

% load experimental design
expdesign0 = load('~/nsd/nsddata/experiments/nsd/nsd_expdesign.mat');
  % expdesign0.subjectim(subjix,ordU) are the 73k IDs for the prepared data that we load in (1 x NIMAGES)

% compute
expdesign0.shared515 = [];  % indices into the 1000 of which ones are special (all 3 shown to all subjects)
for iii=1:1000
  if sum(expdesign0.masterordering(1:750*30)==iii)==3
    expdesign0.shared515 = [expdesign0.shared515 iii];
  end
end

% chunk stuff
numinchunk = 10000;

% knobs
  % [a b e f] for faces (but what about [c g] for faces (maybe?) background??)
  % [a b d e f] for words ([maybe a or f only??] for words??)
manualchoices = {[1 2 5 6] [1 2 4 5 6]};

% calc
valix = find(ordU <= 1000);        % of the images we have, which ones are the shared1000 to use as validation images?
valix2 = ordU(ordU <= 1000);       % of the shared1000 (1:1000), which ones are the ones we have as validation images?
shared515ix = find(ismember(ordU,expdesign0.shared515));  % of the images we have, which ones are the shared515?
shared515ix2 = expdesign0.shared515;                      % of the shared1000 (1:1000), which ones are the shared515?
rep3ix = find(numtrialsperim==3);  % of the images we have, which ones have all 3 trials?

% some resampling stuff
trainixsALL = {};

% idea (AUTOMATED): do NOT use validation images, but use the rest
trainixsALL{1,1} = setdiff(1:length(ordU),valix);  % of the images we have, which ones are the ones to pull for training images?
trainixsALL{1,2} = trainixsALL{1,1}(1:2:end);           % " but odd
trainixsALL{1,3} = trainixsALL{1,1}(2:2:end);           % " but even

% idea (SEMI-AUTOMATED): use only the validation images
trainixsALL{2,1} = valix;                          % of the images we have, which ones are the ones to pull for training images?
trainixsALL{2,2} = trainixsALL{2,1}(1:2:end);           % " but odd
trainixsALL{2,3} = trainixsALL{2,1}(2:2:end);           % " but even

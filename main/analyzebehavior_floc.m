function analyzebehavior_floc(datadir,num_TR,stimulus_delay)

% <num_TR> = 312;         % total number of TRs desired
% <stimulus_delay> = 12;  % delay in seconds at beginning

ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);

% load behavioral files
stimulusfiles = matchfiles([datadir '/mat_files_from_scan/data/*/*run*.mat']); cellfun(@disp,stimulusfiles);
fprintf('=======\n');

% define
design_TR = 1;        % TR in seconds

% process each file
totaldur = [];
hitprop = [];
stimulus = {};
results = [];
for zz=1:length(stimulusfiles)

  % load
  a1 = load(stimulusfiles{zz});
  
  % sanity check
  assert(a1.theSubject.task == 3 && ...
         a1.theSubject.countDown == stimulus_delay && ...
         a1.theSubject.stimSize == 714);

  % record
  totaldur(zz) = a1.theSubject.totalTime;
  hitprop(zz)  = a1.theData.propHit;

  % construct design matrix
  stimulus{zz} = sparse(double(floc_make_design_matrix(stimulusfiles{zz},num_TR,design_TR,stimulus_delay)));

  % prepare nice output
%            keys: {1x600 cell}
%              rt: [1x600 double]
%            resp: [1x600 double]
%           nreps: 20
%            hits: 19
%     falseAlarms: 1
%         propHit: 0.95
   A = a1.theData.nreps;
   B = a1.theData.hits;
   C = a1.theData.falseAlarms;
   D = (B-C)/A*100;
   results(:,zz) = [A B C D];

end

% save
save([ppdir '/behavioralresults_floc.mat'],'totaldur','hitprop','results');
save([ppdir '/designmatrix_floc.mat'],'stimulus');

% make a 3pertrial version
for p=1:length(stimulus)
  assert(sum(flatten(stimulus{p}(1:4:end,:)))==sum(stimulus{p}(:)));
  stimulus{p} = sparse(upsamplematrix(double(stimulus{p}(1:4:end,:)),3,1,[],0));
end
save([ppdir '/designmatrix_3pertrial_floc.mat'],'stimulus');

end

%%%%%%%%%%%%%%%%%%%%%%%

function d_mat = floc_make_design_matrix(run_str, num_TR, design_TR, stimulus_delay)
% FLOC_MAKE_DESIGN_MATRIX makes design matrix out of run string
%
% run_str is the .mat file
% num_TR is the total number of TRs desired
% design_TR is the TR in seconds
% stimulus_delay is number of seconds of delay at the beginning

	tvol = [0:num_TR-1] * design_TR;  % the time corresponding to the data points

	P = load(run_str);
	T = P.theSubject.trials;

	D_stim = zeros(numel(T.block),11);
	D_stim(sub2ind(size(D_stim),1:size(D_stim,1),T.cond+1)) = 1;
	D_stim = D_stim(:,2:end);
	D_tr = interp1(T.onset+stimulus_delay,D_stim,tvol,'linear');
	D_tr(isnan(D_tr)) = 0;

	d_mat = diff([zeros(1,size(D_tr,2)); D_tr])>0;
	    
end

% Expects hh and ii (sorted unique indices) to be defined.
% Creates <somedata> in the workspace.

iimax = max(ii);

somedata = [];
end0 = -Inf;
tic;
for iii=1:length(ii)
  if ii(iii) <= end0
    continue;
  end

  start0 = ii(iii);
  end0 = min(ii(iii)+100,iimax);
  while any(ismember(end0+1:end0+100,ii))
    end0 = min(end0+100,iimax);
  end

  arg1 = [start0 1];
  arg2 = [end0-start0+1 Inf];
  
  fprintf('.');
  temp0 = h5read(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5', ...
                            subjix,datanames{hh}),'/betas',arg1,arg2);
  if isempty(somedata)
    somedata = temp0(ismember(start0:end0,ii),:);
  else
    somedata = cat(1,somedata,temp0(ismember(start0:end0,ii),:));
  end

end
toc;




% SLOW:
if 0 

% h5read parameters
arg1 = [min(ii) 1];
arg2 = [range(ii)+1 Inf];

% load the data
tic;
somedata = h5read(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5', ...
                          subjix,datanames{hh}),'/betas',arg1,arg2);
somedata = somedata(ismember(ii(1):ii(end),ii),:);
toc;

end



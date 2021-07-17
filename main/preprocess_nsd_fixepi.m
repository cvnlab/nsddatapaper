function epis = preprocess_nsd_fixepi(epis,sessix)

switch sessix

case 35  % subj01-nsd02
  epis{2}(:,:,:,128:152) = repmat(epis{2}(:,:,:,128),[1 1 1 152-128+1]);

case 40  % subj08-nsd02

  % 04b, 5, and 6 need to be scaled by 0.78/0.4
  epis{5} = (0.78/0.4) * epis{5};
  epis{6} = (0.78/0.4) * epis{6};
  epis{7} = (0.78/0.4) * epis{7};
  
  % 50 trials is 200 s.  200/1.6 = 125.
  epis{4}(:,:,:,111:125) = repmat(epis{4}(:,:,:,111),[1 1 1 125-111+1]);

case 295  % subj06-nsd31
  epis{1}(:,:,:,181:188) = repmat(epis{1}(:,:,:,181),[1 1 1 188-181+1]);
  epis{2}(:,:,:,181:188) = repmat(epis{2}(:,:,:,181),[1 1 1 188-181+1]);
 
end

% check
epis

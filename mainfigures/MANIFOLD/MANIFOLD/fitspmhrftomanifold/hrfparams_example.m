% load HRF parameters
load('hrfparams.mat','params');  % 20 different HRFs x 7 parameters

% visualize
figure; hold on;
cmap0 = copper(20);
for p=1:size(params,1)
  temp = spm_hrf(0.1,params(p,:));
  temp = temp / max(temp);
  plot(0:0.1:50,temp,'k-','Color',cmap0(p,:));
end
xlabel('Time (s)');

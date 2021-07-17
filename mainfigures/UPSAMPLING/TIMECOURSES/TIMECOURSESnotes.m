subjix = 1;
load(sprintf('~/ext/figurefiles/nsd/onofffirSPECIAL_subj%02d.mat',subjix),'timecourses','R2');

% 1mm: (82,24,73)

offset = 13;
offset = 2;  % just a few
xx = 82+offset:-1:82-offset;

figureprep([0 0 1500 400]);
for ii=1:length(xx)
  subplot(1,length(xx),ii); hold on;
  temp = squish(timecourses(xx(ii),1,73,:,:),4);
  plot(0:30,temp,'-');%,'Color',[1 .7 .7]);
  plot(0:30,mean(temp,2),'k-','LineWidth',2);
  set(errorbar2(0:30,mean(temp,2),std(temp,[],2)/sqrt(40),'v','k-'),'LineWidth',2);
  axis([0 30 -6 6]);
  straightline(0,'h','k-');
  straightline(0:3:30,'v','b-');
end
figurewrite('tcs',[],-2,'~/Dropbox/KKTEMP');

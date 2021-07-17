% orthMID42.png taken from subj01

% generate distribution visualization
n = 100000;
vals = vmrand((360/40/2)/180*pi,3^6,[1 6*n]);  % center in middle of first session
vals2 = rand(1,4*n)*2*pi-pi;
figureprep([100 100 400*.6 250*.6],1); hold on;
h = histogram([vals vals2]/pi*180,5000,'Normalization','pdf');
set(h,'FaceColor','r');
set(h,'EdgeColor','none');
ax = axis;
axis([-180 180 -0.01 ax(4)]);
%straightline(0,'h','k-');
xlabel('Angle (deg)');
ylabel('p(x)');
set(gca,'XTick',-180:90:180);
set(straightline(linspace(-180,180,41),'v','k-',[-.01 .01]),'Color',[.7 .7 .7]);
figurewrite('distribution',[],-1,'~/Desktop');

% visualize special metrics
a1 = load('/research/stimuli/multiclass/prepareexpdesign_nsd/simulations.mat');
pp = 8; qq = 3;
figureprep([100 100 280 220],1); hold on;
h1 = plot(a1.results(pp,qq).novelfrac,'LineWidth',2); 
h2 = plot(a1.results(pp,qq).recentfrac,'LineWidth',2); 
h3 = plot(a1.results(pp,qq).recentfrac2,'LineWidth',2);
axis([0 41 0 100]);
xlabel('Session number');
ylabel('Percent (%)');
legend([h1 h2 h3],{'Out of all trials, percent that are new trials' 'Out of old trials, percent that are easy trials' 'Out of all trials, percent that are easy trials'},'Location','East');
figurewrite('empiricalpercentages',[],-1,'~/Desktop');

% experimental design timecourse of one run
a2 = load('/research/stimuli/multiclass/prepareexpdesign_nsd/nsd_expdesign.mat');
figureprep([100 100 1500 100],1); hold on;
cnt = 0;
xx = [];
yy = [];
for p=1:75
  if a2.stimpattern(1,1,p)==0
    xx = [xx cnt+(1:4)];
    yy = [yy zeros(1,4)];
  else
    xx = [xx cnt+1 cnt+(1:4) cnt+(4)];
    yy = [yy 0 1 1 1 1 0];
  end
  cnt = cnt + 4;
end
xx = [xx cnt+1];
yy = [yy 0];
plot(xx,yy);
ax = axis;
axis([0 302 -2 2]);
figurewrite('timecourserun',[],-1,'~/Desktop');

%%%%%%%%%%%%%%%%%%%%%%%%% determine coco id

aa=h5read('nsd_stimuli.hdf5','/imgBrick');
okok = zeros(64,64,73000);
for p=1:size(aa,4)
  if mod(p,500)==1, fprintf('.'); end
  im0 = permute(aa(:,:,:,p),[3 2 1]);
  im0 = rgb2gray(im0);
  im0 = imresize(im0,[64 64]);
  okok(:,:,p) = im0;
end

trg = imread('~/Dropbox/KKTEMP/11.png');
trg = rgb2gray(trg);
trg = imresize(trg,[64 64]);

rs = calccorrelation(double(squish(okok,2)),double(repmat(trg(:),[1 size(okok,3)])),1);

[~,iix] = max(rs)  % 21254

figure; imagesc(permute(aa(:,:,:,21254),[3 2 1]));

% 21254 corresponds to COCO ID 55402

highlevelvisual ROI..  this was quick and dirty ROI definition.
  roilabels = {'posterior' 'anterior' 'V1' 'V2' 'V3'};

==================== PLOT THE BRAIN MAP FOR ROIs

cd '/research/papers/2022 nsd/figures/SCIENCE.RSA'
vals = cvnloadmgz('*.highlevelvisual.mgz');
Lookup = [];
eo = {'hemibordercolor' 'w' 'rgbnan' 1 'drawscalebar' false};
view0 = {'ventral'         'inflated'                 1  1500    0         [1 1]};
cvnlookup('fsaverage',view0,vals,[0 5],jet,0.5,Lookup,0,eo);
imwrite(rgbimg,'fsaverage ventral inflated highlevelvisual.png');

===================== ADDITIONAL ANALYSES

% separate into 8 subjects
data = np.load('all_sub_modelcomp_corrs.npy')
np.save('all_sub_modelcomp_corrs_subj1.npy',data[0])
np.save('all_sub_modelcomp_corrs_subj2.npy',data[1])
np.save('all_sub_modelcomp_corrs_subj3.npy',data[2])
np.save('all_sub_modelcomp_corrs_subj4.npy',data[3])
np.save('all_sub_modelcomp_corrs_subj5.npy',data[4])
np.save('all_sub_modelcomp_corrs_subj6.npy',data[5])
np.save('all_sub_modelcomp_corrs_subj7.npy',data[6])
np.save('all_sub_modelcomp_corrs_subj8.npy',data[7])

%%%%% COCO category model

% load
files0 = matchfiles('all_sub_modelcomp_corrs_subj*.npy');
data = {};  % each is independent-draws x 5 ROIs (each draw involves 100 images)
for p=1:length(files0)
  data{p} = readNPY(files0{p});
end

% within-subject error
tempmn = [];  % subjects x ROIs
tempse = [];  % subjects x ROIs
for p=1:8
  tempmn(p,:) = mean(data{p},1);
  tempse(p,:) = std(data{p},[],1)/sqrt(size(data{p},1));
end

% plot
ord = [3 4 5 1 2];  % re-ordering to get a nicer order
cmap1 = jet(8);
cmap1(3:5,:) = (cmap1(3:5,:)*5 + [0 0 0])/6;
roilabels = {'pVTC' 'aVTC' 'V1' 'V2' 'V3'};
figureprep([100 100 220 180]); hold on;
for p=1:8
  errorbar3(1:5,tempmn(p,ord),tempse(p,ord),'v',(cmap1(p,:)+2*[1 1 1])/3);
  plot(1:5,tempmn(p,ord),'k-','Color',cmap1(p,:),'LineWidth',2);
end
mn = mean(tempmn,1);
se = std(tempmn,[],1)/sqrt(size(tempmn,1));
errorbar3(1:5,mn(ord),se(ord),'v',[.8 .8 .8]);
plot(1:5,mn(ord),'k-','LineWidth',3);
set(gca,'XTick',1:5,'XTickLabel',roilabels(ord));
ylabel(sprintf('Correlation of category RDM\nwith brain RDM'));
xlim([0.5 5.5]);
figurewrite('cococategorymodel',[],-2,'~/Dropbox/KKTEMP');

%%%%%

% load
data2 = readNPY('all_sub_roi_rdms_correlations.npy');  % 40 x 40
data2(logical(eye(size(data2,1)))) = 1;  % diagonal is actually 1

% define
ord2 = [1:5:40 2:5:40 3:5:40 4:5:40 5:5:40];  % reorder such that outer is ROI, inner is subjects
data3 = data2(ord2,ord2);

% plot
figureprep([100 100 230 180]); hold on;
imagesc(data3);
set(gca,'YDir','reverse');
axis image;
colormap(cmapturbo);
caxis([0 1]);
cb = colorbar;
straightline(0.5:8:40.5,'v','k-');
straightline(0.5:8:40.5,'h','k-');
set(gca,'XTick',(1+8)/2:8:40,'XTickLabel',roilabels(ord));
set(gca,'YTick',(1+8)/2:8:40,'YTickLabel',roilabels(ord));
figurewrite('acrosssubject',[],-2,'~/Dropbox/KKTEMP');

% extract across-subject similarities
ok = [];  % 5 x 5 x boots. the i,j entry indicates average correlation between ROI i and ROI j for DISTINCT subjects
for roi1=1:5
  for roi2=1:5
    temp = data3((roi1-1)*8+(1:8),(roi2-1)*8+(1:8));
    for boot=1:1000
      ix = ceil(rand(1,8)*8);    % bootstrap subjects
      isvalid = logical([8 8]);  % this will indicate valid entries in this bootstrap sample
      isvalid(:) = 0;
      for rowix=1:8
        for colix=1:8
          if ix(rowix) ~= ix(colix)
            isvalid(rowix,colix) = 1;
          end
        end
      end
      temp2 = temp(ix,ix);
      ok(roi1,roi2,boot) = mean(temp2(isvalid));
    end
  end
end

% plot the quantitative summary
figureprep([100 100 220 180]); hold on;
for p=1:5
  mn0 = mean(ok(p,:,:),3);
  se0 = std(ok(p,:,:),[],3);
  color0 = cmaplookup(ord(p),0,5,0,jet(256));
  errorbar3(1:5,mn0,se0,'v',(color0+2*[1 1 1])/3);
  plot(1:5,mn0,'r-','Color',color0,'LineWidth',2);
end
set(gca,'XTick',1:5,'XTickLabel',roilabels(ord));
xlim([0.5 5.5]);
ylabel(sprintf('Average across-subject\ncorrelation of brain RDM'));
figurewrite('acrosssubjectquantitative',[],-2,'~/Dropbox/KKTEMP');

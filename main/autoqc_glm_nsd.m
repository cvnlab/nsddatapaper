function autoqc_glm_nsd(sessix,glmdirname)

% setup
nsdsetup;
datadir = datadirs{sessix};
glmdir = regexprep(datadir,'rawdata',glmdirname);
ppdir = regexprep(datadir,'rawdata','ppdata');
outputdir = sprintf('%s/../%s-%s-glmBASIC',ppdir,nsdidfun(sessix),glmdirname);

% check
isspecial = 0;
if nsdscreeningfun(sessix)
  sessnum = 0;
  isspecial = 1;
elseif ismember(sessix,nsdsyntheticsessions)
  sessnum = totalnsdsessions(nsdsubjixfun(sessix)) + 1;
  isspecial = 1;
elseif ismember(sessix,nsdimagerysessions)
  sessnum = totalnsdsessions(nsdsubjixfun(sessix)) + 2;
  isspecial = 1;
end
if isspecial
  if isequal(glmdirname,'glmdata')
    a1 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGE/mean.nii']);
  elseif isequal(glmdirname,'glmdataLOW')
    a1 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGELOW/mean.nii']);
  elseif isequal(glmdirname,'glmdataSUPER')
    a1 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGESUPER/mean.nii']);
  end
  meanvol = double(a1.img);
else
  sessnum = nsdsessnumfun(sessix);
  a1 = load([glmdir '/GLMdenoise_nsdBASIConoff.mat']);
  meanvol = a1.meanvol;
  %%a2 = load([glmdir '/GLMdenoise_nsdBASICsingletrialfithrf.mat']);
end

% make output dirs
mkdirquiet(outputdir);

% write mean volume
imwrite(uint8(255*makeimagestack(meanvol,1)),gray(256),sprintf('%s/meanvol%02d.png',outputdir,sessnum));

% write mean volume with max orthogonal
mx = max(sizefull(meanvol,3));
mn = mean(meanvol(:));
im1 = placematrix(mn*ones(mx,mx),max(meanvol,[],3),[1 1]);
im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(max(meanvol,[],2),2),1,2,1),[1 1]);
im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(max(meanvol,[],1),2),1,2,1),[1 1]);
im = cat(3,im1,im2,im3);
imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/orthMAX%02d.png',outputdir,sessnum));
im1 = placematrix(mn*ones(mx,mx),mean(meanvol,3),[1 1]);
im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(mean(meanvol,2),2),1,2,1),[1 1]);
im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(mean(meanvol,1),2),1,2,1),[1 1]);
im = cat(3,im1,im2,im3);
imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/orthMEAN%02d.png',outputdir,sessnum));
im1 = placematrix(mn*ones(mx,mx),meanvol(:,:,round(end/2)),[1 1]);
im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(meanvol(:,round(end/2),:),2),1,2,1),[1 1]);
im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(meanvol(round(end/2),:,:),2),1,2,1),[1 1]);
im = cat(3,im1,im2,im3);
imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/orthMID%02d.png',outputdir,sessnum));

if ~isspecial

  %% "onoff.mat"

  % write R2
  imwrite(uint8(255*makeimagestack(signedarraypower(a1.R2/100,0.5),[0 1])),hot(256),sprintf('%s/varexp%02d.png',outputdir,sessnum));

  % write R2 orth
  im1 = placematrix(zeros(mx,mx),signedarraypower(a1.R2(:,:,round(end/2))/100,0.5),[1 1]);
  im2 = placematrix(zeros(mx,mx),signedarraypower(rotatematrix(squish(a1.R2(:,round(end/2),:),2),1,2,1)/100,0.5),[1 1]);
  im3 = placematrix(zeros(mx,mx),signedarraypower(rotatematrix(squish(a1.R2(round(end/2),:,:),2),1,2,1)/100,0.5),[1 1]);
  im = cat(3,im1,im2,im3);
  imwrite(uint8(255*makeimagestack(im,[0 1])),hot(256),sprintf('%s/varexporthMID%02d.png',outputdir,sessnum));

  % write PSC
  imwrite(uint8(255*makeimagestack(a1.modelmd,[-10 10])),cmapsign4(256),sprintf('%s/psc%02d.png',outputdir,sessnum));

  % write grand R2  [NOTE: nanmean]
  files0 = cellfun(@(x) regexprep(x,'rawdata',glmdirname),nsdallsessfun(sessix),'UniformOutput',0);
  a1R2 = nanmean(loadmulti(cellfun(@(x) [x '/GLMdenoise_nsdBASIConoff.mat'],files0,'UniformOutput',0),'R2',4),4);
  imwrite(uint8(255*makeimagestack(signedarraypower(a1R2/100,0.5),[0 1])),hot(256),sprintf('%s/varexpcum.png',outputdir));
  
  % write grand R2 upsampled
  if isequal(glmdirname,'glmdataLOW')
    [xx,yy,zz] = ndgrid(1:1/1.8:size(a1R2,1),1:1/1.8:size(a1R2,2),1:1/1.8:size(a1R2,3));
    a1R2 = reshape(ba_interp3_wrapper(a1R2,[flatten(xx); flatten(yy); flatten(zz)]),size(xx));
    imwrite(uint8(255*makeimagestack(signedarraypower(a1R2/100,0.5),[0 1])),hot(256),sprintf('%s/varexpcumUP.png',outputdir));
  end

end

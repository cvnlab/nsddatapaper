[kendrick@surly ~/nsd/ppdata]$ ls -la *prffloc/GLMdenoise_floc_figures/FinalModel.png 

a = matchfiles('~/nsd/ppdata/*prffloc/GLMdenoi*figures/FinalModel.png');
a = a([1 2 3 4 5 6 7 8 9    11 10   12  14 13]);  % to match order in nsdsetup.m (FINALVERSION HAS RIGHT ORDER)

mkdirquiet('~/Dropbox/KKTEMP/ok/');
for p=1:length(a)
  copyfile(a{p},sprintf('~/Dropbox/KKTEMP/ok/sess%02d.png',p));
end

imcropfiles('*png',[],hot(256))
  % 6th row, 3-8 columns

im = concatimages('*png',num2cell(1:14)');
imwrite(im,hot(256),'final.png');

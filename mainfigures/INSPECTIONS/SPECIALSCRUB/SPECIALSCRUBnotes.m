%%%%%%%%%%%%% pick randomvolumes and write an orthogonal view png

outputdir = '~/nsd/ppdata/specialscrub';
mkdirquiet(outputdir);
fullrec = {};
for subjix=1:8, subjix
  files = matchfiles(sprintf('~/nsd/nsddata_timeseries/ppdata/subj%02d/func1pt8mm/timeseries/timeseries_session*.nii.gz',subjix));
  for rep=1:100, rep
    runix = ceil(rand*length(files));
    a1 = load_untouch_nii(files{runix});
    volix = ceil(rand*size(a1.img,4));
    meanvol = permute(flipdim(flipdim(double(a1.img(:,:,:,volix)),1),2),[2 1 3]);
    
    mx = max(sizefull(meanvol,3));
    mn = mean(meanvol(:));
    im1 = placematrix(mn*ones(mx,mx),meanvol(:,:,round(end/2)),[1 1]);
    im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(meanvol(:,round(end/2),:),2),1,2,1),[1 1]);
    im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(meanvol(round(end/2),:,:),2),1,2,1),[1 1]);
    im = cat(3,im1,im2,im3);
    imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/subj%02d_rep%02d.png',outputdir,subjix,rep));
    
    f = regexp(files{runix},'.+session(\d+).+run(\d+)','tokens');
    fullrec{subjix,rep} = [subjix str2double(f{1}{1}) str2double(f{1}{2}) volix];
  end
end

mkdirquiet('~/Dropbox/KKTEMP/SPECIALSCRUB');
save('~/Dropbox/KKTEMP/SPECIALSCRUB/results.mat','fullrec');

%%%%%%%%%%%%% compile

% this is based on autoqc_nsd_grandfinale.m...

alltempfiles = {};

min0 = -710;

% calc
outputdir = '~/Dropbox/KKTEMP/SPECIALSCRUB/figures';

% make output dir
mkdirquiet(outputdir);

% orthMID
masterfile = {};
for p=1:100
  order = {[1 5] [2 6] [3 7] [4 8]};
  for q=1:8
    a0 = matchfiles(sprintf('~/nsd/ppdata/specialscrub/subj%02d_rep%02d.png',q,p));
    if ~isempty(a0)
      masterfile{q} = a0{1};

      % SPECIAL TITLES!!!
      imtemp = imread(masterfile{q});
      imtemp = imresize(imtemp,2,'nearest');  % notice that we upsample by 2 using nearest
      imtemp = cat(1,255*ones(2,size(imtemp,2)),imtemp);
      imtemp = cat(2,255*ones(size(imtemp,1),2),imtemp);
      im = drawtexts(500,0,0,'FreeMono',0.04,[1 1 1],[0 0 0],sprintf('subj%02d, session %2d, run %02d, vol %03d',fullrec{q,p}(1),fullrec{q,p}(2),fullrec{q,p}(3),fullrec{q,p}(4)),{'FontWeight','bold'});
      im = uint8(255*placematrix(zeros(30,size(imtemp,2)),im));
      imtemp = cat(1,im,imtemp);
      fileX = [tempname '.png'];
      imwrite(imtemp,fileX);
      masterfile{q} = fileX;
      alltempfiles = [alltempfiles {fileX}];
    end
  end
  im0 = concatimages(masterfile,order,0,min0);
  imwrite(uint8(im0),gray(256),sprintf('%s/rep%03d.png',outputdir,p));
end

for zz=1:length(alltempfiles)
  delete(alltempfiles{zz});
end

%%%%%%%%%%%%% postprocess

15 fps, handbrake (HQ 1080, custom 0000 (for cropping)) -> 
  randomscrubbing.mp4

% on surly:

% annotate the individual subjects
for p=1:8
  unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_curvature_subj%02d.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''subj%02d curvature'' -stroke none -annotate +20+20 ''subj%02d curvature'' ~/Dropbox/KKTEMP/APPLE%02d.png',p,p,p,p));
end

% annotate the groupavg result
unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_curvature_groupavg.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''groupavg curvature'' -stroke none -annotate +20+20 ''groupavg curvature'' ~/Dropbox/KKTEMP/APPLE09.png'));

% annotate the fsaverage curvature
unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_curvature_fsaverage.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''fsaverage curvature'' -stroke none -annotate +20+20 ''fsaverage curvature'' ~/Dropbox/KKTEMP/APPLE10.png'));

% make order (this is really wasteful)
ord = [upsamplematrix([1:8 repmat(9,[1 3]) repmat(10,[1 3])],3,2,[],'nearest') [1:8 repmat(9,[1 3*3]) repmat(10,[1 3*3])]];
mkdirquiet('~/Dropbox/KKTEMP/imagelist');
for p=1:length(ord)
  copyfile(sprintf('~/Dropbox/KKTEMP/APPLE%02d.png',ord(p)),sprintf('~/Dropbox/KKTEMP/imagelist/final%03d.png',p));
end

%%%%%%%%%%%%%%%

% on stone:
cd ~/Dropbox/KKTEMP/
ffmpeg -y -framerate 6 -pattern_type glob -i 'imagelist/final*.png' -crf 18 -c:v libx264 -pix_fmt yuv420p /home/stone-ext4/generic/Dropbox/nsddata/inspections/fsaveragecheck.mp4

% convert
handbrake (HQ 720) 
-> fsaveragecheck.mp4

%%%%%%%%%%%%%%% JUNK BELOW

% on stone:  (at fps 2, 14/2*2 = 14 s) (at fps 6, 14/6*6 = 14)
cd ~/Dropbox/KKTEMP/
ffmpeg -y -pattern_type glob -i 'A*.png' -crf 18 -c:v libx264 -pix_fmt yuv420p -vf "zoompan=d=25+'50*eq(in,3)'+'100*eq(in,5)'" movie.mp4

% % on stone:  (at fps 2, 14/2*2 = 14 s) (at fps 6, 14/6*6 = 14)
% cd ~/Dropbox/KKTEMP/
% ffmpeg -y -loop 1 -framerate 2 -pattern_type glob -i 'A*.png' -crf 18 -c:v libx264 -pix_fmt yuv420p -t 14 movie1.mp4
% ffmpeg -y -loop 1 -framerate 6 -pattern_type glob -i 'A*.png' -crf 18 -c:v libx264 -pix_fmt yuv420p -t 14 movie2.mp4

% in bash:
ffmpeg -f concat -safe 0 -i <(find . -name 'movie*.mp4' -printf "file '$PWD/%p'\n") -c copy finalmovie.mp4

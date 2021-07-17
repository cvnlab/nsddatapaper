% -stroke '#FFFFFF
% -strokewidth 2 
% columnsxrows

% on surly

remap = [1 3 5 7 2 4 6 8];

for p=1:8
  unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/subj%02dflat_curvature_subj%02d.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''subj%02d, native surface'' -stroke none -annotate +20+20 ''subj%02d, native surface'' /tmp/kk/A%02d.png',p,p,p,p,remap(p)));
end
unix('montage -background ''#ffffff'' -gravity Center -geometry 1000x+10+10 -tile 2x4 /tmp/kk/A*.png /home/stone-ext4/generic/Dropbox/nsddata/inspections/subjectmontages/native.png');
unix('rm -rf /tmp/kk/A*.png');

for p=1:8
  unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/subj%02dflat_corticalsulc_subj%02d.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''subj%02d, native surface + corticalsulc'' -stroke none -annotate +20+20 ''subj%02d, native surface + corticalsulc'' /tmp/kk/C%02d.png',p,p,p,p,remap(p)));
end
unix('montage -background ''#ffffff'' -gravity Center -geometry 1000x+10+10 -tile 2x4 /tmp/kk/C*.png /home/stone-ext4/generic/Dropbox/nsddata/inspections/subjectmontages/nativewithcorticalsulc.png');
unix('rm -rf /tmp/kk/C*.png');

for p=1:8
  unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_curvature_subj%02d.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''subj%02d, fsaverage surface'' -stroke none -annotate +20+20 ''subj%02d, fsaverage surface'' /tmp/kk/E%02d.png',p,p,p,remap(p)));
end
unix('montage -background ''#ffffff'' -gravity Center -geometry 1000x+10+10 -tile 2x4 /tmp/kk/E*.png /home/stone-ext4/generic/Dropbox/nsddata/inspections/subjectmontages/fsaverage.png');
unix('rm -rf /tmp/kk/E*.png');

for p=1:8
  unix(sprintf('convert /home/stone-ext4/generic/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_corticalsulc_fsaverage.png -gravity northwest -font FreeMono-Bold -fill black -pointsize 120 -stroke ''#eeeeee'' -strokewidth 10 -annotate +20+20 ''fsaverage surface + corticalsulc'' -stroke none -annotate +20+20 ''fsaverage surface + corticalsulc'' /tmp/kk/G%02d.png',remap(p)));
end
unix('montage -background ''#ffffff'' -gravity Center -geometry 1000x+10+10 -tile 2x4 /tmp/kk/G*.png /home/stone-ext4/generic/Dropbox/nsddata/inspections/subjectmontages/fsaveragewithcorticalsulc.png');
unix('rm -rf /tmp/kk/G*.png');


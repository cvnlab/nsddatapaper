function autoqc_acquisition(sessix)

% Check if EPI and fieldmaps have the same SliceLocation

nsdsetup;

rec = []; rec2 = [];

bb = matchfiles([datadirs{sessix} '/dicom/MR*gre*']);
for q=1:length(bb)
  cc = matchfiles([bb{q} '/*.dcm']);
  a = dicominfo(cc{1});
  rec(q) = a.SliceLocation;
end

bb = matchfiles([datadirs{sessix} '/dicom/MR*bold*']);
for q=1:length(bb)
  cc = matchfiles([bb{q} '/*.dcm']);
  a = dicominfo(cc{1});
  rec2(q) = a.SliceLocation;
end

if allzero(std(rec)) && ...
   allzero(std(rec2)) && ...
   allzero((rec(1)-1.8) - (rec2(1)-0.9))
else
  emailme(sprintf('SliceLocation failure for %d',sessix));
end

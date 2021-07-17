function f = SAVEDATA(hh,data,dataname)

% function f = SAVEDATA(hh,data,dataname)
%
% <hh> is 1/2/3/4/5 (lh/rh/subcortex/func1pt8mm/cerebellum)
% <data> is V x QUANTITIES
% <dataname> is a string label
%
% Save the <data> to a "freesurfer" directory location.
% Relies on some variables from the caller workspace (from LOADDATA.m)!!

% define
hemis = {'lh' 'rh'};
fsdir =           evalin('caller','fsdir');
fsdirlabel =      evalin('caller','fsdirlabel');
dsz =             evalin('caller','dsz');
bmii =            evalin('caller','bmii');
subcortexorigin = evalin('caller','subcortexorigin');
subjix          = evalin('caller','subjix');
cdsz            = evalin('caller','cdsz');
cii             = evalin('caller','cii');
cd1             = evalin('caller','cd1');
cd2             = evalin('caller','cd2');
cd3             = evalin('caller','cd3');
matrixsize1mm   = evalin('caller','matrixsize1mm');

% ensure dir exists
mkdirquiet(sprintf('%s/func1pt8mm',fsdirlabel));

% save
if hh <= 2
  nsd_savemgz(data,sprintf('%s/%s.%s.mgz',fsdirlabel,hemis{hh},dataname),fsdir);
elseif hh == 3
%   for qq=1:size(data,2)
%     [~,mn,mx] = robustrange(data(:,qq));
%     data(:,qq) = normalizerange(data(:,qq),mn,mx,mn,mx);
%   end
  nsd_savenifti(reshape(data,dsz(1),dsz(2),dsz(3),[]),[1 1 1], ...
                sprintf('%s/%s.nii.gz',fsdirlabel,dataname),1,subcortexorigin);
elseif hh == 4
  newdata = zeros([prod(size(bmii)) size(data,2)],class(data));
  for ii=1:size(data,2)
    newdata(find(bmii),ii) = data(:,ii);
  end
  nsd_savenifti(reshape(newdata,[size(bmii) size(data,2)]),[1.8 1.8 1.8], ...
                sprintf('%s/func1pt8mm/%s.nii.gz',fsdirlabel,dataname),4/3);
elseif hh == 5
  newdata = NaN*zeros([matrixsize1mm{subjix} size(data,2)],class(data));
  for ii=1:size(data,2)
    temp = NaN*zeros(cdsz);
    temp(cii) = data(:,ii);
    newdata(cd1,cd2,cd3,ii) = temp;
  end
  nsd_savenifti(newdata,[1 1 1], ...
                sprintf('%s/cerebellum/%s.nii.gz',fsdirlabel,dataname),1);
end

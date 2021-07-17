% this is subj08
bet MRAangiogramDATECENSORED_121153TOF3Dmultislabs037a1001.nii.gz angio.nii.gz -f 0.25


cd ~/Dropbox/KKTEMP

%a = load_untouch_nii(gunziptemp('MRAangiogramDATECENSORED_121153TOF3Dmultislabs037a1001.nii.gz'));
a = load_untouch_nii(gunziptemp('angio.nii.gz'));

a = double(a.img);

skip = 3;
skip = 1;
[xx,yy,zz] = ndgrid(1:skip:464,1:skip:512,1:skip:264);

xx = xx - (464/2);
yy = yy - (512/2);
zz = zz - (264/2);

xx = xx * 0.390625;
yy = yy * 0.390625;
zz = zz * 0.5;

C = [xx(:) yy(:) zz(:) ones(numel(xx),1)];

huh = [];
angs = 0:359;

for r=1:2
  for p=1:length(angs), p
    switch r
    case 1
      rot0 = xyzrotate_x(angs(p));
    case 2
      rot0 = xyzrotate_y(angs(p));
    case 3
      rot0 = xyzrotate_z(angs(p));
    end
  
    temp = C * rot0;
  
  %   xx0 = reshape(temp(:,1) + 464/2,size(xx));
  %   yy0 = reshape(temp(:,2) + 512/2,size(yy));
  %   zz0 = reshape(temp(:,3) + 264/2,size(zz));

    data = ba_interp3_wrapper(a,bsxfun(@plus,bsxfun(@rdivide,temp(:,1:3)',[0.390625 0.390625 0.5]'),[464 512 264]'/2),'linear');

    huh(:,:,(r-1)*360+p) = max(reshape(data,size(xx)),[],3);
%    imagesc(huh(:,:,(r-1)*360+p)); drawnow;
  end
end

%%%%%

mx = max(huh(:));
imagesequencetomovie(uint8(255*normalizerange(huh,0,1,mx/10,mx)), ...
  '~/Dropbox/KKTEMP/angio.mov',60);

handbrake -> production standard.

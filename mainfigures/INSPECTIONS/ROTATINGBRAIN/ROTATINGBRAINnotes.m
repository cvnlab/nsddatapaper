%%%%%%%%%%%%%%%%%%%%%%%%%% quick one off to get the values

todo = { ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/%s.%s.mgh']  @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b3nc'} ...
};
hemis = {'lh' 'rh'};

for tt=1:length(todo)
  name =    todo{tt}{1};
  fileloc = todo{tt}{2};
  tfun =    todo{tt}{3};
  rng =     todo{tt}{4};
  cmap =    todo{tt}{5};
  thresh =  todo{tt}{6};
  outname = todo{tt}{7};

  % load data and map to fsaverage
  vals = {};
  for p=1:8
    for hh=1:length(hemis)
      inputfile = matchfiles(sprintf(fileloc,p,hemis{hh},name));
      assert(length(inputfile)==1);
      inputfile = inputfile{1};
      data = cvnloadmgz(inputfile);
      vals{p,hh} = nsd_mapdata(p,[hemis{hh} '.white'],'fsaverage',data);
    end
  end
end

vals = vflatten(mean(tfun(reshape(cell2mat(vals),[],8,2)),2));

save('~/Dropbox/KKTEMP/b3nc.mat','vals');

%%%%%%%%%% interlude to generate nice inflated versions

load('~/Dropbox/KKTEMP/b3nc.mat','vals');  % 328k x 1

extraopts = {'rgbnan',1,'hemibordercolor',[1 1 1],'scalebarcolor','w'};  %% NOTE

cvnlookup('fsaverage',{'medial' 'inflated' 0 1000 0 [1 1]},vals,[0 75],jet(256),15,[],0,extraopts);
imwrite(rgbimg,'~/Dropbox/KKTEMP/b3ncinflatedmedial.png');
  
cvnlookup('fsaverage',{'lateral' 'inflated' 0 1000 0 [1 1]},vals,[0 75],jet(256),15,[],0,extraopts);
imwrite(rgbimg,'~/Dropbox/KKTEMP/b3ncinflatedlateral.png');

%%%%%%%%%%%%%%%%%%%%%%%%%% proceed

islh = 1;
islh = 0;

if islh
  a1 = load('~/ext/anatomicals/fsaverage/lhwhite.mat');
  hemi = 'lh';
else
  a1 = load('~/ext/anatomicals/fsaverage/rhwhite.mat');
  hemi = 'rh';
end

vals = loadmulti('~/Dropbox/KKTEMP/b3nc.mat','vals');
if islh
  val = flatten(vals(1:163842));
else
  val = flatten(vals(163842+1:end));
end

val2 = -restrictrange(a1.curvature',-.5,.5)>0;

rng = [0 75];

viewsurfacedata(sprintf('~/ext/anatomicals/fsaverage/%swhite.mat',hemi), ...
                sprintf('~/ext/anatomicals/fsaverage/%sinflated.mat',hemi), ...   % or lhwhite or lhinflated
                {val2 val},{[-1 2] rng},[],[],{[] val});

% data to jet.  thresh >= 15
% curvature data to gray
% view to corticalsulcfs...
% zoom to 0.8, three single clicks left (scroll left) (LH) or four single clicks right (scroll right) (RH)
% background to white

% PRE FLIGHT TO ENSURE THAT JET ROI is NOT INTERPOLATED (FLAT)
setfigurepos([100 100 800 800]);
test = findobj('type','patch');
set(test(2),'FaceColor','flat');

%%%%%% GENERATE THE IMAGES

if islh
  outdir = 'lhinflated';
else
  outdir = 'rhinflated';
end

% generate the frames
global VS_GUI;
cd ~/Desktop;
dd = 1;
mkdirquiet(outdir);
cnt = 1;

% deal with lh/rh setting
if islh
  sc = 1;   % lh
else
  sc = -1;  % rh
end

% first animation
for p=1:(360/dd)
  handles = guidata(VS_GUI);
  set(handles.rotleft,'Value',dd*sc);
  viewsurfacedata_gui('rotleft_Callback',0,[],handles);
  viewsurfacedata_autosnapshot2; movefile('image000.png',sprintf('%s/frame%03d.png',outdir,cnt));
  cnt = cnt + 1;
end

% roll animation
for p=1:(90/dd)
  handles = guidata(VS_GUI);
  set(handles.rollccw,'Value',dd*sc);
  viewsurfacedata_gui('rollccw_Callback',0,[],handles);
  viewsurfacedata_autosnapshot2; movefile('image000.png',sprintf('%s/frame%03d.png',outdir,cnt));
  cnt = cnt + 1;
end

% another swivel
for p=1:(360/dd)
  handles = guidata(VS_GUI);
  set(handles.rotleft,'Value',dd*sc);
  viewsurfacedata_gui('rotleft_Callback',0,[],handles);
  viewsurfacedata_autosnapshot2; movefile('image000.png',sprintf('%s/frame%03d.png',outdir,cnt));
  cnt = cnt + 1;
end

% roll animation
for p=1:(90/dd)
  handles = guidata(VS_GUI);
  set(handles.rollccw,'Value',-dd*sc);
  viewsurfacedata_gui('rollccw_Callback',0,[],handles);
  viewsurfacedata_autosnapshot2; movefile('image000.png',sprintf('%s/frame%03d.png',outdir,cnt));
  cnt = cnt + 1;
end

%%%%%% STITCH TOGETHER

% Hm, perhaps we could use ImageMagick to do this crop and concatenation in the future...

cd ~/Desktop
mkdirquiet('composite');
blank0 = uint8(zeros(800-2*40,800-2*40,3));
for p=1:900
%  im1 = placematrix(blank0,imread(sprintf('lhwhite/frame%03d.png',p)));
  im2 = placematrix(blank0,imread(sprintf('lhinflated/frame%03d.png',p)));
%  im3 = placematrix(blank0,imread(sprintf('rhwhite/frame%03d.png',p)));
  im4 = placematrix(blank0,imread(sprintf('rhinflated/frame%03d.png',p)));
%  imwrite([im1 im3; im2 im4],sprintf('composite/frame%03d.png',p));
  imwrite([im2 im4],sprintf('composite/frame%03d.png',p));
end

%%%%%% MAKE MOVIE:

QT 7 -> 30 fps -> Handbrake for Fast 720
  -> b3noiseceiling.mp4

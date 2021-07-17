%% %%%%% PHYSIO (just modify/hack physioinspections.m)

% create physio inspection figures
types = {'puls' 'resp'};
for qq=1:length(types)
  typ = types{qq};
  alldata = {};
  for subjix=1
    for sess=25
      for runix=1:14
        file0 = sprintf('~/nsd/nsddata_timeseries/ppdata/subj%02d/physio/physio_session%02d_run%02d_%s.tsv', ...
                        subjix,sess,runix,typ);
        a0 = load(file0);
        alldata{subjix,sess,runix} = a0;
      end

      figureprep([100 100 2000 1600]); hold on;
      for runix=1:14
        subplot(14,1,runix); hold on;
        data0 = alldata{subjix,sess,runix};
        vals = prctile(data0,[25 75]);
        iq = vals(2)-vals(1);
        plot(data0);
        ax = axis;
        if iq~=0
          axis([ax(1:2) vals(1)-2*iq vals(2)+2*iq]);
        end
      end
      figurewrite(sprintf('subj%02d_sess%02d',subjix,sess),[],-1,sprintf('~/Dropbox/KKTEMP/samplephysio%s',typ));

    end
  end
end

% subj 1, sess 25, run 2

%% %%%%% SAMPLE VOLUMES

% copy from inspections/coregistration
% crop to isolate subj06.

%% %%%%% HRT2

itksnap: 

subj08:
HRT2/NSD929_DATECENSORED_174115t2tsecoraniso1pt5mmR2s002a1001.nii.gz
range 0 3113
cursor 292 148 34

% high quality:
a1=load_untouch_nii('NSD929_DATECENSORED_174115t2tsecoraniso1pt5mmR2s002a1001.nii.gz');
imwrite(uint8(255*cmaplookup(rotatematrix(double(a1.img(:,:,34)),1,2,1),0,3113,[],gray(256))),gray(256),'~/Desktop/HRT2.png');

%% %%%%% SWI

subj08
get DATECENSORED-NSD929-nsd06/dicom/MR-SE036-Mag_Images

% raw magnitude data
a = dicomloaddir('MR-SE036-Mag_Images');
  %{384x384x256 double}
  %figure;imagesc(a{1}(:,:,100))

%do min intensity of the raw magnitude data
% 16 slices (at 0.6mm) = 9.6mm
for p=1:100
  imwrite(uint8(255*makeimagestack(min(double(a{1}(:,:,(p-1)*16+(1:16))),[],3),1)),sprintf('~/Dropbox/KKTEMP/image%03d.png',p));
end

% use image009.png

%% %%%%% Diffusion

subj08
nsd/ppdata/DATECENSORED_121153dMRIdir98APs012a001.nii.gz
itksnap.
0 1500
vol 27
72 67 51 

%% %%%%% Angiogram

% see folder

%% %%%%% surfaces

follow analysis_freesurfer.m for subj02
0.70 opacity for surfaceimperfections
click in OFC  (150 190 222)
save screenshot (don't hide anything) (magnification factor 1)
save axial, sagittal, coronal.

%% %%%%% prffloc

pull from the nsddata/experiments and take a screencapture

%% %%%%% nsddatacollection

simply print PDF from nsddatacollection.xlsx
then open in preview and take big screenshot

%% %%%%% eyetracking

from /research/nsd/eyetracking/samples/eyevideo-DATECENSORED-NSD258-nsd30.mp4
20:27
take screenshot

%% %%%%% resting

take screenshot from nsddata/experiments/

%% %%%%% slice

~/nsd/rawdata/DATECENSORED-NSD149-nsd08]$ ls
rsync -av slice.png ~/Dropbox/KKTEMP/

%% %%%%% leaderboard

just copy /Users/kendrick/Dropbox/nsd/autoqc_nsd/leaderboard_forsubjectsB.png

%% %%%%% also

SCREENING/FINALPARTICIPANTSUMMARY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

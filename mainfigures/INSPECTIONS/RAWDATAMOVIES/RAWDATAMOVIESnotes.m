=============== RAW DATA (rawdatamovies)

√ LOOK AT RAW DATA LIKE movie itk-snap of 1 random run. sess10, run06

ls -d *NSD134*nsd10/dicom/*MR*bold*

cp -r DATECENSORED-NSD134-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj01
cp -r DATECENSORED-NSD139-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj02
cp -r DATECENSORED-NSD149-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj03
cp -r DATECENSORED-NSD168-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj04
cp -r DATECENSORED-NSD258-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj05
cp -r DATECENSORED-NSD400-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj06
cp -r DATECENSORED-NSD814-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj07
cp -r DATECENSORED-NSD929-nsd10/dicom/MR-SE019-cmrr_mbep2d_bold_84slices_1pt8mm_1600TR_MB3_R2 ~/Dropbox/KKTEMP/subj08

dcm2nii -> rename the file  ->  itksnap

0 to 1500 for range
use arrow key to do it.
immediately trim the movie...

then convert using handbrake (super hq 720)

============= PP DATA

download from timeseries (1.8mm). same run.
0 to 1500
arrow key
trim
handbrake super hq 720.

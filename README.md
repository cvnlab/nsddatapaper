# nsddatapaper

This folder contains an archive of code used in the manuscript "A massive 7T 
fMRI dataset to bridge cognitive neuroscience and artificial intelligence", 
which is informally known as the Natural Scenes Dataset (NSD) data paper.

================= Part 1 (Author: Kendrick Kay)

"main" - This folder contains core analyses and preparation of the
         NSD data files. This includes data pre-processing, analysis of the
         NSD, fLoc, and pRF experiments, data metric calculations,
         noise ceiling estimation, and ROI definition.
         
"mainfigures" - This folder contains support for several of the figures and
         videos provided in the NSD data paper.

"prfestimation" - This folder contains files for the pRF estimation analysis
                  (results depicted in Extended Data Figure 9).

Note that the code relies on external repositories, including:
* https://github.com/kendrickkay/knkutils/
* https://github.com/kendrickkay/cvncode/
* https://github.com/kendrickkay/GLMdenoise/
* https://github.com/kendrickkay/preprocessfmri/
* https://github.com/kendrickkay/alignvolumedata/
* https://github.com/kendrickkay/viewsurfacedata/
* https://github.com/kendrickkay/analyzePRF/
* https://github.com/kendrickkay/nsdcode/
* https://github.com/nrdg/fracridge/

Some useful notes on the code underlying different figures:

```
Figure 1C: TIMELINEnotes.m
Figure 1D: NSDSUMMARYbehavior.m
Figure 2: ACQUISITIONnotes.m
Figure 3: MANIFOLDnotes.m
Figure 4: SCIENCEMEMORYnotes.m
Figure 5: SCIENCERSAnotes.m
Supplementary Figure 1: METRICS_TSNR.m
Supplementary Figure 2: METRICS_MOTION.m
Supplementary Figure 4: PAIRWISEnotes.m
Supplementary Figure 6: MAPOFSIGNALnotes.m
Supplementary Figure 7: SCIENCEMEMORYnotes.m
Extended Data Figure 1: NSDSUMMARYnotes.m
Extended Data Figure 5: UPSAMPLINGnotes.m
Extended Data Figure 7: ROInotes.m
Extended Data Figure 8: EXPOSITIONnotes.m
Extended Data Figure 9: PARAMETERSnotes.m
Supplementary Videos 1-10: INSPECTIONS
```

================= Part 2 (Author: Matthias Nau)

Eyetracking pre-processing and analysis code is available at main/nsd_et*.m

================= Part 3 (Author: Brad Caron)

Code scripts used to analyze the diffusion data are freely available on GitHub 
via the links provided in the "Cloud processing via brainlife.io" section of
the NSD Data Manual. All code scripts used to generate figures related to the
diffusion processing can be found on GitHub at https://github.com/bacaron/nsd-analyses.

Figures related to the diffusion analysis include
Extended Data Figure 6 and Supplementary Figure 5.

================= Part 4 (Author: Ben Hutchinson)

"nsd_recog" - This folder contains code for the behavioral analysis
              shown in Figure 4A.

================= Part 5 (Author: Ian Charest)

The code for the representational similarity analysis shown in Figure 5
is freely available at https://github.com/Charestlab/nsddatapaper_rsa

================= Part 6 (Author: Ghislain St-Yves)

The code for the neural network analyses shown in Figure 6 and 
Extended Data Figure 10 are available at https://github.com/styvesg/nsd.

manually make T1 T2 EPI loop, 4x.
T2 - SWI loop, 4x.
T1 - TOF loop, 4x.

then:
ffmpeg -framerate 1 -pattern_type glob -i '*.png' -crf 18 -c:v libx264 -pix_fmt yuv420p T1-T2-EPI.mp4
T2-SWI.mp4
T1-TOF.mp4

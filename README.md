# ProAnubisReconstruction
Repository for scripts to analyse data and reconstruct events from ProAnubis at CERN.

Functions necessary contained in ANUBIS_triggered_functions.py ("ANT") and anubisPlotUtils.py. Folder called scripts contains current working build.

Relevant functions are located in ANUBIS_triggered_functions.py

TODO:
- Fix expected muon angular distribution. ✔ Done

- Find eta and Phi angular distributions separtely. ✔ Done
  - As expected, the phi side showed symmetry in its angular distribution and the eta side showed asymmetry. Particularly eta side had a peak at around 50 degrees in one direction (this corresponds to a zenith angle of around 10 degrees).
  - This agrees with the findings from the zenith angle distribution of all the muons. (See angular distributions ppt in "Investigations" folder)
  
-  Heat map of three hour run.

-  Timing between RPCs. Calibrate channels and RPCs for timing diff. (Find TOF From trajectory and find dT distributions).
-  Also, find directionality of muon trajectory (vertical directionality) by considering average hit time for Top and Bottom hit RPCs. Trivial task, but at the momen not so important.
  
-  Find efficiency of RPCs using "Tag and Probe" Method.
  
   https://espace.cern.ch/lip/pub/docs/LIP-STUDENTS-20-18.pdf
   
   https://cms-opendata-workshop.github.io/workshop-lesson-tagandprobe/aio/index.html

   https://iopscience.iop.org/article/10.1088/1748-0221/14/10/C10020/pdf
   
-  Look at new 3 hour runs data sets. Compare with 3Hours_24_3_1 dataset. Find efficiences etc.

- What are the gas mixtures used ?

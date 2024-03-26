# ProANUBIS Reconstruction

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/06402d57-ecf8-4095-9aa8-93e6ff16fdb6)

Repository for scripts to analyse data and reconstruct events from ProAnubis 2024 at CERN.

Functions necessary contained in ANUBIS_triggered_functions.py ("ANT") and anubisPlotUtils.py. Folder called scripts contains current working build.

Relevant functions are located in ANUBIS_triggered_functions.py

Tasklist:

- Cluster events temporally and spatially. ✔ Done
  
- Reconstruct events. ✔ Done
    - Used SVD and minimise residual on possible combos to find best trajectory.
    - Sub-optimal: Still haven't accounted for uncertainties in hit coordinates,
 
- Create event displays, channel hit_maps an 3D plots ✔ Done

- Fix ProANUBIS muon zenith angle distribution. ✔ Done

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

     - Cannot exactly use the same "Tag and Probe" method outilned in the papers above. ProANUBIS triggers a readout whenever any 3 eta channels fire simultaneously within one window, so attempting to use reconstructed muons from these          events to test the efficiencies of the RPCs would greatly bias our efficiencies measurements.
     - Usually with "Tag and Probe" techniques, a decay near resonance is observed. Two decay products are produced, one acts as the "Tag" particle indicating a decay has occured. The expected path of the other decay product, the                "Probe", is then used to measure the efficiencies of detectors being studied.
     - I have another way of doing this... stay tuned.
   
- Look at new 3 hour runs data sets. Compare with 3Hours_24_3_1 dataset. Find efficiences etc.

- What are the gas mixtures used ?

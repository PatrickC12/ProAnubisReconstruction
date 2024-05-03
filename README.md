  # ProANUBIS Reconstruction

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/06402d57-ecf8-4095-9aa8-93e6ff16fdb6)

Repository for scripts to analyse data and reconstruct events from ProAnubis 2024 at CERN.

Mostly all in jupyter notebooks, saves you from having to reread data from HDF5 files.

Functions necessary contained in ANUBIS_triggered_functions.py ("ANT") and anubisPlotUtils.py. Folder called scripts contains current working build.

For example of how to run the code see in scripts folder file : ANUBIS_reconstruction_example_60s.ipynb. This jupyter notebook will run through reconstruction on a small data set.

Reconstruction Algorithm (quick):

1) Cluster events temporally, real particle events will be correlated in time. Should filter out most dark counts.
2) Cluster hits spatially into x and y coordinates on RPCs.
3) Generate all possible combinations of coordiante sets that could reconsctruct an event.
4) Fit each combination with SVD.
5) Find combination that minimises Chi2 value of reconstruction.


Angular Distributions:

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/f32c77a0-e70f-4c06-b33a-0bb9aa48e625)

Tasklist:

- Cluster events temporally and spatially. ✔ Done
  
- Reconstruct events. ✔ Done
    - Used SVD and minimise residual on possible combos to find best trajectory.
 
- Create event displays, channel hit_maps an 3D plots ✔ Done

- Fix ProANUBIS muon zenith angle distribution. ✔ Done

- Find eta and Phi angular distributions separtely. ✔ Done
  - As expected, the phi side showed symmetry in its angular distribution and the eta side showed asymmetry. Particularly eta side had a peak at around 50 degrees in one direction (this corresponds to a zenith angle of around 10 degrees).
  - This agrees with the findings from the zenith angle distribution of all the muons.
  
-  Heat map of three hour run

-  Timing between RPCs. Calibrate channels and RPCs for timing diff. (Find TOF From trajectory and find dT distributions). ✔ Done
-  Also, find directionality of muon trajectory (vertical directionality) by considering average hit time for Top and Bottom hit RPCs. Trivial task, but at the moment not so important. ✔ Done
  
-  Find efficiency of RPCs using "Tag and Probe" Method. ✔ Done
  
   https://espace.cern.ch/lip/pub/docs/LIP-STUDENTS-20-18.pdf
   
   https://cms-opendata-workshop.github.io/workshop-lesson-tagandprobe/aio/index.html

   https://iopscience.iop.org/article/10.1088/1748-0221/14/10/C10020/pdf

     - Cannot exactly use the same "Tag and Probe" method outilned in the papers above. ProANUBIS triggers a readout whenever any 3 eta channels fire simultaneously within a time window (?), so attempting to use reconstructed muons from these events to test the efficiencies of the RPCs would greatly bias our efficiencies measurements.
     - Usually with "Tag and Probe" techniques, a decay near resonance is observed. Two decay products are produced, one acts as the "Tag" particle indicating a decay has occured. The expected path of the other decay product, the "Probe", is then used to measure the efficiencies of detectors being studied.
     - With the current dataset, data is triggered by any group of 3 Eta channels across the 6 RPCs. Therefore, use eta side as trigger and probe if phi side also fires simultaneously. Find efficiency of phi side average is around 0.95.
     - Carried out biased measurement of efficiencies for 3 data sets: Cosmics 24_3_1, Cosmics 03/04/24, ATLAS Luminosity Spike 03/04/24.
     - Attempted to reduce bias as much as possible by requiring all 5 other RPCs (excluding RPC under test) to be in reconstructed event. This reduces the likelihood that the RPC under test was involved in triggering the data readout for an event. 
     - Such efficiences measurements can be done in situ with less bias once LHC is running and muons passing through ATLAS muon spectrometer can be used as Tags to probe ProANUBIS.
     - Furthermore, the new scintillator setup in the lab would allow efficiency variation with operational parameters to be determined.
   
- Looking at angular distributions during ATLAS Luminosity events on 03/04/24 ✔ Done
     - See if timing information can be used to determine direction of trajectory.  ✔ Done
     - Plot timing difference between layers for hit. ✔ Done
 
- Looking at timing distributions for dZ during ATLAS luminosity spike and cosmics only. ✔ Done
    - Compare two and three chamber events ✔ Done
 
    - 
Result of Dr Revering's simulations of cosmic flux in proANUBIS show good agreement with reconstructed events (in terms of distribution form).

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/385e1e88-a6a6-488b-b32f-d08311dfe845)


- Gas mixtures in ProANUBIS: ✔ Done

- Plot Residual Distributions. 

- PLOT RESIDUALS OF RECONSTRUCTED EVENTS

- proANUBIS gas mixture
    - 30% CO2
    - 70% CERN mixture:
      - 94.7% Freon
      - 5% Isobutane
      - 0.3% SF6
ProANUBIS using 5.8 - 6.0 kV
 
- Measure the efficiency of this gas mixture as a function of HV in the lab TODO

# ProANUBIS Reconstruction

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/06402d57-ecf8-4095-9aa8-93e6ff16fdb6)

Repository for scripts to analyse data and reconstruct events from ProAnubis 2024 at CERN.

Mostly all in jupyter notebooks, saves you from having to reread data from HDF5 files.

Functions necessary contained in ANUBIS_triggered_functions.py ("ANT") and anubisPlotUtils.py. Folder called scripts contains current working build. (Note still needs to be cleaned up!)

I will create an example notebook at some point to show how the code can be used. In the meantime, my workings in the various notebooks should be decipherable somewhat to elucidate what is happening. 

## Reconstruction Algorithm (quick):

1) Cluster events temporally, real particle events will be correlated in time. Should filter out most dark counts.
2) Cluster hits spatially into x and y coordinates on RPCs.
3) Generate all possible combinations of coordiante sets that could reconsctruct an event.
4) Fit each combination with SVD (inspired by PCA, https://arxiv.org/abs/1404.1100) 
5) Find combination that minimises Chi2 value of reconstruction.

I will compile a more comprehensive pseudo-code esque description of the reconstruction algorithm at some point.

### Example Result

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/1247454e-79fc-4051-9643-89fa04a117dc)

## Angular Distributions:

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/f32c77a0-e70f-4c06-b33a-0bb9aa48e625)

Project the reconstructed 

## Efficiencies:

The efficiencies of the RPCs can be measured in-Situ if so desired. Note that the data is only recorded when triggered, this will introduce some inherent bias into the efficiency measurements done with the RPCs. In an attempt to minimise this bias, the following in-situ efficiency check algorithm is being utilised:

1) Select a test RPC.
2) Reconstruct tracks using **all 5 other** RPCs. This reduces the chance of the RPC under test being involved in the trigger for the data read.
3) Extrapolate reconstructed track to height of RPC under test and see if the reconstructed track then passes throug the dimensions of the RPC.
4) Check that the RPC under test produces some useable coordinate within some tolerance of where the expected hit is.

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/e3c345aa-7f29-441d-9bad-a1e9f42d4faf)

As is seen from the above plot, the efficiencies of the RPCs plateau within ~ 5cm which corresponds to 2 strip widths in the RPC. Taking the limit of tolerance to larger values, the geometrical constraint on the reconstructed hit is being relaxed, but we are still seeing rather low efficiences! This may be an artifact of the selection constraints for hits in the RPC under test. At some point I will rerun this efficiency with a less constrained test RPC, if the low efficiences are still present then the issue does not lie with the reconstruction but perhaps the RPC or some aspect of the data acquisition (DAQ) system.

Once proANUBIS is synced with the LHC clock, a less biased tag and probe method can be used. This will involve tagging muons passing through the ATLAS detection volume and using these tagged muons to probe the proANUBIS RPC strips. This can achieved using resonance muon, anti-muon pairs (as is done in this paper: https://iopscience.iop.org/article/10.1088/1748-0221/14/10/C10020/pdf), or by using muons identified by the ATLAS detector, as proANUBIS lies outside the main detector volume !

Work has been carried out in the ATLAS project room to upgrade the scintillation counter setup. This has now reached a stage where the scintillation counter can be used to tag muons that can probe RPCs under test. This provides a form of offline efficiency testing of the RPC and will be useful for testing different operational parameter effects on performance:

  - High voltages.
  - Gas mixtures.
  - Threshold Voltages.
  - RPC geometries.

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/08ff0d09-3cf9-437b-8a74-26c6bd90a283)


## Other links:

ANUBIS Twiki: https://twiki.cern.ch/twiki/bin/view/ANUBIS/

## Tasklist:

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

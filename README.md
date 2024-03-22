# ProAnubisReconstruction
Repository for scripts to analyse data and reconstruct events from ProAnubis at CERN.

I will sort of the files soon. But currently the working code is located in TestEnv_Patrick

Relevant functions are located in ANUBIS_triggered_functions.py

TODO:
- Fix expected muon angular distribution. ✔ Done

- Find eta and Phi angular distributions separtely. ✔ Done
  - As expected, the phi side showed symmetry in its angular distribution and the eta side showed asymmetry. Particularly eta side had a peak at around 50 degrees in one direction (this corresponds to a zenith angle of around 10 degrees).
  - This agrees with the findings from the zenith angle distribution of all the muons. (See angular distributions ppt in "Investigations" folder)
  
-  Heat map of three hour run.
  
-  Find efficiency of RPCs using "Tag and Probe" Method.
   https://espace.cern.ch/lip/pub/docs/LIP-STUDENTS-20-18.pdf
   
-  Look at new 3 hour runs data sets. Compare with 3Hours_24_3_1 dataset. Find efficiences etc.

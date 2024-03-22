# ProAnubisReconstruction
Repository for scripts to analyse data and reconstruct events from ProAnubis at CERN.

I will sort of the files soon. But currently the working code is located in TestEnv_Patrick

Relevant functions are located in ANUBIS_triggered_functions.py

TODO:
- Fix expected muon angular distribution. ✔ Done
- Find eta and Phi angular distributions separtely. ✔ Done
  - As expected, the phi side showed symmetry in its angular distribution and the eta side showed asymmetry. Particularly eta side had a peak at around 50 degrees in one direction (this corresponds to a zenith angle of around 10 degrees).
  - This agrees with the findings from the zenith angle distribution of all the muons. (See angular distributions ppt in "Investigations" folder)

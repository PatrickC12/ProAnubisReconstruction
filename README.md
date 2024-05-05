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
4) Fit each combination with Singular Value Decomposition (SVD) (inspired by Principal Component Analysis (PCA), https://arxiv.org/abs/1404.1100) 
5) Find combination that minimises Chi2 (normalised w.r.t. DoF) value of reconstruction.

I will compile a more comprehensive pseudo-code esque description of the reconstruction algorithm at some point.

### Example Result

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/1247454e-79fc-4051-9643-89fa04a117dc)

## Angular Distributions:

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/f32c77a0-e70f-4c06-b33a-0bb9aa48e625)

Project the reconstructed trajectory vector onto the XZ and YZ planes of the proANUBIS setup.

Below is a plot of muon distribution for 3 hours of cosmics. The event reconstruction rate is ~ 0.33, so the counts match the amount you would expect in ~ 1 hour, the black line shows expected distribution of muons through proANUBIS in 1 hour ( from detailed simulation work done by Dr Revering, simulating muon trajectories through Earth and cavern entrance to proANUBIS).

(Filtering here is atleast 3 RPCs, no constraint on cross-chamberdness)

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/0ad34672-e9ee-4d7a-b00c-dc5799accfdb)

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/a6511c0e-33cc-45a1-a2ad-2d7547cd27ce)

The RPCs in the triplet layer of proANUBIS are close together (~12 mm between RPC centres), so the granularity of the detectors have a significant effect on the produced angular distributions.
This causes peaks in the counts at 0 degrees, 32 degrees, 51 degrees and 68 degrees. 

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/acc2c0b7-ccfe-467b-9532-2c61d42a3d58)

The majority of the cosmic events pass through the triplet layer of proANUBIS only, so enforcing cross chamber is a good filter to cut off the cosmic ray muons. The geometry of proANUBIS is such that 60 degrees is the maximum angle (eta) for cross chamber events. The plots below enforce the cross chamber constraint. 

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/99a41b3c-7842-4198-9849-6eaadc469f4a)

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/c00c0e2b-c6c7-4b94-a771-c609224e02e7)

## Efficiencies:

The efficiencies of the RPCs can be measured in-Situ if so desired. Note that the data is only recorded when triggered, this will introduce some inherent bias into the efficiency measurements done with the RPCs. In an attempt to minimise this bias, the following in-situ efficiency check algorithm is being utilised:

1) Select a test RPC.
2) Reconstruct tracks using **all 5 other** RPCs. This reduces the chance of the RPC under test being involved in the trigger for the data read.
3) Extrapolate reconstructed track to height of RPC under test and see if the reconstructed track then passes throug the dimensions of the RPC.
4) Check that the RPC under test produces some useable coordinate within some tolerance of where the expected hit is.

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/e3c345aa-7f29-441d-9bad-a1e9f42d4faf)

As is seen from the above plot, the efficiencies of the RPCs plateau within ~ 5cm which corresponds to 2 strip widths in the RPC. Taking the limit of tolerance to larger values, the geometrical constraint on the reconstructed hit is being relaxed, but we are still seeing rather low efficiences! This may be an artifact of the selection constraints for hits in the RPC under test. 

A more relaxed probe of the test RPCs were carried out, the reconstruction of the tagged muon was still kept strict (all 5 other RPC, 15ns time window and max cluster size 3). No restriction was placed on cluster sizes or the location of a hit, geometrical effects were ignored. Rather a temporal tolerance was placed on the test RPC to ask if it produced a hit (anywhere) within some tolerance of the reconstructed event. The results for the doublet suggest an issue with the TDC chips, indeed the TDC EE04 was entirely resposinble for DoubletHigh phi side and had the largest number of corrupted timestamps of all the TDCs. (Corrupted timestamp here just means a channel recorded a time outside of 0-1250 ns, which is the window in which data is readout on trigger). 

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/2a03720e-b6aa-45af-bbc3-7540832eaa64)

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/80d2e6f7-e6f9-43cf-9a28-7deaa0f5734f)



Once proANUBIS is synced with the LHC clock, a less biased tag and probe method can be used. This will involve tagging muons passing through the ATLAS detection volume and using these tagged muons to probe the proANUBIS RPC strips. This can achieved using resonance muon, anti-muon pairs (as is done in this paper: https://iopscience.iop.org/article/10.1088/1748-0221/14/10/C10020/pdf), or by using muons identified by the ATLAS detector, as proANUBIS lies outside the main detector volume !

Work has been carried out in the ATLAS project room to upgrade the scintillation counter setup. This has now reached a stage where the scintillation counter can be used to tag muons that can probe RPCs under test. This provides a form of offline efficiency testing of the RPC and will be useful for testing different operational parameter effects on performance:

  - High voltages.
  - Gas mixtures.
  - Threshold Voltages.
  - RPC geometries.

![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/08ff0d09-3cf9-437b-8a74-26c6bd90a283)

## Todo:

- The timing of hits needs to be looked into. Current analysis suggests a convolution of a large number of errors, e.g. different FEBs, different TDCs. The analsysis of the timing distributions is really only to first order currently, the first step to improve this would be to include effects of timewalk of signal 

## Other links:

ANUBIS Twiki: https://twiki.cern.ch/twiki/bin/view/ANUBIS/

## Misc:

### proANUBIS location in ATLAS cavern
![image](https://github.com/PatrickC12/ProAnubisReconstruction/assets/123903514/bdbf3f35-e784-494c-bad7-f01dd988a745)

- proANUBIS gas mixture
    - 30% CO2
    - 70% CERN mixture:
      - 94.7% Freon (ionising gas)
      - 5% Isobutane (photon quencher)
      - 0.3% SF6 (electron quenncher)
      - 
ProANUBIS using 5.8 - 6.0 kV HV.
Work being carried out on measurement of RPC efficiency using this (and similar) gas mixture in the ATLAS project room in the Cavendish Lab.


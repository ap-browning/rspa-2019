# rspa-2019
Supplementary code for "Reversible signal transmission in an active mechanical metamaterial" by Browning et al. (2019).

Code was created in `MATLAB R2019a` (Mathworks).

* Run each section of `ExampleUsage.m` to see:
  1. A solution to the discrete model, corresponding to Figure 1e and Figure 1f of the main text
  2. A solution to the discrete model, corresponding to Figure 3b and Figure 3f of the main text. Note that the numerical parameters must be refined as per table S2 of the supporting material document to obtain accurate figures and an accurate estimate of the wavespeed.

* `SpringsDiscrete.m` Solve the continuous metamaterial model described by Eqs 2.1 to Eq 2.3 in the main document.
* `SpringsContinuous.m` Solves the continuous metamaterial model described by Eq 2.6 and Eq 2.7 in the main document.
* `EstimateWavespeed.m` Estimates the transmission speed (wavespeed) using output from `SpringsContinuous.m`.

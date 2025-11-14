Description of scripts / data usage for Glycine project-
Many of these functions require GANNET and/or FID-A, which can be found at: 

https://github.com/markmikkelsen/Gannet

https://github.com/CIC-methods/FID-A

All data files / LCModel scripts can be found at: https://drive.google.com/drive/folders/1agYTJLwY6icosPg6LDW0d-fIPu1lXpgC?usp=sharing

.
.
.
.

---------------------------------------------------------- SCRIPTS ---------------------------------------------------------

SDAT2MAT.m - Runs a loop through Gannet to convert .SDAT to .mat for in vivo spectra. 

RMS_Figure_2.m - Calculates and plots RMS of metabolites (3.55 ppm range) as function of TE

scalemetab.m - function used in RMS_Figure_2 script which simply scales to concentration

Gly_mI_Figure_3.m - Script to plot spectra near 3.55 as a function of [Gly]/[mI] and Linebroadening

GABAsimAreas_Figure_4 - Script to analyze GABA simulated difference spectra as a function of TE.



--------------------------------------------------------  Data (all zip files) -------------------------------------------------------

Gly_mI_sims - .MAT files used in Gly+mI simulations (edit OFF / SUM) TE 64 & 68 ms with variable concentration ratios and exponential LB

GABA_TE_series - .MAT files used for GABA area simulations with TE-modulation and T2-modulation

Gly_mI_TE_series - .MAT files used to analyze RMS contribution in 3.5 - 3.6 ppm range (Glucose, Threonine, Gly, mI) TEs used = 60-88 ms (2 ms step size)

LCModel_Data_and_Code - example control files for edit OFF and SUM (TE 64ms and 68ms), all .RAW files used with anonymous naming, and basis sets used.

Volunteer_Data_Anon - Anonymized Philips (.SDAT, .SPAR) data for all subjects in this study.


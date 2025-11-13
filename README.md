Description of script /simulated data usage for Glycine project-
Some scripts perform almost identically, but may be preferred depending on the objective.
Code and metioned data is in the zip file
Nearly all of functions use Gannet and/or FID-A, which can be found at: 

https://github.com/markmikkelsen/Gannet

https://github.com/CIC-methods/FID-A
.
.
.
.
------------------------------------------------------- GENERAL SCRIPTS ---------------------------------------------------------


WritedotMATjrs_v2.m - Runs loop through Gannet to convert .SDAT to .mat for in vivo spectra, also performs manual phase correction and freq correction. 

WriteRaw4LCMjrs.m - Writes the .mat spectra to LCM .raw file


RMS_simulation_forpaper.m - Calculates and plots RMS of metabolite as function of TE


plotting_coords_v2.m - Plots the .coord file outputs from LCM; used to generate Fig 8 and 5. 


DiffSpecdotMATs_v2.m - Extracts GABA/Glx areas and concentrations; plots results and runs statistics


GABAGlxsimAreas.m - GABA/Glx areas, plotting the multiplet shape and the areas, also has spatial plotting option similar to script below.


generatespatialplots.m - takes FID-A outputs and plots spectra across the 19x19 simulation grid (36mm^2 total area) to visualize transition BW effects


LCMcsvStats.m - Main LCM processing script given CSV outputs; Concentrations, CVs, CRLBs
 

CRLB_comp_GlyvsNoGly & NoGlyAnalysis_v2.m - Stats for Glycine in/exclusion; can be modified for use of removal or inclusion of any metabolite 


Gly_mI_simulations.m - Plots the simulted specrta as a function of Gly/mI concentration and LW. 
 

-------------------------------------------------------- "Short" functions -------------------------------------------------------- 


phase_correction_gui - inputs are spec and ppm axis. This loads an interactive plot where you adjust the (0 order) with the slider, arrows, or typing the desired angle. The phase angle is the output of this script so actual correction requires [.*exp(1i * phase_angle * pi / 180)]. Or just use the phshift function bc I got tired of typing that constantly. 


pshift.m - this is needed for the gui; it simply performs the calculation of the spectra given the phase angle you are inputting in real time to visualize the phase. 


crSNR.m - Takes input spec and calculates SNR from Cr peak. 


freqadjust - Finds theoretical NAA index and shifts the spectrum such that NAA peak aligns with it (requires non-saturated NAA peak)

	
scalemetab.m - Multiplies simulated metabolite by reported / observed in vivo concentration. 


adjustnoise.m - Adds noise to simulated spectra based on PDF of noise in corresponding in vivo spectra, noise factor is modulated by the voxel volume input


plist.m - Creates an array of subjects like: P1, P2, etc... based on input (start and stop) It can be altered for different prefixes.



-------------------------------------------------------- Simulated data (all zip files) -------------------------------------------------------- 

Gly_mI_sims - .MAT files used in Gly+mI simulations (edit OFF / SUM) TE 64 & 68 ms with variable concentration ratios and exponential LB

GABA_TE_series - .MAT files used for GABA area simulations with TE-modulation and T2-modulation

Gly_mI_TE_series - .MAT files used to analyze RMS contribution in 3.4 - 3.7 ppm range (Glucose, Threonine, Gly, mI) TE's used = 60-88 ms (2 ms step-size)


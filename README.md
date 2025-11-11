Description of script usage for Glycine project-
Some scripts perform almost identically, but may be preferred depending on the objective.
Almost all of these functions use Gannet and/or FID-A in one way or another: 
https://github.com/markmikkelsen/Gannet
https://github.com/CIC-methods/FID-A


------------------------------------------------------- GENERAL SCRIPTS ---------------------------------------------------------


WritedotMATjrs_v2.m - Runs loop through Gannet to convert .SDAT to .mat for in vivo spectra, also performs manual phase correction and freq correction. 


WriteRaw4LCMjrs.m - Writes the .mat spectra to LCM .raw file


RMS_simulation_forpaper.m - RMS of metabolite as function of TE


plotting_coords_v2.m - Plots the .coord file outputs from LCM; used to generate Fig 8 and 5. 


DiffSpecdotMATs_v2.m - extracts GABA/Glx areas and concentrations; plots results and runs statistics


GABAGlxsimAreas.m - GABA/Glx areas, plotting the multiplet shape and the areas, also has spatial plotting option similar to script below.


generatespatialplots.m - takes FID-A outputs and plots spectra across the 19x19 simulation grid (36mm^2 total area) to visualize transition BW effects


LCMcsvStats.m - Main LCM processing script given CSV outputs; Concentrations, CVs, CRLBs
 

CRLB_comp_GlyvsNoGly & NoGlyAnalysis_v2.m - stats for Glycine in/exclusion; can be modified for use of removal or inclusion of any metabolite 


Gly_mI_simulations.m - Script that plots the simulted specrta as a function of Gly/mI concentration and LW. Example .MAT outputs are included in "Gly_mI_sims" folder 

 

-------------------------------------------------------- "Short" functions -------------------------------------------------------- 


phase_correction_gui - inputs are spec and ppm axis. This loads an interactive plot where you adjust the (0 order) with the slider, arrows, or typing the desired angle. The phase angle is the output of this script so actual correction requires [.*exp(1i * phase_angle * pi / 180)]. Or just use the phshift function bc I got tired of typing that constantly. 


pshift.m - this is needed for the gui; it simply performs the calculation of the spectra given the phase angle you are inputting in real time to visualize the phase. 


crSNR.m - takes input spec and calculates SNR from Cr peak. 


freqadjust - finds theoretical NAA index and shift spec such that NAA peak aligns with it

	
scalemetab.m - simply multiplies simulated metabolite by reported in vivo concentration. 


plist.m - Creates an array of subjects like: P1, P2, etc... based on input (start and stop) It can be altered for different prefixes






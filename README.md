# radio-imaging
Queries: salmoli.ghosh@gmail.com
Author: Salmoli Ghosh

The uGMRT_band4_polarization_pipeline.py is the pipeline used to image polarized emission from astronomical sources like active galaxies, etc. using uGMRT band 4 (550 - 900 MHz) data.
The uGMRT_band5_polarization_pipeline.py is the pipeline used to image polarized emission from astronomical sources like active galaxies, etc. using uGMRT band 5 (1050 - 1450 MHz) data.
The VLA_polarization_pipeline is the modified version of the EVLA polarization pipeline by Alice Pasetto.

This version of the pipeline is tested for CASA 6.2, and can be used safely for versions above CASA 6.

The files gvfits and listscan should be made executable before running them as commands. ($chmod +x listscan, $chmod +x gvfits)

Continuum Imaging Pipeline modified by Ishwara Chandra in 2018 
Kindly refer to the paper: Ishwara-Chandra et al 2020  

Incorporation of polarization calibration done by Silpa Sasikumar in 2019-2020
Kindly refer to Silpa et al. 2021 (https://ui.adsabs.harvard.edu/abs/2021MNRAS.507..991S/abstract)
 
Pipeline modifications by Janhavi Baghel in 2020-2021
Kindly refer to the script https://github.com/jbaghel/Improved-uGMRT-polarization-pipeline/blob/main/Improved_uGMRT_POL.py 
 
Pipeline modifications by Salmoli Ghosh in 2023. Suggestions taken from CAPTURE pipeline (Kale and Ishwara-Chandra, 2021) and Alice Pasetto's EVLA polarization pipeline. The following changes are made:

(1)Omitting initial clipping on the uncalibrated data

(2)Modifying the functions for polarization calibration

(3)Removing bad channels from the data by splitting (for reducing flagging percentage for band-4)

(4)Baseline dependent flagging for target

(5)Running setjy for the unpolarized calibrator

(6)Clipping high data points from leakage solutions

(7)Optimising parameter values for functions like flagdata, gaincal, polcal, tclean. Keeping mostly default values

(8)Optimising imaging parameters in tclean

(9)ADDING POLARIZATION STEPS FOR uGMRT BAND-5 (1050-1450 MHz), SOLVING EQUATIONS IN LINEAR (X&Y) BASIS



Please inspect your acquired UV data using listobs & plotms in CASA before running any further data reduction step.
Please CHANGE the input parameters as per your data.

In BAND-4, recommended channels are corresponding to ~ 560 MHz to ~ 810 MHz. The sensitivity drops sharply after 810 MHz. In any case, DO NOT use beyond 820 MHz  (--- Silpa S.).
In Band-5, around 1300 MHz channels are usually affected by satellite RFI.

It is highly recommended not to use OQ208 (unpolarized calibrator) to calculate the instrumental leakage for uGMRT since it is a very faint source (~few Jy); therefore  single short scan does not provide sufficient SNR to accurately determine the instrumental polarization.

Also, we do not recommend the use of 3C138 (polarized calibrator) for leakage calibration.
We recommend 3C286 (polarized calibrator) or 3C84/3C147 (unpolarized calibrator) for leakage calibration.



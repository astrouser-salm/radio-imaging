# radio-imaging
Queries: salmoli.ghosh@gmail.com
Author: Salmoli Ghosh

The uGMRT_band4_polarization_pipeline.py is the pipeline used to image polarized emission from astronomical sources like active galaxies, etc. using uGMRT band 4 (550 - 900 MHz) data.


Pipeline originally developed by Russ Taylor in 2011

Modified by Ishwara Chandra in 2018 
Improved version of the polarization pipeline by Silpa Sasikumar
Major improvements in flagging and self-calibration 
Kindly refer to the paper: Ishwara-Chandra et al 2020  

Initial phase only calibration added by Silpa Sasikumar in 2019-2020
Polarization steps added by Silpa Sasikumar in 2019-2020, which involve:
(1) Flagging of polarized and unpolarized calibrators
(2) Polarization calibration (cross-hand delays; leakage terms; R-L polarization angle)
(3) Stokes 'Q' and 'U' imaging
The current version of the pipeline also flags all the four correlations (RR, LL, RL, LR). 
 
Pipeline modifications by Janhavi Baghel in 2020-2021, changes made:
(1)refantmode = 'flex' and only one reference antenna to be specified.
(2)datacoulmn = 'corrected' in the flagging of target field before splitting
(3)datacoulmn = 'data' in tclean during creation of dirty image. Split target file has no datacolumn 'corrected'.
(4)datacoulmn = 'RESIDUAL_DATA' when flagging residuals in self-calibration cycles
(5)Initial flagging of known bad antenna from observer log.
(6)Polarization modelling included within the pipeline for 3C286,3C48,3C138 with the pol_*.txt files.
TEST.FITS is the output of gvfits, which coverts GMRT LTA format to multi-source FITS
 
Pipeline modifications by Salmoli Ghosh in 2023. Suggestions taken from CAPTURE pipeline (Kale and Ishwara-Chandra, 2021) and Alice Pasetto's EVLA polarization pipeline. The following changes are made:
(1)Omitting initial clipping on the uncalibrated data
(2)Modifying the functions for polarization calibration
(3)Removing bad channels from the data by splitting (for reducing flagging percentage for band-4)
(4)Baseline dependent flagging for target
(5)Running setjy for the unpolarized calibrator
(6)Clipping high data points from leakage solutions
(7)Optimising parameter values for functions like flagdata, gaincal, polcal, tclean. Keeping mostly default values

Please CHANGE channels and source fields as per your data.
Also change clip parameters if you have much stronger calibrator and/or target

The parameters below are typical for 550 MHz, 2048 channels at Band-4 (550-750 MHz)
Please change as required for your data.
In BAND-4, recommended channels corresponding to ~ 560 MHz to ~ 810 MHz. The sensitivity drops sharply after 810 MHz.
In any case DO NOT use beyond 820 MHz.
It is highly recommended not to use OQ208 (unpolarized calibrator) to calculate the instrumental leakage for uGMRT since it is a very faint source (~few Jy); therefore  single short scan does not provide sufficient SNR to accurately determine the instrumental polarization.
Also, we do not recommend the use of 3C138 (polarized calibrator) for leakage calibration.
We recommend 3C286 (polarized calibrator) or 3C84/3C147 (unpolarized calibrator) for leakage calibration.



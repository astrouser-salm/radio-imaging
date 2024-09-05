'''
# Pipeline originally developed by Russ Taylor in 2011
# Modified by Ishwara Chandra in 2018 
# Improved version of the polarization pipeline by Silpa Sasikumar, Janhavi Baghel, and Salmoli Ghosh
# Major improvements in flagging and self-calibration 
# Kindly refer to the paper: Ishwara-Chandra et al 2020  
#
# Initial phase only calibration added by Silpa Sasikumar in 2019-2020
# Polarization steps added by Silpa Sasikumar in 2019-2020, which involve:
# (1) Flagging of polarized and unpolarized calibrators
# (2) Polarization calibration (cross-hand delays; leakage terms; R-L polarization angle)
# (3) Stokes 'Q' and 'U' imaging
# The current version of the pipeline also flags all the four correlations (RR, LL, RL, LR). 
# 
# Pipeline modifications by Janhavi Baghel in 2020-2021, changes made:
# (1)refantmode = 'flex' and only one reference antenna to be specified.
# (2)datacoulmn = 'corrected' in the flagging of target field before splitting
# (3)datacoulmn = 'data' in tclean during creation of dirty image. Split target file has no datacolumn 'corrected'.
# (4)datacoulmn = 'RESIDUAL_DATA' when flagging residuals in self-calibration cycles
# (5)Initial flagging of known bad antenna from observer log.
# (6)Polarization modelling included within the pipeline for 3C286,3C48,3C138 with the pol_*.txt files.
# TEST.FITS is the output of gvfits, which coverts GMRT LTA format to multi-source FITS
# 
# Pipeline modifications by Salmoli Ghosh in 2023. Suggestions taken from CAPTURE pipeline (Kale and Ishwara-Chandra, 2021) and Alice Pasetto's EVLA polarization pipeline. The following changes are made:
# (1)Omitting initial clipping on the uncalibrated data
# (2)Modifying the functions for polarization calibration
# (3)Removing bad channels from the data by splitting (for reducing flagging percentage for band-4)
# (4)Baseline dependent flagging for target
# (5)Running setjy for the unpolarized calibrator
# (6)Clipping high data points from leakage solutions
# (7)Optimising parameter values for functions like flagdata, gaincal, polcal, tclean. Keeping mostly default values
#
# Please CHANGE channels and source fields as per your data.
# Also change clip parameters if you have much stronger calibrator and/or target
#
# The parameters below are typical for 550 MHz, 2048 channels at Band-4 (550-750 MHz)
# Please change as required for your data.
# In BAND-4, recommended channels corresponding to ~ 560 MHz to ~ 810 MHz. The sensitivity drops sharply after 810 MHz.
# In any case DO NOT use beyond 820 MHz.
# It is highly recommended not to use OQ208 (unpolarized calibrator) to calculate the instrumental leakage for uGMRT since it is a very faint source (~few Jy); therefore a single short scan does not provide sufficient SNR to accurately determine the instrumental polarization.
# Also, we do not recommend the use of 3C138 (polarized calibrator) for leakage calibration.
# We recommend 3C286 (polarized calibrator) or 3C84/3C147 (unpolarized calibrator) for leakage calibration.
# Queries: salmoli.ghosh@gmail.com
'''
######################################################################################################################################################################################################
#listscan filename.lta

#gvfits filename.log

#open casa

print ("Starting conversion of TEST.FITS to multi.ms")
importgmrt(fitsfile='TEST.FITS', vis='multi.ms')

ms='multi.ms'

print("Listing the observations") 
listobs(vis=ms)


#check the following things in plotms to know how your data looks like

#plot (rr,ll correlations only) amp vs time (without averaging over channels, if reqd. plot for a single channel), amp vs channel for each field iterating over each antenna. Also check the phase vs time/channel for the flux calibrators.

#you can decide the clipping fluxes based on the amplitudes for each field
 

######################################################################################################################################################################################################

print ("Initializing parameters") #Change these according to your data


fluxfield      = '0,3'            # field number of the primary flux calibrator
fluxfieldind=fluxfield.split(',') #this is useful if you have more than one flux calibrators                     ------ Salmoli Ghosh
phasefield = '1'                  # field number of the phase calibrator
polcalib = '3'                    # field number of the polarized calibratorar 1
unpolcalib = '4'                  # field number of the unpolarized calibrator 1
target         = '2'              # If more than one target, use target='2,3,4...', 
badants = 'C02,E02,S03'           #check observer's log for the bad antennas                                     ------ Salmoli Ghosh
badchans ='0:1168~2047'           #persistant RFI after 778 MHz, also sensitivity drops after 810 MHz. Can also consider till 750 MHz.
splitspw = '0:0~1167'             # split channel range if you want to split out only the good channels          ------ Salmoli Ghosh
samptime=2.6                      #sampling time in seconds                                                      ------ Salmoli Ghosh
gainspw = '0:50~1100'             # central ~ 75% good channel range for calibration
specave = 7                       # number of channels to average; suggested post-average BW (approx)
                                 
timeave = '0s'                    # time averaging
gainspw2 = '0:7~150'              # central good channels after split for self-cal

bpassfield     = '3'              # field number of the bandpass calibrator which has been observed almost throughout the observations, add phasecal if strong enough

gaincals       ='0,1,3,4'         # flux and phase calibrators
allcals        ='0,1,3,4'         # all calibrators
kcorrfield     = '3'              # field number of the antenna-based delay calibrator (the primary flux calibrator)

ref_ant='C00'                     # use only one option and refantmode 'flex'. Choose the reference antenna such that it is close to the centre, but not shadowed and shows good response.
                                  # C00, C01 are one of the best options for uGMRT (not much chances of shadowing)------ Salmoli Ghosh			
#
transferfield = '1,4'             #other calibrators except flux calibrators
#Can add more polarization calibrators to compare and verify the polarization calibration ----- Janhavi Baghel


# Change limits as required
clipfluxcal2 =[0.0,70.0]          #for the fluxcal 3C147; atleast twice the expected flux; only to remove high points
clipphasecal2 =[0.0,30.0]         #for phasecal 1150-003; atleast twice the expected flux; only to remove high points
clippolcalib2 =[0.0,60.0]         #for polcal 3C286; atleast twice the expected flux; only to remove high points
clipunpolcalib2 =[0.0,30.0]       #for unpolcal OQ208; atleast twice the expected flux; only to remove high points
clipunpolcalib =[0.0,50.0]        # atleast twice the expected flux; only to remove high points
#clipanofield = [0.0,150.0]       # atleast twice the expected flux; only to remove high points ----- Janhavi Baghel
cliptarget =[0.0,60.0]            # atleast four times the expected flux; only to remove high points
clipresid=[0.0,10.0]              # 10 times the rms for single channel and single baseline
uvracal='>1.0klambda'             # uvrange for gain calibration in main calibration
uvrascal='>0.75klambda'           # uvrange for gain calibration in self calibration
#

# Filenames for initial round of calibration
kcorrfile0 = ms+'.kcal0'
bpassfile0 = ms+'.bcal0'
gainfilep0 =  ms+'.gcalp0'
gainfile0 =  ms+'.gcal0'
fluxfile0 =  ms+'.fluxscale0'
#
# Filenames for final round of calibration
kcorrfile= ms+'.kcal'
bpassfile= ms+'.bcal'
gainfilep =  ms+'.gcalp'
gainfile=  ms+'.gcal'
fluxfile=  ms+'.fluxscale'
#
#Filenames for polarization calibration
#Add as many polarization calibrators you have----- Janhavi Baghel
#For polarized calibrator 1
kcross1 =    ms+'.kcross1'
leakage1 =   ms+'.leakage1'
polang1 =    ms+'.polang1'
#For polarized calibrator 2  #untag these if you have more than one polarized calibrator
#kcross2 =    ms+'.kcross2'
#leakage2 =   ms+'.leakage2'
#polang2 =    ms+'.polang2'
#
#Can also do leakage calibration with unpolarized calibrator----- Janhavi Baghel
unpolleakage1 = ms+'.unpolleakage1'
unpolleakage1a = ms+'.unpolleakage1a'
unpolleakage2 = ms+'.unpolleakage2'

######################################################################################################################################################################################################
#Imaging parameters

pcycles=4              # number of phase-only self-calibration
apcycles=3             # number of amplitude and phase self-calibration
doflag=True            #
sol_int=16.0           # this transaltes to 8, 4, 2 and 1 min for each selfcal loop
apsolint=8.0           # this transaltes to  4, 2 and 1 min for each selfcal loop
startthreshold=0.1     # start threshold to stop flux (mJy, will reduce by count subsequently)
startniter=2500        # start iterations (will double in each phase-selfcal and 4 times in each A&P loop)
imagesize=[4096,4096]  # should cover alteast up to null at lower part of the band. For computational time benefits use a multiple of 2, 3 and 5
cellsize='1arcsec'     # should be atleast 3 pixels in minor axis
wproj=-1               # w projection, default autocalculate

######################################################################################################################################################################################################
#Remove the unusable channels at band-4 from the data                            ------ Salmoli Ghosh
default(split)                                                                
split(vis=ms,outputvis='multi_split.ms',datacolumn='DATA',spw=splitspw)
ms='multi_split.ms'
######################################################################################################################################################################################################




######################################################################################################################################################################################################
#
#Defining functions for polarization calibration  (added/modified by Salmoli Ghosh; taken from Alice Pasetto's EVLA polarization pipeline)
#
######################################################################################################################################################################################################

msmd.open(ms)
fieldnames=msmd.fieldnames()
print ("fieldnames:", fieldnames) #getting the field names from the ms file

print ("Measurement set contains :")
print ("Target : " + fieldnames[int(target)]) 

print ("Flux Calibrator(s):" + fieldnames[int(fluxfieldind[0])],',',fieldnames[int(fluxfieldind[1])]) #check the no. of flux calibrators

print ("Phase Calibrator:" + fieldnames[int(phasefield)])

print ("Polarization Calibrator:" + fieldnames[int(polcalib)])

print ("Polarization Calibrator 2 (unpolarized): " + fieldnames[int(unpolcalib)])

reffreqdef = msmd.reffreq(0)  #getting the reference frequency from the ms file
msmd.done()

reffreqval = (reffreqdef['m0']['value'])
reffrequnit = (reffreqdef['m0']['unit'])
reffreq = str(reffreqval)+reffrequnit

print('reffreq=',reffreq)

import numpy as np
from scipy.optimize import curve_fit
import os


#Reference frequency in GHz
if reffrequnit =="Hz":
	f0=reffreqval /(10**9)
elif reffrequnit =="MHz":
	f0=reffreqval/(10**3)
elif reffrequnit =="GHz":
	f0=reffreqval

#Defining total intensity spectrum of the Flux calibrator
def S(f,Sf,alpha,beta):
        return Sf*(f/f0)**(alpha+beta*np.log10(f/f0))

#Defining fractional polarization spectrum of the Flux calibrator
def PF(f,a,b,c):
        return a+b*((f-f0)/f0)+c*((f-f0)/f0)**2

#Defining polarization angle of the Flux calibrator
def PA(f,a,b):
        return a+b*((f-f0)/f0)
        
        

######################################################################################################################################################################################################
#
#Setting the models 
#
######################################################################################################################################################################################################


#Model for SETJY using the flux calibrator using existing observed values at different frequencies

print("At frequency (in GHz)=",f0)

data = np.loadtxt('Calibrator_data/'+fieldnames[int(polcalib)]+'_2019.txt')
data_index_range=np.where(data[:,0]<3.0)    #change if reffreq lies in different frequency range
data_to_use=data[data_index_range]
popt, pcov = curve_fit(S, data_to_use[:,0], data_to_use[:,1])
print('I = ', popt[0], ' Jy')               #flux density 
print('alpha = ', popt[1])                  #spectral index
print('beta = ', popt[2])
print( 'Covariance : ', pcov)

    
# Stokes I flux density
I = popt[0] 
# Spectral Index and beta coeff
alpha = [popt[1], popt[2]]
    
######################################################################################################################################################################################################
#Model for Polarization fraction: obtaining value at reference frequency

popt, pcov = curve_fit(PF, data_to_use[:,0], data_to_use[:,2])
print("Polarization fraction polynomial : ", popt)
print("Covariance : ",pcov)

    
# Polarization Fraction coeff
polfrac = [*popt]
    
######################################################################################################################################################################################################
#Model for Polarization angle: obtaining value at reference frequency

popt, pcov = curve_fit(PA, data_to_use[:,0], data_to_use[:,3])
print("Polarization angle polynomial : ", popt)
print("Covariance : ", pcov)

    
# Polarization Angle
polangle = [*popt]



######################################################################################################################################################################################################
#
#  Flagging (Mode: manual and quack)              (added by Salmoli Ghosh)
#
######################################################################################################################################################################################################
print ("First Round of Flagging")

default(flagdata)
print("Flagging bad antennas")
flagdata(vis=ms, mode='manual', antenna=badants, flagbackup=True)    

default(flagdata)
print("Flagging the first and last sample from each scan")
flagdata(vis=ms,mode='quack',quackinterval=samptime,quackmode='beg',flagbackup=True)
flagdata(vis=ms,mode='quack',quackinterval=samptime,quackmode='endb',flagbackup=True)

default(flagdata)
print("Flagging the first channel (usually a bad channel)")
flagdata(vis=ms, mode='manual', spw='0:0', flagbackup=True) 
#print("Flagging other bad channels")
#flagdata(vis=ms, mode='manual', spw=badchans, flagbackup=False)  #check if you require the backup

#Check the flagdata percentage at this stage
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="DATA", name=ms+'summary1.split', action="calculate", overwrite=True, writeflags=True) 


'''
######################################################################################################################################################################################################
#
#  Flagging (Mode: clip)                       (Omitted by Salmoli Ghosh as clipping uncalibrated data may be harmful)
#
######################################################################################################################################################################################################

default(flagdata)
#Flag using 'clip' option to remove high points for calibrators
print ("Clipping high and zero values for flux calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=fluxfield, clipminmax=clipfluxcal,
        datacolumn="DATA",clipoutside=True, clipzeros=True, extendpols=False, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
print ("Clipping high and zero values for phase calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=phasefield, clipminmax=clipphasecal,
        datacolumn="DATA",clipoutside=True, clipzeros=True, extendpols=False, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
print ("Clipping high and zero values for polarized calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=polcalib, clipminmax=clippolcalib,
        datacolumn="DATA",clipoutside=True, clipzeros=True, extendpols=False, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
print ("Clipping high and zero values for unpolarized calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=unpolcalib, clipminmax=clipunpolcalib,
        datacolumn="DATA",clipoutside=True, clipzeros=True, extendpols=False, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
print ("Clipping high and zero values for target") 
flagdata(vis=ms,mode="clip", spw='',field=target, clipminmax=cliptarget,
        datacolumn="DATA",clipoutside=True, clipzeros=True, extendpols=False, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
        
        
#Check the flagdata percentage at this stage
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="DATA", name=ms+'summary2.split', action="calculate", overwrite=True, writeflags=True)   

'''
######################################################################################################################################################################################################
#
#  Flagging (Mode: tfcrop and extend)                (Parameters modified by Salmoli Ghosh)
#
######################################################################################################################################################################################################
# After clip, now flag using 'tfcrop' option for calibrator tight flagging
default(flagdata)
print ("Running automated flagging tfcrop to remove RFI from the calibrators") 
flagdata(vis=ms,mode="tfcrop", datacolumn="DATA", field=allcals, ntime="scan",
        timecutoff=5.0, freqcutoff=5.0, timefit="poly",freqfit="poly",flagdimension="freqtime", 
        extendflags=False,growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)
# Now extend the flags (90% more means full flag, change if required)
default(flagdata)
print ("Extending the flags for all correlations for calibrators") 
flagdata(vis=ms,mode="extend",spw='',field=allcals,datacolumn="DATA",
         ntime="scan", extendpols=True,growtime=90.0, growfreq=90.0,growaround=False,
         flagneartime=False, flagnearfreq=False, action="apply", flagbackup=True,overwrite=True, writeflags=True)


# Now flag for target - moderate flagging, more flagging in self-cal cycles
# After clip, now flag using 'tfcrop' option for target
default(flagdata)
print ("Running automated flagging tfcrop to remove RFI from the target") 
flagdata(vis=ms,mode="tfcrop", datacolumn="DATA", field=target, ntime="scan",
        timecutoff=6.0, freqcutoff=6.0, timefit="poly",freqfit="poly",flagdimension="freqtime", 
        extendflags=False, growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)
# Now extend the flags (90% more means full flag, change if required)
default(flagdata)
print ("Extending the flags for all correlations for target") 
flagdata(vis=ms,mode="extend",spw='',field=target,datacolumn="DATA", ntime="scan", extendpols=True,growtime=90.0, growfreq=90.0,growaround=False,
         flagneartime=False, flagnearfreq=False, action="apply", flagbackup=True,overwrite=True, writeflags=True)
# Now summary
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="DATA", name=ms+'summary3.split', action="calculate", overwrite=True, writeflags=True)   


######################################################################################################################################################################################################
#
#  Initial Calibration Round                (Parameters modified by Salmoli Ghosh)
#
######################################################################################################################################################################################################
print ("First round of calibration")

print ("starting initial flux density scaling")
set1=setjy(vis=ms, field = fluxfield, spw = '', scalebychan=True, usescratch=True)
print('setjy output',set1)

# Phase only calibration added - suggested and tested by Silpa Sasikumar
#
default(gaincal)
print (" starting initial phase only gaincal -> %s" % p0)
gaincal(vis=ms,caltable=p0,field=gaincals, spw=gainspw,
         solint="30s",combine="", minsnr=3.0, refant=ref_ant, solnorm=False,
         gaintype="G",calmode="p",append=False, parang=True)

default(gaincal)
print (" starting delay calibration -> %s" % kcorrfile0)
gaincal(vis=ms, caltable = kcorrfile0, field = kcorrfield, spw = '', 
        refant = ref_ant, solnorm = True,  gaintype = 'K', 
        gaintable =[p0], gainfield=gaincals, solint = '1.0min', combine = 'scan', minsnr=3.0,
        parang = True, append = False)

default(bandpass)  
print (" starting bandpass -> %s" % bpassfile0)
bandpass(vis=ms, caltable = bpassfile0, field = bpassfield, spw = '', minsnr=3.0,
         refant = ref_ant, solnorm = True,  solint = 'inf', 
         bandtype = 'B', fillgaps = 0, gaintable = [p0, kcorrfile0], gainfield=[gaincals,kcorrfield], 
         parang = True, append = False)
         
default(gaincal)         
print (" starting gaincal -> %s" % 0)
gaincal(vis=ms, caltable = 0, field = gaincals, spw = gainspw, 
        refant = ref_ant, solint = '1.0min', solnorm = False,  
        gaintype = 'G', combine = '', calmode = 'ap', minsnr=3.0, uvrange=uvracal,
        gaintable = [kcorrfile0,bpassfile0], gainfield = [kcorrfield,bpassfield],
        append = False, parang = True)
        
        
default(fluxscale)
print (" starting fluxscale -> %s" % fluxfile0) 
fluxsc=fluxscale(vis=ms, caltable = 0, reference = [fluxfield], 
          transfer = [transferfield], fluxtable = fluxfile0, 
          listfile = ms+'.fluxscale.txt0',
          append = False)               
print(fluxsc)  


######################################################################################################################################################################################################
#
#  Applying the Initial Calibration                (Parameters modified by Salmoli Ghosh)
#
######################################################################################################################################################################################################

print (" Applying Calibrations:")

#for individual flux calibrators (added by Salmoli Ghosh)
for i in range(0,len(fluxfieldind)):
	print (" applying calibrations: primary calibrator",fieldnames[int(fluxfieldind[i])])
	default(applycal)
	applycal(vis=ms, field = fluxfieldind[i], spw = '', selectdata=False, calwt = False,
    	gaintable = [kcorrfile0,bpassfile0, fluxfile0], interp=['','nearest','linear'],
   	gainfield = [kcorrfield,bpassfield,fluxfieldind[i]], 
    	parang = True)

print (" applying calibrations: phase calibrator",fieldnames[int(phasefield)])
default(applycal)
applycal(vis=ms, field = phasefield, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile0, bpassfile0, fluxfile0], interp=['','nearest','linear'],
    gainfield = [kcorrfield, bpassfield,phasefield],
    parang= True)

#print " applying calibrations: polarized calibrator 2"    #if the polarized calibrator is already a flux calibrator this step is not needed, if not untag this function
#applycal(vis=ms, field = polcalib, spw = flagspw, selectdata = False, calwt = False,
#    gaintable = [kcorrfile0, bpassfile0, fluxfile0],
#    gainfield = [kcorrfield, bpassfield,polcalib2],
#    parang= True)

default(applycal)
print (" applying calibrations: unpolarized calibrator",fieldnames[int(unpolcalib)])
applycal(vis=ms, field = unpolcalib, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile0, bpassfile0, fluxfile0], interp=['','nearest','linear'],
    gainfield = [kcorrfield, bpassfield,unpolcalib],
    parang= True)

#print " applying calibrations: another field"        ---------Janhavi Baghel
#applycal(vis=ms, field = anofield, spw = flagspw, selectdata = False, calwt = False,
#    gaintable = [kcorrfile0, bpassfile0, fluxfile0],
#    gainfield = [kcorrfield, bpassfield, anofield],
#    parang= True)

default(applycal)
print (" applying calibrations: target fields")
applycal(vis=ms, field = target, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile0, bpassfile0, fluxfile0], interp=['','nearest','linear'],
    gainfield = [kcorrfield, bpassfield,phasefield],
    parang= True)




######################################################################################################################################################################################################
#
#  Second Round of Flagging 
#  Flagging (Mode: clip)       
#
######################################################################################################################################################################################################

print ("Second Round of Flagging") 


default(flagdata)
#Flag using 'clip' option to remove high points for calibrators
print ("Clipping high and zero values for flux calibrator" )
flagdata(vis=ms,mode="clip", spw='',field=fluxfieldind[0], clipminmax=clipfluxcal2,
        datacolumn="corrected",clipoutside=True, clipzeros=True, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)

default(flagdata)
print ("Clipping high and zero values for phase calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=phasefield, clipminmax=clipphasecal2,
        datacolumn="corrected",clipoutside=True, clipzeros=True,
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
        
#print "Flagging Step 3/12"                 ---------Janhavi Baghel
#flagdata(vis=ms,mode="clip", spw=flagspw,field=anofield, clipminmax=clipanofield,
#        datacolumn="corrected",clipoutside=True, clipzeros=True, extendpols=False, 
#        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)

default(flagdata)
print ("Clipping high and zero values for polarized calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=polcalib, clipminmax=clippolcalib2,
        datacolumn="corrected",clipoutside=True, clipzeros=True,
       action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)

default(flagdata)
print ("Clipping high and zero values for unpolarized calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=unpolcalib, clipminmax=clipunpolcalib2,
        datacolumn="corrected",clipoutside=True, clipzeros=True, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
        
default(flagdata)
print ("Clipping high and zero values for target") 
flagdata(vis=ms,mode="clip", spw='',field=target, clipminmax=cliptarget,
        datacolumn="corrected",clipoutside=True, clipzeros=True,
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)        
# Now summary
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="corrected", name=ms+'summary4.split', action="calculate", overwrite=True, writeflags=True) 

######################################################################################################################################################################################################
#
#   Flagging the calibrators (Mode: tfcrop & rflag)
#
######################################################################################################################################################################################################


# After clip, now flag using 'tfcrop' option for calibrator tight flagging
default(flagdata)
print ("Running automated flagging tfcrop to remove RFI from the calibrators") 
flagdata(vis=ms,mode="tfcrop", datacolumn="residual", field=gaincals, ntime="scan",
        timecutoff=4.0, freqcutoff=4.0, timefit="line",freqfit="line",flagdimension="freqtime", 
        extendflags=False,growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)

# Now flag using 'rflag' option for calibrator tight flagging
default(flagdata)
print ("Running automated flagging rflag to remove RFI from the calibrators") 
flagdata(vis=ms,mode="rflag",datacolumn="residual",field=gaincals, extendflags=False,
        timedevscale=5.0,freqdevscale=5.0,spectralmax=1000.0, growaround=False,
        flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)


# Now extend the flags (90% more means full flag, change if required)
print ("Extending the flags for all correlations for calibrators") 
flagdata(vis=ms,mode="extend",spw='',field=gaincals,datacolumn="residual",
         ntime="scan", extendpols=True, growtime=90.0, growfreq=90.0,growaround=False,
         flagneartime=False, flagnearfreq=False, action="apply", flagbackup=True,overwrite=True, writeflags=True)

# Now summary
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="residual", name=ms+'summary5.split', action="calculate", overwrite=True, writeflags=True) 

######################################################################################################################################################################################################

#   Flagging the target (Mode: tfcrop & rflag)

######################################################################################################################################################################################################   

#Baseline-dependent rflag for target           (added by Salmoli Ghosh, idea taken from CAPTURE pipeline) 
#Less flagging for shorter baselines (bad RFI condition), normal flagging for longer baselines (normal RFI condition). Otherwise shorter baselines will be completely flagged
msmd.open(ms)
antnames=[]
for i in range(0,30,1):                   #GMRT has 30 antennas
	antnames.append(msmd.antennanames(i)[0])
msmd.done()
baselines=[]
for j in range(0,len(antnames)):
	for k in range(0,len(antnames)):
		if k>j:
			baselines.append(antnames[j]+'&'+antnames[k])
ccant=[]
carmant=[]
for i in range(0,len(baselines)):
	if baselines[i].count('C')==2:
		ccant.append(baselines[i])
	else:
		carmant.append(baselines[i])

shortbls=[] 
shortbls.append(str('; '.join(ccant)))
longbls=[]
longbls.append(str('; '.join(carmant)))

# After clip, now flag using 'tfcrop' option for target
default(flagdata)
print ("Running automated flagging tfcrop to remove RFI from the target") 
flagdata(vis=ms,mode="tfcrop", datacolumn="corrected", field=target, ntime="scan", antenna=shortbls[0],
        timecutoff=8.0, freqcutoff=8.0, timefit="poly",freqfit="poly",flagdimension="freqtime", 
        extendflags=False,growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)
flagdata(vis=ms,mode="tfcrop", datacolumn="corrected", field=target, ntime="scan", antenna=longbls[0],
        timecutoff=6.0, freqcutoff=5.0, timefit="poly",freqfit="poly",flagdimension="freqtime", 
        extendflags=False,growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)

# Now flag using 'rflag' option for target
default(flagdata)
print ("Running automated flagging rflag to remove RFI from the target") 
flagdata(vis=ms,mode="rflag",datacolumn="corrected",field=target, antenna=shortbls[0], extendflags=False,
        timedevscale=8.0,freqdevscale=5.0,spectralmax=1000.0, growaround=False,
        flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
flagdata(vis=ms,mode="rflag",datacolumn="corrected",field=target, antenna=longbls[0], extendflags=False,
        timedevscale=5.0,freqdevscale=5.0,spectralmax=1000.0, growaround=False,
        flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
# Now summary
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="corrected", name=ms+'summary6.split', action="calculate", overwrite=True, writeflags=True) 
#

######################################################################################################################################################################################################
#
print (" Deleting existing model column")
clearcal(ms)

######################################################################################################################################################################################################
#
######################################################################################################################################################################################################
#
#  Second Calibration Round                (Parameters modified by Salmoli Ghosh)
#
######################################################################################################################################################################################################
print ("Calibrating measurement set %s" % ms)

print (" starting initial flux density scaling")
set2=setjy(vis=ms, field = fluxfield, spw = '', scalebychan=True)
print(set2)

default(gaincal)
print (" starting initial phase only gaincal -> %s" % gainfilep)
gaincal(vis=ms,caltable=gainfilep,field=gaincals,spw=gainspw,intent="",solint="30s",
         refant=ref_ant, refantmode='strict', minsnr=2.0,solnorm=False,
         gaintype="G",calmode="p",append=False,gaintable=[],parang=True)

default(gaincal)
print (" starting initial gaincal -> %s" % kcorrfile)
gaincal(vis=ms, caltable = kcorrfile, field = kcorrfield, spw = '', 
        refant = ref_ant, refantmode='strict', solnorm = True,  gaintype = 'K', 
        gaintable =[gainfilep], gainfield=gaincals, solint = '10min', minsnr=2.0,
        parang = True, append = False)
  
default(bandpass)
print (" starting bandpass -> %s" % bpassfile)
bandpass(vis=ms, caltable = bpassfile, field = bpassfield, spw = '', minsnr=3.0,
         refant = ref_ant,solnorm = True,  solint = 'inf', combine='scan',
         bandtype = 'B', fillgaps = 0, gaintable = [gainfilep, kcorrfile], gainfield=[gaincals, kcorrfield], 
         parang = True, append = False)
     
default(gaincal)
print (" starting gaincal -> %s" % gainfile)
gaincal(vis=ms, caltable = gainfile, field = gaincals, spw = gainspw, 
        refant = ref_ant, refantmode='strict', solint = '1.0min', solnorm = False,  
        gaintype = 'G', calmode = 'ap', minsnr=2.0, uvrange=uvracal,
        gaintable = [kcorrfile,bpassfile], gainfield = [kcorrfield,bpassfield],
        append = False, parang = True)

default(fluxscale)
print (" starting fluxscale -> %s" % fluxfile)
fluxsc2=fluxscale(vis=ms, caltable = gainfile, reference = [fluxfield], 
          transfer = [transferfield], fluxtable = fluxfile, 
          listfile = ms+'.fluxscale.txt2',
          append = False) 
print(fluxsc2) 
######################################################################################################################################################################################################
#
#  Polarization Calibration   
#          
######################################################################################################################################################################################################
fluxdensity_UNPOL=0.311265592978078          #from fluxscale task (added by Salmoli ghosh)
reffreq_UNPOL='6.62402344e+08Hz'
spix_UNPOL=[0., 0., 0.]
######################################################################################################################################################################################################
print (" starting Polarization calibration -> ")
#
print (" setting the polarization calibrator models")
set3=setjy(vis=ms, field = polcalib, spw = '', scalebychan=True, standard='manual', 
      fluxdensity=[I,0,0,0], spix=alpha, reffreq=reffreq, polindex=polfrac, polangle=polangle, usescratch=True)
set4=setjy(vis=ms, field = unpolcalib, spw = '', scalebychan=True, standard='manual', 
      fluxdensity=[fluxdensity_UNPOL,0,0,0], spix=spix_UNPOL, reffreq=reffreq_UNPOL, polindex=[], polangle=[], usescratch=True)
  
print('setjy output for polarized calibrator', set3)
print('setjy output for unpolarized calibrator', set4)

print (" starting cross-hand delay calibration -> %s" % kcross1)
gaincal(vis=ms, caltable = kcross1, field = polcalib, spw = '', 
        refant = ref_ant, solint = 'inf', gaintype = 'KCROSS', combine = 'scan', calmode='ap',append=False,
        gaintable = [kcorrfile, bpassfile, gainfile], gainfield = [kcorrfield,bpassfield,polcalib],
        parang = True) 
#print " starting cross-hand delay calibration -> %s" % kcross2
#gaincal(vis=ms, caltable = kcross2, field = polcalib2, spw = flagspw, 
#        refant = refant, solint = 'inf', gaintype = 'KCROSS', combine = 'scan',
#        gaintable = [kcorrfile, bpassfile, gainfile], gainfield = [kcorrfield,bpassfield,polcalib2],
#        parang = True) 

######################################################################################################################################################################################################
#Analyse the various tables and choose the correct ones to apply

kcross = kcross1 #or kcross2
kcrosscalib = polcalib #or polcalib2
######################################################################################################################################################################################################

print (" starting leakage calibration -> %s" % leakage1)
polcal(vis=ms, caltable = leakage1, field = polcalib, spw = '', 
       refant = ref_ant, solint = 'inf', poltype = 'Df+QU', combine = 'scan',
       gaintable = [kcorrfile, bpassfile, gainfile, kcross], gainfield = [kcorrfield, bpassfield, polcalib, kcrosscalib],append=False)
#print " starting leakage calibration -> %s" % leakage2
#polcal(vis=ms, caltable = leakage2, field = polcalib2, spw = flagspw, 
#       refant = refant, solint = 'inf', poltype = 'Df+QU', combine = 'scan',
#       gaintable = [kcorrfile, bpassfile, gainfile, kcross], gainfield = [kcorrfield, bpassfield, polcalib2, kcrosscalib])

print (" starting leakage calibration :unpolarized -> %s" % unpolleakage1)
polcal(vis=ms, caltable = unpolleakage1, field = unpolcalib, spw = '', refant = ref_ant, solint = 'inf', poltype = 'Df', combine = 'scan', gaintable = [kcorrfile, bpassfile, gainfile, kcross], gainfield = [kcorrfield,bpassfield,unpolcalib,kcrosscalib])

print (" starting leakage calibration :unpolarized -> %s" % unpolleakage2)
polcal(vis=ms, caltable = unpolleakage2, field = unpolcalib2, spw = flagspw, refant = refant, solint = 'inf', poltype = 'Df', combine = 'scan', gaintable = [kcorrfile, bpassfile, gainfile, kcross], gainfield = [kcorrfield,bpassfield,unpolcalib2,kcrosscalib])

######################################################################################################################################################################################################
#
#  Open the leakage plots to check if there is any abnormality
#  Flag the high data point (say leakages above 45%)  (added by Salmoli Ghosh)  
#          
######################################################################################################################################################################################################
default(flagdata)
flagdata(vis=leakage1, mode='clip', correlation='ABS_ALL', clipminmax=[0.0, 0.25], datacolumn='CPARAM', clipoutside=True, action='apply', flagbackup=True, savepars=False)

flagdata(vis=unpolleakage1, mode='clip', correlation='ABS_ALL', clipminmax=[0.0, 0.45], datacolumn='CPARAM', clipoutside=True, action='apply', flagbackup=False, savepars=False)
flagdata(vis=unpolleakage2, mode='clip', correlation='ABS_ALL', clipminmax=[0.0, 0.45], datacolumn='CPARAM', clipoutside=True, action='apply', flagbackup=False, savepars=False)

######################################################################################################################################################################################################

leakage = leakage1#leakage1 or leakage2 or unpolleakage1
leakagecalib = polcalib#polcalib1 or polcalib2 or unpolcalib1

######################################################################################################################################################################################################


print (" starting polarization angle calibration -> %s" % polang1)
polcal(vis=ms, caltable = polang1, spw='', field = polcalib, refant = ref_ant, solint = 'inf', poltype = 'Xf',combine = 'scan', gaintable = [kcorrfile, bpassfile, gainfile, kcross, leakage], 
        gainfield = [kcorrfield,bpassfield,polcalib,kcrosscalib,leakagecalib],append=False)
#print " starting polarization angle calibration -> %s" % polang2
#polcal(vis=ms, caltable = polang2, field = polcalib2, refant = refant, solint = 'inf', poltype = 'Xf',combine = 'scan', gaintable = [kcorrfile, bpassfile, gainfile, kcross, leakage], 
#        gainfield = [kcorrfield,bpassfield,polcalib2,kcrosscalib,leakagecalib])

######################################################################################################################################################################################################
polang = polang1 #or polang2
polangcalib = polcalib #or polcalib2
######################################################################################################################################################################################################
#
#Applying the final calibration
#
######################################################################################################################################################################################################
#
print  (" Applying Calibrations:" )
for i in range(0,len(fluxfieldind)):
	print (" applying calibrations: primary calibrator",fieldnames[int(fluxfieldind[i])])
	applycal(vis=ms, field = fluxfieldind[i], spw = '', selectdata=False, calwt = False,
    	gaintable = [kcorrfile,bpassfile, fluxfile, kcross, leakage, polang], interp=['','nearest','nearest','','',''],
    	gainfield = [kcorrfield,bpassfield,fluxfieldind[i], kcrosscalib, leakagecalib, polangcalib],
    	parang = True)

print (" applying calibrations: phase calibrators")
applycal(vis=ms, field = phasefield, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile, bpassfile, fluxfile, kcross, leakage, polang], interp=['','nearest','nearest','','',''],
    gainfield = [kcorrfield, bpassfield,phasefield, kcrosscalib, leakagecalib, polangcalib],
    parang= True)

#print " applying calibrations: polarized calibrator"
#applycal(vis=ms, field = polcalib2, spw = flagspw, selectdata = False, calwt = False,
#    gaintable = [kcorrfile, bpassfile, fluxfile, kcross, leakage, polang],
#    gainfield = [kcorrfield, bpassfield,polcalib2, kcrosscalib, leakagecalib, polangcalib],
#    parang= True)

print (" applying calibrations: unpolarized calibrator")
applycal(vis=ms, field = unpolcalib, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile, bpassfile, fluxfile, kcross, leakage, polang], interp=['','nearest','nearest','','',''],
    gainfield = [kcorrfield, bpassfield,unpolcalib, kcrosscalib, leakagecalib, polangcalib],
    parang= True)

print (" applying calibrations: target fields")
applycal(vis=ms, field = target, spw = '', selectdata = False, calwt = False,
    gaintable = [kcorrfile, bpassfile, fluxfile, kcross, leakage, polang], interp=['','nearest','linear','','',''],
    gainfield = [kcorrfield, bpassfield,phasefield, kcrosscalib, leakagecalib, polangcalib],
    parang= True)
    
#print " applying calibrations: Another field"
#applycal(vis=ms, field = anofield, spw = flagspw, selectdata = False, calwt = False,
#    gaintable = [kcorrfile, bpassfile, fluxfile, kcross, leakage, polang],
#    gainfield = [kcorrfield, bpassfield,anofield, kcrosscalib, leakagecalib, polangcalib],
#    parang= True)


######################################################################################################################################################################################################
#
#  Some more flagging on the target; baseline dependent flagging for mode=rflag
#
######################################################################################################################################################################################################
print ("Flagging Target")
# Now flag for target - with tfcrop
default(flagdata)
print ("Running automated flagging tfcrop to remove RFI from the target") 
flagdata(vis=ms,mode="tfcrop",field=target, datacolumn="corrected", ntime="300s",
        timecutoff=8.0, freqcutoff=8.0, timefit="line",freqfit="line",flagdimension="freqtime", 
        extendflags=False,usewindowstats='sum',growaround=False,
        action="apply", flagbackup=True,overwrite=True, writeflags=True)
# Now flag using 'rflag' option for target
default(flagdata)
print ("Running automated flagging rflag to remove RFI from the target") 
flagdata(vis=ms,mode="rflag",field=target,datacolumn="corrected",antenna=shortbls[0], ntime='scan',extendflags=False,
        timedevscale=6.0,freqdevscale=6.0,spectralmax=1e6, growaround=False,
        flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
flagdata(vis=ms,mode="rflag",field=target,datacolumn="corrected",antenna=longbls[0], ntime='scan', extendflags=False,
        timedevscale=5.0,freqdevscale=5.0,spectralmax=1e6, growaround=False,
        flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
# Now summary
default(flagdata)
print("Flag summary till now")
flagdata(vis=ms,mode='summary',datacolumn="corrected", name=ms+'summary7.split', action="calculate", overwrite=True, writeflags=True) 

######################################################################################################################################################################################################
#
# Splitting the sources
#
######################################################################################################################################################################################################
print ("Splitting target field")
split(vis=ms, outputvis = fieldnames[int(target)]+'.ms', datacolumn='corrected', 
          field = target, spw = '', keepflags=False, width = specave)

split(vis=ms, outputvis = fieldnames[int(polcalib)]+'.ms', datacolumn='corrected', 
          field = polcalib, spw = '', keepflags=False, width = specave)

split(vis=ms, outputvis = fieldnames[int(unpolcalib)]+'.ms', datacolumn='corrected', 
          field = unpolcalib, spw = '', keepflags=False, width = specave)
          
split(vis=ms, outputvis = fieldnames[int(phasefield)]+'.ms', datacolumn='corrected', 
          field = phasefield, spw = '', keepflags=False, width = specave)


######################################################################################################################################################################################################
#
# Imaging  (parameters modified by Salmoli Ghosh)
#
######################################################################################################################################################################################################



ms=fieldnames[int(target)]+'.ms'


	
print ("\n Cleaning up. Starting imaging...")
#start self-calibration cycles    
count=1
scmode='p'
#
     
print ("Prepaing dirty image")
default(tclean)
tclean(vis=ms,
       datacolumn='data',
       imagename=ms+'.'+scmode+str(count-1),
       specmode='mfs',
       nchan=-1,
       nterms=2,
       niter=int(startniter),
       gain=0.1,
       gridder='widefield', 
       wprojplanes=-1,
       pblimit=-1.0,
       deconvolver='mtmfs',
       interactive=True,
       imsize=imagesize, cell=cellsize,
       stokes='I',
       projection="SIN",
       threshold=str(startthreshold/count)+'mJy',
       weighting='briggs',robust=0.5, 
       usemask='user',
       savemodel='modelcolumn')

exportfits(imagename=ms+'.'+scmode+str(count-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'.fits')
os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-1)+".ms")
print ("Made : " +scmode+str(count-1))

'''
#for creating dirty Q,U images

print ("Prepaing dirty Q image") 
 
tclean(vis=ms, 
	datacolumn='data', 
	imagename=ms+'.'+scmode+str(count-1)+'_Q', 
	specmode='mfs', 
	nchan=-1, 
	nterms=2, 
	niter=int(startniter), 
	gain=0.1,
	gridder='widefield',  
	wprojplanes=-1, 
	pblimit=-1.0, 
	deconvolver='mtmfs', 
	interactive=False, 
	imsize=imagesize, cell=cellsize, 
	stokes='Q', 
	projection="SIN", 
	threshold=str(startthreshold/count)+'mJy', 
	sidelobethreshold=2.0, 
  	growiterations=75, 
	weighting='briggs',robust=0.5,  
	usemask='user',mask=ms+'.'+scmode+str(count-1)+'.mask', 
	savemodel=None) 
      
exportfits(imagename=ms+'.'+scmode+str(count-1)+'_Q.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'_Q.fits') 

print ("Prepaing dirty U image") 
tclean(vis=ms, 
	datacolumn='data', 
	imagename=ms+'.'+scmode+str(count-1)+'_U', 
	specmode='mfs', 
	nchan=-1, 
	nterms=2, 
	niter=int(startniter), 
	gain=0.1,
	gridder='widefield',  
	wprojplanes=-1, 
	pblimit=-1.0, 
	deconvolver='mtmfs', 
	interactive=False, 
	imsize=imagesize, cell=cellsize, 
	stokes='U', 
	projection="SIN",  
	threshold=str(startthreshold/count)+'mJy', 
	sidelobethreshold=2.0, 
	growiterations=75, 
	weighting='briggs',robust=0.5,  
	usemask='user',mask=ms+'.'+scmode+str(count-1)+'.mask', 
	savemodel=None)
      
exportfits(imagename=ms+'.'+scmode+str(count-1)+'_U.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'_U.fits') 

tclean(vis=ms, 
	datacolumn='data', 
	imagename=ms+'.'+scmode+str(count-1)+'_V', 
	specmode='mfs', 
	nchan=-1, 
	nterms=2, 
	niter=int(startniter), 
	gain=0.1, 
	smallscalebias=0.6,
	gridder='widefield',  
	wprojplanes=-1, 
	pblimit=-1.0, 
	deconvolver='mtmfs', 
	interactive=False, 
	imsize=imagesize, cell=cellsize, 
	stokes='V', 
	projection="SIN", 
	threshold=str(startthreshold/count)+'mJy', 
	sidelobethreshold=2.0, 
  	growiterations=75, 
	weighting='briggs',robust=0,  
	usemask='user',mask=ms+'.'+scmode+str(count-1)+'.mask', 
	savemodel=None) 
      
exportfits(imagename=ms+'.'+scmode+str(count-1)+'_V.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'_V.fits') 


print ("Made : Dirty Q, U image " +scmode+str(count-1))
'''
#start self-calibration cycles  
print ("Starting self-calibration, going to phase only calibration Cycle")
casalog.post("Staring self-calibration, going to phase only calibration Cycle")
count=2
for j in range(pcycles):  
  scmode='p'
  if(doflag==True and count>=2):
     print ("Began flagging :"+scmode+str(count-1))
     default(flagdata)
     flagdata(vis=ms,mode="clip", clipminmax=clipresid,
              datacolumn="RESIDUAL_DATA",clipoutside=True, clipzeros=True, extendpols=False, 
              action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="rflag",datacolumn="RESIDUAL_DATA",antenna=shortbls[0], ntime='scan',extendflags=False,
              timedevscale=6.0,freqdevscale=6.0,spectralmax=1e6, growaround=False,
              flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="rflag",datacolumn="RESIDUAL_DATA",antenna=longbls[0], ntime='scan', extendflags=False,
              timedevscale=5.0,freqdevscale=5.0,spectralmax=1e6, growaround=False,
              flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="summary",datacolumn="RESIDUAL_DATA", extendflags=False, 
              name=ms+'temp.summary'+scmode+str(count-1), action="apply", flagbackup=True,overwrite=True, writeflags=True)
#
  print ("Began doing self-cal on :"+scmode+str(count-1))
  default(gaincal)
  gaincal(vis=ms,caltable=ms+'.'+scmode+str(count-1),selectdata=False,solint=str(sol_int/2**(count-1))+'min',refant=ref_ant,refantmode="strict",
        minblperant=6, spw=gainspw2,minsnr=3.0,solnorm=True,gaintype="G",calmode=scmode,append=False, uvrange=uvrascal, parang=False)
# 
  print ("Began processing :"+scmode+str(count-1))
  default(applycal)
  applycal(vis=ms, selectdata=False,gaintable=ms+'.'+scmode+str(count-1), parang=False,calwt=False,applymode="calflag",flagbackup=True)  
#
  default(tclean)
  tclean(vis=ms,
       datacolumn='corrected',
       imagename=ms+'.'+scmode+str(count-1),
       specmode='mfs',
       nchan=-1,
       nterms=2,
       niter=int(startniter*2**count),
       cycleniter=int(startniter*2**(count-1)),
       gain=0.1,
       smallscalebias=0.6,
       gridder='widefield', 
       wprojplanes=-1,
       pblimit=-1.0,
       deconvolver='mtmfs',
       interactive=True,
       imsize=imagesize, cell=cellsize,
       stokes='I',
       projection="SIN",
       threshold=str(startthreshold/count)+'mJy',
       weighting='briggs',robust=0.5, 
       usemask='user',mask=ms+'.'+scmode+str(count-2)+'.mask',
       savemodel='modelcolumn')

  exportfits(imagename=ms+'.'+scmode+str(count-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'.fits')
  os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-1)+".ms")
  print ("Made : " +scmode+str(count-1))
  
  count = count + 1

#
print ("Completed phase only self-calibration, going to A&P calibration Cycle")
casalog.post("Completed phase only self-calibration, going to A&P calibration Cycle")
#
count=count
for j in range(apcycles):  
  scmode='ap'
  if j>=3:
   sfactor=4
  else:
   sfactor=count-pcycles-1 
#
  if(doflag==True):
     print ("Began flagging :"+scmode+str(count-pcycles-1))
     default(flagdata)
     flagdata(vis=ms,mode="clip", spw="",field='', clipminmax=clipresid,
              datacolumn="RESIDUAL_DATA",clipoutside=True, clipzeros=True, extendpols=False, 
              action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="rflag",datacolumn="RESIDUAL_DATA",antenna=shortbls[0], ntime='scan',extendflags=False,
              timedevscale=6.0,freqdevscale=6.0,spectralmax=1e6, growaround=False,
              flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="rflag",datacolumn="RESIDUAL_DATA",antenna=longbls[0], ntime='scan', extendflags=False,
              timedevscale=5.0,freqdevscale=5.0,spectralmax=1e6, growaround=False,
              flagneartime=False,flagnearfreq=False,action="apply",flagbackup=True,overwrite=True, writeflags=True)
     flagdata(vis=ms,mode="summary",datacolumn="RESIDUAL_DATA", extendflags=False, 
              name=ms+'temp.summary'+scmode+str(count-1), action="apply", flagbackup=True,overwrite=True, writeflags=True)
#
  print ("Began doing self-cal on :"+scmode+str(count-pcycles-1))
  gaincal(vis=ms,caltable=ms+'.'+scmode+str(count-pcycles-1),selectdata=False,solint=str(apsolint/2**sfactor)+'min',refant=ref_ant,
        refantmode="strict",spw=gainspw2,minblperant=6, minsnr=3.0,solnorm=True,gaintype="G",calmode=scmode,append=False, parang=False)
#
  applycal(vis=ms, selectdata=False,gaintable=ms+'.'+scmode+str(count-pcycles-1), parang=False,calwt=False,applymode="calflag",flagbackup=True) 

#
  print ("Began processing :"+scmode+str(count-pcycles-1))
  if count-pcycles-1==1:
  	maskfile=ms+'.p'+str(pcycles)+'.mask'
  else:
  	maskfile=ms+'.ap'+str(count-pcycles-2)+'.mask'
  tclean(vis=ms,
       datacolumn='corrected',
       imagename=ms+'.'+scmode+str(count-pcycles-1),
       specmode='mfs',
       nchan=-1,
       nterms=2,
       niter=int(startniter*2**count),
       cycleniter=int(startniter*2**(count-1)),
       gain=0.1,
       gridder='widefield', 
       wprojplanes=-1,
       pblimit=-1.0,
       deconvolver='mtmfs',
       interactive=True,
       imsize=imagesize, cell=cellsize,
       stokes='I',
       projection="SIN",
       threshold=str(startthreshold/count)+'mJy',
       weighting='briggs',robust=0.5, 
       usemask='user',mask=maskfile,
       savemodel='modelcolumn')
  
  exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'.fits')
  os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-pcycles-1)+".ms")
  print ("Made : " +scmode+str(count-pcycles-1))
  
  count = count + 1
#
print ("Completed processing AP self-calibrations\n")
casalog.post("Completed processing A&P self-calibrations")

######################################################################################################################################################################################################
# Creating the Q, U, V images from the self-calibrated ms file


count=count-1

print ("Creating Stokes Q image")
tclean(vis=ms,
	selectdata=True,
	datacolumn='corrected', 
	imagename=ms+'.'+scmode+str(count-pcycles-1)+'_Q',
	imsize=imagesize,
	cell=cellsize,
       stokes="Q",
       projection="SIN",
       specmode="mfs",
       reffreq="",
       nchan=-1,
       gridder="widefield",
       wprojplanes=wproj,
       pblimit=-1,
       gain=0.1,
       deconvolver="mtmfs",
       pbcor=False,
       weighting="briggs",robust=0.5,
       niter=int(startniter*2**count),
       threshold="0.005mJy",
       cycleniter=-1,
       interactive=False,usemask="user",mask=ms+'.'+scmode+str(count-pcycles-1)+'.mask',
       savemodel="none")

exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'_Q'+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'_Q'+'.fits' )

print ("Creating Stokes U image")

tclean(vis=ms,
	selectdata=True,
	datacolumn='corrected', 
	imagename=ms+'.'+scmode+str(count-pcycles-1)+'_U',
	imsize=imagesize,
	cell=cellsize,
       stokes="U",
       projection="SIN",
       specmode="mfs",
       reffreq="",
       nchan=-1,
       gridder="widefield",
       wprojplanes=wproj,
       pblimit=-1,
       deconvolver="mtmfs",
       pbcor=False,
       gain=0.1,
       weighting="briggs",robust=0.5,
       niter=int(startniter*2**count),
       threshold="0.005mJy",
       cycleniter=-1,
       interactive=False,usemask="user",mask=ms+'.'+scmode+str(count-pcycles-1)+'.mask',
       savemodel="none")

exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'_U'+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'_U'+'.fits' )


print ("Done")



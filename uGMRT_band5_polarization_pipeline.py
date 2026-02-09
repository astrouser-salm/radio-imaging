'''
# Pipeline originally developed by Russ Taylor in 2011
# Modified by Ishwara Chandra in 2018 
# This is an improved version of the polarization pipeline by Silpa Sasikumar
# Major improvements in flagging and self-calibration 
# Kindly refer to the paper: Ishwara-Chandra et al 2020  
# NOT tested for uGMRT band-2 (150 MHz) 
# Queries: jbaghel@ncra.tifr.res.in
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
# Please CHANGE channels and source fields as per your data.
# Also change clip parameters if you have much stronger calibrator and/or target
#
#
# The parameters below are typical for 550 MHz, 2048 channels at Band-4 (550-750 MHz)
# Please change as required for your data.
# In BAND-4, recommended channels corresponding to ~ 560 MHz to ~ 810 MHz. The sensitivity drops sharply after 810 MHz.
# In any case DO NOT use beyond 820 MHz.
# It is highly recommended not to use OQ208 (unpolarized calibrator) to calculate the instrumental leakage for uGMRT since it is a very faint source (~few Jy); therefore a single short scan does not provide sufficient SNR to accurately determine the instrumental polarization.
# Also, we do not recommend the use of 3C138 (polarized calibrator) for leakage calibration.
# We recommend 3C286 (polarized calibrator) or 3C84 (unpolarized calibrator) for leakage calibration.
'''
####################################################################################################################################
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

#you can decide the clipping fluxes on the basis of the amplitudes for each field
 

####################################################################################################################################


print ("Initializing parameters") #Change these according to your data


####################################################################################################################################

fluxfield      = '0,3'            # field number of the primary flux calibrator
#If any of your fluxfield is absent in the catalog, better not to add it as a fluxfield, CASA solves it incorrectly 
#always check fluxscale results
#verified by Salmoli

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
specave = 10                       # number of channels to average; suggested post-average BW (approx)
                                 
timeave = '0s'                    # time averaging
gainspw2 = '0:7~150'              # central good channels after split for self-cal

bpassfield     = '3'              # field number of the bandpass calibrator which has been observed atleast twice during the observations, add phasecal only if it is very strong 
#tested by Salmoli Ghosh

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
####################################################################################################################################
#Imaging parameters

pcycles=4              # number of phase-only self-calibration
apcycles=3             # number of amplitude and phase self-calibration
doflag=True            #
sol_int=16.0            # this transaltes to 8, 4, 2 and 1 min for each selfcal loop
apsolint=8.0          # this transaltes to  4, 2 and 1 min for each selfcal loop
startthreshold=0.1     # start threshold to stop flux (mJy, will reduce by count subsequently)
startniter=2500        # start iterations (will double in each phase-selfcal and 4 times in each A&P loop)
imagesize=[7200,7200]  # should cover alteast up to null at lower part of the band
cellsize='0.5arcsec'   # should be atleast 3 pixels in minor axis
wproj=-1               # w projection, default autocalculate
eachQUV=False	       # Create Stokes Q and U images for each self-cal iteration. ----- Janhavi Baghel
dirtyQUV=True           # Create Stokes Q and U images for dirty image. ----- Janhavi Baghel
createV=True		 # Create Stokes V image ----- Janhavi Baghel

####################################################################################################################################
#default(split)
#split(vis=ms,outputvis='multi_split.ms',datacolumn='DATA',spw=splitspw)
#ms='multi_split.ms'
####################################################################################################################################
#For polarization calibration; ----- Janhavi Baghel
#
##################################################################################################################
#Defining functions
##################################################################################################################

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
print (" starting initial phase only gaincal -> %s" % gainfilep0)
gaincal(vis=ms,caltable=gainfilep0,field=gaincals, spw=gainspw,
         solint="30s",combine="", minsnr=3.0, refant=ref_ant, solnorm=False,
         gaintype="G",calmode="p",append=False, parang=True)

default(gaincal)
print (" starting delay calibration -> %s" % kcorrfile0)
gaincal(vis=ms, caltable = kcorrfile0, field = kcorrfield, spw = '', 
        refant = ref_ant, solnorm = True,  gaintype = 'K', 
        gaintable =[gainfilep0], gainfield=gaincals, solint = '1.0min', combine = 'scan', minsnr=3.0,
        parang = True, append = False)

default(bandpass)  
print (" starting bandpass -> %s" % bpassfile0)
bandpass(vis=ms, caltable = bpassfile0, field = bpassfield, spw = '', minsnr=3.0,
         refant = ref_ant, solnorm = True,  solint = 'inf', 
         bandtype = 'B', fillgaps = 0, gaintable = [gainfilep0, kcorrfile0], gainfield=[gaincals,kcorrfield], 
         parang = True, append = False)
         
default(gaincal)         
print (" starting gaincal -> %s" % gainfile0)
gaincal(vis=ms, caltable = gainfile0, field = gaincals, spw = gainspw, 
        refant = ref_ant, solint = '1.0min', solnorm = False,  
        gaintype = 'G', combine = '', calmode = 'ap', minsnr=3.0, uvrange=uvracal,
        gaintable = [kcorrfile0,bpassfile0], gainfield = [kcorrfield,bpassfield],
        append = False, parang = True)
        
        
default(fluxscale)
print (" starting fluxscale -> %s" % fluxfile0) 
fluxsc=fluxscale(vis=ms, caltable = gainfile0, reference = [fluxfield], 
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
'''
default(flagdata)
print ("Clipping high and zero values for unpolarized calibrator") 
flagdata(vis=ms,mode="clip", spw='',field=unpolcalib, clipminmax=clipunpolcalib2,
        datacolumn="corrected",clipoutside=True, clipzeros=True, 
        action="apply",flagbackup=True, savepars=False, overwrite=True, writeflags=True)
'''       
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



######################################################################################################################################################################################################
#
print (" Starting Polarization Calibration")
#Polarization calibration steps are for solving equations in linear basis (X and Y)
#The following steps have been adapted from ALMA polarization calibration procedure (checked by Salmoli Ghosh, Silpa S.; J. Baghel tested the other method)
#Includes outputs for reference
######################################################################################################################################################################################################
kcross = kcross1
kcrosscalib = polcalib


from casarecipes.almapolhelpers import *
'''
applycal(vis=ms, field = phasefield, spw = '', selectdata = False, calwt = False,
gaintable = [kcorrfile, bpassfile, fluxfile], interp=['','nearest','nearest'], gainfield = [kcorrfield, bpassfield,phasefield],parang= True)
'''     
gcalpol=ms+'.gcalpol'
gaincal(vis=ms,caltable=gcalpol,field=polcalib,solint='int',smodel=[1,0,0,0], gaintype='G',
gaintable=[bpassfile], refant=ref_ant, refantmode='strict', interp='nearest',parang=True)
      
     
qu=qufromgain(gcalpol)
'''
Latitude =  19.1000701545
Found as many as 4 fields.
Can't discern an ALMA bandname from: none
Found as many as 1 spws.
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=0.0deg) Gx/Gy= 1.31574458345 Q= 0.0735198593522 U= -0.0659361906026 P= 0.0987560172868 X= -20.943650689
For field id =  3  there are  1 good spws.
Spw mean: Fld= 3 Q= 0.0735198593522 U= -0.0659361906026 (rms= 0.0 0.0 ) P= 0.0987560172868 X= -20.943650689
'''
'''
Latitude =  19.10007015450946
Found as many as 4 fields.
Can't discern an ALMA bandname from: none
Found as many as 1 spws.
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
/data/salmoli/Downloads/casa-6.5.4-9-pipeline-2023.1.0.125/lib/py/lib/python3.8/site-packages/casarecipes/almapolhelpers.py:168: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
  fit=pl.lstsq(A,pl.square(ratio))
/data/salmoli/Downloads/casa-6.5.4-9-pipeline-2023.1.0.125/lib/py/lib/python3.8/site-packages/casarecipes/almapolhelpers.py:174: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
  fit=pl.lstsq(A,pl.square(rsum))
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 1 Spw= 0 (B=none, PA offset=0.0deg) Gx/Gy= -2.5360252695541585 Q= 4.361403508875255 U= -2.5388642870464326 P= 5.046550548173394 X= -15.102290381137898
For field id =  1  there are  1 good spws.
Spw mean: Fld= 1 Q= 4.361403508875255 U= -2.5388642870464326 (rms= 0.0 0.0 ) P= 5.046550548173394 X= -15.102290381137898

'''
gaincal(vis=ms, caltable=kcross1, selectdata=True, gaintype='KCROSS', field=kcrosscalib, solint='inf',refant=ref_ant,  refantmode='strict', smodel=[1,0,1,0], gaintable=[bpassfile,gcalpol], interp=['nearest','linear'])
'''      
applycal(vis=ms, field=polcalib, calwt=True, gaintable=[bpassfile,gcalpol,kcross1], interp=['nearest','linear', 'nearest'])
'''      
gaincal(vis=ms,caltable=ms+'.XY0amb', field=polcalib, gaintype='XYf+QU', solint='inf', combine='scan', preavg=300, refant=ref_ant, refantmode='strict', smodel=[1,0,1,0], gaintable=[bpassfile,gcalpol,kcross1], interp=['nearest','linear','nearest'])
     
'''     
Spw = 0 (ich=1024/2048): 
 X-Y phase = -54.5277668048 deg.
 Fractional Poln: Q = -0.290090739727, U = 0.0144347893074; P = 0.290449659346, X = 88.5756682681deg.
 Net (over baselines) instrumental polarization: 0.114551083591
''' 
xy0amb=ms+'.XY0amb'

S=xyamb(xytab=xy0amb,qu=list(qu.values())[0],xyout=ms+'.XY0')
'''
Expected QU =  (4.361403508875255, -2.5388642870464326)
Spw = 0: Found QU = [-0.29009074  0.01443479]
   ...CONVERTING X-Y phase from 23.572392498095944 to -156.42760750190408 deg
Ambiguity resolved (spw mean): Q= 0.29009073972702026 U= -0.014434789307415485 (rms= 0.0 0.0 ) P= 0.2904496521218769 X= -1.424331780744425
Returning the following Stokes vector: [1.0, 0.29009073972702026, -0.014434789307415485, 0.0]
'''
xy0=ms+'.XY0'

'''
applycal(vis=ms, field=polcalib, calwt=[True,True,False,False], gaintable=[bpassfile,gcalpol,kcross1,xy0], interp=['nearest','linear','nearest','nearest'])
'''

gcalpol2=ms+'.gcalpol2'

gaincal(vis=ms, caltable=gcalpol2, field=polcalib, solint='int', refant=ref_ant, refantmode='strict', smodel=S, gaintable=[bpassfile],interp=['nearest'], parang=True)
      
qu1=qufromgain(gcalpol2)
'''
Latitude =  19.1000701545
Found as many as 4 fields.
Can't discern an ALMA bandname from: none
Found as many as 1 spws.
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=0.0deg) Gx/Gy= 1.2736415716 Q= -0.0585953563219 U= -0.0827969321728 P= 0.101433464693 X= -62.643513808
For field id =  3  there are  1 good spws.
Spw mean: Fld= 3 Q= -0.0585953563219 U= -0.0827969321728 (rms= 0.0 0.0 ) P= 0.101433464693 X= -62.643513808
Out[58]: {3: (-0.058595356321934175, -0.082796932172799736)}

Latitude =  19.10007015450946
Found as many as 4 fields.
Can't discern an ALMA bandname from: none
Found as many as 1 spws.
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
/data/salmoli/Downloads/casa-6.5.4-9-pipeline-2023.1.0.125/lib/py/lib/python3.8/site-packages/casarecipes/almapolhelpers.py:168: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
  fit=pl.lstsq(A,pl.square(ratio))
/data/salmoli/Downloads/casa-6.5.4-9-pipeline-2023.1.0.125/lib/py/lib/python3.8/site-packages/casarecipes/almapolhelpers.py:174: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
  fit=pl.lstsq(A,pl.square(rsum))
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 1 Spw= 0 (B=none, PA offset=0.0deg) Gx/Gy= -8.359516454790715 Q= 2.673681871991404 U= -1.4111989618780396 P= 3.0232527619471576 X= -13.912794750857126
For field id =  1  there are  1 good spws.
Spw mean: Fld= 1 Q= 2.673681871991404 U= -1.4111989618780396 (rms= 0.0 0.0 ) P= 3.0232527619471576 X= -13.912794750857126
'''
polcal(vis=ms, caltable=ms+'.Df0', field=polcalib, solint='inf',combine='scan', preavg=300, poltype='Dflls', refant='', smodel=S, gaintable=[bpassfile,gcalpol2,kcross1,xy0], gainfield=['', '', '', ''], interp=['nearest','linear','nearest','nearest']
      
df0=ms+'.Df0'
      
Dgen(dtab=df0,dout=ms+'.Df0gen')

df0gen=ms+'.Df0gen'

applycal(vis=ms, field=polcalib, calwt=[False,True,True,False,False,False], gaintable=[kcorrfile,bpassfile,gcalpol2,kcross1,xy0,df0gen], interp=['','nearest','linear','linear','nearest','nearest'], gainfield=[kcorrfield,bpassfield, '', '','', '',], parang=True, applymode='calonly')
'''
tclean(vis=ms,
      ...:       field='3',
      ...:       imagename='3C286.noDterm.Stokes.clean',
      ...:       cell=['0.5arcsec'],
      ...:       imsize=[1000,1000],
      ...:       stokes='IQUV',
      ...:       deconvolver='clarkstokes',
      ...:       interactive=True, 
      ...:       weighting='briggs',
      ...:       robust=0.5,
      ...:       niter=1000)
      
df0gen=ms+'.Df0gen'
      ...: applycal(vis=ms, 
      ...:          field='3', 
      ...:          calwt=[True,True,False,False,False],
      ...:          gaintable=[bpassfile,gcalpol2,kcross1,xy0,df0gen],
      ...:          interp=['nearest','linear', 'linear','nearest', 'nearest'],  
      ...:          gainfield=['', '','', '', ''],
      ...:          parang=True)
'''     
######################################################################################################################################################################################################
#
#Imaging polarized calibrator after polarization calibration
#
######################################################################################################################################################################################################


tclean(vis=ms,
	field=polcalib,
	imagename=fieldnames[int(polcalib)]+'.withDterm.Stokes.clean',
  	cell=['0.1arcsec'],
	imsize=[250,250],
	stokes='IQUV',
	deconvolver='hogbom',
	weighting='briggs',
	robust=0.5,
	interactive=True, niter=10000)
	
######################################################################################################################################################################################################
#
#Applying the final calibration to other sources and target
#
######################################################################################################################################################################################################
      
applycal(vis=ms, field=unpolcalib, calwt=[False,True,False,False,False,False], gaintable = [kcorrfile,bpassfile,gainfile, kcross1,xy0,df0gen], interp=['','nearest','linear','nearest','nearest','nearest'], gainfield=[kcorrfield,bpassfield,unpolcal,'', '', ''], parang=True,applymode='calonly')

applycal(vis=ms, field=phasefield,target, calwt=[False,True,False,False,False,False], gaintable = [kcorrfile,bpassfile,fluxfile, kcross1,xy0,df0gen], interp=['','nearest','linear','nearest','nearest','nearest'], gainfield=[kcorrfield,bpassfield,phasefield,'', '', ''], parang=True, applymode='calonly')

######################################################################################################################################################################################################
#
# Splitting the sources
#
######################################################################################################################################################################################################
      
split(vis=ms, outputvis = fieldnames[int(target)]+'.ms', datacolumn='corrected', field = target, spw = '', keepflags=False, width = specave)
      
split(vis=ms, outputvis = fieldnames[int(fluxfieldind[0])]+'.ms', datacolumn='corrected', field = unpolcalib, spw = '', keepflags=False, width = specave)

split(vis=ms, outputvis = fieldnames[int(polcalib)]+'.ms', datacolumn='corrected', field = polcalib, spw = '', keepflags=False, width = specave)

split(vis=ms, outputvis = fieldnames[int(phasefield)]+'.ms', datacolumn='corrected', field = phasefield, spw = '', keepflags=False, width = specave)

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
       deconvolver='mtmfs', scales=[0,5,15,30],
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
	deconvolver='mtmfs', scales=[0,5,15,30],
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
	deconvolver='mtmfs', scales=[0,5,15,30],
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
	deconvolver='mtmfs', scales=[0,5,15,30],
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
        minblperant=6, spw=gainspw2,minsnr=3.0,solnorm=False,gaintype="G",calmode=scmode,append=False, uvrange=uvrascal, parang=False)
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
       cycleniter=1000,
       gain=0.1,
       smallscalebias=0.6,
       gridder='widefield', 
       wprojplanes=-1,
       pblimit=-1.0,
       deconvolver='mtmfs', scales=[0,5,15,30],
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
       cycleniter=2000,
       gain=0.1,
       gridder='widefield', 
       wprojplanes=-1,
       pblimit=-1.0,
       deconvolver='mtmfs', scales=[0,5,15,30],
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
       deconvolver="mtmfs", scales=[0,5,15,30],
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
       deconvolver="mtmfs", scales=[0,5,15,30],
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


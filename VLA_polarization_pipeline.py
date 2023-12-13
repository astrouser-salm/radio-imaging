#This is a modified version of EVLA polarization pipeline by Alice Pasetto

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

data='19B-198.sb37296228.eb37511549.58797.518064895834'
ms=data+'.ms'


fluxcal='3C138' # calibrator used to determine solutions for stokeI, cross hand delay and RL pol angle
fluxfield='0'
polcalib='3C138'
polfield='0'
target='XXXX'
targetfield='2'
phasecal='XXXX'
phasefield='3'
unpolcal='3C84'
unpolfield='1'
rfreq='9.75278139576GHz'        #check from casa log
reffreq=9.75278139576
fluxdensity_UNPOL=31.17500435967818  #from the weblog, fluxboot task
spix_UNPOL=[-0.17523711119766397,-0.88798061264207295] #from the weblog, fluxboot task
inifreq=8   #initial frequency #check from listobs  (for giving frequency range for interpolation)
endfreq=12   #final frequency 
ref_ant='ea12' #consider an antenna close to the centre but not the central one to avoid shadowing


##################################################################################################################
#
#  FLAGGING
#
##################################################################################################################
# Use plotms to plot amp vs spw to for different fields, flag the spws giving too high values
# Check for RL, LR correlations, if the amplitudes are too high, flag them


##################################################################################################################
#Defining functions
##################################################################################################################
#Defining total intensity spectrum of the Flux calibrator
def S(f,Sf,alpha,beta):
        return Sf*(f/reffreq)**(alpha+beta*np.log10(f/reffreq))

#Defining fractional polarization spectrum of the Flux calibrator
def PF(f,a,b,c):
        return a+b*((f-reffreq)/reffreq)+c*((f-reffreq)/reffreq)**2

#Defining polarization angle of the Flux calibrator
def PA(f,a,b):
        return a+b*((f-reffreq)/reffreq)


#If you have uncalibrated data (not calibrated using VLA pipeline while retrieving data) comment out the following chunk
'''
##################################################################################################################
NRAO pipeline for basic calibration
##################################################################################################################

context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')

hifv_importdata(vis=[mySDM], createmms='automatic',\
                 asis='Receiver CalAtmosphere', ocorr_mode='co',\
                 nocopy=False, overwrite=False)
hifv_hanning(pipelinemode="automatic")
hifv_flagdata(hm_tbuff='1.5int', fracspw=0.01, intents='*POINTING*,*FOCUS*,\
              *ATMOSPHERE*,*SIDEBAND_RATIO*, *UNKNOWN*, *SYSTEM_CONFIGURATION*,\
              *UNSPECIFIED#UNSPECIFIED*', template=True, filetemplate='17B-294-629_AliceFlags.txt')
hifv_vlasetjy(pipelinemode="automatic")
hifv_priorcals(pipelinemode="automatic")
hifv_testBPdcals(pipelinemode="automatic")
hifv_checkflag(checkflagmode='bpd-vla')
hifv_semiFinalBPdcals(pipelinemode="automatic")
hifv_checkflag(checkflagmode='allcals-vla')
hifv_solint(pipelinemode="automatic")
hifv_fluxboot(pipelinemode="automatic")
hifv_finalcals(pipelinemode="automatic")
hifv_applycals(pipelinemode="automatic")
hifv_checkflag(checkflagmode='target-vla')
hifv_targetflag(intents='*TARGET*')
hifv_statwt(datacolumn='corrected')
'''

##################################################################################################################
#Setting the models 
##################################################################################################################


#Model for SETJY using the flux calibrator using existing observed values at different frequencies

print("At frequency (in GHz)=",reffreq)

data = np.loadtxt('/data/kharb/salmoli/Calibrator_data/'+fluxcal+'_2019.txt')
data_index_range=np.where((data[:,0]>inifreq-3.0) & (data[:,0]<endfreq+3.0)) 
data_to_use=data[data_index_range]
popt, pcov = curve_fit(S, data_to_use[:,0], data_to_use[:,1])
print('I', popt[0], ' Jy') #flux density 
print('alpha', popt[1])    #spectral index
print('beta', popt[2])
print( 'Covariance')
print(pcov)
    
# Stokes I flux density
I = popt[0] 
# Spectral Index and beta coeff
alpha = [popt[1], popt[2]]
    
##################################################################################################################
#Model for Polarization fraction: obtaining value at reference frequency

popt, pcov = curve_fit(PF, data_to_use[:,0], data_to_use[:,2])
print("Polfrac Polynomial: ", popt)
print("Covariance")
print(pcov)
    
# Polarization Fraction coeff
polfrac = [*popt]
    
##################################################################################################################
#Model for Polarization angle: obtaining value at reference frequency

popt, pcov = curve_fit(PA, data_to_use[:,0], data_to_use[:,3])
print("Polangle Polynomial: ", popt)
print("Covariance")
print(pcov)
    
# Polarization Angle
polangle = [*popt]
    

##################################################################################################################
#
#####   SETJY -- FLUX CALIBRATOR   
#
##################################################################################################################
print (" starting Polarization calibration -> ")
#
print (" setting the polarization calibrator models for polarized calibrator-> %s" % polcalib)
 
setjy(vis=ms,field=fluxfield,spw='',selectdata=False,timerange='',scan='',intent='', 
observation='',scalebychan=True,standard='manual',model='',modimage='',listmodels=False, 
fluxdensity=[I,0,0,0],spix=alpha,reffreq=rfreq,polindex=polfrac,polangle=polangle,
rotmeas=0,fluxdict={},useephemdir=False,interpolation='nearest',usescratch=True,ismms=False)
    
##################################################################################################################
#
#####   SETJY -- UNPOLARIZED CALIBRATOR   
#
##################################################################################################################

print (" setting the polarization calibrator models for unpolarized calibrator-> %s" % unpolcal) 
   
setjy(vis=ms,field=unpolfield,spw='',selectdata=False,timerange="",scan="",intent="",
observation="",scalebychan=True,standard="manual",model="",modimage="",listmodels=False,
fluxdensity=[fluxdensity_UNPOL , 0, 0, 0],spix=spix_UNPOL,reffreq=rfreq,polindex=[],polangle=[],
rotmeas=0,fluxdict={},useephemdir=False,interpolation="nearest",usescratch=True,ismms=False)



##################################################################################################################    
#
####    Defining gtable and gfield
#
##################################################################################################################


gtable=[ms+'.hifv_priorcals.s5'+'_2.gc.tbl', 
ms+'.hifv_priorcals.s5'+'_3.opac.tbl', 
ms+'.hifv_priorcals.s5'+'_4.rq.tbl', 
ms+'.hifv_priorcals.s5'+'_6.ants.tbl',
ms+'.hifv_finalcals.s14'+'_2.finaldelay.tbl', 
ms+'.hifv_finalcals.s14'+'_4.finalBPcal.tbl', 
ms+'.hifv_finalcals.s14'+'_5.averagephasegain.tbl', 
ms+'.hifv_finalcals.s14'+'_7.finalampgaincal.tbl', 
ms+'.hifv_finalcals.s14'+'_8.finalphasegaincal.tbl']
    
    
interpol=['', '', '', '', '', 'linear,linearflag', '', '', '']
gfield=['']*(len(gtable)-2)
gfieldpol = gfield+[polfield,polfield]
gfieldunpol = gfield+[unpolfield,unpolfield]
gfieldphase = gfield+[phasefield,phasefield]


##################################################################################################################
#
#####   Solving for CROSS-HAND DELAY (use flux calibrator)   
#
##################################################################################################################
#To solve using Multiband Delay use combine="scan,spw"
#To solve for Singleband Delay (solving for each spw), use Combine="scan
    
##################################################################################################################

print (" starting cross-hand delay calibration -> %s" % polcalib)
gaincal(vis=ms,
        caltable=ms+'.kcross',
        field=polfield,
        spw='',
        refant=ref_ant,
        gaintype="KCROSS",
        solint="inf",
        combine="scan",
        calmode="ap",
        append=False,
        gaintable=gtable,
        gainfield=gfieldpol,
        interp=interpol,
        parang=True)
   
gtable.append(ms+'.kcross')
interpol.append('')
gfieldpol.append('')  
gfieldunpol.append('') 
gfieldphase.append('')

##################################################################################################################
#
#    Solving for the LEAKAGE TERMS   (both polarized and unpolarized calibrator can be used)
#
# for using polarized calib with unknown polarization, more than 2 scans are required
##################################################################################################################
    
print ("starting leakage calibration -> %s" % polcalib)
polcal(vis=ms,
           caltable=ms+'.polleakage',
           field=polfield,
           spw='',
           refant=ref_ant,
           poltype='Df+QU',
           solint='inf',
           combine='scan',
           gaintable=gtable,
           gainfield=gfieldpol,
           interp=interpol,
           append=False)
           
print ("starting leakage calibration -> %s" % unpolcal)
polcal(vis=ms,
           caltable=ms+'.unpolleakage',
           field=unpolfield,
           spw='',
           refant=ref_ant,
           poltype='Df',
           solint='inf',
           combine='scan',
           gaintable=gtable,
           gainfield=gfieldunpol,
           interp=interpol,
           append=False)


#stop the script here
#use plotms to plot amp vs frequency for each antenna with colourising the different correlations
#flag if the leakage values are too high or go back and check the issue


#if you are using polarized calibrator for leakage calib
gtable.append(ms+'.polleakage')    
gfieldpol.append('')  
gfieldunpol.append('') 
gfieldphase.append('')   

#if you are using unpolarized calibrator for leakage calib 
gtable.append(ms+'.unpolleakage')    
gfieldpol.append(unpolfield)  
gfieldunpol.append(unpolfield) 
gfieldphase.append(unpolfield)
    
##################################################################################################################
#
# Solving for the R-L POLARIZATION ANGLE (using polarization calibrator with known EVPA)
#
##################################################################################################################
   
print (" starting polarization angle calibration -> %s" % polcalib)
polcal(vis=ms,
           caltable=ms+'.polang',
           spw='',
           field=polfield,
           solint='inf',
           combine='scan',
           poltype='Xf',
           refant = ref_ant,
           gaintable=gtable,
           gainfield=gfieldpol,
           interp=interpol,
           append=False)
    
gtable.append(ms+'.polang')
interpol.append('')
gfieldpol.append('')  
gfieldunpol.append('') 
gfieldphase.append('')     


##################################################################################################################
#  END POLARIZATION CALIBRATION
##################################################################################################################
#
# APPLYING THE CALIBRATION
#
##################################################################################################################

print  (" Applying Calibrations:" )

print (" applying calibrations: primary calibrator")
applycal(vis = ms,
             field=fluxfield,
             gainfield=gfieldpol, 
             flagbackup=True,
             interp=interpol,
             gaintable=gtable,
             spw='', 
             calwt=False, 
             applymode='calflagstrict', 
             parang=True)  
             
print (" applying calibrations: unpolarized calibrator")
applycal(vis = ms,
             field=unpolfield,
             gainfield=gfieldunpol, 
             flagbackup=True,
             interp=interpol,
             gaintable=gtable,
             spw='', 
             calwt=False, 
             applymode='calflagstrict', 
             parang=True) 
             
print (" applying calibrations: phase calibrator")
applycal(vis = ms,
             field=phasefield,
             gainfield=gfieldphase, 
             flagbackup=True,
             interp=interpol,
             gaintable=gtable,
             spw='', 
             calwt=False, 
             applymode='calflagstrict', 
             parang=True) 
          
print (" applying calibrations: target fields")
applycal(vis = ms,
             field=targetfield,
             gainfield=gfieldphase, 
             flagbackup=True,
             interp=interpol,
             gaintable=gtable,
             spw='', 
             calwt=False, 
             applymode='calflagstrict', 
             parang=True)  
             
             
             
             
##################################################################################################################
#
# SPLITTING THE .MS FILE INTO DIFFERENT .MS FILES FOR DIFFERENT SOURCES
#
##################################################################################################################
             
specave=8 #VLA has 64 channels for each spw; remember about bandwidth smearing!             
split(vis=ms, outputvis = fluxcal+'.ms', datacolumn='corrected', 
          field = fluxfield, spw = '', keepflags=False, width = specave)  
split(vis=ms, outputvis = target+'.ms', datacolumn='corrected', 
          field = targetfield, spw = '', keepflags=False, width = specave)
split(vis=ms, outputvis = unpolcal+'.ms', datacolumn='corrected', 
          field = unpolfield, spw = '', keepflags=False, width = specave)   


####################################################################################################################################
#
# IMAGING
#
####################################################################################################################################
pcycles=3              # number of phase-only self-calibration
apcycles=2             # number of amplitude and phase self-calibration            
solint=8.0            # this transaltes to 8, 4, 2 min for each selfcal loop
apsolint=4.0          # this transaltes to  8, 4 min for each selfcal loop
startthreshold=0.1     # start threshold to stop flux (mJy, will reduce by count subsequently)
startniter=2500        # start iterations (will increase in further cycles)
imagesize=[4096,4096]  # should cover alteast up to null at lower part of the band
cellsize='1arcsec'     # should be atleast 3 pixels in minor axis
wproj=-1               # w projection, default autocalculate
uvrascal='>0.75klambda'
####################################################################################################################################

ms=target+'.ms'

	
print ("\n Cleaning up. Starting imaging...")
#start self-calibration cycles    
count=1
scmode='p'


print ("Preparing dirty image")

tclean(vis=ms,
       datacolumn='data',
       imagename=ms+'.'+scmode+str(count-1),
       specmode='mfs',
       nchan=-1,
       nterms=2,
       niter=int(startniter*2**count),
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
       sidelobethreshold=2.0,
       growiterations=75,
       weighting='briggs',robust=0, 
       usemask='user',
       savemodel='modelcolumn')


exportfits(imagename=ms+'.'+scmode+str(count-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'.fits')
os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-1)+".ms")
print ("Made : Dirty image " +scmode+str(count-1))

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
	cyclefactor=1.3, 
	threshold=str(startthreshold/count)+'mJy', 
	sidelobethreshold=2.0, 
  	growiterations=75, 
	weighting='briggs',robust=0,  
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
	cyclefactor=1.3, 
	threshold=str(startthreshold/count)+'mJy', 
	sidelobethreshold=2.0, 
	growiterations=75, 
	weighting='briggs',robust=0,  
	usemask='user',mask=ms+'.'+scmode+str(count-1)+'.mask', 
	savemodel=None)
      
exportfits(imagename=ms+'.'+scmode+str(count-1)+'_U.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'_U.fits') 

print ("Made : Dirty Q, U image " +scmode+str(count-1))
'''

#start self-calibration cycles  
print ("Starting self-calibration, going to phase only calibration Cycle")

count=2
for j in range(pcycles):
  scmode='p'
  
  print ("Began doing self-cal on :"+scmode+str(count-1))
  gaincal(vis=ms,
  	caltable=ms+'.'+scmode+str(count-1),
  	selectdata=False,
  	solint=str(solint/2**(count-1))+'min',
  	refant=ref_ant,
  	refantmode="strict",
        minblperant=6, 
        minsnr=3.0,
        gaintype="G",
        calmode=scmode,
        append=False, 
        parang=False)
# 
  print ("Began processing :"+scmode+str(count-1))
  applycal(vis=ms, 
  	selectdata=False,
  	gaintable=ms+'.'+scmode+str(count-1), 
  	parang=False,calwt=False,
  	applymode="calflag",
  	flagbackup=True)  
#
  tclean(vis=ms,
       datacolumn='corrected',
       imagename=ms+'.'+scmode+str(count-1),
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
       cyclefactor=1.3,
       sidelobethreshold=1.0,
       growiterations=75,
       weighting='briggs',robust=0, 
       usemask='user',mask=ms+'.'+scmode+str(count-2)+'.mask',
       savemodel='modelcolumn')
 
  exportfits(imagename=ms+'.'+scmode+str(count-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-1)+'.fits')
  os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-1)+".ms")
  print ("Made : " +scmode+str(count-1))
 
  count = count + 1

print ("Completed phase only self-calibration, going to A&P calibration Cycle")
#


count=count
for j in range(apcycles):  
  scmode='ap'
  if j>=3:
   sfactor=4
  else:
   sfactor=count-pcycles-1   
#
  print ("Began doing self-cal on :"+scmode+str(count-pcycles-1))
  gaincal(vis=ms,
  	caltable=ms+'.'+scmode+str(count-pcycles-1),
  	selectdata=False,
  	solint=str(apsolint/2**sfactor)+'min',
  	refant=ref_ant,
  	refantmode="strict", 
  	solnorm=True,
        minsnr=3.0,
        gaintype="G",
        calmode=scmode,
        append=False, 
        parang=False)
 
  applycal(vis=ms, 
  	selectdata=False,
  	gaintable=ms+'.'+scmode+str(count-pcycles-1), 
  	parang=False,
  	calwt=False,
  	applymode="calflag",
  	flagbackup=True)  
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
       cyclefactor=1.3,
       sidelobethreshold=1.0,
       growiterations=75,
       weighting='briggs',robust=0, 
       usemask='user',mask=maskfile,
       savemodel='modelcolumn')

  exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'.fits')
  os.system("cp -r "+ms+" "+ms+"."+scmode+str(count-pcycles-1)+".ms")
  print ("Made : " +scmode+str(count-pcycles-1))
  
  count = count + 1
#

print ("Completed processing AP self-calibrations\n")


####################################################################################################################################
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
       deconvolver="mtmfs",
       pbcor=False,
       weighting="briggs",robust=0,
       niter=int(startniter*2**count),gain=0.1,
       threshold="0.005mJy",
       cycleniter=-1,cyclefactor=1.3,
       interactive=False,usemask="user",mask=ms+'.'+scmode+str(count-pcycles-1)+'.mask',
       pbmask=0.0,
       sidelobethreshold=1.0,
       growiterations=75,
       restart=True,
       savemodel="none",calcres=True,
       calcpsf=True,parallel=False)

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
       weighting="briggs",robust=0,
       niter=int(startniter*2**count),gain=0.1,
       threshold="0.005mJy",
       cycleniter=-1,cyclefactor=1.3,
       interactive=False,usemask="user",mask=ms+'.'+scmode+str(count-pcycles-1)+'.mask',
       pbmask=0.0,
       sidelobethreshold=1.0,
       growiterations=75,
       restart=True,
       savemodel="none",calcres=True,
       calcpsf=True,parallel=False)

exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'_U'+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'_U'+'.fits' )

print ("Creating Stokes V image")

tclean(vis=ms,
	selectdata=True,
	datacolumn='corrected', 
	imagename=ms+'.'+scmode+str(count-pcycles-1)+'_V',
	imsize=imagesize,
	cell=cellsize,
       stokes="V",
       projection="SIN",
       specmode="mfs",
       reffreq="",
       nchan=-1,
       gridder="widefield",
       wprojplanes=wproj,
       pblimit=-1,
       deconvolver="mtmfs",
       pbcor=False,
       weighting="briggs",robust=0,
       niter=int(startniter*2**count),gain=0.1,
       threshold="0.005mJy",
       cycleniter=-1,cyclefactor=1.3,
       interactive=False,usemask="user",mask=ms+'.'+scmode+str(count-pcycles-1)+'.mask',
       pbmask=0.0,
       sidelobethreshold=1.0,
       growiterations=75,
       restart=True,
       savemodel="none",calcres=True,
       calcpsf=True,parallel=False)

exportfits(imagename=ms+'.'+scmode+str(count-pcycles-1)+'_V'+'.image.tt0', fitsimage=ms+'.'+scmode+str(count-pcycles-1)+'_V'+'.fits' )



print ("Done")

#The following steps are for imaging using VLA pipeline
'''

#Recipes for making images using NRAO pipeline
hifv_plotsummary(pipelinemode="automatic")
#hif_makeimlist(intent='PHASE,BANDPASS', specmode='cont')
#hif_makeimages(hm_masking='centralregion')
#hifv_exportdata(pipelinemode="automatic")

h_save()

###########################################
#  END NRAO PIPELINE
###########################################
'''
                              

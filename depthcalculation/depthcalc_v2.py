#This is the workable version of j160116depth.py
#latest update: June 17 2020
#KMJ
########################

import numpy as np
import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,random,time
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import scipy.stats as ss
from astropy.modeling import models, fitting
import sys
from array import array
import pylab as pl

###########read in contents of file and set up output files
infile=sys.argv[1]
fitsfilename='../HSTdata/uds_xmmlss0218/uds-100mas-f160w_drz_sci.fits'
valscaled=np.loadtxt(infile,dtype=float,delimiter='\t',usecols=2,comments='#')
hdu1=fits.open(fitsfilename)
hdr1=hdu1[0].header
#print(hdr1['PHOTPLAM'])
tmpfile=infile
tmplist=tmpfile.split('/')
tmpname=tmplist[len(tmplist)-1]
outname=tmpname.replace('.dat','')
#print(outname)


#########calculate average flux and stdev
valsin=valscaled
stdev_empty=np.std(valsin)
avg_emptyflux=np.average(valsin)
minval=np.amin(valsin)
maxval=np.amax(valsin)


#########set up parameters
temp_photflam=hdr1['PHOTFLAM']
temp_photplam=hdr1['PHOTPLAM']
temp_photfnu=hdr1['PHOTFNU']
#temp_exptime=hdr1['EXPTIME']

abmagzpt=-2.5*np.log10(temp_photflam) - 21.10 - 5.0*np.log10(temp_photplam) + 18.692
#experimental_abzpt=abmagzpt+2.5*np.log10(temp_exptime)
abmag=-2.5*np.log10(5.0*stdev_empty)+abmagzpt
print(abmag)



#for F814W
#photflam=7.039170025e-20
#photplam= 8.0449922e3
#exptime = 10000.0
#fl=photflam*photplam*photplam
#fv=fl*3.34*(10.0**4.0)*stdev_empty
#fv2=photplam*photplam/(2.998*(10.0**18.0)) # lambda squared/c in angstroms/s so that angstrom *angstrom /angstrom/s = Ang*s= Ang/hz
#then ang/hz * erg/cm2/A/electron * electron/s = erg/cm2/s/A *A/hz = erg/cm2/s/hz, perfectz
#fv3= photflam*5.0*stdev_empty*fv2
#mabcgs=-2.5*np.log10(fv3)-48.6
# according to SE for dummies, IF the image input is the sum of N frames, the magnitude zeropoint has to then include +2.5np.log10(exptime)
# now, our units are in counts/second, buuuut I think our astrodrizzled images are technically summed. so let me try it
#mabcgs_extra = mabcgs - 2.5*np.log10(exptime)
#
#print(mabcgs_extra)
#
#
##print(fv3)
#print('stdev==',stdev_empty)
#print('mabcgs=======',mabcgs)
#STMAG = -2.5*np.log10(5.0*stdev_empty*photflam) - 21.10
#ABMAG = STMAG - 5.0*np.log10(photplam) + 18.692
#print('stmag======',STMAG)
#print('abmagfromstmag =======',ABMAG)
#mabjy=-2.5*np.log10(5.0*fv)+8.90
#abmagnaut=-2.5*np.log10(fv)
#abmagnormzero=abmagnaut-48.6
#abmagzero=-2.5*np.log10(photflam) -21.10 -5.0*np.log10(photplam) +18.692
#abmagfull=abmagnaut-abmagzero

#print('mabjy==================',mabjy)
#print('abmagnaut============',abmagnaut)
#print('abmagnormzero===============',abmagnormzero)
#print('abmagzero==================',abmagzero)
#print('abmagfull==============',abmagfull)
#print('naut,normzero,zeropt,combined')

#magnitude calculationsA
#all these below ain't workin
#mag1=-2.5*np.log10(stdev_empty*photflam)-21.10
#mag3=-2.5*np.log10(3*stdev_empty*photflam)-21.10
#mag5=-2.5*np.log10(5*stdev_empty*photflam)-21.10
#print(mag5)
###below: my attempts to convert things to ABmagnitude. not apparently successful
#`part1=-2.5*np.log10(photflam)
#part2=-5.0*np.log10(photplam)
#abmag_zpt=part1+part2+18.692-21.10
#abmag_zpt=-2.5*np.log10(photflam) -21.10 -5*np.log10(photplam)+18.692
#print(abmag_zpt)


#convert to magnitude
#abmag_zpt is what you subtract off of the flux I think.
#mag_depth=-2.5log10(flux)-abmag_zpt?
#questions: is that + or -abmag_zpt? abmagzpt is negative
#does flux need to be in a different unit? I think it's relative Jy, but it might need to be in erg/s/cm2/Hz or A since photflam is in that format
#photometric zeropoint represents the magnitude of a star-like object that produces 1 count per second in a given aperture and system.i

#print(infile)
#print('1sig,3sig,5sig magnitude/depth')
#mag1out=-2.5*np.log10(stdev_empty)+abmag_zpt
#mag3out=-2.5*np.log10(3*stdev_empty)+abmag_zpt
#mag5out=-2.5*np.log10(5*stdev_empty)+abmag_zpt
#print(mag5out)

hdu1.close()
#yes, actually need this hdu1.close() line above

import os
import sys
import numpy as np
import scipy as sp
import astropy
from astropy.io import fits

filename=sys.argv[1]
whtfile=sys.argv[2]
hdu1=astropy.io.fits.open(filename)
hdu2=astropy.io.fits.open(whtfile)


#print(hdu1.info())
#print(hdu2.info())

#a quick investigation of headers and what info they contain
#hdr=hdu1[1].header
#print(list(hdr.keys()))
#it looks like the header for [1] is not as useful as the header for [0]. [0] is the science info and [1] is... the header info? so the header of the header info describes the format that the header itself has, which isn't as useful as describing the values of the header itself. 
#hdr2=hdu1[0].header
#print(list(hdr2.keys()))

data1=hdu1[0].data
data2=hdu2[0].data
print(data2.shape)
print(data1.shape)

outdata=data1*np.sqrt(data2)
hdout=fits.PrimaryHDU(outdata)
tmpfile=filename
tmplist=tmpfile.split('/')
tmpname=tmplist[len(tmplist)-1]
outname=tmpname.replace('.fits','')
newname=outname+'_norm.fits'
print(newname)
hdout.writeto(newname)

data3=hdu2[0].data

i1=np.where(data3>2)
i0=np.where(data3<2)

data3[i1]=1
data3[i0]=0

hdo2=fits.PrimaryHDU(data3)
newname2=outname+'_mask.fits'
hdo2.writeto(newname2)

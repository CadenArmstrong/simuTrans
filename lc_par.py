#!/usr/bin/env python
import os
import sys
import numpy as np
#import matplotlib
#from matplotlib import pyplot as plt
class lightcurve(object):

    def __init__(self,jd,mag,err=None,description="",cadence=0.02044):
        self.jd=jd
        self.mag=mag
        self.description=description
        self.cadence=cadence
        if err is None:
            if np.abs(cadence-1./60./24.)<1.e-7:
                self.err=np.zeros(len(jd))+2.7e-4
            else:
                self.err=np.zeros(len(jd))+2.7e-4/(cadence*60.*24.)**0.5
        else:
            self.err=err
        return

    def plot(self,tofile=False):
        if not tofile:
            plt.plot(self.jd,self.mag-np.median(self.mag),'.')
            plt.show()
        else:
            if self.description=="":
                raise IOError,"no output given for light curve figures"
            plt.plot(self.jd,self.mag-np.median(self.mag),'.')
            plt.savefig(os.path.splitext(self.description)[0]+'.png')

def read_lc(options):
    if not options.lc.inpath=="":
        if not os.path.exists(options.lc.inpath):
            raise IOError, "the input path for light curve file %s does not exist."
    try:
        if not os.path.exists(options.lc.inpath+'/'+options.lc.infile):
            raise IOError, "the light curve file %s does not exist."

        #tobedone, how to figure out the cadence
        #let's go with whitespace lcs first
        jd,mag=np.loadtxt(options.lc.inpath+'/'+options.lc.infile,usecols=(options.lc.coljd-1,options.lc.colmag-1),unpack=True)
        lcdata= [lightcurve(jd,mag,description=options.lc.inpath+'/'+options.lc.infile,cadence=options.lc.cadence)]
    except AttributeError:
        if not os.path.exists(options.lc.inpath+'/'+options.lc.inlist):
            raise IOError, "the light curve list %s does not exist."
        lcdata=[]
        #assuming the inlist only gives the file name, and files are under inpath
        fin=open(options.lc.inpath+'/'+options.lc.inlist,'r')
        for line in fin.readlines():
            if line.startswith("#"):
                continue
            infile,cadence=line.rstrip().split()
            infile=options.lc.inpath+'/'+infile
            if not os.path.exists(infile):
                raise IOError, "the light curve file %s does not exist." % infile

            jd,mag=np.loadtxt(infile,usecols=(options.lc.coljd-1,options.lc.colmag-1),unpack=True)
            lcdata.append(lightcurve(jd,mag,description=infile,cadence=float(cadence)))
    return lcdata


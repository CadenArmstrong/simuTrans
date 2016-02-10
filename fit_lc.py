import scipy as sp
import numpy as np
import cfg_parse as cfg
import cmd_parse as cmd
from lc_par import read_lc
import os
import sys
import emcee
from parameter import parameter
from Emceewrapper import mcmc_engine
from simplemodel import SimpleModel as tmodel 
#import matplotlib 
#from matplotlib import pyplot as plt

class Params():
    def __init__(self,options):
        self.paramarr=options.params
        vararr=np.zeros(len(self.paramarr)) 
        keyarr=[]
        self.lenfix=0
        self.transitmodel=tmodel()
        for i in xrange(len(self.paramarr)):
            vararr[i]=self.paramarr[i].val
            keyarr.append(self.paramarr[i].name)
            #print self.paramarr[i]
            if self.paramarr[i].fitflag==0:
                self.lenfix+=1
        self.paradic=dict(zip(keyarr,vararr))
        self.checkparam()
        
        return
    def __str__(self):
        string=""
        for key,value in self.paradic.iteritems():
            string+="%s:%f\n" % (key,value)
        return string

    def checkparam(self):
        #check and fill in all the default parameters

        #need to modify the value depend on q1 and q2
        if 'u1' not in self.paradic:
            self.paradic['u1']= self.paradic['q1']
        if 'u2' not in self.paradic:
            self.paradic['u2']= self.paradic['q2']
        #if 'b' not in self.paradic:
        #    self.paradic['b']=self.paradic[]
        if 'star_gridsize' not in self.paradic:
            self.paradic['star_gridsize']=5000
        if 'planet_gridsize' not in self.paradic:
            self.paradic['planet_gridsize']=200
        if 'gd_beta' not in self.paradic:
            self.paradic['gd_beta']=0
        if 'star_f' not in self.paradic:
            self.paradic['star_f']=0
        if 'planet_f' not in self.paradic:
            self.paradic['planet_f']=0
        if 'e' not in self.paradic:
            self.paradic['e']=0
        return

    def get_freeparams(self):
        return self.paramarr[self.lenfix:]

    def updatedic(self,freeparamarr):
        for i in xrange(self.lenfix,len(self.paramarr)):
            self.paradic[self.paramarr[i].name]=freeparamarr[i-self.lenfix]
        return 

    def model(self,cadence):
        self.transitmodel.SetupStar(np.array([self.paradic['star_gridsize'],self.paradic['u1'],self.paradic['u2'],self.paradic['gd_beta'],self.paradic['star_f']]))
        self.transitmodel.SetupPlanet(np.array([self.paradic['planet_gridsize'],self.paradic['b'],self.paradic['Rratio'],1./self.paradic['sma'],self.paradic['planet_f'],self.paradic['e']]))
        phase=self.cal_phase(cadence)
        #phase=np.arcsin((np.arange(50)-25.)/25.*0.75/5000.)
        print phase
        model_lc=np.zeros(len(phase))
        #print type(phase),type(model_lc)
        self.transitmodel.RelativeFlux(phase,model_lc)
        return model_lc

    def cal_phase(self,cadence):
        #calculate the phase of the planet orbit from the cadence
        period=self.paradic['P']
        epoch=self.paradic['T0']
        phase=np.pi*((cadence-epoch)/period-np.round((cadence-epoch)/period))
        return phase 

    def check_init(self,lcdata):
        model_lc=self.model(lcdata[0].jd)
        #phase=np.arcsin((np.arange(50)-25.)/25.*0.75/5000.)
        for i in xrange(len(lcdata[0].jd)):
            print lcdata[0].jd[i],model_lc[i]
        #plt.plot(lcdata[0].jd,lcdata[0].mag,'.')
        #plt.plot(lcdata[0].jd,1-model_lc+np.median(lcdata[0].mag),'+')
        #plt.show()
        return 
    def lc_chisq(self,freeparamarr,lcdata):
        lci = lcdata.copy()
        cadence = find_cadence(lci)
        self.updatedic(freeparamarr)
        model_lc=self.model(cadence)
        x0 = [median(lcdata[:,1])]
        def minfunc(x0):
            flux_ii = lcdata[:,1] + x0[0]
            return (flux_ii-model_lc)/lcdata[:,2]

        x0 = optimize.leastsq(minfunc,x0)
        lc[:,1] += x0[0]

        diff = lcdata[:,1] - model_lc
        chisq =  sum((diff/lcdata[:,2])**2)

        if std(model_lc) == 0:
            chisq = nan ### if it doesn't transit, return -1*inf
        return chisq

def main():
    options=cmd.fitlc_parse()
    cfg.fitlc_parse(options)
    #print options
    fitparams=Params(options)
    MC=mcmc_engine(options)
    print fitparams
    #return
    lcdata=read_lc(options)
    #lcdata[0].plot()
    fitparams.check_init(lcdata)
    return
    MC.run_mcmc(fitparams,lcdata)
if __name__=='__main__':
    main()

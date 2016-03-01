#!/usr/bin/env python
import scipy as sp
import numpy as np
from scipy import optimize
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
import time
class Params():
    def __init__(self,options):
        self.paramarr=options.params
        vararr=np.zeros(len(self.paramarr)) 
        keyarr=[]
        self.lenfree=0
        self.transitmodel=tmodel()
        self.qflag=False
        for i in xrange(len(self.paramarr)):
            vararr[i]=i
            keyarr.append(self.paramarr[i].name)
            #print self.paramarr[i]
            if self.paramarr[i].fitflag==1:
                self.lenfree+=1
            if self.paramarr[i].name=='q1':
                self.qflag=True
        self.paradic=dict(zip(keyarr,vararr))
        self.checkparam()
        
        return
    def __str__(self):
        string=""
        for key,value in self.paradic.iteritems():
            #print key,value
            string+="%s:%f\n" % (key,self.paramarr[int(value)].val)
        return string

    def readpara(self,name):
        #print self.paradic[name]
        return self.paramarr[int(self.paradic[name])]

    def checkparam(self):
        #check and fill in all the default parameters

        #need to modify the value depend on q1 and q2
        if 'u1' not in self.paradic:
            self.paramarr.append(self.readpara('q1'))
            self.paramarr[-1].name='u1'
            self.paradic['u1']=len(self.paramarr)-1 
        if 'u2' not in self.paradic:
            self.paramarr.append(self.readpara('q2'))
            self.paramarr[-1].name='u2'
            self.paradic['u2']=len(self.paramarr)-1 
        #if 'b' not in self.paradic:
        #    self.paradic['b']=self.paradic[]
        if 'star_gridsize' not in self.paradic:
            self.paramarr.append(parameter(1000,0,0,'star_gridsize'))
            self.paradic['star_gridsize']=len(self.paramarr)-1
        if 'planet_gridsize' not in self.paradic:
            self.paramarr.append(parameter(200,0,0,'planet_gridsize'))
            self.paradic['planet_gridsize']=len(self.paramarr)-1
        if 'gd_beta' not in self.paradic:
            self.paramarr.append(parameter(0,0,0,'gd_beta'))
            self.paradic['gd_beta']=len(self.paramarr)-1
        if 'star_f' not in self.paradic:
            self.paramarr.append(parameter(0,0,0,'star_f'))
            self.paradic['star_f']=len(self.paramarr)-1
        if 'planet_f' not in self.paradic:
            self.paramarr.append(parameter(0,0,0,'planet_f'))
            self.paradic['planet_f']=len(self.paramarr)-1
        if 'e' not in self.paradic:
            self.paramarr.append(parameter(0,0,0,'e'))
            self.paradic['e']=len(self.paramarr)-1
        return

    def get_freeparams(self):
        return self.paramarr[:self.lenfree]

    def updatedic(self,freeparamarr):
        #qflag=False
        for i in xrange(self.lenfree):
            self.paramarr[i].val=freeparamarr[i]
        #    if self.paramarr[i].name=='q1':
        #        qflag=True
        if self.qflag:
            self.readpara('u1').val=self.readpara('q1').val
            self.readpara('u2').val=self.readpara('q2').val
        return 

    def model(self,cadence):
        phase=self.cal_phase(cadence)
        #phase=np.arcsin((np.arange(50)-25.)/25.*0.75/5000.)
        #print phase
        model_lc=np.zeros(len(phase))
        #print type(phase),type(model_lc)
        self.transitmodel.RelativeFlux(phase,model_lc)
        return model_lc

    def cal_phase(self,cadence):
        #calculate the phase of the planet orbit from the cadence
        period=self.readpara('P').val
        epoch=self.readpara('T0').val
#phase=np.pi*((cadence-epoch)/period-np.round((cadence-epoch)/period))
        phase=np.pi*2.*((cadence-epoch)/period-np.round((cadence-epoch)/period))
        return phase 

    def check_init(self,lcdata):
        self.update()
        for i in xrange(len(lcdata)):
            model_lc=self.model(lcdata[0].jd)
            #for l in xrange(len(lcdata[0].jd)):
            #    print lcdata[0].jd[i],model_lc[i]
            #plt.plot(lcdata[0].jd,lcdata[0].mag,'.')
            #plt.plot(lcdata[0].jd,1-model_lc+np.median(lcdata[0].mag),'+')
            #plt.show()
        return

    def update(self):

        self.transitmodel.SetupStar(np.array([self.readpara('star_gridsize').val,self.readpara('u1').val,self.readpara('u2').val,self.readpara('gd_beta').val,self.readpara('star_f').val]))
        self.transitmodel.SetupPlanet(np.array([self.readpara('planet_gridsize').val,self.readpara('b').val,self.readpara('Rratio').val,1./self.readpara('sma').val,self.readpara('planet_f').val,self.readpara('e').val]))
        return

    def checkbound(self):
        for i in xrange(self.lenfree):
            if self.paramarr[i].val>self.paramarr[i].xmax or self.paramarr[i].val<self.paramarr[i].xmin:
                return False
        return True

    def lc_chisq(self,freeparamarr,lcdata):
        #lci = lcdata.copy()
        #cadence = find_cadence(lci)
        self.updatedic(freeparamarr)
        if not self.checkbound():
            return -np.inf
        chisq=0
        for i in xrange(len(lcdata)):
            self.update()
            model_lc=self.model(lcdata[i].jd)
            x0 = [np.median(lcdata[i].mag)]
            
            def minfunc(x0):
                flux_ii = lcdata[i].mag + x0[0]
                
                return (flux_ii-model_lc)/lcdata[i].err**2.

            x0 = optimize.leastsq(minfunc,x0)

            if np.std(model_lc) == 0:
                chisq = -np.inf ### if it doesn't transit, return -1*inf
                return chisq
            diff = lcdata[i].mag+x0[0] - model_lc
            chisq_i=sum((diff/lcdata[i].err)**2)
            chisq +=chisq_i 

        print self.readpara('u1').val,self.readpara('u2').val,self.readpara('b').val,self.readpara('sma').val,chisq
        return -chisq

def main():
    options=cmd.fitlc_parse()
    cfg.fitlc_parse(options)
    #print options
    fitparams=Params(options)
    MC=mcmc_engine(options)
    #print fitparams
    #return
    lcdata=read_lc(options)
    #lcdata[0].plot()
    starttime=time.time()
    fitparams.check_init(lcdata)
    print time.time()-starttime
    #print "before del"
    #del fitparams.transitmodel
    #print "after del"
    #print "end of check_init"
    return
    MC.run_mcmc(fitparams,lcdata)
    return
if __name__=='__main__':
    main()

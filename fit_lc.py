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

class Params():
    def __init__(self,options):
        #self.model=
        self.paramarr=options.params
        vararr=np.zeros(len(self.paramarr)) 
        keyarr=[]
        self.lenfix=0
        for i in xrange(len(self.paramarr)):
            vararr[i]=self.paramarr[i].val
            keyarr.append(self.paramarr[i].name)
            if self.paramarr[i].fitflag==0:
                self.lenfix+=1
        self.paradic=dict(zip(keyarr,vararr))
        return

    def get_freeparams(self):
        return self.paramarr[self.lenfix:]

    def updatedic(self,freeparamarr):
        for i in xrange(self.lenfix,len(self.paramarr)):
            self.paradic[self.paramarr[i].name]=freeparamarr[i-self.lenfix]
        return 

    def lc_chisq(self,freeparamarr,lcdata):
        lci = lcdata.copy()
        cadence = find_cadence(lci)
        self.updatedic(freeparamarr)
        model_lc=self.model(self.paradic,cadence)
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
    return
    lcdata=read_lc(options)
    MC.run_mcmc(fitparams,lcdata)

if __name__=='__main__':
    main()

import scipy as sp
import numpy as np
import cfg_parse as cfg
import cmd_parse as cmd
from lc_par import read_lc
import os
import sys
import emcee

class parameter():
    def __init__(self,val,upper,lower,name,fitflag=0):
        self.val=val
        self.upper=upper
        self.lower=lower
        self.name=name
        #fitflag=0 indicate fixed parameter
        self.fitflag=fitflag 
        return
    def __str__(self):
        return "%s:%f + %f - %f" % (self.name,self.val,self.upper-self.val,self.val=self.lower)

    def __call__(self):
        return self.val

class Params():
    def __init__(self,options):
        #self.model=
        #self.paramarr=
        return

    def get_freeparams(self):
        free_params=[]
        for i in xrange(len(self.paramarr)):
            if self.paramarr[i].fitflag==1:
                free_params.append(self.paramarr)
        return free_params

    def lc_chisq(self,lcdata):
        lci = lcdata.copy()
        cadence = find_cadence(lci)
        model_lc=self.model(self.paramarr,cadence)
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

class mcmc_engine():
    def __init__(self,options):
        self.nwalkers=options.nwalkers
        self.niter=options.niter
        self.nburn=options.nburn
        self.nthreads=options.nthreads
        self.output=options.output
        return

    def run_mcmc(self,fitparams,lcdata):
        free_params=fitparams.get_freeparams()
        ndim = len(free_params)
        #initial position
        p0 = []
        for i in range(self.nwalkers):
            pi = []
            for j in range(len(free_params)):
                print free_params[j]
                pi_i = random.normal(free_params[j].var,0.01*(free_params[j].upper-free_params[j].lower))
                while pi_i > free_params[j].upper or pi_i < free_params[j].lower:
                    pi_i = random.normal(free_params[j].var,0.01*(free_params[j].upper-free_params[j].lower))

            
                pi.append(pi_i)
            p0.append(pi)

        #burn

        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, fitparams.lc_chisq, args=[lcdata],threads=nthreads)

        pos, prob, state = sampler.run_mcmc(p0, self.nburn)


        #real iteration
        sampler.reset()
    
        master_pos = []

        for result in sampler.sample(pos, iterations=self.nmcmc, storechain=False):
            position,probability = result[0],result[1]

            for i in range(len(position)):
                if functions.isnan(probability[i]):
                    master_pos.append(list(position[i])+[probability[i]])

        print "iteration finished"
        master_pos = np.array(master_pos)
        np.savetxt(self.output,master_pos,fmt="%.10f")



        return

def main():
    options=cmd.fitlc_parse()
    cfg.fitlc_parse(options)
    fitparams=Params(options)
    MC=mcmc_engin(options)
    lcdata=read_lc(options)
    MC.run_mcmc(fitparams,lcdata)

if __name__=='__main__':
    main()

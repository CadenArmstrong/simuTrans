#!/usr/bin/env python
import numpy as np
from numpy import random
import parameter
import emcee
class mcmc_engine():
    def __init__(self,options):
        self.nwalkers=options.mcmc.nwalkers
        self.niter=options.mcmc.niter
        self.nburn=options.mcmc.nburn
        self.nthreads=options.mcmc.nthreads
        self.output=options.mcmc.output
        return

    def run_mcmc(self,fitparams,lcdata):
        free_params=fitparams.get_freeparams()
        ndim = len(free_params)
        #initial position
        p0 = []
        for i in range(self.nwalkers):
            pi = []
            for j in range(len(free_params)):
                #print free_params[j]
                pi_i = random.normal(free_params[j].val,0.01*(free_params[j].upper-free_params[j].lower))
                while pi_i > free_params[j].upper or pi_i < free_params[j].lower:
                    pi_i = random.normal(free_params[j].val,0.01*(free_params[j].upper-free_params[j].lower))

            
                pi.append(pi_i)
            p0.append(pi)

        #burn

        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, fitparams.lc_chisq, args=[lcdata],threads=self.nthreads)

        pos, prob, state = sampler.run_mcmc(p0, self.nburn)


        #real iteration
        sampler.reset()
    
        #master_pos = []
        f=open(self.output,"a")
        #f.close()
        for result in sampler.sample(pos, iterations=self.niter, storechain=False):
            position,probability = result[0],result[1]
            #f=open(self.output,"a")
            for i in range(position.shape[0]):
                if np.isnan(probability[i]):
                    continue
                #f.write("{0:4d} {1:s} {1:s}\n".format(i," ".join(str(position[i])),probability[i]))
                f.write(str(i)) 
                for j in xrange(len(position[i])):
                    f.write(" %f" % position[i][j])
                f.write(" %f\n" % probability[i])
                #f.write("{0:4d} {1:s} \n".format(i," ".join(list(str(position[i])))))
        f.close()
        print "iteration finished"



        return



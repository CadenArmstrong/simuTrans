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



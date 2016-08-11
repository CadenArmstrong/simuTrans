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
import ZeipelModel
import LaraModel
from simplemodel import SimpleModel as tmodel 
import time
import copy


class PickalableSWIG:
    def __setstate__(self, state):
        self.__init__(*state['args'])

    def __getstate__(self):
        return {'args': self.args}

class PickalableC(tmodel, PickalableSWIG):

    def __init__(self, *args):
        self.args = args
        tmodel.__init__(self,*args)




class Params(object):
    def __init__(self,options):
        self.paramarr=options.params
        vararr=np.zeros(len(self.paramarr)) 
        keyarr=[]
        self.lenfree=0
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
        self.requiredpara={'star_gridsize':1000,'u1':0.0,'u2':0.0,'gd_beta':0.0,'star_f':0.0,'phi':0.0,'Mstar':1.0,'Rstar':1.0,'Prot':8.4,'planet_gridsize':200,'b':0,'Rratio':0.1,'sma':0.03,'lambda':0.0,'e':0.0,'planet_f':0.0,'gd_flag':1,'P':3.0,'T0':0.0,'b2':0} 
        self.checkparam()
        self.transitmodel=PickalableC(int(self.readpara('star_gridsize').val), int(self.readpara('planet_gridsize').val))
        return
    def __str__(self):
        string=""
        for key,value in self.paradic.iteritems():
            string+="%s\n" % (self.paramarr[int(value)])
        return string
    def __call__(self, freeparamarr,lcdata):

        return self.lc_chisq(freeparamarr,lcdata)
    def almosteq(self,a,b,tol=1.e-7):
        
        return np.abs(a-b)<tol

    def readpara(self,name):
        #print self.paradic[name]
        return self.paramarr[int(self.paradic[name])]

    def checkparam(self):
        #check and fill in all the default parameters

        #need to modify the value depend on q1 and q2
        if 'u1' not in self.paradic or 'u2' not in self.paradic:
            if 'q1' not in self.paradic or 'q2' not in self.paradic:

                self.paramarr.append(parameter(self.requiredpara['u1'],0.,0.,'u1'))
                self.paradic['u1']=len(self.paramarr)-1
                self.paramarr.append(parameter(self.requiredpara['u2'],0.,0.,'u2'))
                self.paradic['u2']=len(self.paramarr)-1
            else:
                u1,u2=self.cal_LD_qtou(self.readpara('q1').val,self.readpara('q2').val)
                self.paramarr.append(parameter(u1,0.,0.,'u1'))
                self.paradic['u1']=len(self.paramarr)-1 
                self.paramarr.append(parameter(u2,0.,0.,'u2'))
                self.paradic['u2']=len(self.paramarr)-1
        if 'b2' not in self.paradic:
            if 'b' not in self.paradic:
                self.paramarr.append(parameter(self.requiredpara['b'],0.,0.,'b'))
                self.paramarr.append(parameter(self.requiredpara['b2'],0.,0.,'b2'))
            else:
                para_b=copy.deepcopy(self.readpara('b'))
                para_b2=parameter(para_b.val**2.,para_b.upper**2.,para_b.lower**2.,'b2',fitflag=para_b.fitflag,xmin=0,xmax=1)
                if para_b.fitflag==1:
                    #print self.paradic['b']
                    self.paramarr[int(self.paradic['b'])]=para_b2
                    self.paradic['b2']=self.paradic['b']
                    #self.paramarr.append(para_b)
                    #self.paradic['b']=len(self.paramarr)-1
                    #self.readpara('b').fitflag=0
        for key,value in self.requiredpara.iteritems():
            if key not in self.paradic:
                self.paramarr.append(parameter(value,0.,0.,key))
                self.paradic[key]=len(self.paramarr)-1
        return

    def cal_LD_qtou(self,q1,q2):
        #kipping 2013, eq 15,16
        u1=2.*np.sqrt(q1)*q2
        u2=np.sqrt(q1)*(1-2*q2)
        return [u1,u2]

    def get_freeparams(self):
        return self.paramarr[:self.lenfree]

    def updatedic(self,freeparamarr):
        #qflag=False
        for i in xrange(self.lenfree):
            self.paramarr[i].val=freeparamarr[i]
        #    if self.paramarr[i].name=='q1':
        #        qflag=True
        if self.qflag:
            u1,u2=self.cal_LD_qtou(self.readpara('q1').val,self.readpara('q2').val)
            self.readpara('u1').val=u1
            self.readpara('u2').val=u2
        return 
    def getshortcadence(self,jd,cadence,Nresample=5):
        tmax=np.max(jd)+cadence/2.
        tmin=np.min(jd)-cadence/2.
        ncadence=(tmax-tmin)/(cadence/Nresample)
        shortcadence=tmin+np.arange(ncadence)*(cadence/Nresample) 
        return shortcadence
    def getlc(self,jd,cadence,shortcadence,model_sc,Nresample=5):
        model_lc=np.zeros(len(jd))
        t0=shortcadence[0]
        length=len(shortcadence)
        index1=(((jd+0.5*cadence)-t0)/(cadence/Nresample)).astype(int)
        index2=(((jd-0.5*cadence)-t0)/(cadence/Nresample)).astype(int)
        np.where(index2<0,index2,0)
        np.where(index1>(length-1),index1,length-1)
        for i in xrange(len(model_lc)): 
            model_lc[i]=np.mean(model_sc[index2[i]:(index1[i]+1)])
        return model_lc
    def model(self,jd,cadence=1./60./24.):
        if self.almosteq(cadence,1./60./24):
            shortcadence=jd
        else:
            shortcadence=self.getshortcadence(jd,cadence)
        phase=self.cal_phase(shortcadence)
        #phase=np.arcsin((np.arange(50)-25.)/25.*0.75/5000.)
        #print phase

        model_sc=np.zeros(len(phase))
        #print type(phase),type(model_lc)
        self.transitmodel.RelativeFlux(phase,model_sc)
        if self.almosteq(cadence,1./60./24):
            return model_sc
        else:
            model_lc=self.getlc(jd,cadence,shortcadence,model_sc)
            return model_lc
    def cal_phase(self,jd):
        #calculate the phase of the planet orbit from the cadence
        period=self.readpara('P').val
        epoch=self.readpara('T0').val
#phase=np.pi*((cadence-epoch)/period-np.round((cadence-epoch)/period))
        phase=np.pi*2.*((jd-epoch)/period-np.round((jd-epoch)/period))
        return phase 

    def check_init(self,lcdata):
        self.update()
        for i in xrange(len(lcdata)):
            model_lc=self.model(lcdata[i].jd,lcdata[i].cadence)
            #for l in xrange(len(lcdata[i].jd)):
            #    print lcdata[i].jd[l],model_lc[l]
            try:
                fig=plt.figure()
                ax=fig.add_subplot(111)
                #ax.plot(lcdata[i].jd,lcdata[i].mag,'.')
                ax.plot(lcdata[i].jd,model_lc,'.')
                #ax.plot(lcdata[i].jd,1-model_lc+np.median(lcdata[i].mag)-lcdata[i].mag,'+')
                #ax.plot(lcdata[i].jd,model_lc-1+np.median(lcdata[i].mag)-lcdata[i].mag,'+')
                #ax.plot(lcdata[i].jd,1-model_lc+np.median(lcdata[i].mag),'+')
                #ax.plot(lcdata[i].jd,model_lc-1+np.median(lcdata[i].mag),'+')
                #phase=self.cal_phase(lcdata[i].jd)
                #ax.plot(lcdata[i].jd,1-model_lc,'+')
                #ax.plot(phase,1-model_lc,'+')
                #phase=self.cal_phase(lcdata[i].jd)/2./np.pi*self.readpara('P').val*3600.*24.
                #ax.plot(phase,model_lc,'+')
                #ax.set_xlim([-8000,8000])
                y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
                ax.yaxis.set_major_formatter(y_formatter)
                plt.show()
                model_lc=1.-model_lc
                x0 = [np.median(lcdata[i].mag)]
                
                def minfunc(x0):
                    flux_ii = lcdata[i].mag + x0[0]
                    
                    return ((flux_ii-model_lc)/lcdata[i].err)**2.

                x0 = optimize.leastsq(minfunc,x0)

                if np.std(model_lc) == 0:
                    chisq = -np.inf ### if it doesn't transit, return -1*inf
                    return chisq
                diff = lcdata[i].mag+x0[0] - model_lc
                chisq_i=sum((diff/lcdata[i].err)**2)
                print  chisq_i, np.mean(lcdata[i].err),np.mean(diff),np.std(diff)


            except NameError:
                #raise
                continue
        return

    def update(self):
        #Prot in hours, Mstar in Msun, Rstar in Rsun, graivity constant G is already in the coeff
        groteq=(2.*np.pi)**2*0.19567/self.readpara('Prot').val**2./self.readpara('Mstar').val*self.readpara('Rstar').val**3.       
        #groteq=1.9567/(self.readpara('Prot').val)**2./self.readpara('Mstar').val*self.readpara('Rstar').val**3.       
        #print groteq
        #return
        #print np.array([self.readpara('star_gridsize').val,self.readpara('u1').val,self.readpara('u2').val,self.readpara('star_f').val,self.readpara('phi').val*np.pi/180.,groteq,self.readpara('gd_beta').val,self.readpara('gd_flag').val])
        self.transitmodel.SetupStar(np.array([self.readpara('star_gridsize').val,self.readpara('u1').val,self.readpara('u2').val,self.readpara('star_f').val,self.readpara('phi').val*np.pi/180.,groteq,self.readpara('gd_beta').val,self.readpara('gd_flag').val]))
        self.transitmodel.SetupPlanet(np.array([self.readpara('planet_gridsize').val,np.sqrt(self.readpara('b2').val),self.readpara('Rratio').val,1./self.readpara('sma').val,self.readpara('lambda').val*np.pi/180.,self.readpara('e').val, self.readpara('planet_f').val]))
        return

    def checkbound(self):
        for i in xrange(self.lenfree):
            if self.paramarr[i].name=='phi' or self.paramarr[i].name=='lambda':
                if self.paramarr[i].val>self.paramarr[i].xmax:
                    self.paramarr[i].val=self.paramarr[i].val-180.*int((self.paramarr[i].val-self.paramarr[i].xmin)/180.)
                if self.paramarr[i].val<self.paramarr[i].xmin:
                    self.paramarr[i].val=self.paramarr[i].val+180.*int((-self.paramarr[i].val+self.paramarr[i].xmax)/180.)
                continue
            if self.paramarr[i].val>self.paramarr[i].xmax or self.paramarr[i].val<self.paramarr[i].xmin:
                return False
        return True

    def lc_chisq(self,freeparamarr,lcdata):
        #lci = lcdata.copy()
        #cadence = find_cadence(lci)
        self.updatedic(freeparamarr)
        print freeparamarr
        if not self.checkbound():
            return -np.inf
        chisq=0
        for i in xrange(len(lcdata)):
            self.update()
            model_lc=1.-self.model(lcdata[i].jd,lcdata[i].cadence)
            x0 = [np.median(lcdata[i].mag)]
            
            def minfunc(x0):
                flux_ii = lcdata[i].mag + x0[0]
                
                return ((flux_ii-model_lc)/lcdata[i].err)**2.

            x0 = optimize.leastsq(minfunc,x0)

            if np.std(model_lc) == 0:
                chisq = -np.inf ### if it doesn't transit, return -1*inf
                return chisq
            diff = lcdata[i].mag+x0[0] - model_lc
            chisq_i=sum((diff/lcdata[i].err)**2)
            chisq +=chisq_i 

        #print 'u1,u2,b2,sma=',self.readpara('u1').val,self.readpara('u2').val,self.readpara('b2').val,self.readpara('sma').val,chisq
        return -chisq

def main():
    options=cmd.fitlc_parse()
    if options.plot:
        import matplotlib 
        from matplotlib import pyplot as plt
        global matplotlib
        global plt
    cfg.fitlc_parse(options)
    #print options
    #return
    fitparams=Params(options)

    MC=mcmc_engine(options)
    print fitparams
    #return
    lcdata=read_lc(options)
    #lcdata[0].plot()
    starttime=time.time()
    fitparams.check_init(lcdata)
    print '#',time.time()-starttime
    #print "before del"
    #del fitparams.transitmodel
    #print "after del"
    #print "end of check_init"
    #return
    if not options.plot:
        MC.run_mcmc(fitparams,lcdata)
    return
if __name__=='__main__':
    main()

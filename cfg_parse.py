#!/usr/bin/python
#cfg_parse.py use to parse all the configure parameters from .cfg file for 
#functions related to run lc function. 

import ConfigParser
import cmd_parse as cmdp
from parameter import parameter
def set_parse(infile):
    p = ConfigParser.RawConfigParser()
  
    p.add_section('Section LC')
    p.set('Section LC','infile','') 
    p.set('Section LC','inpath','') 
    p.set('Section LC','coljd',1) 
    p.set('Section LC','colmag',7) 
    p.set('Section LC','cadence',0.02044) 
    p.set('Section LC','inlist','')

    p.add_section('Section MCMC')
    p.set('Section MCMC','nwalkers',1)
    p.set('Section MCMC','nburn',50)
    p.set('Section MCMC','niter',1000)
    p.set('Section MCMC','nthreads',1)
    p.set('Section MCMC','output','')

    p.add_section('Section FixedParameters')
    p.set('Section FixedParameters','fixparams','P,T0,Rratio') 
    p.set('Section FixedParameters','P','100.0') 
    p.set('Section FixedParameters','T0','891.4') 
    p.set('Section FixedParameters','Rratio','0.1') 
    
    p.add_section('Section FreeParameters')
    p.set('Section FreeParameters','freeparams','sma,q1,q2,b') 
    p.set('Section FreeParameters','sma','0.1,0.02,0.02') 
    p.set('Section FreeParameters','q1','0.2,0.1,0.1') 
    p.set('Section FreeParameters','q2','0.2,0.1,0.1') 
    p.set('Section FreeParameters','b','0.6,0.1,0.1') 
    
    p.add_section('This is a configure file for xxx package')
    
    with open (infile,'wb') as configfile:
	p.write(configfile)	
    return 


def fitlc_parse(options):
    p = ConfigParser.RawConfigParser()
    try:
        p.read(options.cfgfile)
    except IOError:
        print "configure file missing" 
        raise
    #parse lc format
    try:
        setattr(options.lc,'infile',p.get('Section LC',"infile"))
        setattr(options.lc,'cadence',p.get('Section LC',"cadence"))
    except ConfigParser.MissingSectionHeaderError:	
        raise 'Error: Section LC missing, excute set_parse to see example.cfg'
    except ConfigParser.NoOptionError:
        try:
            setattr(options.lc,'inlist',p.get('Section LC',"inlist"))
        except ConfigParser.NoOptionError:
            raise 'Error: Can not find neither a light curve file and its cadence or a list of light curve files for fitting purpose'

    try:
        setattr(options.lc,'inpath',p.get('Section LC',"inpath"))
    except ConfigParser.NoOptionError:
        print 'Warning: Can not find the input path of light curves, use current directory by default'
        setattr(options.lc,'inpath',"")
    try:
        setattr(options.lc,'coljd',int(p.get('Section LC',"coljd")))
    except ConfigParser.NoOptionError:
        print 'Warning: Can not find the input column of jds in light curves, use column 1 by default'
        setattr(options.lc,'coljd',1)
    try:
        setattr(options.lc,'colmag',int(p.get('Section LC',"colmag")))
    except ConfigParser.NoOptionError:
        print 'Warning: Can not find the input column of jds in light curves, use column 2 by default'
        setattr(options.lc,'colmag',2)

    #parse fitting options
    setattr(options.mcmc,'nwalkers',int(p.get('Section MCMC',"nwalkers")))
    setattr(options.mcmc,'nburn',int(p.get('Section MCMC',"nburn")))
    setattr(options.mcmc,'niter',int(p.get('Section MCMC',"niter")))
    setattr(options.mcmc,'nthreads',int(p.get('Section MCMC',"nthreads")))
    setattr(options.mcmc,'output',p.get('Section MCMC',"output"))

    #parse parameter options
    freeparams=p.get('Section FreeParameters','freeparams').split(',')
    for i in xrange(len(freeparams)):
        try:
            freeargs=p.get('Section FreeParameters',freeparams[i]).split(',')
            if len(freeargs)==3:
                val,upper,lower=freeargs
                options.params.append(parameter(float(val),float(upper),float(lower),freeparams[i],fitflag=1))
            elif len(freeargs)==5:
                val,upper,lower,xmin,xmax=freeargs
                options.params.append(parameter(float(val),float(upper),float(lower),freeparams[i],fitflag=1,xmin=xmin,xmax=xmax))
            else:
                raise ValueError
        except ValueError:
            raise ValueError,'Free Parameters %s need to be given as val,upper,lower (,xmin,xmax) format in configure file' % freeparams[i]

    fixparams=p.get('Section FixedParameters','fixparams').split(',')
    for i in xrange(len(fixparams)):
        options.params.append(parameter(float(p.get('Section FixedParameters',fixparams[i])),0,0,fixparams[i],fitflag=0))


    return



if __name__=='__main__':
	options=cmdp.runcfg()
	if(options.eflag):	
		infile=options.infile
		set_parse(infile)
		


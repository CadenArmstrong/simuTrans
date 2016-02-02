#!/usr/bin/python
#cfg_parse.py use to parse all the configure parameters from .cfg file for 
#functions related to run lc function. 

import ConfigParser
import cmd_parse as cmdp
from fit_lc import parameter
def set_parse(infile):
    p = ConfigParser.RawConfigParser()
  
    p.add_section('Section LC')
    p.set('Section LC','infile','') 
    p.set('Section LC','inpath','') 
    p.set('Section LC','coljd',1) 
    p.set('Section LC','colmag',7) 
    p.set('Section LC','inlist','')

    p.add_section('Section MCMC')
    p.set('Section MCMC','nwalker','')
    p.set('Section MCMC','nburn','')
    p.set('Section MCMC','niter','')

    p.add_section('Section FixedParameters')
    p.set('Section FixedParameters','fixparams','P,T0,Rratio') 
    p.set('Section FixedParameters','P','100.0') 
    p.set('Section FixedParameters','T0','891.4') 
    p.set('Section FixedParameters','Rratio','0.1') 
    
    p.add_section('Section FreeParameters')
    p.set('Section FreeParameters','freeparams','sma,q1,q2,b') 
    p.set('Section FreeParameters','sma','0.1,0.02,0.02') 
    p.set('Section FreeParameters','q1','0.2,0.1,0.1') 
    p.set('Section FreeParameters','q1','0.2,0.1,0.1') 
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
    except ConfigParser.MissingSectionHeaderError:	
        raise 'Error: Section LC missing, excute set_parse to see example.cfg'
    except ConfigParser.NoOptionError:
        try:
            setattr(options.lc,'inlist',p.get('Section LC',"inlist"))
        except ConfigParser.NoOptionError:
            raise 'Error: Can not find neither a light curve file or a list of light curve files for fitting purpose'

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
    setattr(options.mcmc,'nwalker',int(p.get('Section MCMC',"nwalker")))
    setattr(options.mcmc,'nburn',int(p.get('Section MCMC',"nburn")))
    setattr(options.mcmc,'niter',int(p.get('Section MCMC',"niter")))

    #parse parameter options
    fixparams=p.get('Section FixedParameters','fixparams').split(',')
    for i in xrange(len(fixparams)):
        options.params.append(parameter(float(p.get('Section FixedParameters',fixparams[i])),0,0,fixparams[i],0))

    freeparams=p.get('Section FreeParameters','freeparams').split(',')
    for i in xrange(len(freeparams)):
        try:
            val,upper,lower=p.get('Section FreeParameters',freeparams[i]).split(',')
        except ValueError:
            raise 'Free Parameters need to be given as val,upper,lower format in configure file'
        options.params.append(parameter(float(val),float(upper),float(lower),freeparams[i],1))
    return



if __name__=='__main__':
	options=cmdp.runcfg()
	if(options.eflag):	
		infile=options.infile
		set_parse(infile)
		


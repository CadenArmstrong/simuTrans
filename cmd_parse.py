#!/usr/bin/python
#lc_parse.py use to parse all the command lines for functions realted to run lc function. 
import optparse
from fancyfont import colors
import os
tc=colors()

def runcfg():
    p = optparse.OptionParser()
    p.add_option('--eflag','-e',default=False,action='store_true',help='Creat example.cfg for all the parameter settings')
    p.add_option('--infile','-i',default='example.cfg',help='the file to write cfg parameters')
    options,arguments=p.parse_args()
    return options

class fitlc_parse(object):
    def __init__(self):
        p=optparse.OptionParser()
        p.add_option('--config','-c',default='example.cfg',help='the configuration file for the program')
        p.add_option('--plot','-p',default=False,action='store_true',help='enable plotting')
        p.add_option('--infile','-i',default='',help='overwrite the input file from configure file')
        p.add_option('--outfile','-o',default='',help='overwrite the output file from configure file')
        options,arguments=p.parse_args()
        #self.cfgfile="example.cfg"
        self.infile=options.infile
        self.output=options.outfile
        self.cfgfile=options.config
        self.plot=options.plot
        if not os.path.exists(self.cfgfile):
            raise IOError, "configuration file %s does not exist" % self.cfgfile
        self.mcmc=lambda:None
        self.lc=lambda:None
        self.params=[]
        return

    def __str__(self):
        strparse=tc.red("Fitting has set up as the following\n")
        strparse+=tc.blue("MCMC Params:\n")
        for key,value in self.mcmc.__dict__.iteritems(): 
            strparse+=(key+"="+str(value)+"\n") 
        strparse+=tc.blue("LC Params:\n")
        for key,value in self.lc.__dict__.iteritems(): 
            strparse+=(key+"="+str(value)+"\n") 
        strparse+=tc.blue("Model Params:\n")
        for i in xrange(len(self.params)):
            strparse+=str(self.params[i])+"\n"
        return strparse

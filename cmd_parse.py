#!/usr/bin/python
#lc_parse.py use to parse all the command lines for functions realted to run lc function. 
import optparse
from fancyfont import colors

tc=colors()

def runcfg():
    p = optparse.OptionParser()
    p.add_option('--eflag','-e',default=False,action='store_true',help='Creat example.cfg for all the parameter settings')
    p.add_option('--infile','-i',default='example.cfg',help='the file to write cfg parameters')
    options,arguments=p.parse_args()
    return options

class fitlc_parse(object):
    def __init__(self):
        self.cfgfile="example.cfg"
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

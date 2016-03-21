#!/usr/bin/env python
import numpy as np

class parameter():
    def __init__(self,val,upper,lower,name,fitflag=0,xmin=-np.inf,xmax=np.inf):
        self.val=val
        self.upper=upper
        self.lower=lower
        self.name=name
        self.fitflag=fitflag 
        self.xmin=-np.inf
        self.xmax=np.inf
        return
    def __str__(self):
        if self.fitflag==0:
            return "#%s:%f" % (self.name,self.val)
        else:
            return "#%s:%f + %f - %f >%f <%f" % (self.name,self.val,self.upper-self.val,self.val-self.lower,self.xmin,self.xmax)


    def __call__(self):
        return self.val


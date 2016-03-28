#!/usr/bin/env python
import numpy as np

class parameter():
    def __init__(self,val,upper,lower,name,fitflag=0,xmin=-np.inf,xmax=np.inf):
        self.val=float(val)
        self.upper=float(upper)
        self.lower=float(lower)
        if fitflag:
            if self.upper<self.val or self.lower>self.val:
                raise ValueError, 'The upper (lower) error bound for the fitted variable %s is smaller (larger) than the median value' % name
            if float(xmax)<self.val or float(xmin)>self.val:
                raise ValueError, 'The upper %f (lower %f) hard limit for the fitted variable %s is smaller (larger) than the median value %f' % (float(xmin),float(xmax),name,float(val))
        self.name=name
        self.fitflag=fitflag 
        self.xmin=float(xmin)
        self.xmax=float(xmax)
        return
    def __str__(self):
        if self.fitflag==0:
            return "#%s:%f" % (self.name,self.val)
        else:
            return "#%s:%f + %f - %f >%f <%f" % (self.name,self.val,self.upper-self.val,self.val-self.lower,self.xmin,self.xmax)


    def __call__(self):
        return self.val


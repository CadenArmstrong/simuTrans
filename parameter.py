class parameter():
    def __init__(self,val,upper,lower,name,fitflag=0):
        self.val=val
        self.upper=upper
        self.lower=lower
        self.name=name
        self.fitflag=fitflag 

        return
    def __str__(self):
        if self.fitflag==0:
            return "%s:%f" % (self.name,self.val)
        else:
            return "%s:%f + %f - %f" % (self.name,self.val,self.upper-self.val,self.val-self.lower)


    def __call__(self):
        return self.val


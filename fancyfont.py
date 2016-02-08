
class colors(object):
    def __init__(self):
        self.HEADER = '\033[95m'
        self.OKBLUE = '\033[94m'
        self.OKGREEN = '\033[92m'
        self.WARNING = '\033[93m'
        self.FAIL = '\033[91m'
        self.ENDC = '\033[0m'
        self.BOLD = "\033[1m"
        self.UNDERLINE = '\033[4m'
    def red(self,msg):
        return self.FAIL + msg + self.ENDC
    
    def green(self,msg):
        return self.OKGREEN + msg + self.ENDC
    
    def blue(self,msg):
        return self.OKBLUE + msg + self.ENDC
    
    def yellow(self,msg):
        return self.WARNING + msg + self.ENDC
    def magenta(self,msg):
        return self.HEADER+msg+self.ENDC
    def bold(self,msg):
        return self.BOLD+msg+self.ENDC
    def underline(self,msg):
        return self.UNDERLINE+msg+self.ENDC

if __name__=='__main__':
    t=colors()
    print t.magenta('m') 
    print t.bold('k') 
    print t.underline('underline') 
    print t.red('r')
    print t.green('g')
    print t.blue('b')
    print t.yellow('y')



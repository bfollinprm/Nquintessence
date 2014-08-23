from cosmoslik import SlikPlugin
from numpy import hstack, zeros, arange, pi, inf, nan

class cmb_lik_try(SlikPlugin):
    """
    
    """
    
    def __init__(self):
        
        super(cmb_lik_try,self).__init__()#**nuisance)
        
      
    def __call__(self, **kwargs):
        cmb = {}
        for key, value in kwargs.iteritems():
            cmb[key] = value
	
        #for key, value in cmb.iteritems():
        #    print "cmb[",key,"]=",value        

        #nuisance = [self[k] for k in self.clik.get_extra_parameter_names()]
        lnl = 1        
        #lnl = -self.clik(hstack([cl,nuisance]))[0]
        if lnl==0: return inf
        else: return lnl
            
    
#def tocl(dl): 
#    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
    

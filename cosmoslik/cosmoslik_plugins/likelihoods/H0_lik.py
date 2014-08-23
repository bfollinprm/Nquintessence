from cosmoslik import SlikPlugin
from numpy import dot

class H0_lik(SlikPlugin):
    """
    get likelihood for H0 given model predictions
    UNITS ARE Mpc^-1 to match classy output. 
    """
    
    def __init__(self):
        """
        initialize with H0. units are eV
        """
        super(H0_lik,self).__init__()
        self.H0_dat = 2.46e-4 ##measured value from Reiss 2011
        self.sig_H0 = 8.0e-6           
        
    def __call__(self, H_0=1.5e33):
        """
        returns ln(likelihood) given H0 for the particular model.
        H should be given in eV
        """     
	#print "H0c2=",((H_0-self.H0_dat)*(1.0/self.sig_H0))**2
        lnl = 0.5*(((H_0-self.H0_dat)*(1.0/self.sig_H0))**2)
        return lnl
            
    

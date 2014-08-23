from cosmoslik import SlikPlugin
from numpy import dot

class bao236_lik(SlikPlugin):
    """
    get likelihood for bao236 given model predictions
    can use any units so long as they are consistent 
    """
    
    def __init__(self):
        """
        initialize with bao data for 2.36. can use any units so long as they are consistent
        """
        super(bao236_lik,self).__init__()
        self.c_over_Hrd = 9.0 ## c/(H(2.36)*r_d) measured value from BAO at z = 2.36
        self.sig_c_over_Hrd =  0.3## uncertainty  
        self.D_A_over_rd = 10.8## D_A/r_d measured value from BAO at z = 2.36
        self.sig_D_A_over_rd = 0.4
        
        
    def __call__(self, **kwargs):
        """
        returns ln(likelihood) given H(2.36) and D_A(2.36) for the particular model.
        can use any units so long as they are consistent.
        """     
        pred = { 'H_236':70.0,
                 'D_A_236' : 2000.0,
                 'r_d' : 149.28,
                 'c' : 2.9979e5
        }
        for key, value in kwargs.iteritems():
            pred[key] = value  
	#print "pred['D_A_236']'",pred['D_A_236'], "---measured=",self.D_A_236_dat
	#print "H236c2 =",((pred['H_236']-self.H_236_dat)*(1.0/self.sig_H_236))**2,"DA236c2 =", ((pred['D_A_236']-self.D_A_236_dat)*(1.0/self.sig_D_A_236))**2
        self.c_over_Hrd_theo = pred['c']/pred['H_236']/pred['r_d']
        self.D_A_over_rd_theo = pred['D_A_236']/pred['r_d']   	
        lnl = 0.5*(((self.c_over_Hrd_theo-self.c_over_Hrd)*(1.0/self.sig_c_over_Hrd))**2 + ((self.D_A_over_rd_theo-self.D_A_over_rd)*(1.0/self.sig_D_A_over_rd))**2)
        return lnl
            
    

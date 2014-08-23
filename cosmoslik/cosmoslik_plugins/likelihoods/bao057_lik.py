from cosmoslik import SlikPlugin
from numpy import array, dot

class bao057_lik(SlikPlugin):
    """
    get likelihood for bao057 given model predictions
    """
    
    def __init__(self):
        """
        initialize with bao data for 0.57
        """
        super(bao057_lik,self).__init__()
        self.DV_over_rd =  13.88    ##UNITLESS measured from BAO at z = 0.57 using rd from CAMB rd = 149.28Mpc
        self.F_AP = 0.683 ##UNITLESS measured from BAO at z = 0.57
        self.Cinverse = array([[31.0326, 77.6828], [77.6828, 2687.372]])  ##covariance matrix for above two points
        self.rd_camb =  149.28# in Mpc, in eV it's 2.3343e31 
        
        
    def __call__(self, **kwargs):
        """
        returns ln(likelihood) given H(0.57) D_A(0.57) and r_d for the particular model. 
        Units can be anything so long as they are consistent.
        """
        pred = { 'H_057':70.0,
                 'D_A_057' : 2000.0,
                 'r_d' : 149.28,
                 'c' :2.9979e5
        }
        for key, value in kwargs.iteritems():
            pred[key] = value     
        self.DV_over_rd_theo = ((pred['c']*(0.57)*((1.0+0.57)**2)*(pred['D_A_057']**2)*(1.0/pred['H_057']))**(1.0/3.0))/pred['r_d']
        self.F_AP_theo = ((1.0+0.57)*pred['D_A_057']*pred['H_057'])/pred['c']
        DV_diff = self.DV_over_rd_theo - self.DV_over_rd
        F_AP_diff = self.F_AP_theo - self.F_AP
	#print "057c2=",dot([DV_diff,F_AP_diff],dot(self.Cinverse,[DV_diff,F_AP_diff]))
        lnl = 0.5*dot([DV_diff,F_AP_diff],dot(self.Cinverse,[DV_diff,F_AP_diff]))
        return lnl
            
    

from cosmoslik import SlikPlugin
from numpy import array, dot

class bao032_lik(SlikPlugin):
    """
    get likelihood for bao at z=0.10 given model predictions
    """
    
    def __init__(self):
        """
        initialize with bao data for 0.32
        """
        super(bao032_lik,self).__init__()
        self.rd_camb =  149.28# in Mpc, in eV it's 2.3343e31 
        self.DV_over_rd = 8.467
        self.sig_DV_over_rd = 0.167
        
    def __call__(self, **kwargs):
        """
        returns ln(likelihood) given H(0.32) D_A(0.32) and r_d for the particular model. 
        Units can be anything so long as they are consistent.
        """
        pred = { 'D_A_032' : 2000.0,
                 'r_d' : 149.28,
                 'H_032' : 70,
                 'c' :2.9979e5
        }
        for key, value in kwargs.iteritems():
            pred[key] = value     
        self.DV_over_rd_theo = ((pred['c']*(0.32)*((1.0+0.32)**2)*(pred['D_A_032']**2)*(1.0/pred['H_032']))**(1.0/3.0))/pred['r_d']
        DV_diff = self.DV_over_rd_theo - self.DV_over_rd
	#print "057c2=",dot([DV_diff,F_AP_diff],dot(self.Cinverse,[DV_diff,F_AP_diff]))
        lnl = 0.5*DV_diff/self.sig_DV_over_rd
        return lnl
            
    

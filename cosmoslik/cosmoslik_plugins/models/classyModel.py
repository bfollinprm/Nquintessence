from cosmoslik import SlikPlugin
from classy import Class
from numpy import arange, zeros, exp, pi

class classyModel(SlikPlugin):
    """
    Runs the classy code (for LCDM for now) and saves the theory 
    predictions needed for H0_lik, bao057_lik, bao236_lik, and clik_like. 
    """
    def __init__(self):
        super(classyModel,self).__init__()
        


    def __call__(self,
                 **kwargs):
        names = {
            'A_s' : 'As',
            'n_s' : 'ns',
            'omega_b' : 'ombh2',
            'omega_cdm' : 'omch2',
            'tau_reio' : 'tau',
            'H0' : 'H0',
            'N_ur' : 'massive_neutrinos'
        }
        self.pars = {
            'A_s' : 2e-9,
            'n_s' : 0.6,
            'omega_b' : 0.05,
            'omega_cdm' : 0.2,
            'tau_reio' : 0.1,
            'H0' : 67.0,
            'N_ur' : 3.6,
        }

        for key, value in names.iteritems():
            for key2, value2 in kwargs.iteritems():
                if (names[key]==key2):
                    self.pars[key] = value2


 

        self.model= Class()
        self.model.set(**self.pars)
        self.model.set(output = 'tCl, lCl, pCl')
        self.model.set({'lensing':'yes', 'l_max_scalars':4000})
        print "before1 compute"
        self.model.compute()
        print "before bao"
        self.bao057 = { 'H_057' : self.model.Hubble(0.57),
                        'D_A_057' : self.model.angular_distance(0.57),
                        'c' : 1.0,
                        'r_d' : (self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']
                      }
        self.bao236 = { 'H_236' : self.model.Hubble(2.36),
                        'D_A_236' : self.model.angular_distance(2.36),
                        'c' : 1.0,
                        'r_d' : (self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']
                      }
        self.bao032 = { 'D_A_032' : self.model.angular_distance(0.32),
                        'H_032' : self.model.Hubble(0.32),
                        'c' : 1.0,
                        'r_d' : (self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']
                      }
        
        self.H0 = self.model.Hubble(0.0)
	blank = zeros(4001)
        ell = arange(4001)
        #make this in the right units!
        #l_max_scalars should be adjustable??   
        self.cmb = { 'cl_TT' : (self.model.lensed_cl(4000)['tt'])*ell*(ell+1)*((2.726e6)**2)/2/pi#,
#                     'cl_EE' : (self.model.lensed_cl(4000)['ee'])*ell*(ell+1)*((2.726e6)**2),
#                     'cl_BB' : (self.model.lensed_cl(4000)['bb'])*ell*(ell+1)*((2.726e6)**2),
#                     'cl_TE' : (self.model.lensed_cl(4000)['te'])*ell*(ell+1)*((2.726e6)**2),
#                     'cl_PP' : (self.model.lensed_cl(4000)['pp'])*ell*(ell+1)*((2.726e6)**2), 
#                     'cl_TP' : (self.model.lensed_cl(4000)['tp'])*ell*(ell+1)*((2.726e6)**2),
#                     'cl_EB' : blank,
#                     'cl_TB' : blank
                   }
        return

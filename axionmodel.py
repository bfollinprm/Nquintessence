from cosmoslik import param_shortcut, lsum, get_plugin, SlikDict, SlikPlugin, Slik
from numpy import identity, exp, inf, arange, hstack, loadtxt, zeros, ones
import sys
        
param = param_shortcut('start','scale')


class main(SlikPlugin):
    def __init__(self):
        super(SlikPlugin,self).__init__()
        
        self.cosmo = get_plugin('models.cosmology')(
                 omch2 = param(0.12,min=0),
                 theta = param(0.0104),  
                 ombh2 = param(0.0225, min=0),
                 tau = param(0.09,min=0,gaussian_prior=(.085,.015)),
                 logA = param(3.17805383035),
                 ns = param(0.96),
                 m_chi_over_H0 = param(2.0,scale=0.001,min=0.0),
                 chi0 = param(0.2,scale=0.001,min=0.0) 
                 
        )
 
        self.get_predictions = get_plugin('models.classyModel')()
 
        self.bao057 = get_plugin('likelihoods.bao057_lik')()
 
        self.bao236 = get_plugin('likelihoods.bao236_lik')()
       
        self.bao032 = get_plugin('likelihoods.bao032_lik')()
 
        self.H0 = get_plugin('likelihoods.H0_lik')()
 	
        self.camspec = get_plugin('likelihoods.clik_like')(
            clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/CAMspec_v6.2TN_2013_02_26_dist.clik',
            A_ps_100=param(150,min=0),
            A_ps_143=param(60,min=0),
            A_ps_217=param(60,min=0),
            A_cib_143=param(10,min=0),
            A_cib_217=param(40,min=0),
            A_sz=param(5,scale=1,range=(0,20)),
            r_ps=param(0.7,range=(0,1)),
            r_cib=param(0.7,range=(0,1)),
            n_Dl_cib=param(0.8,scale=0.2,gaussian_prior=(0.8,0.2)),
            cal_100=param(1,scale=0.001),
            cal_217=param(1,scale=0.001),
            xi_sz_cib=param(0.5,range=(-1,1),scale=0.2),
            A_ksz=param(1,range=(0,5)),
            Bm_1_1=param(0,gaussian_prior=(0,1),scale=1)
        )
        
        self.hubble_theta = get_plugin('models.hubble_theta')()
#        print self.cosmo
#	self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        
        self.priors = get_plugin('likelihoods.priors')(self)


        self.sampler = get_plugin('samplers.metropolis_hastings')(
              self,
              num_samples=10000,
              output_file=sys.argv[1]+'.chain',
              proposal_cov='proposal.covmat',
              proposal_scale=1
        )
 
 
    def __call__(self):
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        self.cosmo['As'] = 1e-10*exp(self.cosmo['logA'])

        self.predictions = self.get_predictions(**self.cosmo)
	#print "b57c2 =", 2*self.bao057(**self.get_predictions.bao057), "-- b236c2 =", 2*self.bao236(**self.get_predictions.bao236), "-- H0chi2", 2*self.H0(H_0=self.get_predictions.H0)
       
        return lsum(lambda: self.priors(self),
                    lambda: self.bao057(**self.get_predictions.bao057),
                    lambda: self.bao236(**self.get_predictions.bao236),
                    lambda: self.bao032(**self.get_predictions.bao032),
                    lambda: self.H0(H_0=self.get_predictions.H0),
                    lambda:  self.camspec(self.get_predictions.cmb))


if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass

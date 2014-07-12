

import sys
sys.path.append('../julia_code/cosmoslik')
from numpy import identity, exp, inf, hstack, load, zeros, arange, pi, sqrt, loadtxt, save, array
from numpy.random import normal
from get_background import get_hubble_rate

def lnl(cosmo):


	BAO = {
						'redshifts':[0.35, 0.57, 2.36], 
						'Da_rs':[6.875, 9.191, 10.8], 
						'Da_err': [0.246, 0.294, 0.4], 
						'H_rs':[12895, 14231, 3e5/9], 
						'H_err': [1070,1195,1e30]
					}


	params = dict()
	for key in cosmo:
		params[key] = cosmo[key]

	params['z0'] = 10
	params['step_density'] = 5e2


	angular_diameter, hubble = get_hubble_rate(**params)

		


	### Likelihoods
	chi2 = 0
	##CMB prior
	chi2 += (hubble(0)*2.998e5 - 73.8)**2/2.4**2
	for i,redshift in enumerate(BAO['redshifts']):
		chi2+= (BAO['Da_rs'][i] - angular_diameter(redshift) / cosmo['r_s'])**2/BAO['Da_err'][i]**2
		chi2+= (BAO['H_rs'][i] - hubble(redshift) * cosmo['r_s'] * 2.998e5)**2/BAO['H_err'][i]**2
		chi2+= (9.0 - 1/hubble(2.36)/cosmo['r_s'])**2/.3**2
	chi2 += (cosmo['omega_m'] - 0.1423)**2/0.0029**2


	lnl = chi2/2.0
	print 'lnl is %5f'%lnl
	return lnl




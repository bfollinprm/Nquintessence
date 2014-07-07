from classy import Class
from cosmoslik_plugins.models import cosmology
from evolve_background import get_hubble_rate




def get_model(**CMB_parameters):
	'''
	Required arguments for this to run successfully:
	'A_s',
	'H0',
	'N_ur',
	'omega_b',
	'omega_cdm',
	'tau_reio',
	'n_s',
	'step_density',
	'z0', (value of z where the inital conditions for the chi field are set)
	'm_chi_over_H0'
	'chi0'
	'chidot0'
	'''
	model = Class()
	valid_params = [  ###Can find more in class explanatory.ini file.  These are the only ones we need, though.
					'A_s',
					'H0',
					'N_ur',
					'omega_b',
					'omega_cdm',
					'tau_reio',
					'n_s'
					]
	params = dict()
	for key, value in CMB_parameters.iteritems():
		if key in valid params:
			params[key] = value

	model.set({'output': 'lCl, tCl', 'lensing':'yes', 'l_max_scalars':4000})
	model.set(params)

	### Compute Background quantities

	model.compute(lvl = ['background'])

	### Sub out the hubble module for one that includes the background effects of the chi field.
	oldHubble = model.Hubble
	output, newHubble = get_hubble_rate(CMB_parameters)

	def Hubble(z):
		try:
			newHubble(z)
		else:
			oldHubble(z) ##Revert to old history for z > z0

	cmb.Hubble = Hubble

	###this changes what H0 is
	CMB_parameters['H0'] = hubble_rate(0) * 1e5

	cmb.compute()
	_ell = arange(3001)
	clTT = cmb.lensed_cls(3000)['tt'] * _ell * (_ell + 1) * (1.0e6*cmb.Tcmb())**2


	

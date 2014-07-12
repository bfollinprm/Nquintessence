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
	tmp_model = Class()
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

	tmp_model.set({'output': 'lCl, tCl', 'lensing':'yes', 'l_max_scalars':4000})
	tmp_model.set(params)

	### Compute Background quantities

	tmp_model.compute(lvl = ['background'])

	### Sub out the hubble module for one that includes the background effects of the chi field.
	oldHubble = tmp_model.Hubble
	oldDa = tmp_model.angular_diameter_distance
	del tmp_model
	output, newHubble = get_hubble_rate(CMB_parameters)
	redshift = linspace(0, output['z'].max(), output['z'].size)
	newDa = cumsum(array([quad(lambda z: 1/newHubble(z) for i in arange(redshift.size - 1)]))* 1.0/(1+redshift)

	class NewClass(Class):
		def Hubble(self, z):
			try:
				newHubble(z)
			else:
				oldHubble(z) ##Revert to old history for z > z0
		def angular_diameter_distance(self, z):
			try:
				newDa(z)
			else:
				oldDa(z)

	model = NewClass()
	model.set({'output': 'lCl, tCl', 'lensing':'yes', 'l_max_scalars':4000})
	model.set(params)
	model.compute()
	###this changes what H0 is
	CMB_parameters['H0_true'] = model.Hubble(0) * 1e5
	derived_parameters = ['age', 'z_reio', 'z_rec','rs_rec','100*theta_s', 'YHe', 'sigma8', 'rs_drag']
	for key in derived_parameters:
	CMB_parameters[key] = model.get_current_derived_parameters(derived_parameters)
	CMB_parameters['age'] = model.age()
	_ell = arange(3001)
	CMB_parameters['clTT'] = model.lensed_cls(3000)['tt'] * _ell * (_ell + 1) * (1.0e6*model.Tcmb())**2

	###### ADD MORE REQUIRED NUMBERS HERE, e.g. below ######
	CMB_parameters['H057'] = model.Hubble(0.57) * 1e5
	CMB_parameters['DA057'] = model.angular_diameter_distance(0.57) ## In Mpc

	return CMB_parameters

	

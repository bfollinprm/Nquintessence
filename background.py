try:
	from numpy import pi, sqrt, append, array, arange, cumsum, load
except:
	print 'cant find numpy'
try:
	from scipy.integrate import ode
	from scipy.interpolate import interp1d
except:
	print 'cant find scipy'


class background():
	def __init__(self, *args, **kwargs):
		##unpack the args tuple
		h = .7
		omega_m = .12
		omega_rad = 1e-5
		chi0 = .2
		mchi_over_H0 = 2
		self.c = 2.998e5 #km/s
		self.kb = 8.6173e-5 #eV/K
		_H0 = h * 1e2
		self.cosmo = {
						'H0_guess' : _H0/self.c,
						'Omega_m': omega_m/h**2,
						'Omega_rad': omega_rad/h**2,
						#'omega_ncdm': omega_ncdm/h**2,  ## To be implemented
						'chi0': chi0,
						'mchi_over_H0': mchi_over_H0 
					 }
		self.cosmo['Omega_L'] = 1 - self.cosmo['Omega_m'] - self.cosmo['Omega_rad']
		try:
			self.z0 = z0
		except:
			self.z0 = (self.cosmo['mchi_over_H0'])**(2.0/3.0)*21.54 -1 ## This should put us deep enough in matter domination that rho_chi << rho_dm

		try:
			self.step = step
		except:
			self.step = 0.005 ## per hubble time
	def hubble(self,z):
		if (self.cosmo['chi0'] == 0) or (z > self.z0):
			hubble = self.cosmo['H0_guess'] * sqrt(self.cosmo['Omega_m'] * (1+z)**3 + self.cosmo['Omega_rad'] * (1+z)**4 + self.cosmo['Omega_L'])
			return hubble
		else:
			if hasattr(self, 'hubble_func') == True:
				return self.hubble_func(z).tolist()
			else:
				##We need to do the integration
				self.hubble_func = self.calc_hubble_rate()
				return self.hubble_func(z).tolist()

	def calc_hubble_rate(self):
		y0 = array([self.z0, self.cosmo['chi0'], 0])
		r = self._integrator().set_initial_value(y0)
		zlist = [self.z0]
		self.H = self._hofz(y0)
		self.dH = self._dhofz(y0)
		hlist = [self.H]

		z = self.z0
		while r.successful() and r.y[0] > 0:
			r.integrate(r.t + self.step/max(self.H, 1.4*self.cosmo['mchi_over_H0']*self.cosmo['H0_guess']))  ## Sample at the fastest of two scales; hubble time or oscillation time w = sqrt(2m^2)
			self.H = self._hofz(r.y)
			self.dH = self._dhofz(r.y) 
			zlist.append(r.y[0])
			hlist.append(self.H)

		hubble_array = array(hlist)[::-1]
		redshift_array = array(zlist)[::-1]
		##Interpolate
		hubble_interpolated_func = interp1d(redshift_array, hubble_array, kind ='linear', copy = False)#, assume_sorted = True)

		return hubble_interpolated_func




	def _integrator(self):
		integrator = ode(self._chi_field_eom)#, jac = self._chi_field_jacobian)
		integrator.set_integrator('lsoda')#, with_jacobian = True)
		#integrator.set_f_params(self)
		#integrator.set_jac_params(self)
		return integrator


	def _hofz(self,x):
		H = self.cosmo['H0_guess'] * sqrt(
						self.cosmo['Omega_m'] * ((1+x[0]))**3  	
						+ self.cosmo['Omega_rad'] * ((1+x[0]))**4 	
						+ self.cosmo['Omega_L']
						+ self.cosmo['mchi_over_H0']**2 * x[1]**2/3.0 + x[2]**2/6.0/self.cosmo['H0_guess']**2 
						)
		return H

	def _dhofz(self, x):
		dH = self.cosmo['H0_guess'] / sqrt(
						self.cosmo['Omega_m'] * ((1+x[0]))**3  	
						+ self.cosmo['Omega_rad'] * ((1+x[0]))**4 	
						+ self.cosmo['Omega_L']
						+ self.cosmo['mchi_over_H0']**2 * x[1]**2/3.0 + x[2]**2/6.0/self.cosmo['H0_guess']**2 
						)/2
		return dH


	def _chi_field_eom(self, t, x):
		dy =array([
				 -( 1 + x[0])*self.H,
				  x[2],
				 -( 2* (self.cosmo['mchi_over_H0'] * self.cosmo['H0_guess'])**2 * x[1] + 3 * self.H * x[2])
			])
		return dy

	def _chi_field_jacobian(self, t, x):
		jac = 0
		return jac

	def _check_for_hubble(self):

				##Check temporary numpy file to see if the cosmology agrees
		try:				
			cosmo_test = load('tmp/cosmology.npz')
			for key in cosmo_test.keys():
				if (cosmo_test[key] != self.cosmo[key]):
					return False
				return True
		except:
			return False


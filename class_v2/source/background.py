
from numpy import pi, sqrt, append, array, arange, cumsum, load
#from numpy import array
from scipy.integrate import ode, quad
from scipy.interpolate import interp1d



class background():
	def __init__(self, *args, **kwargs):
		##unpack the args tuple
		h = args[0]
		omega_m = args[1]
		omega_rad = args[2]
		chi0 = args[3]
		mchi_over_H0 = args[4]
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


	def density(self,*args):
		if args:
			try:
				z = args[0]
			except:
				z = args
		else:
			z = 0
		if (self.cosmo['chi0'] == 0) or (z > self.z0):
			density = self.cosmo['H0_guess']**2 * (self.cosmo['Omega_m'] * (1+z)**3 + self.cosmo['Omega_rad'] * (1+z)**4 + self.cosmo['Omega_L'])
			return float(density)
		else:
			if hasattr(self, 'density_interpolated_func') == True:
				density =  self.density_interpolated_func(z).tolist()
				return float(density)
			else:
				##We need to do the integration
				self.get_pressure_density_funcs()
				density = self.density_interpolated_func(z).tolist()
				return float(density)


	def angular_diameter_distance(self, z):
		Da = 1.0/(1+z) * quad(lambda x: 1.0/self.hubble(x), 0, z)[0]
		return Da

	def get_pressure_density_funcs(self):
		y0 = array([self.z0, self.cosmo['chi0'], 0])
		r = self._integrator().set_initial_value(y0)
		zlist = [self.z0]
		self.H = self._hofz(y0)		
		pressure_list = [self._Pofz(y0)]
		density_list = [self.H**2]
		z = self.z0
		while r.successful() and r.y[0] > 0:
			r.integrate(r.t + self.step/max(self.H, 1.4*self.cosmo['mchi_over_H0']*self.cosmo['H0_guess']))  ## Sample at the fastest of two scales; hubble time or oscillation time w = sqrt(2m^2)
			self.H = self._hofz(r.y)
			self.P = self._Pofz(r.y)
			self.rho = (self.H)**2
			zlist.append(r.y[0])
			pressure_list.append(self.P)
			density_list.append(self.rho)

		density_array = array(density_list)[::-1]
		pressure_array = array(pressure_list)[::-1]
		redshift_array = array(zlist)[::-1]
		##Interpolate
		self.pressure_interpolated_func = interp1d(redshift_array, pressure_array, kind ='linear', copy = False)#, assume_sorted = True)
		self.density_interpolated_func = interp1d(redshift_array, density_array, kind ='linear', copy = False)#, assume_sorted = True)

		return 0




	def _integrator(self):
		integrator = ode(self._chi_field_eom)
		integrator.set_integrator('lsoda')
		return integrator

	def _Pofz(self, x):
		Pressure = self.cosmo['H0_guess']**2 * (
												self.cosmo['Omega_rad']/3.0 * (1 + x[0])**4
												- self.cosmo['Omega_L']
												+ x[2]**2/6.0/self.cosmo['H0_guess']**2
												- self.cosmo['mchi_over_H0'] ** 2 *x[1]**2/3.0
												)
		return Pressure

	def _hofz(self,x):
		H = self.cosmo['H0_guess'] * sqrt(
						self.cosmo['Omega_m'] * ((1+x[0]))**3  	
						+ self.cosmo['Omega_rad'] * ((1+x[0]))**4 	
						+ self.cosmo['Omega_L']
						+ self.cosmo['mchi_over_H0']**2 * x[1]**2/3.0 + x[2]**2/6.0/self.cosmo['H0_guess']**2 
						)
		return H


	def _chi_field_eom(self, t, x):
		dy =array([
				 -( 1 + x[0])*self.H,
				  x[2],
				 -( 2* (self.cosmo['mchi_over_H0'] * self.cosmo['H0_guess'])**2 * x[1] + 3 * self.H * x[2])
			])
		return dy


	def pressure(self,*args):
		if args:
			try:
				z = args[0]
			except:
				z = args
		else:
			z = 0
		if (self.cosmo['chi0'] == 0) or (z > self.z0):
			pressure = self.cosmo['H0_guess']**2  * (self.cosmo['Omega_rad']/3.0 * (1+z)**4 - self.cosmo['Omega_L'])
			return float(pressure)
		else:
			if hasattr(self, 'pressure_interpolated_func') == True:
				pressure =  self.pressure_interpolated_func(z).tolist()
				return float(pressure)
			else:
				##We need to do the integration
				self.get_pressure_density_funcs()
				pressure = self.pressure_interpolated_func(z).tolist()
				return float(pressure)


def get_background_density(*args, **kwargs):
	background_class_instance = background(*args)
	density = background_class_instance.density
	return density

def get_background_pressure(*args, **kwargs):
	background_class_instance = background(*args)
	pressure = background_class_instance.pressure
	return pressure

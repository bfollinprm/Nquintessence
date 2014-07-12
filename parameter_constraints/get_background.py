from numpy import pi, sqrt, vstack, array, arange, cumsum, roll, insert, linspace
from numpy.polynomial import Polynomial, Legendre
from scipy.integrate import ode, quad, romberg
from scipy.interpolate import UnivariateSpline, interp1d
import itertools
from pylab import *

def get_hubble_rate(**kwargs):
		## define default cosmology initial conditions.  values depend on when you're setting intiial conditions.
	c = 2.998e5 #m/s
	k_b = 8.6173e-5 #eV/K
	params = {
			'H0' : 70,  ##What H0 would be if there were no chi field.
			'omega_m':0.14,
			'omega_rad': 0,
			'step_density': 1e2,
			'z0':30, ##redshift of your initial conditions (reheating)??  Shouldn't matter much as long as it's far into radiation era, and t(z0) << period of DM oscillator
			'm_chi_over_H0': 3,
			'chi0': 1,##value of chi field in units of Mp
			'chidot0':0  ## In units of Mp/Mpc
			}
	for key, value in kwargs.iteritems():
		params[key] = value

	## Set initial conditions
	initial_conditions = get_initial_conditions(params)
	del params

	initial_conditions['Vchi'] = Polynomial([0,0, (initial_conditions['m_chi_over_H0'] * initial_conditions['H0'])**2])

	output = evolve_background(initial_conditions)   ##Get the background
	hubble_func = interp1d(
							output['z'],
							output['H'],
							kind = 'linear',
							bounds_error = False,
							fill_value = output['H'].min()
							)


					# UnivariateSpline(
					# 				output['z'], 
					# 				output['H'], 
					# 			 	bbox = [0, initial_conditions['z0']], 
					# 			 	k = 1,
					# 			 	s = 0
					# 			  )
	redshift = linspace(0, initial_conditions['z'], 200)

	angular_diameter = cumsum(array([
	 								quad(
	 									lambda x: 1/hubble_func(x),redshift[i], redshift[i+1])[0]
	 								for i in arange(redshift.size - 1)]
	 								)) / (1+redshift[1:])
	

	angular_distance_func = interp1d(
									redshift,
									insert(angular_diameter, 0, 0),
									kind = 'linear'

									)

	rho_chi = interp1d(
						output['z'],
						initial_conditions['Vchi'](output['chi']) + output['chidot']**2/2,
						kind = 'linear',
						)

	# loglog(1/(1+output['z']), initial_conditions['Vchi'](output['chi'] + output['chidot']**2/2))
	# loglog(1/(1+output['z']), initial_conditions['rho_matter'] * ((1+output['z'])/(1+initial_conditions['z']))**3 )
	# loglog(1/(1+output['z']), initial_conditions['Lambda'] + output['z']*0)

	# show()
	
	return rho_chi, angular_distance_func, hubble_func

def get_initial_conditions(params):
	##first, find the Omega's to find Omega_L
	Omega_m = params['omega_m']*(1e4)/params['H0']**2
	Omega_ur = params['omega_rad'] *(1e4)/params['H0']**2
	Omega_L = 1 - Omega_m - Omega_ur
	initial_conditions = dict()
	##Now, turn them into physical densities at redshift z0 (units of Mp**2)
	###Pass all the other values into the inital_conditions dict
	for key, value in params.iteritems():
		initial_conditions[key] = value
	z = params['z0'] 
	initial_conditions['z'] = z
	H0 = params['H0'] / 2.998e5
	initial_conditions['H0'] = H0
	initial_conditions['rho_matter'] = Omega_m * (1+z)**3 * H0**2
	initial_conditions['rho_ur'] = Omega_ur * (1+z)**4 * H0 ** 2
	initial_conditions['Lambda'] = Omega_L * H0 ** 2

	return initial_conditions


def get_H_of_z(z, chi, chidot, params):

	H_squared = 1.0 *   (
						params['rho_matter'] * ((1+z)/(1+params['z']))**3  	##Baryons
						+ params['rho_ur'] * ((1+z)/(1+params['z']))**4 	##Neutrinos
						+ params['Lambda']
						+ potential_chi(chi,params) + chidot**2/2  ##Phi-field (DE)
						)
	return sqrt(H_squared)


def potential_chi(chi, global_vars):
	potential = global_vars['Vchi']
	return potential(chi)


def vprime_chi(chi, global_vars):
	vprime = global_vars['Vchi'].deriv()
	return vprime(chi)

def evolve_background(params):

	#y = array([
	#	z,
	#	phi,
	#	phidot,
	#	chi,
	#	chidot
	#	])###
	t0 = 0 
	chidot0 = 0

	##Recursiveley find initial field velocities consistent with slow-roll \ddot{\phi} ~ 0
	integrator = ode(equation_of_motion)
	integrator.set_integrator('dopri5')  ## Choose a method
	integrator.set_initial_value(array([		## Set the initial conditions
										params['z'], 
										params['chi0'], 
										params['chidot0']
											]), t0)
	H_local = get_H_of_z(integrator.y[0], integrator.y[1],integrator.y[2], params)
	cosmic_history = [[integrator.t], [H_local], integrator.y.tolist()]  ##array to save values in for later:  time, H, z, phi, phidot, chi, and chidot.
	cosmic_history = [num for elem in cosmic_history for num in elem] ##flatten, important for vstack

	#print cosmic_history
	integrator.set_f_params(params)  ##pass non-evolving parameters to function 'equation_of_motion'
	count = 0

	while integrator.successful() and integrator.y[0] >= 0:  ## Integrate
		count = count + 1
		integrator.integrate(integrator.t + 1.0/params['step_density']/H_local)  ## Time step between saves goes as 1/H(z), with 100 saves per hubble time.  Is this what we want???
		H_local = get_H_of_z(integrator.y[0], integrator.y[1],integrator.y[2], params)
		newrow = [[integrator.t], [H_local], integrator.y.tolist()]
		newrow = [num for elem in newrow for num in elem]

		cosmic_history = vstack([cosmic_history, newrow])  ##add a new row with updated values for parameters of interest.

	history = {} ###repackage history as a dictionary to make it nicer to look things up later.
	history['t'] = cosmic_history[:,0]
	history['H'] = cosmic_history[:,1]  
	history['z'] = cosmic_history[:,2]
	history['chi'] = cosmic_history[:,3]
	history['chidot'] = cosmic_history[:,4]  ##This should have units of Mp/Mpc.  I'm not positive--though I'm pretty sure--it does.

	return history

def equation_of_motion(t, y, params):
	H = get_H_of_z(y[0], y[1], y[2], params)

	dy =  array([
				-(1 + y[0]) * H,		##dz =  -(1 + z) * H
				y[2], 							 							##dchi = chidot
				-( 3 * H *y[2] + vprime_chi(y[1], params)),				##dchidot = -3 H chidot - V'(chi)
				])
	return dy


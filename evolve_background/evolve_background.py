from numpy import pi, sqrt, vstack, array, arange, cumsum, roll, insert
from numpy.polynomial import Polynomial, Legendre
from scipy.integrate import ode, quad, romberg
from scipy.interpolate import interp1d
import itertools
Mp = 2.435e27 ##eV
k_B = 8.617332478e-5 ##eV/K
H0 = 1.494e-24 ##70 km/s/Mpc in eV
rho_crit = 3*(H0*Mp)**2
z0 = 5
def run(**kwargs):
		## define default cosmology initial conditions.  values depend on when you're setting intiial conditions.
	m_chi = 50 * H0
	m_phi = 1.0/20.0 * H0
	global_vars = {
			'verbose': True,
			'step_density': 1e3,
			'backwards_evolve': False, ##Flag to say whether these are initial conditions today, or far into the past.
			'z0':300, ##redshift of your initial conditions (reheating)??  Shouldn't matter much as long as it's far into radiation era, and t(z0) << period of DM oscillator
			'z_end':0, ##Redshift where you end the integration.  
			'rho_b0':.02 * sqrt(3) * (Mp*H0)**2 * (1+z0)**3,##Baryons.  Some thought goes into what these should be.  What we know is that \Omega_b h^2 = 0.022, but we don't know what h = H_0/100 is in this model.
						##Should be enough information, though, to work it out in terms of the other numbers here...not quite seeing how right now though.
			'rho_c0':.28 * sqrt(3) * (Mp * H0)**2 * (1+z0)**3, ##Dark Matter (non oscillating).  Get this from \Omega_c h^2 = 0.120
			'rho_nu0':7.0/8.0 * 3.046 * (4.0/11.0)**(4/3) * pi**2 * (k_B * 2.73)**4/15.0 * (1+z0)**4, ##Neutrinos.  Get this from \Omega_{\nu} h^2 = \frac{7}{8} N_{\nu} \left(\frac{4}{11}\right)^{4/3} \Omega_{\gamma} h^2, where N_{\nu} = 3.046 including 1 loop.
			'rho_gamma0':pi**2 * (k_B * 2.73)**4/15.0 * (1+z0)**4,##This is just the energy density of a black body spectrum with temperature T_{CMB}.  Easy to figure out, but will depend on units, 
						 ##and where you're setting your initial conditions (when is reheating?).
			'V_phi':Polynomial([0, 0, m_phi**2/2, 0]), ##Can add more terms in the effective field theory here.  Polynomial expands in powers of phi.
			'V_chi':Polynomial([0, 0, m_chi**2/2.0, 0]),  ##Also can add here.   For true cosine potentials, switch to Legendre expansion
			'phi0': 0.6 * sqrt(3.0/2) * (Mp*H0)/(m_phi), ##value of phi field at initial conditions (reheating)?
			'chi0': 10 * sqrt(3.0/2) * (Mp*H0)/(m_chi)##value of chi field.  Code assumes initial velocities are given by slow roll; chidot0 = -V'/3H, etc.
			}
	for key, value in kwargs.iteritems():
		global_vars[key] = value



	output = evolve_background(global_vars)   ##Get the background
	#### Derived Parameters ####
	global_vars['H_of_z_interpolated'] = interp1d(output['z'][::-1], output['H'][::-1], bounds_error = True)
	if global_vars['verbose'] == True: print 'evaluating derived parameters'
	differential_DA = array([quad(one_over_H_of_z, a[0], a[1], args = (global_vars,))[0] for a in zip(output['z'][::-1], roll(output['z'][::-1],-1)[:-1])])
	output['D_A'] = insert(cumsum(differential_DA), 0, 0.0)[::-1]
	output['Vphi'] = array([potential_phi(b, global_vars) for b in output['phi']])
	output['Vchi'] = array([potential_chi(b, global_vars) for b in output['chi']])
	output['wchi'] = (-output['Vchi'] + output['chidot']**2/2)/(output['chidot']**2/2 + output['Vchi'])
	output['wphi'] = (-output['Vphi'] + output['phidot']**2/2)/(output['phidot']**2/2 + output['Vphi'])
	output['wtot'] = ((-output['Vchi'] + output['chidot']**2/2 -output['Vphi'] + output['phidot']**2/2) / 
							(output['chidot']**2/2 + output['Vchi'] + output['phidot']**2/2 + output['Vphi'] ))
	output['wdark'] = ((-output['Vchi'] + output['chidot']**2/2 -output['Vphi'] + output['phidot']**2/2) / 
							(output['chidot']**2/2 + output['Vchi'] + output['phidot']**2/2 + output['Vphi'] + global_vars['rho_c0'] * (1+output['z'])**3/(1+z0)**3))

	#print output['D_A']
	#output['D_A_star'] = array([quad(one_over_H_of_z, b, 1100, args = (global_vars,))[0] for b in output['z']])
	return output

def one_over_H_of_z(z, global_vars):
	H_of_z = global_vars['H_of_z_interpolated']
	return 1.0/H_of_z(z)


def get_H_of_z(z, phi, phidot, chi, chidot, global_vars):

	H_squared = 1.0/(3 * Mp**2) * (global_vars['rho_b0'] * ((1+z)/(1+global_vars['z0']))**3  	##Baryons
						+ global_vars['rho_c0'] * (1+z)**3 /(1+global_vars['z0'])**3   	##(Non-axion) Dark Matter
						+ global_vars['rho_nu0'] * (1+z)**4 /(1+global_vars['z0'])**4  	##Neutrinos
						+ global_vars['rho_gamma0'] * (1+z)**4 /(1+global_vars['z0'])**4  ##Photons
						+ potential_phi(phi,global_vars) + phidot**2/2  ##Phi-field (DE)
						+ potential_chi(chi,global_vars) + chidot**2/2) ##Chi-field (slipping axion)
	#print 'contributions to H^2 = %g from'%H_squared
	#cont = (1.0/Mp**2/3 * (global_vars['rho_b0'] * (1+z)**3 + global_vars['rho_c0'] * (1+z)**3))
	#print '...matter: %g eV'%cont
	#cont = (1.0/Mp**2/3 * (global_vars['rho_gamma0'] * (1+z)**4 + global_vars['rho_nu0'] * (1+z)**4))
	#print '...radiation: %g eV'%cont
	#cont = 1.0/Mp**2/3 *(potential_phi(phi,global_vars) + phidot**2/2)
	#print '...phi: %g eV'%cont
	#cont = 1.0/Mp**2/3 *(potential_chi(chi,global_vars) + chidot**2/2)
	#print '...chi: %g eV'%cont
	return sqrt(H_squared)


def potential_chi(chi, global_vars):
	potential = global_vars['V_chi']
	return potential(chi)

def potential_phi(phi, global_vars):
	potential = global_vars['V_phi']
	return potential(phi)


def vprime_phi(phi, global_vars):
	vprime = global_vars['V_phi'].deriv()
	return vprime(phi)

def vprime_chi(chi, global_vars):
	vprime = global_vars['V_chi'].deriv()
	return vprime(chi)

def evolve_background(global_vars):

	#y = array([
	#	z,
	#	phi,
	#	phidot,
	#	chi,
	#	chidot
	#	])###
	t0 = 0 

	##Recursiveley find initial field velocities consistent with slow-roll \ddot{\phi} ~ 0
	if global_vars['verbose'] == True: print 'finding consistent initial conditions'
	phidot0 = 0#get_initial_phi('phi', global_vars)
	chidot0 = 0#get_initial_phi('chi', global_vars)
	if global_vars['verbose'] == True: 
		print 'phidot is initially %g'%phidot0
		print 'chidot is initially %g'%chidot0

	integrator = ode(equation_of_motion)
	integrator.set_integrator('dopri5')  ## Choose a method
	integrator.set_initial_value(array([		## Set the initial conditions
										global_vars['z0'], 
										global_vars['phi0'], 
										phidot0,
										global_vars['chi0'], 
										chidot0
											]), t0)
	H_local = get_H_of_z(integrator.y[0], integrator.y[1],integrator.y[2],integrator.y[3],integrator.y[4], global_vars)
	cosmic_history = [[integrator.t], [H_local], integrator.y.tolist()]  ##array to save values in for later:  time, H, z, phi, phidot, chi, and chidot.
	cosmic_history = [num for elem in cosmic_history for num in elem] ##flatten, important for vstack

	#print cosmic_history
	integrator.set_f_params(global_vars)  ##pass non-evolving parameters to function 'equation_of_motion'
	count = 0
	if global_vars['backwards_evolve'] ==True:  
		while integrator.successful() and integrator.y[0] < global_vars['z_end']:  ## Integrate
			count = count + 1
			integrator.integrate(integrator.t + 1.0/global_vars['step_density']/H_local)  ## Time step between saves goes as 1/H(z)
			H_local = get_H_of_z(integrator.y[0], integrator.y[1],integrator.y[2],integrator.y[3],integrator.y[4], global_vars
)
			newrow = [[integrator.t], [H_local], integrator.y.tolist()]
			newrow = [num for elem in newrow for num in elem]
			if (global_vars['verbose'] == True) and (count % global_vars['step_density'] == 0):
				print 'evolved back to z = %g'%integrator.y[0]
				print 'H(z) is %g'%H_local
			cosmic_history = vstack([cosmic_history, newrow])  ##add a new row with updated values for parameters of interest.

	else:
		while integrator.successful() and integrator.y[0] > global_vars['z_end']:  ## Integrate
			count = count + 1
			integrator.integrate(integrator.t + 1.0/global_vars['step_density']/H_local)  ## Time step between saves goes as 1/H(z), with 100 saves per hubble time.  Is this what we want???
			H_local = get_H_of_z(integrator.y[0], integrator.y[1],integrator.y[2],integrator.y[3],integrator.y[4], global_vars
)
			newrow = [[integrator.t], [H_local], integrator.y.tolist()]
			newrow = [num for elem in newrow for num in elem]
			if (global_vars['verbose'] == True) and (count % global_vars['setp_density'] == 0):
				print 'evolved forwards to z = %g'%integrator.y[0]
				print 'H(z) is %g'%H_local
			cosmic_history = vstack([cosmic_history, newrow])  ##add a new row with updated values for parameters of interest.

	history = {} ###repackage history as a dictionary to make it nicer to look things up later.
	history['t'] = -cosmic_history[:,0] ## This is probably a useless variable.  I don't know its units, for sure.
	history['H'] = cosmic_history[:,1]  
	history['z'] = cosmic_history[:,2]
	history['phi'] = cosmic_history[:,3]
	history['phidot'] = cosmic_history[:,4]
	history['chi'] = cosmic_history[:,5]
	history['chidot'] = cosmic_history[:,6]

	return history

def equation_of_motion(t, y, global_vars):
	H = get_H_of_z(y[0], y[1], y[2], y[3], y[4], global_vars)

	dy =  (2 * (global_vars['backwards_evolve']) - 1) * array([
				(1 + y[0]) * H,		##dz = (1 + z) * H
				-y[2], 							 							##dphi = phidot
				( 3 * H *y[2] + vprime_phi(y[1], global_vars)),				##dphidot = -3 H phidot - V'(phi)
				-y[4], 							 							##dchi = chidot
				( 3 * H *y[4] + vprime_chi(y[3], global_vars))				##dchidot = -3 H chidot - V'(chi)
				])
	return dy


def get_initial_phi(field, global_vars):
	guess = 1e50
	guess_new = 0
	while abs(guess - guess_new) > 1e-10:
		guess = guess_new
		if field == 'phi': guess_new = (2 * (global_vars['backwards_evolve']) - 1) * vprime_phi(global_vars['phi0'], global_vars)/(3 * get_H_of_z(0, global_vars['phi0'],0,global_vars['chi0'],0, global_vars))
		if field == 'chi': guess_new = (2 * (global_vars['backwards_evolve']) - 1) * vprime_chi(global_vars['chi0'], global_vars)/(3 * get_H_of_z(0, global_vars['phi0'],0,global_vars['chi0'],0, global_vars))
	return guess_new
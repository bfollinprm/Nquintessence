from pylab import *
from scipy.interpolate import interp1d

import evolve_background
Mp = 2.435e27 ##eV
k_B = 8.617332478e-5 ##eV/K
H0 = 1.494e-24 ##70 km/s/Mpc in eV
rho_crit = 3*(H0*Mp)**2
m_phi = 1.0/20.0 * H0
output = evolve_background.run()
nochi = evolve_background.run(**{'chi0':0, 'phi0': 0.7 * sqrt(3.0/2) * (Mp*H0)/(m_phi), 'rho_c0':.28 * sqrt(3) * (Mp * H0)**2})


print 'plotting'

## Plot w(z)
wchi = output['wchi']#(-output['Vchi'] + output['chidot']**2/2)/(output['chidot']**2/2 + output['Vchi'])
wphi = output['wphi']#(-output['Vphi'] + output['phidot']**2/2)/(output['phidot']**2/2 + output['Vphi'])
wtot = output['wtot']
z = output['z']
figure()
loglog(z, wchi+1, label = '$1 + w_{\chi}$')
loglog(z, wphi+1, label = '$1 + w_{\phi}$')
loglog(z, wtot+1, label = '$1 + w_{tot}$')
loglog(z, z*0+1, linestyle='dashed', label = '$1 + w_{dm}$')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('w_of_z.png')


## Plot V(z)

Vchi = output['Vchi']
Vphi = output['Vphi']
MH = 2.435e27 * output['H']

figure()
loglog(z, Vchi, label = '$V_{\chi}$')
loglog(z, Vphi, label = '$V_{\phi}$')
loglog(z, MH**2, label = '$(M_{pl}H(z))^2$')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('v_of_phi.png')



## Plot \dot{\phi}^2/2
phidot = output['phidot']
chidot = output['chidot']

figure()
loglog(z, chidot**2/2.0, label = '$\dot{\chi}^2/2$')
loglog(z, phidot**2/2.0, label = '$\dot{\phi}^2/2$')
loglog(z, MH**2, label = '$(M_{pl}H(z))^2$')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('kinetic.png')


## Plot field values

phi = abs(output['phi'])
chi = abs(output['chi'])

figure()
loglog(z, chi, label = '${\chi}$')
loglog(z, phi, label = '${\phi}$')
#loglog(z, output['H'], label = '$H(z)$')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('fields.png')


## Plot D_A(z)

D_A = output['D_A']
D_A2 = nochi['D_A']
z2 = nochi['z']
D_A_func = interp1d(z[1:].tolist(),D_A.tolist(),bounds_error = False, fill_value = 0)
D_A2_func = interp1d(z2[1:].tolist(),D_A2.tolist(),bounds_error = False, fill_value = 0)
z_reg = linspace(0,5,1000)
D_A = array([D_A_func(a) for a in z_reg])
D_A2 = array([D_A2_func(a) for a in z_reg])
figure()
plot(z_reg,D_A/D_A2, label = '$D_A/D_A^{F}$')
#loglog(z_reg, D_A2, label = 'time')
#loglog(z, 1.0/output['H'], label = 'dD_A/dz')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('D_A.png')


## Plot H(z)
figure()
H_func = interp1d(z.tolist(), output['H'].tolist())
H2_func = interp1d(z2.tolist(), nochi['H'].tolist())
H = array([H_func(a) for a in z_reg])
H2 = array([H2_func(a) for a in z_reg])
plot(z_reg, H/H2, label = '$H(z)/H(z)^F$')
legend()
xlim(z.min(), z.max() * 1.2)
xlabel('z', fontsize = 20)
savefig('Hubble.png')






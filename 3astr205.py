
"""
Need to find at least velocities in orbit, masses, radii,
type of orbit, centre of mass velocity, and semi major axes.

"""


import numpy as np
import scipy.optimize as optim
import matplotlib.pyplot as plt

c = 3*(10**8)  # m/s
p = 700816.5792  # in seconds
G = 6.67*(10**-11)  # gravitational constant in kg, km, s
m_sun = 1.989*(10**30)
s_b_const = 5.67*(10**-8.0)


def open_file(filen):
    return eval(open(filen).read())

# read in the data and make it into arrays from json format

veldata = open_file("vel.json")
v_s = np.asarray(veldata["v_s"])
v_p = np.asarray(veldata["v_p"])
days = veldata["days"]
phase = np.asarray(veldata["phase"])
photdata = open_file("phot.json")
v_mag = photdata["vmag"]
phase_phot = np.asarray(photdata["phase"])

# these are estimates of the parameters on the sit in order to bring the plots
# down to be centered around zero
# these get used in the sin function for the
guess_freq = 1/(2 * np.pi)
guess_amplitude = 60.0
guess_phase = 0.0
guess_offset = 80.0
p_0 = [guess_freq, guess_amplitude, guess_phase, guess_offset]

# this is the function and plotting for the curve fit
# the fitter is a sin function, and we use the parameters of the optimized function
# this will give us the plot later
# also tells us that e = o since the isn function fits nearly perfectly


def sin(x, freq, amplitude, phase, offset):
    return np.sin(x * freq + phase) * amplitude + offset


fit_1 = optim.curve_fit(sin, phase, v_s, p_0)
fit_2 = optim.curve_fit(sin, phase, v_p, p_0)
print fit_1
print fit_2
print(fit_1[0][1])
print(fit_2[0][1])

v_max_1 = fit_1[0][1]
v_max_2 = fit_2[0][1]

x_val = np.linspace(-0.5, 0.5, 1000)  # the linspace creates the x_val for the curve fit
plt.scatter(phase, v_s, color='#ff0000')
plt.scatter(phase, v_p)
plt.plot(x_val, sin(x_val, 6.27977076e+00, 6.25933287e+01, -4.17810864e-03, 7.88068450e+01))
plt.plot(x_val, sin(x_val, -6.28907633e+00, 5.94690758e+01, -3.96979814e-03, 7.87244524e+01))
plt.xlabel("Phase")
plt.ylabel("Velocities in km/s")
plt.title("Velocity Curve")
plt.show()

vel_sys = fit_1[0][3]  # velocity of the centre of mass of the system
print vel_sys

# now we use the ratios of semi major axis and kepler's 3rd to find the masses


def semi_major(v_max):
    return (v_max*p)/(2*np.pi)

semi_maj_1 = semi_major(fit_1[0][1])
semi_maj_2 = semi_major(fit_2[0][1])
a_ratio = semi_maj_1/semi_maj_2
print semi_maj_1
print semi_maj_2

mtot = p*((v_max_1 * 1000 + v_max_2 * 1000)**3) / (2*np.pi*G)
print mtot

m_1 = mtot/(1+a_ratio)
m_2 = mtot - m_1
print m_1
print m_2

# next we plot the photometry and the phase
# we can then use the difference in the minors to find the radii

plt.scatter(phase_phot, v_mag)
plt.gca().invert_yaxis()
plt.ylabel("Magnitude in V band")
plt.xlabel("Phase")
plt.title("Magnitudes in the V band VS Phase")
plt.show()

# assume that secondary minimum is caused by the primary passing in front, since it it less massive and
# therefore most likely likely to be less bright as they are both main sequence stars

total_eclipse_time = p*0.035
total_rad = (((v_max_2 + v_max_1)*total_eclipse_time)/2)
eclip_time_sec = p*0.014
sec_rad = (((v_max_2 + v_max_1)*eclip_time_sec)/2)
prim_rad = total_rad - sec_rad
print total_rad
print sec_rad
print prim_rad

# now that we have masses and radii, we can find the density

density_1 = m_1/((4/3)*np.pi*(prim_rad*1000)**3)
density_2 = m_2/((4/3)*np.pi*(sec_rad*1000)**3)
print density_1
print density_2

# using the mass luminosity relationship, we are able to find the luminosities
# as a fraction of solar luminosities
# we can the use the luminosities to find the temperatures of the stars
# and then we can use all we know about them (mass, radius, temp) to classify them


lum_1 = (m_1/m_sun)**3.5
lum_2 = (m_2/m_sun)**3.5
print lum_1
print lum_2

# this is the equation relating luminosity and temp rearranged for temp

temp_1 = ((lum_1*3.828*(10**26))/(4*np.pi*((prim_rad*1000)**2)*s_b_const))**0.25
temp_2 = (lum_2*(3.828*(10**26))/(4*np.pi*((sec_rad*1000)**2)*s_b_const))**0.25
print temp_1
print temp_2

# Definitions of physical constants and a few atomic-related quantitites

import numpy as np

Gamma583 = 2*np.pi * 180e3 # Gamma of the 583 nm transition (in Hz)
Gamma841 = 2*np.pi * 8e3 # Gamma of the 841 nm transition (in Hz)
Gamma626 = 2*np.pi * 135e3 # Gamma of the 626 nm transition (Dy) (in Hz)

lambd583 = 583e-9
lambd841 = 841e-9
lambd626 = 626e-9
lambd631 = 631e-9

amu = 1.66e-27 # Atomic ma87ss unit --> kg
m = 166 * amu # Erbium166 isotope mass in kg
# m = 87 * amu # Rb
kB = 1.38e-23 # Boltzmann constant, m^2 kg s^-2 K^-1
e0 = 8.854e-12 # electric permittivity, C^2 N^-1 m^-2
c = 299792458 # speed of light, m/s
g = 9.81 # gravitational acceleration, m/s^2
hbar = 1.05457e-34 # reduced Planck's constant
lambd_trap = 488e-9

electronCharge = 1.60217663e-19
a0 = 5.29e-11 
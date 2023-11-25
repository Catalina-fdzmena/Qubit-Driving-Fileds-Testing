# -*- coding: utf-8 -*-
"""
Research: Qbit control for Quantum Computing 

Families of Exactly Solvable One-Qubit Time-dependent Hamiltonians for a given function Ω(t),
which is chosen in such a way that the solutions to the parametric oscillator equation are
known. 

This algorithm discuss the case for which [Ln R_0]'= x + iδ, λ = δ and R_0 = -ig.
Ω(t) = constant.

The Case studied in this algorithm is Ω(t)= 0

Made on: 20/10/2023 by: Andrea Catalina Fernandez Mena
Modified on: 12/11/2023 by: Andrea Catalina Fernandez Mena & Félix Abdiel Rodriguez Rios.

"""
# Libraries
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.optimize import minimize

# Test parameters for which the periodicity condition holds

#g
g1 = 0.8

#δ
delta1 = 3

#kappa
kappa1 = 1

#Delta
Delta1 = 10

# Test parameters for which the periodicity condition does not hold

Omega0 = 1

#Omega1

Omega1_1 = Delta1/2  

# Re[[Ln R_0]']
x1 = float(input("Please enter an integer for the real part of the logarithmic derivative of R_0:"))

# Calculation function of η(t) for a)
def eta(t_values):
    dt = 0.001
    integral = np.zeros_like(t_values)

    for i, t in enumerate(t_values):
        s_values = np.linspace(0, t, int(t / dt))
        integrand = 1 / ((np.cos(Omega1_1 * s_values)) ** 2 + (kappa1 ** 2) * (np.sin(Omega1_1 * s_values))** 2)
        integral[i] = np.trapz(integrand, dx=dt)

    return integral

# R(t) function for a)
def R(t_values):
    eta_t = eta(t_values)
    return ((-1j * g1)/((Omega0**2)*(t_values**2)-(x1/2)*t_values+1))*np.exp(((x1 + 1j * delta1)/Omega0)*np.arctan(Omega0*t_values))

# P(t) Population inversion
def P(time):
    return ((-(x1**2+2*1j*x1*delta1-delta1**2)/4-g1**2)*(time**2)+1)/((Omega0**2)*(time**2)-(x1/2)*time+1)
            
# Time values
t_values = np.linspace(0, 10, 1000)
time = np.linspace(0,10*np.pi, 1000)

# Calculate ReR(t) and ImR(t)
ReR = np.real(R(t_values))
ImR = np.imag(R(t_values))

# Calculate R(t) 
result_values = R(t_values)

# Calculate P(t)
population_values = P(time)


# Graph for real vs imaginary parts of result_values
plt.subplot(121)
plt.plot(np.real(result_values), np.imag(result_values), color='Purple')
plt.title('Driving Field R(t)')
plt.xlabel('Re R(t)')
plt.ylabel('Im R(t)')

# Population graph
plt.subplot(122)
plt.plot(time, population_values, color="pink")
plt.title('Population Inversion vs Time')
plt.xlabel('Time')
plt.ylabel('P(t)')
plt.xticks([0, 5 * math.pi, 10 * math.pi],
           ['0', '5π', '10π'])

plt.tight_layout()  # Adjust layout to prevent overlapping

plt.show()
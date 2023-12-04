# -*- coding: utf-8 -*-
"""
Research: Qbit control for Quantum Computing 

Families of Exactly Solvable One-Qubit Time-dependent Hamiltonians for a given function Ω(t),
which is chosen in such a way that the solutions to the parametric oscillator equation are
known. 

This algorithm discuss the case for which [Ln R_0]'= x + iδ, λ = δ and R_0 = -ig.
Ω(t) = constant.

The Case studied in this algorithm is Ω(t)= Ω_1

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
g1 = math.sqrt(np.pi)
#g1 = math.sqrt(5)  
g2 = math.sqrt(160)

#δ
delta1 = 10
#delta1 = 4
delta2 = 6
delta3 = 2*math.sqrt(math.pi)

#kappa
kappa1 = 3
#kappa1 = 0.6
kappa2 = 3.1
kappa3 = 0.8
kappa4 = 2.5

#Delta
Delta1 = 6
#Delta1 = 10
Delta2 = 1.9
Delta3 = 32.5
Delta4 = 10.4

# Test parameters for which the periodicity condition does not hold

Omega0 = 1

#g
g3 = math.sqrt(math.pi)

#δ
delta3 = 2*math.sqrt(math.pi)

#kappa
kappa5 = 0.6
kappa6 = 2.5

#Delta
Delta5 = (10*kappa5*math.pi)/Omega0
Delta6 = (10*kappa6*math.pi)/Omega0

#Omega1

Omega1_1 = Delta1/2  
#Omega1_2 = Delta2/2
#Omega1_3 = Delta3/2
#Omega1_4 = Delta4/2  
#Omega1_5 = Delta5/2
#Omega1_6 = Delta6/2

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

# Calculation function of η(t) for  b)
# def eta(t_values):
#     dt = 0.001
#     integral = np.zeros_like(t_values)

#     for i, t in enumerate(t_values):
#         s_values = np.linspace(0, t, int(t / dt))
#         integrand = 1 / ((np.cos(Omega1_2 * s_values)) ** 2 + (kappa2 ** 2) * (np.sin(Omega1_2 * s_values))** 2)
#         integral[i] = np.trapz(integrand, dx=dt)

#     return integral

# Calculation function of η(t) for  c)
# def eta(t_values):
#     dt = 0.001
#     integral = np.zeros_like(t_values)

#     for i, t in enumerate(t_values):
#         s_values = np.linspace(0, t, int(t / dt))
#         integrand = 1 / ((np.cos(Omega1_3 * s_values)) ** 2 + (kappa3 ** 2) * (np.sin(Omega1_3 * s_values))** 2)
#         integral[i] = np.trapz(integrand, dx=dt)

#     return integral

# Calculation function of η(t) for  d)
# def eta(t_values):
#     dt = 0.001
#     integral = np.zeros_like(t_values)

#     for i, t in enumerate(t_values):
#         s_values = np.linspace(0, t, int(t / dt))
#         integrand = 1 / ((np.cos(Omega1_4 * s_values)) ** 2 + (kappa4 ** 2) * (np.sin(Omega1_4 * s_values))** 2)
#         integral[i] = np.trapz(integrand, dx=dt)

#     return integral

# Calculation function of η(t) for  a) (Periodicity condition not fulfilled)
# def eta(t_values):
#     dt = 0.001
#     integral = np.zeros_like(t_values)

#     for i, t in enumerate(t_values):
#         s_values = np.linspace(0, t, int(t / dt))
#         integrand = 1 / ((np.cos(Omega1_5 * s_values)) ** 2 + (kappa5 ** 2) * (np.sin(Omega1_5 * s_values))** 2)
#         integral[i] = np.trapz(integrand, dx=dt)

#     return integral

# Calculation function of η(t) for  b) (Periodicity condition not fulfilled)
# def eta(t_values):
#     dt = 0.001
#     integral = np.zeros_like(t_values)

#     for i, t in enumerate(t_values):
#         s_values = np.linspace(0, t, int(t / dt))
#         integrand = 1 / ((np.cos(Omega1_6 * s_values)) ** 2 + (kappa6 ** 2) * (np.sin(Omega1_6 * s_values))** 2)
#         integral[i] = np.trapz(integrand, dx=dt)

#     return integral

# R(t) function for a)
def R(t_values):
    eta_t = eta(t_values)
    return (-1j * g1 * np.exp((x1 + 1j * delta1) * eta_t)) / ((np.cos(Omega1_1 * t_values) ** 2) + kappa1 ** 2 * np.sin(Omega1_1 * t_values) ** 2)

# R(t) function for b)
# def R(t_values):
#     eta_t = eta(t_values)
#     return (-1j * g1 * np.exp((x1 + 1j * delta2) * eta_t)) / ((np.cos(Omega1_2 * t_values) ** 2) + kappa2 ** 2 * np.sin(Omega1_2 * t_values) ** 2)

# R(t) function for c)
# def R(t_values):
#     eta_t = eta(t_values)
#     return (-1j * g2 * np.exp((x1 + 1j * delta3) * eta_t)) / ((np.cos(Omega1_3 * t_values) ** 2) + kappa3 ** 2 * np.sin(Omega1_3 * t_values) ** 2)

# R(t) function for d)
# def R(t_values):
#     eta_t = eta(t_values)
#     return (-1j * g2 * np.exp((x1 + 1j * delta4) * eta_t)) / ((np.cos(Omega1_4 * t_values) ** 2) + kappa4 ** 2 * np.sin(Omega1_4 * t_values) ** 2)

# R(t) function for a) (Periodicity condition not fulfilled)
# def R(t_values):
#     eta_t = eta(t_values)
#     return (-1j * g3 * np.exp((x1 + 1j * delta5) * eta_t)) / ((np.cos(Omega1_5 * t_values) ** 2) + kappa5 ** 2 * np.sin(Omega1_5 * t_values) ** 2)

# R(t) function for b) (Periodicity condition not fulfilled)
# def R(t_values):
#     eta_t = eta(t_values)
#     return (-1j * g3 * np.exp((x1 + 1j * delta6) * eta_t)) / ((np.cos(Omega1_6 * t_values) ** 2) + kappa6 ** 2 * np.sin(Omega1_6 * t_values) ** 2)

# P(t) Population inversion
def P(time):
    #test 03_12_2023
    return ((4*(Omega1_1**2)*(np.cos(Omega1_1*time)**2)+4*Omega1_1*x1*np.cos(Omega1_1*time)*np.sin(Omega1_1*time)+(x1**2+delta1**2-4*g1**2)*np.sin(Omega1_1*time)**2)/(((4*(Omega1_1**2))*((np.cos(Omega1_1*time)**2)+(kappa1**2)*(np.sin(Omega1_1*time)**2)))))   
    #calculo original 
    #return ((4*(Omega1_1**2)*(np.cos(Omega1_1*time)**2)+(((delta1**2)-(4*(g1**2)))*(np.sin(Omega1_1*time)**2)))/(((4*(Omega1_1**2))*((np.cos(Omega1_1*time)**2)+(kappa1**2)*(np.sin(Omega1_1*time)**2)))))
    #Calculo modificado
    #return ((4*(Omega1_1**2)*(np.cos(Omega1_1*time)**2)+((-(x1**2+ 2*1j*x1*delta1 - delta1**2)-(4*(g1**2)))*(np.sin(Omega1_1*time)**2)))/(((4*(Omega1_1**2))*((np.cos(Omega1_1*time)**2)+(kappa1**2)*(np.sin(Omega1_1*time)**2)))))
        
# Time values
t_values = np.linspace(0, 10, 1000)
time = np.linspace(0,3*np.pi/Omega1_1, 1000)

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

# Blue circle of radius |g|
theta = np.linspace(0, 2 * np.pi, 100)
radius = np.abs(g1)
x = radius * np.cos(theta)
y = radius * np.sin(theta)
#plt.plot(y, x, color='blue', alpha=0.5)

# Gray circle of radius |g|/κ^2
radius = np.abs(g1) / (kappa1 ** 2)
x = radius * np.cos(theta)
y = radius * np.sin(theta)
#plt.plot(y, x, color='green', alpha=0.5)

# Population graph
plt.subplot(122)
plt.plot(time, population_values, color="pink")
plt.title('Population Inversion vs Time')
plt.xlabel('Time')
plt.ylabel('P(t)')
plt.xticks([0, math.pi/Omega1_1, 2 * math.pi/Omega1_1, 3 * math.pi/Omega1_1],
           ['0', 'π/Ω', '2π/Ω', '3π/Ω'])

plt.tight_layout()  # Adjust layout to prevent overlapping

# plt.subplot(133)
# plt.plot(t_values,result_values)

plt.show()

#Periodicity condition
p = 2
v = (Delta1/Omega1_1)+((delta1+x1)*p)/Omega0
print(v)


"""
Research: Qbit control for Quantum Computing 
Driving Fields first tests for case Ω(t) = Ω1: A Precessing Field with Oscillating Amplitude
Made on: 20/10/2023 by: Andrea Catalina Fernandez Mena
"""
# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Parameters
g = 5
delta = 4
kappa = 0.6
Omega1 = 10  #Example ∆
t_max = 10    #max time por integration

# Calculation function of η(t)
def eta(t_values):
    dt = 0.001  # Paso de integración
    integral = np.zeros_like(t_values)
    
    for i, t in enumerate(t_values):
        s_values = np.linspace(0, t, int(t / dt))
        integrand = 1 / (np.cos(Omega1 * s_values) ** 2 + kappa ** 2 * np.sin(Omega1 * s_values) ** 2)
        integral[i] = np.trapz(integrand, dx=dt)
    
    return integral

# R(t) function
def R(t_values):
    eta_t = eta(t_values)
    return -1j * g * np.exp(1j * delta * eta_t) / (np.cos(Omega1 * t_values) ** 2) + kappa ** 2 * np.sin(Omega1 * t_values) ** 2

# Time values
t_values = np.linspace(0, t_max, 1000)

# Define real and imaginary values
ReR = np.real(R(t_values))
ImR = np.imag(R(t_values))

# Graphics
plt.figure(figsize=(10, 6))
plt.plot(t_values, ReR, 'x', label='ReR(t)')
plt.plot(t_values, ImR, 'y', label='ImR(t)')
plt.xlabel('Tiempo')
plt.ylabel('ReR(t) and ImR(t)')
plt.legend()
plt.grid(True)
plt.title('Comportamiento de ReR(t) e ImR(t)')
plt.show()
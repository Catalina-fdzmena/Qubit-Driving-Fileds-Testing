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
Omega1 = 10  # Valor de ejemplo para ∆
t_max = 10    # Valor de ejemplo para el tiempo máximo

# Calculation function of η(t)
def eta(t):
    integral = 0
    dt = 0.001  # Paso de integración
    for s in np.arange(0, t, dt):
        integral += dt / (np.cos(Omega1 * s) ** 2 + kappa ** 2 * np.sin(Omega1 * s) ** 2)
    return integral

# R(t) function
def R(t):
    eta_t = eta(t)
    return -1j * g * np.exp(1j * delta * eta_t) / (np.cos(Omega1 * t) ** 2) + kappa ** 2 * np.sin(Omega1 * t) ** 2

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

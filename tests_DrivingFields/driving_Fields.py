"""
Research: Qbit control for Quantum Computing 
Driving Fields first tests for case Ω(t) = Ω1: A Precessing Field with Oscillating Amplitude
Made on: 20/10/2023 by: Andrea Catalina Fernandez Mena
"""
# Import libraries
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Parameters
g = math.sqrt(5)  #periodicity condition
delta = 4 #Delta imaginario δ
kappa = 0.6   #Tests parameters 
Omega1 = 5  
t_max = 10  #Tests parameters 

# Calculation function of η(t)
def eta(t_values):
    dt = 0.01
    integral = np.zeros_like(t_values)

    for i, t in enumerate(t_values):
        s_values = np.linspace(0, t, int(t / dt))
        integrand = 1 / (np.cos(Omega1 * s_values) ** 2 + kappa ** 2 * np.sin(Omega1 * s_values) ** 2)
        integral[i] = np.trapz(integrand, dx=dt)

    return integral

# R(t) function
def R(t_values):
    eta_t = eta(t_values)
    return (-1j * g * np.exp(1j * delta * eta_t)) / ((np.cos(Omega1 * t_values) ** 2) + kappa ** 2 * np.sin(Omega1 * t_values) ** 2)

# Time values
t_values = np.linspace(0, t_max, 1000)

# Calculate ReR(t) and ImR(t)
ReR = np.real(R(t_values))
ImR = np.imag(R(t_values))

# 3D Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot ReR(t) and ImR(t)
ax.plot(t_values, ReR, ImR, color='red', label='ReR and ImR')

# Circles
if kappa < 1:
    # Blue circle of radius |g|
    theta = np.linspace(0, 2 * np.pi, 100)
    radius = np.abs(g)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.zeros_like(x)
    ax.plot(z, y, x, color='blue', alpha=0.5)

    #Gray circle of radius |g|/κ^2
    radius = np.abs(g) / (kappa ** 2)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    ax.plot(z, y, x,  color='gray', alpha=0.5)

ax.set_xlabel('Time')
ax.set_ylabel('Re R(t)')
ax.set_zlabel('Im R(t)')
plt.title('3D Plot of Re R(t) and Im R(t) vs. Time')

plt.show()

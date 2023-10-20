#Andrea Catalina Fernandez Mena  A01197705
#Actividad Qbits Aleatorios

#Inciso D)
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#------------------------------------------------------------
n_rand = 1

theta_rand = np.arccos(2 * np.random.rand(n_rand) - 1)
phi_rand = np.random.uniform(0, 2 * np.pi, n_rand)

state_0 = np.array([1, 0])
state_0 = state_0.reshape(-1, 1)

state_1 = np.array([0, 1])
state_1 = state_1.reshape(-1, 1)

psi_rand = np.cos(theta_rand / 2) * state_0 + np.sin(theta_rand / 2) * np.exp(1j * phi_rand) * state_1

n_x_rand = np.sin(theta_rand) * np.cos(phi_rand)
n_y_rand = np.sin(theta_rand) * np.sin(phi_rand)
n_z_rand = np.cos(theta_rand)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_surface(x, y, z, color='r', alpha=0.1)

ax.plot([-1, 1], [0, 0], [0, 0], linewidth=1, color='black')
ax.plot([0, 0], [-1, 1], [0, 0], linewidth=1, color='black')
ax.plot([0, 0], [0, 0], [-1, 1], linewidth=1, color='black')

ax.text(0, 0, 1.12, r'$\left| 0 \right>$', fontsize=15, ha='center', va='center')
ax.text(1.12, 0, 0, r'$\left| + \right>$', fontsize=15, ha='center', va='center')
ax.text(-1.12, 0, 0, r'$\left| - \right>$', fontsize=15, ha='center', va='center')
ax.text(0, 0, -1.12, r'$\left| 1 \right>$', fontsize=15, ha='center', va='center')
ax.text(0, 1.12, 0, r'$\left| i \right>$', fontsize=15, ha='center', va='center')
ax.text(0, -1.12, 0, r'$\left| -i \right>$', fontsize=15, ha='center', va='center')

inBS_rand = np.sqrt(n_x_rand**2 + n_y_rand**2 + n_z_rand**2) <= 1
ax.scatter(n_x_rand[inBS_rand], n_y_rand[inBS_rand], n_z_rand[inBS_rand], c='red', marker='.', s=10)

ax.set_title('Ensemble of 1000 Random Qubits in Qubit Projective Space - A01197705')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

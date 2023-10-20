#Andrea Catalina Fernández Mena A01197705

#Bloch Sphere Dynamics Animation

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation

# Condiciones iniciales para esfera de Bloch (test clauses)
delta = 1
theta = np.pi / 4
phi = 0
a = np.cos(theta)
b = np.exp(1j * phi) * np.sin(theta)

# Guardar funcionalidades de Inciso D para coordenadas dentro de funciones lambda
U = lambda t: np.array([[np.cos((delta * t) / 2) - 1j * np.sin((delta * t) / 2), 0],[0, np.cos((delta * t) / 2) + 1j * np.sin((delta * t) / 2)]])
psi = lambda uT: np.dot(uT, np.array([a, b]))
bloch = lambda a, b: [2*np.real(np.conj(a)*b), 2*np.imag(np.conj(a)*b), np.abs(a)**2 - np.abs(b)**2]

#Incluir funciones lamda dentro de un método funcion 
def bloch_Psi_Plot(t):
  uT = U(t)
  a, b = psi(uT)
  n = bloch(a, b)
  return n

# Subplots para proyecciones de  X-Y, X-Z, & Y-Z 
fig, axes = plt.subplots(1, 3, figsize=(8, 4))
plt.suptitle('Bloch Sphere Dynamics Animation - A01197705', fontsize=13)
plt.subplots_adjust(wspace=0.5)

#Actualización de frames para posicionamiento
def update(frame):
    t = frame / ratio_frames
    n_vector = bloch_Psi_Plot(t)
    projections = [(axes[0], "<X>: (+ -)", "<Y>: (± i)", [n_vector[0], n_vector[1]]),
                   (axes[1], "<X>: (+ -)", "<Z>: (0 1)", [n_vector[0], n_vector[2]]),
                   (axes[2], "<Y>: (± i)", "<Z>: (0 1)", [n_vector[1], n_vector[2]])]

    for ax, x_label, y_label, vector in projections:
        ax.clear()
        ax.add_patch(Circle((0, 0), 1, fill=True, color='lightblue', alpha=0.7))
        ax.add_patch(Circle((0, 0), 1, fill=False, color='b', alpha=1))
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='-', linewidth=0.5, alpha=0.5)
        ax.quiver(0, 0, vector[0], vector[1], angles='xy', scale_units='xy', scale=1, color='k')

#Animacion
ratio_frames = 10  # Adjust as needed
animation_duration_s = 10
frames_total = animation_duration_s * ratio_frames
interval_frame = 1000 / ratio_frames

animation = FuncAnimation(fig, update, frames=frames_total, interval=interval_frame)

plt.show()  

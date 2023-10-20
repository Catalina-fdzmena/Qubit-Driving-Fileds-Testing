#Inciso D) Apply the operator U to the state (0), and find the Bloch vector coordinates

#Andrea Catalina Fernández Mena

import numpy as np

#Valores iniciales 

time = 0 
theta = np.pi/4
phi = 0
delta = 1 

#Representar valores iniciales del ket

a = np.cos(theta/ 2 )
b=  np.exp(1j * phi) * np.sin(theta/2)  #1j representa valores imaginarios

#Definir U(t)
Ut = np.array([ [(np.cos(delta*time)/ 2) - 1j*np.sin((delta*time)/ 2) , 0] , [0 ,(np.cos(delta*time)/ 2) + 1j*np.sin((delta*time)/ 2) ] ])

#Calcula el estado time-dependent del Qbit
   #Ket puede representarse como una matriz
   #Np.dot sirve para hacer dot producto multiplicacion de matrices
psi_t = np.dot(Ut, np.array([a,b]))

print("Time dependent qubit state: ", psi_t)

#Calculate the Bloch Vector Component
nX = 2* np.real(np.conj(psi_t[0]) * psi_t[1])  #posicion 0 es A y posicion 1 es B 
nY = 2* np.imag(np.conj(psi_t[0]) * psi_t[1])  #posicion 0 es A y posicion 1 es B 
nZ = np.abs(psi_t[0])**2 - np.abs(psi_t[1])**2

print("Bloch vector coordinates")
print("Nx = ", nX)
print("Ny = ", nY)
print("Nz = ", nZ)



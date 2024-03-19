import math as mt
import matplotlib.pyplot as plt
import numpy as np


def coeff(nx,ny,L,h):
    
    k = int(2*L/h)
    
    U = np.zeros((k,k))
    x = np.arange(-L,L,h)
    y = np.arange(-L,L,h)
    
    for xi in range(k):
        for yi in range(k):
            
            U[xi,yi] = np.cos(nx*np.pi*x[xi]/L)*np.cos(ny*np.pi*y[yi]/L) + np.cos(ny*np.pi*x[xi]/L)*np.cos(nx*np.pi*y[yi]/L)
            
            #this is the coefficient of the Acos(wt) term in the equation for resonant modes
            
            #so where this is zero in the mesh, there is a node
    
    return U


nx = 5
ny = 3

#these are the mode numbers in the x and y directions
#varying these delivers different resonant modes

L = 0.5

#L is distance from centre to edge of plate

h = 0.005
# for this i just used the same step size h for x and y

colourscale = 0.0001
#i adjusted this to make the nodes clearer
#i.e nodes are at colour boundaries



U = coeff(nx,ny,L,h)

x = np.arange(-L,L,h)
y = np.arange(-L,L,h)

plt.imshow(U, vmin=-colourscale, vmax=colourscale, cmap='coolwarm', extent=[-L,L,-L,L])

plt.xlabel('x')
plt.ylabel('y')



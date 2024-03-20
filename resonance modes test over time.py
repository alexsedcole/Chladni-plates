import math as mt
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def coeff(nx,ny,L,h,A,tspan,dt,f):
    
    nt = int(tspan/dt)
    k = int(2*L/h)
    
    U = np.zeros((nt,k,k))
    x = np.arange(-L,L,h)
    y = np.arange(-L,L,h)
    t = np.arange(0,tspan,dt)
  
    w = 2*np.pi*f
    
    for ti in range(nt):
        for xi in range(k):
            for yi in range(k):
                
                U[ti,xi,yi] = A*np.cos(w*t[ti]) * ( np.cos(nx*np.pi*x[xi]/L)*np.cos(ny*np.pi*y[yi]/L) + np.cos(ny*np.pi*x[xi]/L)*np.cos(nx*np.pi*y[yi]/L) )
            
                #this is the coefficient of the Acos(wt) term in the equation for resonant modes
            
                #so where this is zero in the mesh, there is a node
    
    return U


nx = 3
ny = 1

#these are the mode numbers in the x and y directions
#varying these delivers different resonant modes

L = 0.5

#L is distance from centre to edge of plate

h = 0.01
# for this i just used the same step size h for x and y

colourscale = 0.0001
#i adjusted this to make the nodes clearer
#i.e nodes are at colour boundaries

A = 1

tspan = 10
dt = 0.1

f = 0.5

U = coeff(nx,ny,L,h,A,tspan,dt,f)

x = np.arange(-L,L,h)
y = np.arange(-L,L,h)
t = np.arange(0,tspan,dt)

fig = plt.figure()
ax = plt.axes(projection='3d')

def update_plot(frame):
    ax.clear()
    ax.plot_surface(x, y, U[frame], cmap='viridis')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('U')
    ax.set_title('Resonance Modes Over Time')
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-A, A)

ani = animation.FuncAnimation(fig, update_plot, frames=len(t), interval=100)
plt.show()




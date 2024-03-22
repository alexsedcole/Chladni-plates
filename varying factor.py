# In this section I am importing all the libraries I will need
import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import keyboard  # just so erronous plots can be stopped more QUICKLy


# Resources
# https://uwaterloo.ca/computational-mathematics/sites/ca.computational-mathematics/files/uploads/files/linqi_shao_report_pdf.pdf#page16
# http://spiff.rit.edu/classes/phys283/lectures/two_d/two_d.html
# https://www.comsol.com/blogs/how-do-chladni-plates-make-it-possible-to-visualize-sound/
# https://thelig.ht/chladni/


class Plate:
    def __init__(self, L, h, k, t_end,freq,c,amp,m,n,stty,a,pv):
        self.L = L
        self.h = h
        self.k = k
       
        self.t_end = t_end
        self.c = c
        self.amp = amp
        self.m = m
        self.n = n
        self.stty = stty
        #self.freq = freq
        
        self.a = a
        
        self.pv = pv
        
        self.freq = a*((pi*c)/(L)) * sqrt((m**2 + n**2))
        
        
        self.nx = int(L / h) + 1
        self.ny = self.nx
        self.nt = int(t_end / k) + 1

    def mesh(self):
        x = np.arange(-self.L/2,self.L/2+self.h,self.h)  #meshgrid for surface
        y = np.arange(-self.L/2,self.L/2+self.h,self.h)

        Xg, Yg = np.meshgrid(x, y)
        return (Xg, Yg)
    
    
    def analytical(self):
        
        Ua = np.zeros((self.nx,self.ny))
        xa = np.arange(-self.L/2,self.L/2+self.h,self.h)  #meshgrid for surface
        ya = np.arange(-self.L/2,self.L/2+self.h,self.h)
        
        for xi in range(self.nx):
            for yi in range(self.ny):
                
                Ua[xi,yi] = np.cos(2*self.m*np.pi*xa[xi]/L)*np.cos(2*self.n*np.pi*ya[yi]/L) + np.cos(2*self.n*np.pi*xa[xi]/L)*np.cos(2*self.m*np.pi*ya[yi]/L)
                
                #this is the coefficient of the Acos(wt) term in the equation for resonant modes
                
                #so where this is zero in the mesh, there is a node
        
        nodeplot = np.zeros((self.nx,self.ny))
        
        #new matrix which is 1 where nodes exist in Ua
        
        for xi in range(self.nx):
            for yi in range(self.ny):
                if -self.stty < Ua[xi,yi] <self.stty:
                    nodeplot[xi,yi] = 1
        
        return nodeplot
    
    
    def solve_matrix(self):
        # In this section I am defining arrays I would need (if needed)

        U = np.zeros((self.nt, self.nx, self.ny))

        # In this section I am setting the boundary conditions/initial values

        # first version: square box bounding values

        U[0,:,:] = 0   # no initial displacement everywhere - boundary value in time
    
        U[1,:,:] = 0   # need another boundary condition for accel (should improve later!)
        #U[:,:,0], U[:,:,-1] = 0, 0 # no displacement at the boundaries always
        #U[:,0,:], U[:,-1,:] = 0, 0



        
        U[:,int(self.nx/2),int(self.ny/2)] = [self.amp*sin((k)*self.freq*i) for i in range(self.nt)]# + [0 for i in range(nt-20)]   # oscillating point at the centre
        # In this section I am implementing the numerical method

        for t in range(2, self.nt):
            for x in range(1, self.nx-1):#[i for i in range(1, nx - 1) if i != int(nx/2)]:
                for y in range(1, self.ny-1):#[i for i in range(1, ny - 1) if i != int(ny/2)]:
                    if x != int(self.nx/2) or y != int(self.ny/2):# central area version: sqrt((x-nx/2)**2+(y-ny/2)**2) >= 1:
                        uxx = (1/self.h**2) * (U[t-1,x+1,y] - 2*U[t-1,x,y] + U[t-1,x-1,y])
                        uyy = (1/self.h**2) * (U[t-1,x,y+1] - 2*U[t-1,x,y] + U[t-1,x,y-1])
                        U[t,x,y] = 2*U[t-1,x,y] - U[t-2,x,y] + self.k**2*self.c**2*uxx + self.k**2*self.c**2*uyy
                        #U[t,x,y] = (t*c/(hx*hy)) * (U[t-1,x-1,y] + U[t-1,x+1,y] + U[t-1,x,y-1] + U[t-1,x,y+1] - 4*U[t-1,x,y]) + 2*U[t-1,x,y] - U[t-2,x,y]
            
            U[t,0,:] = U[t,1,:]
            U[t,-1,:] = U[t,-2,:]
            U[t,:,0] = U[t,:,1]
            U[t,:,-1] = U[t,:,-2]
        
        return U
    
    def max_amplitude(self):
        U = self.solve_matrix()
        max_amp = np.max(np.abs(U))
        return max_amp
    
    
    
    def visualise(self):
        # In this section I am showing the results

        U = self.solve_matrix()
        Xg = self.mesh()[0]
        Yg = self.mesh()[1]

        

        '''
        for i in range(int(t_end / k) + 1):
            plt.imshow(U[i], interpolation='bilinear', norm=norm)
            plt.colorbar()
            plt.pause(0.000001)
            plt.clf()
            if keyboard.is_pressed('e'):  # if key 'q' is pressed 
                    print('Abort!')
                    break  # finishing the loop
        '''


        #'''
        fig = plt.figure()
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('U')


        ax2 = fig.add_subplot(222)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        
        

        ax3 = fig.add_subplot(223)
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')
        

        ax1.set_box_aspect([1, 1, 1])  # set aspect ratio of 3D plot
        ax2.set_aspect('equal')  # set aspect ratio of contour plot
        ax3.set_aspect('equal')  # set aspect ratio of colour plot

        colour = 'coolwarm'

 
        vmin = -self.amp
        vmax = self.amp
        
        norm = matplotlib.colors.Normalize(vmin=-self.amp, vmax=self.amp)
        
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=colour), ax=ax3, orientation='vertical', label='Displacement (m)')
      

        #modified this to only plot 'pv' proportion of the time
        
        
        for i in range(int(self.pv*(self.t_end / self.k)) ,int(self.t_end / self.k) + 1):
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax1.plot_surface(Xg, Yg, U[i], cmap=colour, vmin=vmin, vmax=vmax)
            ax1.set_zlim([-2, 2])  # set the z-axis limits
            
            #contour plot of 0 displacement
    
            ax2.contour(Xg, Yg, U[i], levels=[0], cmap=colour, vmin=vmin, vmax=vmax)
            ax2.set_title('Zero displacement')
            
            ax3.imshow(U[i], interpolation='bilinear', norm=norm, extent=[0, self.L, 0, self.L], origin='lower', cmap=colour)
  
            
            


            plt.pause(0.000001)

            if keyboard.is_pressed('e'):  # if key 'e' is pressed 
                plt.close()
                print('Abort!')
                break  # finishing the loop
        
        plt.show()

        #'''


L = 1
h = 0.02
k = 0.02
t_end =10
freq = 1
#in rad/s
cs = 0.5
amp = 1


#####################

m = 2

n = 1

#####################


pv = 0.8
#this is how far through the timespan the visualisation is started
#just to speed up the process since we only care about resonance for now
#this means we can increase the timespan (closer to steady state) without having to watch for ages

a = 2.139704
#correction factor - just for experimenting relationship between m,n and frequency


#as we use a mesh for the analytical solution, the values aren't exactly zero at the nodes
#so this value is used to determine what is considered zero
stty = 0.00001

cspan = 1


test_plate = Plate(L,h,k,t_end,freq,cs,amp,m,n,stty,a,pv)


test_plate.visualise()

print('factor =',a)

print('w =', test_plate.freq,'rad/s')

print('Max amplitide =',test_plate.max_amplitude())


#plt.imshow(test_plate.analytical(), interpolation='bilinear', extent=[-L/2, L/2, -L/2, L/2], origin='lower',vmin=-cspan, vmax=cspan)


'''

#I am iterating throug a window of correction factors
#and seeing which ones lead to the greatest peak amplitude
#i.e resonance


da = 0.000001
amin = 2.139695
amax = 2.139705

a = np.arange(amin,amax,da)

for ai in range(len(a)):
    
    test_plate = Plate(L,h,k,t_end,freq,cs,amp,m,n,stty,a[ai],pv)
    print('factor =',a[ai])
    print('w =', test_plate.freq,'rad/s')
    print('Max amplitide =',test_plate.max_amplitude())
    print('-----------------')

#found that a peak amplitude for T = 10s, m=1,n=1 exists at a = 2.1397 ish
#for some reason this maximum factor value seems to change when i change the timespan
#and is not always the same for different m,n values
#but it delivers pretty close to the expected resonance shapes

'''











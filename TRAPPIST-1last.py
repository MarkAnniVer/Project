import numpy as np
import matplotlib.pyplot as plt
import math

deg2rad = np.pi/180.0

N_bodies = 8
dt = 0.001 # stepsize in days
N_steps = 365000

G = 8.89042e-10 # gravity constant in AU^3/M_Earth/day^2

# masses in Earth mass, mass[0] is for the central star
m = [29600, 1.02, 1.16, 0.297, 0.772, 0.934, 1.148, 0.331]

#major semiaxes
R = [0.0, 0.01150,  0.01576, 0.02219, 0.02916, 0.03836, 0.0467, 0.0617]

#orbital periods
P = [0.0, 1.51087637, 2.42180746, 4.049959, 6.099043, 9.205585, 12.354473, 18.767953]

#initial orbital phases
phase = [0.0, 0.000, 217.470, 300.378, 142.558, 323.471, 269.932, 42.487]

# initialize arrays for planets positions
x = np.tile(0.0, (N_bodies, N_steps))
y = np.tile(0.0, (N_bodies, N_steps))

# initialize arrays for planets velocities
vx = np.tile(0.0, (N_bodies, N_steps))
vy = np.tile(0.0, (N_bodies, N_steps))

# initialize array for energy
E = np.tile(0.0, N_steps - 1)

#initialize initial energy
E0 = 0.0

# set planets initial positions and velocities
for i in range(1, N_bodies):    
    x[i, 0] = R[i]*np.cos(phase[i]*deg2rad)
    y[i, 0] = R[i]*np.sin(phase[i]*deg2rad)
    vx[i, 0] = -(2.0*np.pi*R[i]/P[i])*np.sin(phase[i]*deg2rad)
    vy[i, 0] = (2.0*np.pi*R[i]/P[i])*np.cos(phase[i]*deg2rad)
    
#count full initial energy
for i in range(0, N_bodies):
    E0 =E0+ m[i]*(vx[i,0]**2+vy[i,0]**2)/2.0
    for j in range(0, N_bodies):
        if j == i: continue
        E0=E0-G*m[j]*m[i]/(2.0*math.sqrt((x[j,0]-x[i,0])**2+(y[i,0]-y[j,0])**2))
            
# simulate TRAPPIST-1 evolution
for n in range (0, N_steps - 1): # n - step number
    for i in range (0, N_bodies):    
        dvx = 0.0
        dvy = 0.0
        for k in range(0, N_bodies):
            if i == k: continue            
            dx = x[i, n] - x[k, n]
            dy = y[i, n] - y[k, n]
            dr = (dx*dx + dy*dy)**0.5
            dvx = dvx - G*m[k]/dr**2*dx/dr
            dvy = dvy - G*m[k]/dr**2*dy/dr
            
        x[i, n + 1] = x[i, n] + vx[i, n]*dt 
        y[i, n + 1] = y[i, n] + vy[i, n]*dt    
        vx[i, n + 1] = vx[i, n] + dvx*dt
        vy[i, n + 1] = vy[i, n] + dvy*dt
    
    for i in range(0, N_bodies):
        E[n] =E[n]+ m[i]*(vx[i,n]**2+vy[i,n]**2)/2.0/E0
        for j in range(0, N_bodies):
            if j == i: continue
            E[n]=E[n]-G*m[j]*m[i]/(2.0*math.sqrt((x[j,n]-x[i,n])**2+(y[i,n]-y[j,n])**2)) /E0

# trajectories
plt.xlabel('AU')
plt.ylabel('AU')
plt.title('TRAPPIST-1')
plt.plot(x[0,:], y[0,:],'.')
for i in range(1, N_bodies): 
    plt.plot(x[i,:], y[i,:])
plt.show()

#energy
plt.xlabel('Steps')
plt.ylabel('E/E0')
plt.title('Energy')
plt.plot(E)
plt.show()

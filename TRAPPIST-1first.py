import numpy as np
import matplotlib.pyplot as plt
import math

deg2rad = np.pi/180.0

N_bodies = 8
dt = 0.001 # stepsize in days
N_steps = 30000

G = 8.89042e-10

m = [29600, 1.02, 1.16, 0.297, 0.772, 0.934, 1.148, 0.331]

R = [0.0, 0.01150,  0.01576, 0.02219, 0.02916, 0.03836, 0.0467, 0.0617]

P = [0.0, 1.51087637, 2.42180746, 4.049959, 6.099043, 9.205585, 12.354473, 18.767953]

phase = [0.0, 0.000, 217.470, 300.378, 142.558, 323.471, 269.932, 42.487]

x = np.tile(0.0, (N_bodies, N_steps))
y = np.tile(0.0, (N_bodies, N_steps))

vx = np.tile(0.0, (N_bodies, N_steps))
vy = np.tile(0.0, (N_bodies, N_steps))

for i in range(1, N_bodies):    
    x[i, 0] = R[i]*np.cos(phase[i]*deg2rad)
    y[i, 0] = R[i]*np.sin(phase[i]*deg2rad)
    vx[i, 0] = -(2.0*np.pi*R[i]/P[i])*np.sin(phase[i]*deg2rad)
    vy[i, 0] = (2.0*np.pi*R[i]/P[i])*np.cos(phase[i]*deg2rad)
    
for n in range (0, N_steps - 1): # n - number
    for i in range (0, N_bodies):    
        dvx = 0.0
        dvy = 0.0
        for k in range(0, N_bodies):
            if i == k: continue            
            dx = x[i, n] - x[k, n]
            dy = y[i, n] - y[k, n]
            dr = (dx*dx + dy*dy)**1.5
            dvx = dvx - G*m[k]*dx*dt/dr
            dvy = dvy - G*m[k]*dy*dt/dr
            
        x[i, n + 1] = x[i, n] + vx[i, n]*dt 
        y[i, n + 1] = y[i, n] + vy[i, n]*dt    
        vx[i, n + 1] = vx[i, n] + dvx
        vy[i, n + 1] = vy[i, n] + dvy


plt.xlabel('AU')
plt.ylabel('AU')
plt.title('TRAPPIST-1')
plt.plot(x[0,:], y[0,:],'.')
for i in range(1, N_bodies):
    plt.plot(x[i,:], y[i,:])
plt.show()

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(y,t,L):
    dydt = np.sqrt(((y**4)/(L**2))+(2*(y**3)/(L**2))+(2*y)-(L**2))
    return dydt
# initial condition
y0 =1
# time points
t = np.linspace(0,5,1000)
# solve ODE
L=-1
y = odeint(model,y0,t,args=(L,))
#print results
print(t)
print(y)
#definition of ground in coordinates x an y
plt.figure(figsize=(10,10))
# plot results
plt.plot(t,y,label=f'solución con L={L}')
plt.rcParams['legend.fontsize'] = 10
plt.legend(loc="upper right")
plt.xlabel(r'$\phi$')
plt.ylabel(r'$r(\phi)$')
#Grid
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
#definition of ground in coordinates r and phi
plt.figure(figsize=(10,10))
# plot results
plt.plot(y*np.cos(t),y*np.sin(t))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
#Grid
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
#Show the graph's
plt.show()

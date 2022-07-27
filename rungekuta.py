import numpy as np
import matplotlib.pyplot as plt
#parametros para el metodo de runge kutta.
n=1000
a=0
b=0.1
h=(b-a)/n
#Definicion de funciónes.
def f(r):
    L=1.5
    return np.sqrt(((1/(L**2))*((r**4)+(2*(r**3))))+(2*r)-(L**2))
#variables:
r=np.zeros(n)
t=np.zeros(n)
#Condiciones iniciales:
r[0]=1
t[0]=0
#aplicación del método:
for i in np.arange(1,n):
  t[i]=t[i-1]+h
  k1=h*f(r[i-1])
  k2=h*f(r[i-1]+(k1/2))
  k3=h*f(r[i-1]+(k2/2))
  k4=h*f(r[i-1]+k3)
  r[i]=r[i-1]+((1/6)*(k1+2*k2+2*k3+k4))
plt.figure()
plt.plot(t,r)
plt.savefig('plot1.png')
plt.figure()
plt.plot(r*np.cos(t),r*np.sin(t))
plt.savefig('plot2.png')
plt.show()

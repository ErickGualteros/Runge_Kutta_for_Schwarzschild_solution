import numpy as np
import matplotlib.pyplot as plt
#parametros para el metodo de runge kutta.
n=100000
a=0
b=90
h=(b-a)/n
#Definicion de funciónes.
def f(x,L):
  return 1-x+3*x**2/L**2
def g(w):
    return w
#Aplicación del método
#variables:
x=np.zeros(n)
w=np.zeros(n)
t=np.zeros(n)
#Condiciones iniciales:
L=4.5
x[0]=L**2/9
w[0]=0.1
t[0]=0
#aplicación del método:
for i in np.arange(1,n):
    t[i]=t[i-1]+h
    k1=h*f(x[i-1],L)
    k2=h*f(x[i-1]+(k1/2),L)
    k3=h*f(x[i-1]+(k2/2),L)
    k4=h*f(x[i-1]+k3,L)
    w[i]=w[i-1]+(k1+2*k2+2*k3+k4)/6
    l1=h*g(w[i-1])
    l2=h*g(w[i-1]+(l1/2))
    l3=h*g(w[i-1]+(l2/2))
    l4=h*g(w[i-1]+l3)
    x[i]=x[i-1]+(l1+2*l2+2*l3+l4)/6
#grafica de las soluciones 
plt.figure(figsize=(10,7))
plt.title('Solución Rungekutta de orden 4 'r'$\omega(\phi)$'' y 'r'$x(\phi)$'' Vs 'r'$\phi$')
plt.plot(t,x,label=r'$x(\phi)$')
plt.plot(t,w,label=r'$\omega(\phi)$')
plt.rcParams['legend.fontsize'] = 10
plt.legend(loc="upper right")
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\omega(\phi)$'' y 'r'$x(\phi)$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.savefig('plot3.png')
#grafica de las orbitas
r=L**2/x
plt.figure(figsize=(10,10))
plt.title('Orbita Obtenida')
plt.plot(r*np.cos(t),r*np.sin(t),label='Orbita')
plt.plot(2*np.cos(t),2*np.sin(t),label=r'$r=2GM$')
plt.scatter(0,0, color='red',label='Masa Central')
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.rcParams['legend.fontsize'] = 10
plt.legend(loc="upper right")
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.savefig('plot4.png')

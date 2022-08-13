import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
from pylab import *
#parametros para el metodo de runge kutta.
n=100000
a=0
b=10
h=(b-a)/n
#Definicion de funciónes.
def df(x):
    return x
def domega(omega,r,rho):
    return (2*omega*rho)/(r*(r-2))
def drho(r,rho,omega,theta,mu,nu):
    return -((rho**2)/(r*(r-2)))+((r-2)*((mu**2)+(((np.sin(theta))**2)*nu**2)-((omega**2)/(r**3))))
def dmu(theta,nu,rho,r):
    return (np.sin(theta)*np.cos(theta)*nu**2)-(2*rho/r)
def dnu(nu,theta,mu,rho,r):
    f=np.cos(theta)/np.sin(theta)
    return -(2*f*nu*mu)-(2*rho/r)
#variables
r=np.zeros(n)
t=np.zeros(n)
theta=np.zeros(n)
phi=np.zeros(n)
omega=np.zeros(n)
rho=np.zeros(n)
mu=np.zeros(n)
nu=np.zeros(n)
tau=np.zeros(n)
#condiciones iniciales
r[0]=0.4
theta[0]=np.pi/2
phi[0]=np.pi/2
omega[0]=0.4
rho[0]=0
mu[0]=0.4
nu[0]=0.4
#aplicación del método
for i in np.arange(1,n):
    tau[i]=tau[i-1]+h
    #primeras constantes
    k11=h*domega(omega[i-1],r[i-1],rho[i-1])
    k21=h*drho(r[i-1],rho[i-1],omega[i-1],theta[i-1],mu[i-1],nu[i-1])
    k31=h*dmu(theta[i-1],nu[i-1],rho[i-1],r[i-1])
    k41=h*dnu(nu[i-1],theta[i-1],mu[i-1],rho[i-1],r[i-1])
    k51=h*df(omega[i-1])
    k61=h*df(rho[i-1])
    k71=h*df(mu[i-1])
    k81=h*df(nu[i-1])
    #segundas constantes
    k12=h*domega(omega[i-1]+k11/2,r[i-1]+k61/2,rho[i-1]+k21/2)
    k22=h*drho(r[i-1]+k61/2,rho[i-1]+k21/2,omega[i-1]+k11/2,theta[i-1]+k71/2,mu[i-1]+k31/2,nu[i-1]+k41/2)
    k32=h*dmu(theta[i-1]+k71/2,nu[i-1]+k41/2,rho[i-1]+k21/2,r[i-1]+k61/2)
    k42=h*dnu(nu[i-1]+k41/2,theta[i-1]+k71/2,mu[i-1]+k31/2,rho[i-1]+k21/2,r[i-1]+k61/2)
    k52=h*df(omega[i-1]+k11/2)
    k62=h*df(rho[i-1]+k21/2)
    k72=h*df(mu[i-1]+k31/2)
    k82=h*df(nu[i-1]+k41/2)
    #terceras constantes
    k13=h*domega(omega[i-1]+k12/2,r[i-1]+k62/2,rho[i-1]+k22/2)
    k23=h*drho(r[i-1]+k62/2,rho[i-1]+k22/2,omega[i-1]+k12/2,theta[i-1]+k72/2,mu[i-1]+k32/2,nu[i-1]+k42/2)
    k33=h*dmu(theta[i-1]+k72/2,nu[i-1]+k42/2,rho[i-1]+k22/2,r[i-1]+k62/2)
    k43=h*dnu(nu[i-1]+k42/2,theta[i-1]+k72/2,mu[i-1]+k32/2,rho[i-1]+k22/2,r[i-1]+k62/2)
    k53=h*df(omega[i-1]+k12/2)
    k63=h*df(rho[i-1]+k22/2)
    k73=h*df(mu[i-1]+k32/2)
    k83=h*df(nu[i-1]+k42/2)
    #cuartas constantes
    k14=h*domega(omega[i-1]+k13,r[i-1]+k63,rho[i-1]+k23)
    k24=h*drho(r[i-1]+k63,rho[i-1]+k23,omega[i-1]+k13,theta[i-1]+k73,mu[i-1]+k33,nu[i-1]+k43)
    k34=h*dmu(theta[i-1]+k73,nu[i-1]+k43,rho[i-1]+k23,r[i-1]+k63)
    k44=h*dnu(nu[i-1]+k43,theta[i-1]+k73,mu[i-1]+k33,rho[i-1]+k23,r[i-1]+k63)
    k54=h*df(omega[i-1]+k13)
    k64=h*df(rho[i-1]+k23)
    k74=h*df(mu[i-1]+k33)
    k84=h*df(nu[i-1]+k43)
    #final
    r[i]=r[i-1]+(k61+2*k62+2*k63+k64)/6
    t[i]=t[i-1]+(k51+2*k52+2*k53+k54)/6
    theta[i]=theta[i-1]+(k71+2*k72+2*k73+k74)/6
    phi[i]=phi[i-1]+(k81+2*k82+2*k83+k84)/6
    omega[i]=omega[i-1]+(k11+2*k12+2*k13+k14)/6
    rho[i]=rho[i-1]+(k21+2*k22+2*k23+k24)/6
    mu[i]=mu[i-1]+(k31+2*k32+2*k33+k34)/6
    nu[i]=nu[i-1]+(k41+2*k42+2*k43+k44)/6
#grafica 1
plt.figure(figsize=(7,7))
plt.plot(tau,r)
plt.title(r'$r$'' Vs 'r'$\tau$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$r$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot1.png', dpi=400)
#grafica 2
plt.figure(figsize=(7,7))
plt.plot(tau,t)
plt.title(r'$t$'' Vs 'r'$\tau$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$t$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot2.png', dpi=400)
#grafica 3
plt.figure(figsize=(7,7))
plt.plot(tau,theta)
plt.title(r'$\theta$'' Vs 'r'$\tau$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\theta$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot3.png', dpi=400)
#grafica 4
plt.figure(figsize=(7,7))
plt.plot(tau,phi)
plt.title(r'$\phi$'' Vs 'r'$\tau$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\phi$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot4.png', dpi=400)
#grafica 5
plt.figure(figsize=(7,7))
plt.plot(phi,r)
plt.title(r'$r$'' Vs 'r'$\phi$')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$r$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot5.png', dpi=400)
#grafica 6
plt.figure(figsize=(7,7))
plt.plot(theta,r)
plt.title(r'$r$'' Vs 'r'$\theta$')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$r$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot6.png', dpi=400)
#grafica 7
plt.figure(figsize=(7,7))
plt.plot(2*np.cos(tau),2*np.sin(tau),label=r'$r=\frac{2GM}{c^{2}}$')
plt.plot(r*np.cos(theta),r*np.sin(theta),label='orbita')
plt.scatter(0,0,label='masa central', color='red')
plt.title('Orbita con: 'r'$r$'' Vs 'r'$\theta$')
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.rcParams['legend.fontsize']= 10
plt.legend(loc="upper right")
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot7.png', dpi=400)
#grafica 8
plt.figure(figsize=(7,7))
plt.plot(2*np.cos(tau),2*np.sin(tau),label=r'$r=\frac{2GM}{c^{2}}$')
plt.plot(r*np.cos(phi),r*np.sin(phi),label='orbita')
plt.scatter(0,0,label='masa central', color='red')
plt.title('Orbita con: 'r'$r$'' Vs 'r'$\phi$')
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.rcParams['legend.fontsize']= 10
plt.legend(loc="upper right")
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot8.png', dpi=400)
#grafica 9
plt.figure(figsize=(7,7))
plt.plot(t,r)
plt.title(r'$r$'' Vs 'r'$t$')
plt.xlabel(r'$t$')
plt.ylabel(r'$r$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot9.png', dpi=400)
#grafica 10
plt.figure(figsize=(7,7))
plt.plot(t,theta)
plt.title(r'$\theta$'' Vs 'r'$t$')
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot10.png', dpi=400)
#grafica 10
plt.figure(figsize=(7,7))
plt.plot(t,phi)
plt.title(r'$\phi$'' Vs 'r'$t$')
plt.xlabel(r'$t$')
plt.ylabel(r'$\phi$')
plt.grid(visible=True, which='major', color=(0,0,0), linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color=(0.2,0.2,0.2), linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.axvline(0,linewidth=1, color=(1,0,0))
plt.axhline(0,linewidth=1, color=(1,0,0))
plt.savefig('plot11.png', dpi=400)
# Figura en 3D
fig = plt.figure(figsize=(7,7))
# Agrrgamos un plano 3D
ax1 = fig.subplot(projection='3d')
ax1.set_title('Superficie Esférica')
#crear la orbita en 3d
A = np.zeros(n)
B = np.zeros(n)
C = np.zeros(n)
A[0]=r[0]*np.sin(theta[0])*np.cos(phi[0])
B[0]=r[0]*np.sin(theta[0])*np.sin(phi[0])
C[0]=r[0]*np.cos(theta[0])
for i in np.arange(1,n):
    A[i]=r[i]*np.sin(theta[i])*np.cos(phi[i])
    B[i]=r[i]*np.sin(theta[i])*np.sin(phi[i])
    C[i]=r[i]*np.cos(theta[i])
#grafica en 3d
ax1.plot(A,B,C,label='Orbita')
#Nombre de los ejes
ax1.set_xlabel('Eje X')
ax1.set_ylabel('Eje Y')
ax1.set_zlabel('Eje Z')
ax1.set_title('Orbita')
def sphere(r):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    d = r * np.outer(np.cos(u), np.sin(v))
    e = r * np.outer(np.sin(u), np.sin(v))
    f = r * np.outer(np.ones(np.size(u)), np.cos(v))
    return d,e,f
d,e,f = sphere(2)
ax1.plot_surface(d, e, f, rstride=3, cstride=3, color=(0,0.8,0.8,0.1))
ax1.legend()
ax1.view_init(60,35)
fig
plt.savefig('plotorbita.png', dpi=400)
plt.show()
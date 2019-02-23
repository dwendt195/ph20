import sys
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

x0=float(sys.argv[1])
v0=float(sys.argv[2])
n=int(sys.argv[3])
h=float(sys.argv[4])
name=sys.argv[5]

#mass on a spring: d2x/dt2=-k/m x

def euler(x0,v0,n,h,name):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=v[i]-h*x[i]
        x[i+1]=x[i]+h*v[i]
    fig=plt.figure()
    plt.plot(t,x,label='Position (x)')
    plt.plot(t,v,label='Velocity (v)')
    fig.savefig(name)
   
def eulererror(x0,v0,n,h,ax):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=v[i]-h*x[i]
        x[i+1]=x[i]+h*v[i]
    xerr=np.cos(t)-x
    verr=-np.sin(t)-v
    
    ax[0].plot(t,xerr,label='h = '+str(h))
    ax[1].plot(t,verr,label='h = '+str(h))
    
def eulerenergy(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=v[i]-h*x[i]
        x[i+1]=x[i]+h*v[i]
    E = x*x+v*v
    
    plt.plot(t,E)
    
def eulerimplicit(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=(v[i]-h*x[i])/(1+h*h)
        x[i+1]=(x[i]+h*v[i])/(1+h*h)
    fig=plt.plot(t,x,label='Position (x)')
    fig=plt.plot(t,v,label='Velocity (v)')
    
def eulererrorimp(x0,v0,n,h,ax):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=(v[i]-h*x[i])/(1+h*h)
        x[i+1]=(x[i]+h*v[i])/(1+h*h)
    xerr=np.cos(t)-x
    verr=-np.sin(t)-v
    
    ax[0].plot(t,xerr,label='h = '+str(h))
    ax[1].plot(t,verr,label='h = '+str(h))
    
def eulerenergyimp(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=(v[i]-h*x[i])/(1+h*h)
        x[i+1]=(x[i]+h*v[i])/(1+h*h)
    E = x*x+v*v
    
    plt.plot(t,E)
    
def eulerphaseexp(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=v[i]-h*x[i]
        x[i+1]=x[i]+h*v[i]
    fig=plt.plot(x,v,label='Explicit')
    fig=plt.plot(np.cos(t),-np.sin(t),label='Analytic')
    
def eulerphaseimp(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        v[i+1]=(v[i]-h*x[i])/(1+h*h)
        x[i+1]=(x[i]+h*v[i])/(1+h*h)
    fig=plt.plot(x,v,label='Implicit')
    
def eulerphasesymp(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        x[i+1]=x[i]+h*v[i]
        v[i+1]=v[i]-h*x[i+1]
    fig=plt.plot(x,v,label='Symplectic')
    
def eulerenergysymp(x0,v0,n,h,fig):
    t = np.zeros(n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0
    for i in range (n-1):
        t[i+1]=h*(i+1)
        x[i+1]=x[i]+h*v[i]
        v[i+1]=v[i]-h*x[i+1]
    E = x*x+v*v
    fig=plt.plot(t,E,label='Symplectic')

euler(x0,v0,n,h,name)


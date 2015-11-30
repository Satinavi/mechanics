# Hello!
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import ode
import numpy as np
from scipy.integrate import odeint
import pylab
import numpy
import sys

#lobal parameters are here
h = -(1.5*(3.0)**(1.0/3.0)+0.1)

def CreateBoundary(x, h):
        x_ret = []
        y_ret = []
        for q in x:
            if (h - 3.0/2.0*q*q > 0):
                if(1/(h-3.0/2.0*q*q)**2-q**2 > 0):
                    x_ret.append(q)
                    y_ret.append((1.0/(h-3.0/2.0*q*q)**2-q**2)**(0.5))
                    x_ret.append(q)
                    y_ret.append(-(1/(h-3.0/2.0*q*q)**2-q**2)**(0.5))
        return x_ret, y_ret

def CreateBoundaryNew(h, n):
    coeff = [1.5,0.0, -h, 1.0]
    a = np.roots(coeff)
    x = numpy.linspace(-max(a)*0.999999,max(a)*0.999999,n)
    x_ret, y_ret = CreateBoundary(x,h)
    return x_ret, y_ret
    
def Energy(q1,p1,q2,p2,h):
    ret = h + 0.5*(p1**2 + p2**2) - 1.5*q1**2 - (q1**2 + q2**2)**(-0.5)
    return ret

def NewEnergy(xi, pxi, eta, peta):
    ret = 0.5*(pxi*pxi + peta*peta) - (4 + 4*(xi*xi + eta*eta)*(h) + 6*(xi*xi - eta*eta)*(xi*xi - eta*eta)*(xi*xi + eta*eta))
    return ret

def DirTrans(xi, eta):
    return xi*xi - eta*eta, 2*xi*eta

def InvTrans(x, y):
    a = 0.5*(x + (x*x + y*y)**0.5)
    b = 0.5*(-x + (x*x + y*y)**0.5)
    if(y<=0.0):
        xi = -a**0.5
        eta = b**0.5
    if(y>0.0):
        xi = -a**0.5
        eta = -b**0.5
    return xi, eta

def f(s, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    f1 = 4.0*(x*x + y*y)
    f2 = px
#    f3 = 8.0*(x*x + y*y)*py + 8.0*x*h + 24.0*x*x*x
    f3 = 8.0*(x*x + y*y)*py + 8.0*x*h + 12*x*(x*x-y*y)**2.0 + 24.0*x*(x**4-y**4)
    f4 = py
#    f5 = -8.0*(x*x+y*y)*px + 8.0*y*h - 24.0*y*y*y
    f5 = -8.0*(x*x+y*y)*px + 8.0*y*h + 12*y*(x*x-y*y)**2.0 - 24.0*y*(x**4-y**4)
    
    return [f1, f2, f3, f4, f5]



def Df(s, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    return [[0,8.0,0,8.0,0],
            [0,0,1.0,0,0],
            [0,16.0*py*x+8.0*h+180*x**4-12*y**4-72*x**2*y**2,0,16*y*py - 48*(x*y**3+x**3*y),8*(x*x+y*y)],
            [0,0,0,0,1],
            [0,-16*x*px - 48*(y*x**3+x*y**3),-8*(x*x+y*y),-16*y*px+8*h-12*x**4+180*y**4-72*x**2*y**2,0]]

#TotalEnergy =  1.5*(3)**(1.0/3.0)+1
TotalTimeSteps = 10000.0
TotalTime = 1000.0

x_plot, y_plot = CreateBoundaryNew(-h,100)

MaxNumberOfMoons = 100000.0

fig = plt.figure(figsize=(7, 5))

plt.axis('equal')
x_plot_temp, y_plot_temp = CreateBoundaryNew(-h, MaxNumberOfMoons)

fig.tight_layout()


#print(Energy(x_plot[50],0.0,y_plot[50],0.0,-h))
#print(Energy(x_plot[150],0.0,y_plot[150],0.0,-h))
#print(Energy(x_plot[250],0.0,y_plot[250],0.0,-h))

my_plotx = []
my_ploty = []

for i in range(len(x_plot)):
    g = InvTrans(x_plot[i], y_plot[i])
    my_plotx.append(g[0])
    my_ploty.append(g[1])

plotx = []
ploty = []

for i in range(len(my_plotx)):
    g = DirTrans(my_plotx[i], my_ploty[i])
    plotx.append(g[0])
    ploty.append(g[1])


t0 = 0.0
tt = numpy.linspace(0,TotalTime,TotalTimeSteps)

for j in range(len(my_plotx)):   
    init = 0.0, my_plotx[j], 0.0, my_ploty[j], 0.0
    sol11 = numpy.zeros((TotalTimeSteps+1, 5))
    r = ode(f, Df).set_integrator('dop853')
    r.set_initial_value(init, t0)
    t1 = TotalTime
    dt = t1/TotalTimeSteps  
    i = 0
    while r.successful() and r.t <= TotalTime:
        a = r.integrate(r.t+dt)
        sol11[i,0]=a[0]
        sol11[i,1]=a[1]
        sol11[i,2]=a[2]
        sol11[i,3]=a[3]
        sol11[i,4]=a[4]
        i=i+1    
    print(NewEnergy(sol11[0,1],sol11[0,2],sol11[0,3],sol11[0,4]), NewEnergy(sol11[i-5,1],sol11[i-5,2],sol11[i-5,3],sol11[i-5,4]))
    print(sol11[i-5,1],sol11[i-5,2],sol11[i-5,3],sol11[i-5,4])
    sol1new = []
    sol2new = []  
    for i in range(len(sol11[:,1])-5):
        g = DirTrans(sol11[i,1], sol11[i,3])
        sol1new.append(g[0])
        sol2new.append(g[1])
    plt.plot(sol1new, sol2new, marker='', linestyle='-', color='r')
plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)
plt.show()  

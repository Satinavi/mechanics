import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import ode
import numpy as np
from scipy.integrate import odeint
import pylab
import numpy
import sys

#lobal parameters are here
h = -(1.5*(3.0)**(1.0/3.0))

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
    x = numpy.linspace(-max(a)*0.999,max(a)*0.999,n)
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

def f(var, s):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    f1 = 4.0*(x*x + y*y)
    f2 = px
#    f3 = 8.0*(x*x + y*y)*py + 8.0*x*h + 24.0*x*x*x
    f3 = 8.0*(x*x + y*y)*py + 8.0*x*h + 12*x*(x*x-y*y)**2.0 + 24*x*(x**4-y**4)
    f4 = py
#    f5 = -8.0*(x*x+y*y)*px + 8.0*y*h - 24.0*y*y*y
    f5 = -8.0*(x*x+y*y)*px + 8.0*y*h + 12*y*(x*x-y*y)**2.0 - 24*y*(x**4-y**4)
    
    return f1, f2, f3, f4, f5



def Df(t, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    return [[0,8,0,8,0],
            [0,0,1,0,0],
            [0,16*py*x+8*h+48*x*x,0,16*y*py,8*(x*x+y*y)],
            [0,0,0,0,1],
            [0,-16*x*px,-8*(x*x+y*y),-16*y*px+8*h-48*y*y,0]]

#TotalEnergy =  1.5*(3)**(1.0/3.0)+1
TotalTimeSteps = 100000.0
TotalTime = 20.0

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
    
init = 0.0, my_plotx[75], 0.0, my_ploty[75], 0.0
print(init)
#print("+++")
#print(NewEnergy(my_plotx[50],0.0,my_ploty[50],0.0))
#print(NewEnergy(my_plotx[150],0.0,my_ploty[150],0.0))
#print(NewEnergy(my_plotx[250],0.0,my_ploty[250],0.0))
#print(init)

t0 = 0.0
tt = numpy.linspace(0,TotalTime,TotalTimeSteps)
sol11 = odeint(f, init, tt)

#print(NewEnergy(sol11[10,1],sol11[10,2],sol11[10,3],sol11[10,4]))
#print(NewEnergy(sol11[100,1],sol11[100,2],sol11[100,3],sol11[100,4]))
#plt.plot(plotx, ploty, marker='.', linestyle='', color='black')
plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)
sol1new = []
sol2new = []

for i in range(len(sol11[:,1])):
    g = DirTrans(sol11[i,1], sol11[i,3])
    sol1new.append(g[0])
    sol2new.append(g[1])    

plt.plot(sol1new, sol2new, marker='', linestyle='-', color='r')
plt.show()  
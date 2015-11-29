import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import ode
import numpy as np
from scipy.integrate import odeint
import pylab
import numpy

#lobal parameters are here
h = 1.0

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
    
def f(s, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    f1 = 4*(x*x + y*y)
    f2 = px
    f3 = 8*(x*x + y*y)*py + 8*x*h + 24*x*x*x
    f4 = py
    f5 = -8*(x*x+y*y)*px + 8*y*h - 24*y*y*y
    return [f1, f2, f3, f4, f5]

def Df(t, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
        
    
    q1 = var[0]
    p1 = var[1]
    q2 = var[2]
    p2 = var[3]
    return [[0,1,0,0],
            [3.0 - 1.0/((q1**2 + q2**2)**(1.5)) + 3.0*q1**2/((q1**2 + q2**2)**(2.5)),0,3.0*q1*q2/((q1**2 + q2**2)**(2.5)),2.0],
            [0,0,0,1],
            [3*q1*q2/((q1**2 + q2**2)**(2.5)), -2.0, -1.0/((q1**2 + q2**2)**(1.5)) + 3*q1**2/((q1**2 + q2**2)**(2.5)),0.0]]

TotalEnergy =  1.5*(3)**(1.0/3.0)+3
TotalTimeSteps = 100000.0
TotalTime = 100.0

x_plot, y_plot = CreateBoundaryNew(1.5*(3.0)**(1.0/3.0),100)

init = [x_plot[3], 0.0, y_plot[3], 0.0]

print(init)

t0 = 0.0

#sol1 = odeint(f, init, t, Dfun=Df)

sol1 = numpy.zeros((TotalTimeSteps+1, 4))

r = ode(f, Df).set_integrator('dop853')
r.set_initial_value(init, t0)
t1 = TotalTime
dt = t1/TotalTimeSteps

i = 0
while r.successful() and r.t <= TotalTime:
    a = r.integrate(r.t+dt)
    sol1[i,0]=a[0]
    sol1[i,1]=a[1]
    sol1[i,2]=a[2]
    sol1[i,3]=a[3]
    i=i+1
    
#print(sol1)

#solx = sol1[1000000-1]
#print(Energy(solx[0], solx[1], solx[2], solx[3], 1.5*(3.0)**(1.0/3.0)))

min = 100
imin = -1
for i in range(int(TotalTimeSteps)):
    if(i>TotalTimeSteps/2):
        if(sol1[i,1]**2 + sol1[i,3]**2 < min):
            min = sol1[i,1]**2 + sol1[i,3]**2
            imin = i

print("%10.7e"% sol1[imin,1])
print("%10.7e"% sol1[imin,3])
print("%10.7e"% min)
print(imin)

print("%10.7e"% sol1[TotalTimeSteps,0])
print("%10.7e"% sol1[TotalTimeSteps,2])

#plt.plot(x_plot, y_plot, marker='.', linestyle='', color='r')
#plt.show()

#init = x_plot[250], 0.0, y_plot[250]-0.04, 0.0

#sol2 = odeint(f, init, t)

#plt.plot(x_plot, y_plot, marker='', linestyle='-', color='black')

plt.plot(sol1[0:TotalTimeSteps,0], sol1[0:TotalTimeSteps,2], marker='', linestyle='-', color='black')
#plt.plot(sol2[:,0], sol2[:,2], marker='.', linestyle='', color='b')
ax = plt.gca()

ax.set_axis_bgcolor('white')
plt.axis('equal')
plt.show()

#print(sol)

#levels = numpy.linspace(0,3.0,10) 

#X,Y = np.meshgrid(x,y)

#H = 3/2*X**2 + (X**2+Y**2)**(-0.5)

#plt.figure()
#CS = plt.contour(X, Y, H, levels, colors='k')
#plt.show()
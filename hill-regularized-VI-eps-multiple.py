# Hello!
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import ode
import numpy as np
from scipy.integrate import odeint
import pylab
import numpy
import sys
import math

#lobal parameters are here
h = -(1.5*(3.0)**(1.0/3.0)+1.0)

scale = 1.1

ang_glob = []

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

def CreateAngleBoundary(nn, add_ang):
    x_ret = []
    y_ret = []
    start_ang = math.pi*0.0 # rad, to avoid critical cases
    for i in range(int(nn)):
        cur_ang = i*math.pi/nn+start_ang
        if(add_ang):
            ang_glob.append(cur_ang)
        r_sol = SolveCubic(1.5*(math.cos(cur_ang))**2, h)
        x_sol = r_sol*math.cos(cur_ang)
        y_sol = r_sol*math.sin(cur_ang)
        x_ret.append(x_sol)
        y_ret.append(y_sol)
    return x_ret, y_ret
        
        
def SolveCubic(a,b):
    eps = 0.000000000000001
    doit = True
    r_start = 0.0
    r_cur = 0.0
    r_step = 0.00001
    j = 0
    while doit:
        r_cur = r_start + j*r_step
        if(a*r_cur**3  + b*r_cur + 1.0 < 0):
            doit = False
        j=j+1
    x1 = r_cur - r_step
    x2 = r_cur
    err = 1000.0
    #print("OK", a*x1**3  + b*x1 + 1.0, a*x2**3  + b*x2 + 1.0)
    while abs(err) > eps:
        mid = (x1+x2)*0.5
        err = a*(mid)**3 + b*mid + 1.0
        if(err > 0.0):
            x1 = mid
        else:
            x2 = mid
    #print(mid)
    return mid
    

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
    f3 = 8.0*(x*x + y*y)*scale*py + 8.0*x*h + 12*x*(x*x-y*y)**2.0 + 24.0*x*(x**4-y**4)
    f4 = py
#    f5 = -8.0*(x*x+y*y)*px + 8.0*y*h - 24.0*y*y*y
    f5 = -8.0*(x*x+y*y)*scale*px + 8.0*y*h + 12*y*(x*x-y*y)**2.0 - 24.0*y*(x**4-y**4)
    
    return [f1, f2, f3, f4, f5]



def Df(s, var):
    t = var[0]
    x = var[1]
    px = var[2]
    y = var[3]
    py = var[4]
    return [[0,8.0,0,8.0,0],
            [0,0,1.0,0,0],
            [0,scale*16.0*py*x+8.0*h+180*x**4-12*y**4-72*x**2*y**2,0,scale*16*y*py - 48*(x*y**3+x**3*y),8*(x*x+y*y)],
            [0,0,0,0,1],
            [0,-scale*16*x*px - 48*(y*x**3+x*y**3),-8*(x*x+y*y),-scale*16*y*px+8*h-12*x**4+180*y**4-72*x**2*y**2,0]]

#TotalEnergy =  1.5*(3)**(1.0/3.0)+1
TotalTimeSteps = 100000.0
TotalTime = 1000.0

MaxNumberOfMoons = 50.0

x_plot, y_plot = CreateAngleBoundary(MaxNumberOfMoons, True)
x_plot_temp, y_plot_temp = CreateAngleBoundary(MaxNumberOfMoons, False)

fig = plt.figure(figsize=(7, 5))

plt.axis('equal')

fig.tight_layout()

# 1.4470751940907589e-05 - 10   -1.0480505352461478e-13
# 1.23841256186787e-05 - 50 -2.5757174171303632e-13
# 1.2384125619881664e-05 - 100 1.1235457009206584e-13
# 1.2364317562106447e-05 - 1000 -5.6266102888002933e-13
# 1.2364072819163203e-05 - 5000 -2.6237900740966325e-12
# 1.2364072819163203e-05 - 10000 -5.4388715753361794e-12

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

nearest = []

for j in range(len(my_plotx)):
    init = 0.0, my_plotx[j], 0.0, my_ploty[j], 0.0
    sol11 = numpy.zeros((int(TotalTimeSteps)+1, 5))
    r = ode(f, Df).set_integrator('dop853')
    r.set_initial_value(init, t0)
    t1 = TotalTime
    dt = t1/TotalTimeSteps  
    i = 0
    fine_sol0 = []
    fine_sol1 = []
    fine_sol2 = []
    fine_sol3 = []
    fine_sol4 = []
    while r.successful() and r.t <= TotalTime:
        a = r.integrate(r.t+dt)
        fine_sol0.append(a[0])
        fine_sol1.append(a[1])
        fine_sol2.append(a[2])
        fine_sol3.append(a[3])
        fine_sol4.append(a[4])
        xx, yy = DirTrans(fine_sol1[i], fine_sol3[i])
        if(abs(xx) + abs(yy) < 0.01):
            dt = t1/TotalTimeSteps/100.0
        else:
            dt = t1/TotalTimeSteps
        i=i+1    
    print(j, len(my_plotx), NewEnergy(fine_sol1[0],fine_sol2[0],fine_sol3[0],fine_sol4[0]), NewEnergy(fine_sol1[i-5],fine_sol2[i-5],fine_sol3[i-5],fine_sol4[i-5]))
    #print(j, len(my_plotx), sol11[i-5,1],sol11[i-5,2],sol11[i-5,3],sol11[i-5,4])
    sol1new = []
    sol2new = []
    sol1newinv = []
    sol2newinv = []
    near = 1000.0
    for i in range(len(fine_sol1)-5):
        g = DirTrans(fine_sol1[i], fine_sol3[i])
        sol1new.append(g[0])
        sol2new.append(g[1])
        sol1newinv.append(-g[0])
        sol2newinv.append(-g[1])
        if((g[0]*g[0] + g[1]*g[1])**0.5 < near):
            near = (g[0]*g[0] + g[1]*g[1])**0.5
    nearest.append(near)
    plt.plot(sol1new, sol2new, marker='', linestyle='-', color='r')
    plt.plot(sol1newinv, sol2newinv, marker='', linestyle='-', color='r')
plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)
print(nearest)
plt.show()
#plt.plot(ang_glob, nearest, marker='', linestyle='-', color='r')
#plt.show()  

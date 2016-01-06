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
import time

#lobal parameters are here
r_00 = 0.9
omega=1.0
kappa = 1.0

h = -2.0

ee = 0.2

mu = 0.33

scale = 1.1

ang_glob = []

r_0 = 1.6

r_1 = 1.65

#def Beta(r):
#    if(r<=1):
#        return 0.0
#    if(r>1 and r <2.0 ):
#        return omega*(3*(r-1)**2-2*(r-1)**3)
#    if(r>=2.0):
#        return omega
#
#def DBetaDr(r):
#    if(r<=1):
#        return 0.0
#    if(r>1 and r <2.0 ):
#        return omega*6*(r-1)*(2-r)
#    if(r>=2.0):
#        return 0.0
#    return 0.0

def Beta(r):
    if(r<=0):
        return -omega
    if(r<=1):
        return (3*r**2-2*r**3-1)*omega
    if(r>1):
        return 0.0

def DBetaDr(r):
    if(r<=0):
        return 0.0
    if(r<=1):
        return (6*r-6*r**2)*omega
    if(r>1):
        return 0.0


def Sgm(r):
    if(r<r_0):
        return -1.0
    if(r>r_1):
        return 0.0
    if( r>=r_0 and r<=r_1 ):
        return -(1.0 - (6*((r-r_0)/(r_1-r_0))**5-15*((r-r_0)/(r_1-r_0))**4+10*((r-r_0)/(r_1-r_0))**3))

def DSgmDr(r):
    if(r<r_0):
        return 0.0
    if(r>r_1):
        return 0.0
    if( r>=r_0 and r<=r_1 ):
        return (30*((r-r_0)/(r_1-r_0))**4-60*((r-r_0)/(r_1-r_0))**3+30*((r-r_0)/(r_1-r_0))**2) 


def CreateBoundary(N):
    r_ret = []
    phi_ret = []
    for i in range(int(N)):
        r_ret.append(r_0)
        phi_ret.append(i*2*math.pi/N)
    return r_ret, phi_ret    

def ConvertToXY(r_in, phi_in):
    x_out = []
    y_out = []
    for i in range(len(r_in)):
        x_out.append(r_in[i]*math.cos(phi_in[i]))
        y_out.append(r_in[i]*math.sin(phi_in[i]))
    return x_out, y_out

#def f2(var, t):
#    r = var[0]
#    dr = var[1]
#    beta0 = Beta(r_0)
#    f1 = dr
#    f2 = (r*Beta(r)*Beta(r)+beta0**2*r_0**4/r**3 - 2*Beta(r)*beta0*r_0**2/r)+r-1/r + (-Beta(r)+beta0*r_0**2/r**2)*(2*r*Beta(r)+r**2*DBetaDr(r))
#    return f1, f2

def f22(var, t):
    r = var[0]
    dr = var[1]
    sgm0 = Sgm(r_00)
    f1 = dr
    f2 = (r*Sgm(r)*Sgm(r)+sgm0**2*r_00**4/r**3 - 2*Sgm(r)*sgm0*r_0**2/r)+r-1/r/r + (-Sgm(r)+sgm0*r_00**2/r**2)*(2*r*Sgm(r)+r**2*DSgmDr(r))
    return f1, f2

def f33(var, t):
    r = var[0]
    dr = var[1]
    sgm0 = Sgm(r_00)
    f1 = dr
    f2 = r*Sgm(r)*Sgm(r)+r-r**2*Sgm(r)*DSgmDr(r)-2*r*Sgm(r)*Sgm(r)-1/r/r
    return f1, f2

def f11(var, t):
    r = var[0]
    dr = var[1]
    f1 = dr
    f2 = r*Beta(r)*Beta(r)-r-r**2*Beta(r)*DBetaDr(r)-2*r*Beta(r)*Beta(r)
    return f1, f2

def f(var, t):
    r = var[0]
    pr = var[1]
    ph = var[2]
    pph = var[3]
    f1 = pr
    f2 = r*pph**2 - kappa*r + 2*r*Beta(r)*pph + DBetaDr(r)*r**2*pph
    f3 = pph
    f4 = -2.0/r*pr*pph - DBetaDr(r)*pr - 2.0*Beta(r)*pr/r    
    return f1, f2, f3, f4




print("HELLO!")
TotalTimeSteps = 100000.0
TotalTime = 5
MaxNumberOfMoons = 10.0


#fig = plt.figure(figsize=(7, 5))

#plt.axis('equal')
#
#fig.tight_layout()
#
#r = numpy.linspace(-1.0,2.0,TotalTimeSteps)
#fr = []
#
#for j in range(len(r)):
#    fr.append(Beta(r[j]))
#
#plt.plot(r, fr, marker='', linestyle='-', color='r',linewidth=3)  
#
#plt.show()
##
#sys.exit(0)
#
#
#print((r[j]*Beta(r[j])*Beta(r[j])+beta0**2*r_0**4/r[j]**3 - 2*Beta(r[j])*beta0*r_0**2/r[j])-r[j] + (-Beta(r[j])+beta0*r_0**2/r[j]**2)*(2*r[j]*Beta(r[j])+r[j]**2*DBetaDr(r[j])))
#
#fig = plt.figure(figsize=(7, 5))
#
#plt.axis('equal')
#
#fig.tight_layout()
#
#r = numpy.linspace(0.000000000,r_0,TotalTimeSteps)
#fr = []
#beta0 = Beta(r_0)
#print((r[0]*Beta(r[0])*Beta(r[0])+beta0**2*r_0**4/r[0]**3 - 2*Beta(r[0])*beta0*r_0**2/r[0])-r[0] + (-Beta(r[0])+beta0*r_0**2/r[0]**2)*(2*r[0]*Beta(r[0])+r[0]**2*DBetaDr(r[0])))
#
#for j in range(len(r)):
#    fr.append((r[j]*Beta(r[j])*Beta(r[j])+beta0**2*r_0**4/r[j]**3 - 2*Beta(r[j])*beta0*r_0**2/r[j])-r[j] + (-Beta(r[j])+beta0*r_0**2/r[j]**2)*(2*r[j]*Beta(r[j])+r[j]**2*DBetaDr(r[j])))
#
#plt.plot(r, fr, marker='', linestyle='-', color='r')  
#
#plt.show()
##
#sys.exit(0)


print("HELLO!")

r_plot_temp, phi_plot_temp = CreateBoundary(MaxNumberOfMoons)

fig = plt.figure(figsize=(7, 5))

plt.axis('equal')



fig.tight_layout()


init = r_00, 0.0

t = numpy.linspace(0,TotalTime,TotalTimeSteps)

sol = odeint(f22, init, t,atol=1.0e-13, rtol=1.0e-13)
r_out = sol[:,0]
dr_out = sol[:,1] 
plt.plot(r_out, dr_out, marker='', linestyle='-', color='r')  

plt.show()
#
sys.exit(0)

x_plot_temp, y_plot_temp = ConvertToXY(r_plot_temp, phi_plot_temp)

plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)



for j in range(len(r_plot_temp)):    
    init = r_plot_temp[j], 0.0, phi_plot_temp[j], 0.0
    t = numpy.linspace(0,TotalTime,TotalTimeSteps)
    sol = odeint(f, init, t)
    rr_out = sol[:,0]
    phi_out = sol[:,2] 
    x_out = []
    y_out = []
    x_out, y_out  = ConvertToXY(rr_out, phi_out)
    plt.plot(x_out, y_out, marker='', linestyle='-', color='r')   
    
plt.show()

sys.exit(0)


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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import odeint
from matplotlib import gridspec
import sys
import math
import pylab
import numpy

# Parameters: A, B, C, a, b, c, h
# We assume r_s to be (0, 0, c)
# and A > B > C (almost Lagrange case)

# Here we define parameters
A = 2.0
B = 3.0
C = 1.0
a = 0.0
b = 0.0
c = 1.0

p0 = 0
q0 = 0
r0 = 0.0

SolNum = 100
TotEnerg = -0.8
TotalTimeSteps = 1000000.0
TotalTime = 100.0

# Create initial condtions for a given h
def CreateInitCond(r_s, h, n):
    gamma_1_ret = []
    gamma_2_ret = []
    gamma_3_ret = []
    alpha = math.acos(h/r_s[2])
    for i in range(n):
        gamma_1_ret.append(math.sin(alpha)*math.cos(i*2.0*math.pi/n))
        gamma_2_ret.append(math.sin(alpha)*math.sin(i*2.0*math.pi/n))
        gamma_3_ret.append(h/r_s[2])
    return gamma_1_ret, gamma_2_ret, gamma_3_ret

def CreateInitCond2(r_s, h, n):
    gamma_1_ret = []
    gamma_2_ret = []
    gamma_3_ret = []
    alpha = math.acos((h-0.5*(A*p0*p0 + B*q0*q0 + C*r0*r0))/r_s[2])
    for i in range(n):
        gamma_1_ret.append(math.sin(alpha)*math.cos(i*2.0*math.pi/n))
        gamma_2_ret.append(math.sin(alpha)*math.sin(i*2.0*math.pi/n))
        gamma_3_ret.append((h-0.5*(A*p0*p0 + B*q0*q0 + C*r0*r0))/r_s[2])
    return gamma_1_ret, gamma_2_ret, gamma_3_ret

# (x,y,z) -> (xi, eta)
def StereographicProj(x, y, z):
    xi=[]
    eta=[]
    for i in range(len(x)):
        xi.append(2.0*x[i]/(1.0-z[i]))
        eta.append(2.0*y[i]/(1.0-z[i]))
    return xi, eta

def f(var, t):
    p = var[0]
    q = var[1]
    r = var[2]
    gamma_1 = var[3]
    gamma_2 = var[4]
    gamma_3 = var[5]
    f1 = 1.0/A*(-(C-B)*q*r + c*gamma_2 - b*gamma_3)
    f2 = 1.0/B*(-(A-C)*r*p + a*gamma_3 - c*gamma_1)
    f3 = 1.0/C*(-(B-A)*p*q + b*gamma_1 - a*gamma_2)
    f4 = -(q*gamma_3 - r*gamma_2)
    f5 = -(r*gamma_1 - p*gamma_3)
    f6 = -(p*gamma_2 - q*gamma_1)
    return f1, f2, f3, f4, f5, f6

def CreateBoundary(x, h_in):
        x_ret = []
        y_ret = []
        h=-1.0*h_in
        print(h)
        for q in x:
            if (h - 3.0/2.0*q*q > 0):
                if(1/(h-3.0/2.0*q*q)**2-q**2 > 0):
                    x_ret.append(q)
                    y_ret.append((1.0/(h-3.0/2.0*q*q)**2-q**2)**(0.5))
                    x_ret.append(q)
                    y_ret.append(-(1/(h-3.0/2.0*q*q)**2-q**2)**(0.5))
        return x_ret, y_ret

def Energy(q1,p1,q2,p2,h_in):
    h=-h_in
    ret = h + 0.5*(p1**2 + p2**2) - 1.5*q1**2 - (q1**2 + q2**2)**(-0.5)
    return ret

def CreateBoundaryNew(h, n):
    coeff = [1.5,0.0, h, 1.0]
    a = np.roots(coeff)
    x = numpy.linspace(-max(a)*0.999,max(a)*0.999,int(n))
    x_ret, y_ret = CreateBoundary(x,h)
    return x_ret, y_ret

#def f(var, t):
#    q1 = var[0]
#    p1 = var[1]
#    q2 = var[2]
#    p2 = var[3]
#    f1 = p1
#    f3 = p2
#    f2 = 2*p2 + 3*q1 - q1/((q1**2 + q2**2)**(1.5))
#    f4 = -2*p1 - q2/((q1**2 + q2**2)**(1.5))
#    return f1, f2, f3, f4

r_s = [a, b, c]

gamma_12, gamma_22, gamma_32 = CreateInitCond(r_s, TotEnerg, SolNum)
gamma_1, gamma_2, gamma_3 = CreateInitCond2(r_s, TotEnerg, SolNum)
xi, eta = StereographicProj(gamma_12, gamma_22, gamma_32)

fig = plt.figure(figsize=(7, 5))
plt.axis('equal')
t = numpy.linspace(0,TotalTime,TotalTimeSteps)
#init = 0.0, 0.0, 0.0, gamma_1[0], gamma_2[0], gamma_3[0]
init = p0, q0, r0, gamma_1[5], gamma_2[5], gamma_3[5]
sol = odeint(f, init, t)

x_out = []
y_out = []
z_out = []
x_out = sol[:,3]
y_out = sol[:,4]
z_out = sol[:,5]

#print(gamma_1[50]*gamma_1[50] + gamma_2[50]*gamma_2[50] + gamma_3[50]*gamma_3[50])
#
#print(x_out[5000]*x_out[5000] + y_out[5000]*y_out[5000] + z_out[5000]*z_out[5000])
#
#print(x_out[0]*a + y_out[0]*b + z_out[0]*c)
#print(x_out[1]*a + y_out[1]*b + z_out[1]*c)
#print(x_out[1000]*a + y_out[1000]*b + z_out[1000]*c)

xi_out = []
eta_out = []
xi_out, eta_out = StereographicProj(x_out, y_out, z_out)
      
plt.plot(xi_out, eta_out, marker='', linestyle='-', color='r')
plt.scatter(xi,eta, s=5, facecolors='red', edgecolors='none', alpha=1)
#plt.scatter(x_out,y_out, s=5, facecolors='red', edgecolors='none', alpha=1)

#plt.show()
#
##Picture a sphere
#X_sp = []
#Y_sp = []
#Z_sp = []
#theta = numpy.linspace(0,math.pi,10)
#phi = numpy.linspace(0,2.0*math.pi,20)
#for tt in range(len(theta)):
#    x_sp = []
#    y_sp = []
#    z_sp = []
#    for pp in range(len(phi)):
#        x_sp.append(math.sin(theta[tt])*math.cos(phi[pp]))
#        y_sp.append(math.sin(theta[tt])*math.sin(phi[pp]))
#        z_sp.append(math.cos(theta[tt]))
#    X_sp.append(x_sp)
#    Y_sp.append(y_sp)
#    Z_sp.append(z_sp)
#
#ax = fig.gca(projection='3d')
#ax.axis('equal')
#for i in range(10):
#    ax.plot(X_sp[i], Y_sp[i], Z_sp[i], marker=' ', linestyle='-', color='b')    
#
#
#ax.plot(x_out, y_out, z_out, marker='', linestyle='-', color='r')
#plt.show()

sys.exit(0)

TotalEnergy =  -1.5*(3)**(1.0/3.0)-0.1
TotalTimeSteps = 1000000.0
TotalTime = 400.0
Accel = 5.0
MaxNumberOfMoons = 100000.0

fig = plt.figure(figsize=(7, 5))

plt.axis('equal')
x_plot_temp, y_plot_temp = CreateBoundaryNew(TotalEnergy, MaxNumberOfMoons)

fig.tight_layout()
#plt.show()

#sys.exit(0)

MaxNumberOfMoons = 100.0
x_plot, y_plot = CreateBoundaryNew(TotalEnergy, MaxNumberOfMoons)

t = numpy.linspace(0,TotalTime,TotalTimeSteps)

x_out = numpy.zeros((len(x_plot), len(t)))
y_out = numpy.zeros((len(y_plot), len(t)))

#for i in range(len(x_plot)):

i=75;
print(x_plot[50], y_plot[50])
init = x_plot[i], 0.0, y_plot[i], 0.0
sol1 = odeint(f, init, t)
x_out[i] = sol1[:,0]
y_out[i] = sol1[:,2]
plt.plot(sol1[:,0], sol1[:,2], marker='', linestyle='-', color='r')

plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)
plt.show()
    

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
import time

# Parameters: A, B, C, a, b, c, h
# We assume r_s to be (0, 0, c)
# and A > B > C (almost Lagrange case)

# Here we define parameters
A = 1.0
B = A
C = 3.0
a = 0.5
b = 0.0
c = 1.0

k = 1.5

p0 = 0
q0 = 0
r0 = 0.0

SolNum = 100
TotEnerg = 1.0
TotalTimeSteps = 10000.0
TotalTime = 20.0

def VLevels(theta, phi):
    ret = 0.5*k*k/(A*math.sin(theta)*math.sin(theta) + C*math.cos(theta)*math.cos(theta))
    ret = ret + a*math.sin(theta)*math.sin(phi) + b*math.sin(theta)*math.cos(phi) + c*math.cos(theta)
    return ret
 
def SolveEq(phi, a1, b1, h1):
    tol = 0.000001
    fun = 1.0 + h1
    while (abs(fun-h1) > tol):
        c1 = 0.5*(b1+a1)
        if((VLevels(c1, phi)-h1)*(VLevels(a1,phi)-h1) > 0):
#            print(VLevels(a1, phi)-h1, VLevels(c1,phi)-h1)            
            a1 = c1
#            time.sleep(1)
        else:
#            print(VLevels(c1,phi)-h1, VLevels(b1, phi)-h1)
            b1 = c1
#            time.sleep(1)
        fun = VLevels(c1,phi)
#    print(fun)
    return c1

def CalcEn(theta, pt, phi, pp):
    ret = 0.0
#    print("GGG")
    div = A*math.sin(theta)*math.sin(theta) + C*math.cos(theta)*math.cos(theta)
#    print("GGG")
    R2 = 0.5*A*pt*pt + 0.5*pp*pp*A*C*math.sin(theta)*math.sin(theta)/div
#    print("GGG")    
    R0 = -0.5*k*k/div - a*math.sin(theta)*math.sin(phi)-b*math.sin(theta)*math.cos(phi)-c*math.cos(theta)
#    print("GGG")
    ret = R2-R0
    print(ret)
    return ret
 
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

def TransIntoGamma(phi, theta):
    return math.sin(theta)*math.sin(phi), math.sin(theta)*math.cos(phi), math.cos(theta)

# (x,y,z) -> (xi, eta)
def StereographicProj(x, y, z):
    xi=[]
    eta=[]
    for i in range(len(x)):
        xi.append(2.0*x[i]/(1.0-z[i]))
        eta.append(2.0*y[i]/(1.0-z[i]))
    return xi, eta

def vec(var,t):
    theta = var[0]
    pt = var[1]
    phi = var[2]
    pp = var[3]
    
    vec1 = pt
    vec3 = pp
    
    div = A*math.sin(theta)*math.sin(theta) + C*math.cos(theta)*math.cos(theta)
    
    vec2_part1 = -(k - C*pp*math.cos(theta))*C*pp*math.sin(theta)/div
    vec2_part2 = (A-C)*math.sin(theta)*math.cos(theta)*(k - C*pp*math.cos(theta))*(k - C*pp*math.cos(theta))/div/div
    vec2_part3 = -a*math.cos(theta)*math.sin(phi)-b*math.cos(theta)*math.cos(phi)+c*math.sin(theta)
    vec2 = (vec2_part1 + vec2_part2 + vec2_part3)/A
    
    div2 = (C - C*C*math.cos(theta)*math.cos(theta)/div)   
    
    vec4_part1 = -(2*C*C*pp*pt*math.sin(theta)*math.cos(theta)-C*k*pt*math.sin(theta))/div
    vec4_part2 =2.0*pt*(k-C*pp*math.cos(theta))*C*math.cos(theta)
    vec4_part2 = vec4_part2*(A-C)*math.sin(theta)*math.cos(theta)/div/div
    vec4_part3 = -a*math.sin(theta)*math.cos(phi)+b*math.sin(theta)*math.sin(phi)

    vec4 = (vec4_part1 + vec4_part2 + vec4_part3)/div2
    
    return vec1, vec2, vec3, vec4

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

theta = []
phi = numpy.linspace(0,2*math.pi,200)

for i in range(len(phi)):
    theta.append(SolveEq(phi[i], 0.0, math.pi, TotEnerg))

#print(CalcEn(theta[50], 0.0, phi[50], 0.0))
#print("AAA")

fig = plt.figure(figsize=(7, 5))
plt.axis('equal')

out1 = []
out2 = []
out3 = []

for i in range(len(phi)):
    g_ret = TransIntoGamma(phi[i], theta[i])
    out1.append(g_ret[0])
    out2.append(g_ret[1])
    out3.append(g_ret[2])
    
xi, eta = StereographicProj(out1, out2, out3)

for j in range(len(theta)):    
    init = theta[j], 0.0, phi[j], 0.0
    t = numpy.linspace(0,TotalTime,TotalTimeSteps)
    sol = odeint(vec, init, t)
    x_out = []
    y_out = []
    x_out = sol[:,0]
    y_out = sol[:,2]
    out1a = []
    out2a = []
    out3a = []
    
    for i in range(len(x_out)):
        g_ret = TransIntoGamma(y_out[i], x_out[i])
        out1a.append(g_ret[0])
        out2a.append(g_ret[1])
        out3a.append(g_ret[2])
    
    xi_out = []
    eta_out = []
    xi_out, eta_out = StereographicProj(out1a, out2a, out3a)   
    plt.plot(xi_out, eta_out, marker='', linestyle='-', color='r')   
    
#print(CalcEn(sol[0,0],sol[0,1],sol[0,2], sol[0,3]))
#print("JJJ")
#print(CalcEn(sol[10,0],sol[10,1],sol[10,2], sol[10,3]))
#print("HHH")
#print(CalcEn(sol[5001,0],sol[5001,1],sol[5001,2], sol[5001,3]))



#plt.plot(t, x_out, marker='', linestyle='-', color='r')
#
#
#plt.show()
#
#sys.exit(0)

#print(x_out[10], theta[5], y_out[10], phi[5])




plt.plot(xi, eta, marker='', linestyle='-', color='r')



plt.show()

sys.exit(0)
















V = numpy.zeros((len(phi), len(theta)))
for i in range(len(theta)):
    for j in range(len(phi)):
        V[j][i] = (VLevels(theta[i], phi[j]))



plt.figure()
CS = plt.contourf(theta, phi, V, 100)
plt.show()


sys.exit(0)

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
    

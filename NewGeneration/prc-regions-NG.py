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
h = -2

ee = 0.2

mu = 0.33

scale = 1.1

ang_glob = []

r_0 = 0.35

r_1 = 0.37

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
    ret = 0.5*(pxi*pxi + peta*peta) - (4.0*(xi*xi + eta*eta)*(h) + 2.0*(xi*xi+eta*eta)*((xi*xi-eta*eta - mu)**2+4.0*xi*xi*eta*eta) + 4.0*(1.0-mu) + 4.0*(xi*xi+eta*eta)*mu/((xi**2-eta**2-1)**2+4.0*xi*xi*eta*eta)**0.5)
    return ret

def DirTrans(xi, eta):
    return xi*xi - eta*eta - mu, 2*xi*eta

def InvTrans(x, y):
    x = x + mu
    a = 0.5*(x + (x*x + y*y)**0.5)
    b = 0.5*(-x + (x*x + y*y)**0.5)
    if(y<=0.0):
        xi = -a**0.5
        eta = b**0.5
    if(y>0.0):
        xi = -a**0.5
        eta = -b**0.5
    return xi, eta

def f(var, t):
    s = var[0]
    xi = var[1]
    px = var[2]
    eta = var[3]
    py = var[4]
    f0 = 4*(xi*xi+eta*eta)
    f1 = px
    f2 = 8*(xi*xi+eta*eta)*py + 8.0*h*xi + 8.0*mu*xi*(4.0*eta**2*xi**2 + (-eta**2 + xi**2 - 1)**2)**(-0.5) + mu*(4.0*eta**2 + 4.0*xi**2)*(-4.0*eta**2*xi - 2.0*xi*(-eta**2 + xi**2 - 1))*(4.0*eta**2*xi**2 + (-eta**2 + xi**2 - 1)**2)**(-1.5) + 4.0*xi*(4.0*eta**2*xi**2 + (-eta**2 - mu + xi**2)**2) + (2.0*eta**2 + 2.0*xi**2)*(8.0*eta**2*xi + 4*xi*(-eta**2 - mu + xi**2))
    f3 = py
    f4 = -8*(xi*xi+eta*eta)*px+ 8.0*eta*h + 8.0*eta*mu*(4.0*eta**2*xi**2 + (-eta**2 + xi**2 - 1)**2)**(-0.5) + 4.0*eta*(4.0*eta**2*xi**2 + (-eta**2 - mu + xi**2)**2) + mu*(4.0*eta**2 + 4.0*xi**2)*(-4.0*eta*xi**2 + 2.0*eta*(-eta**2 + xi**2 - 1))*(4.0*eta**2*xi**2 + (-eta**2 + xi**2 - 1)**2)**(-1.5) + (2.0*eta**2 + 2.0*xi**2)*(8.0*eta*xi**2 - 4*eta*(-eta**2 - mu + xi**2))
    return f0, f1, f2, f3, f4



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


def CreateAngleBoundary(nn, x0, y0, step):
    x_ret = []
    y_ret = []
    start_ang = 0.0
    for i in range(int(nn)):
        find_a_root = True
        x_cur = x0
        y_cur = y0
        cur_ang = i*2*math.pi/nn+start_ang
        while find_a_root:
            print(TotEnergy(0,0,x_cur,y_cur), TotEnergy(0,0,x_cur+step*math.cos(cur_ang),y_cur+step*math.sin(cur_ang)), h)
            #time.sleep(1)
            if( TotEnergy(0,0,x_cur,y_cur) <= h and TotEnergy(0,0,x_cur+step*math.cos(cur_ang),y_cur+step*math.sin(cur_ang)) > h):
                print("TRYING TO FIND A ROOT")
                root = FineRoot(x_cur,y_cur,x_cur+step*math.cos(cur_ang),y_cur+step*math.sin(cur_ang))
                x_ret.append(root[0])
                y_ret.append(root[1])
                find_a_root = False
            x_cur = x_cur+step*math.cos(cur_ang)
            y_cur = y_cur+step*math.sin(cur_ang)
    return x_ret, y_ret


def FineRoot(x1,y1,x2,y2):
    epsil = 0.000000000000001
    x_mid = 0.5*(x1+x2)
    y_mid = 0.5*(y1+y2)
    print(x_mid, y_mid, TotEnergy(0,0,x1,y1), TotEnergy(0,0,x2,y2))
    while (abs(TotEnergy(0,0,x_mid,y_mid)-h)>epsil):
#        print(abs(TotEnergy(0,0,x_mid,y_mid)), epsil)
        if( TotEnergy(0,0,x_mid,y_mid) < h):
            x1 = 0.5*(x1+x2)
            y1 = 0.5*(y1+y2)
#            print("A")
        else:
            x2 = 0.5*(x1+x2)
            y2 = 0.5*(y1+y2)
#            print("B")
        x_mid = 0.5*(x1+x2)
        y_mid = 0.5*(y1+y2)
    print(x_mid, y_mid, abs(TotEnergy(0,0,x_mid,y_mid)))
    return x_mid, y_mid

def TotEnergy(px, py, x, y):
    T = (px**2 + py**2)*0.5
    V = (x**2 + y**2)*0.5 + (1-mu)/((x+mu)**2+y**2)**0.5+(mu)/((x-1+mu)**2+y**2)**0.5
    return T-V

def Sgm(r):
    if(r<r_0):
        return -1.0
    if(r>r_1):
        return 0.0
    if( r>=r_0 and r<=r_1 ):
        return -(1.0 - (6*((r-r_0)/(r_1-r_0))**5-15*((r-r_0)/(r_1-r_0))**4+10*((r-r_0)/(r_1-r_0))**3))

def L0(x,y):
    return 0.5*(x**2 + y**2) + (1-mu)/((x+mu)**2+y**2)**0.5 + mu/((x-1+mu)**2+y**2)**0.5

def L1(x,y, px, py):
    return (x*py - y*px)*Sgm(((x+mu)**2+y**2)**0.5)

def L2(px,py):
    return 0.5*(px**2+py**2)

def FindBhPlus(xBh, yBh, NN):
    MM = 10
    retX = []
    retY = []
    for iii in range(len(xBh)):
        xx = xBh[iii]
        yy0 = yBh[iii]
        for jjj in range(NN):
            addPoint = True
            yy = yy0 + jjj*(-2*yy0/NN)
            for kkk in range(MM):
                nu = kkk*math.pi*2/MM
                l1 = L1(xx,yy,math.cos(nu),math.sin(nu))
                l0 = L0(xx,yy)
                l2 = L2(math.cos(nu),math.sin(nu))
                if(4*(l0+h)*l2<=l1*l1):
                    addPoint = False
            if(addPoint):
                retX.append(xx)
                retY.append(yy)
    return retX, retY
            

print("HELLO!")
TotalTimeSteps = 10000.0
TotalTime = 0.33
MaxNumberOfMoons = 800.0

#delta = 0.033
#x = np.arange(-2.0, 2.0, delta)
#y = np.arange(-2.0, 2.0, delta)
#V = numpy.zeros((len(x), len(y)))
#for xx in range(len(x)):
#    for yy in range(len(y)):
#        V[xx][yy] = TotEnergy(0,0,x[xx],y[yy])
#
#plt.figure()
#CS = plt.contourf(x, y, V, 50)
#cbar = plt.colorbar(CS)
#plt.show()
#
#sys.exit(0)


print("HELLO!")

x_plot_temp, y_plot_temp = CreateAngleBoundary(MaxNumberOfMoons, 0.0, 0.001, 0.05)


BhX, BhY = FindBhPlus(x_plot_temp, y_plot_temp, 1000)




xi_plot_temp = []
eta_plot_temp = []

#for j in range(len(x_plot_temp)):
#    xe = InvTrans(x_plot_temp[j], y_plot_temp[j])
#    xi_plot_temp.append(xe[0])
#    eta_plot_temp.append(xe[1])

fig = plt.figure(figsize=(7, 5))

plt.axis('equal')

fig.tight_layout()

plt.scatter(x_plot_temp,y_plot_temp, s=5, facecolors='red', edgecolors='none', alpha=1)

plt.scatter(BhX,BhY, s=5, facecolors='navy', edgecolors='none', alpha=1)

fff = InvTrans(0.123, -0.7)

print(DirTrans(fff[0], fff[1]))



for j in range(len(x_plot_temp)):    
    print(j, len(x_plot_temp))
    init = 0.0, xi_plot_temp[j], 0.0, eta_plot_temp[j], 0.0
    t = numpy.linspace(0,TotalTime,TotalTimeSteps)
    sol = odeint(f, init, t)
    xi_out = sol[:,1]
    eta_out = sol[:,3]
    print(NewEnergy(sol[TotalTimeSteps-5,1], sol[TotalTimeSteps-5,2], sol[TotalTimeSteps-5,3], sol[TotalTimeSteps-5,4]))
    print(NewEnergy(sol[5,1], sol[5,2], sol[5,3], sol[5,4]))    
    x_out = []
    y_out = []
    for jj in range(len(t)):
        x_out.append(DirTrans(xi_out[jj], eta_out[jj])[0])
        y_out.append(DirTrans(xi_out[jj], eta_out[jj])[1])
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

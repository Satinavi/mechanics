import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi
sqrt = np.sqrt
cos = np.cos
sin = np.sin

def f(x, t):
    q1 = x[0]
    p1 = x[1]
    q2 = x[2]
    p2 = x[3]

    f1 = p1
    f2 = 3*q1 + 2*p2 + q1/((q1**2 + q2**2)**(3/2))
    f3 = p2
    f4 =  -2*p1 + q2/((q1**2 + q2**2)**(3/2))
    
    return [f1, f2, f3, f4]

t = np.linspace(0, 50, 500)
x0 = [1.0, 0, 1,0, 0]
xx = integrate.odeint(f, x0, t)

#print(t)

#print(xx)

#print(xx[:,0])

#print(xx[:,2])

#sol = np.zeros((20,2))
#sol[:,0] = xx[:,0]
#sol[:,1] = xx[:,2]


plt.figure(1)
plt.plot(xx[:,0],xx[:,2],'.')
plt.show()


"""
Created on Thu Feb 21 20:24:15 2019
@author: Subhajit Pramanik
"""

def vectorfield(w, t, p):
    """
Defines the differential equations for the coupled spring-mass system.
Arguments:
w : vector of the state variables:
w = [x1,y1,x2,y2]
t : time
p : vector of the parameters
p = [m1,m2,k1,k2,L1,L2,b1,b2]
"""
x1, y1, x2, y2 = w
m1, m2, k1, k2, L1, L2, b1, b2 = p
#def f -> (x1',y1',x2',y2'):
f = [y1,
     (-b1 * y1 - k1 * (x1 - L1) + k2 * (x2 - x1 - L2)) / m1,
     y2,(-b2 * y2 - k2 * (x2 - x1 - L2)) / m2]
return f

from scipy.integrate import odeint

# Parameter values
# Masses:
m1 = 10.0
m2 = 15.0

# Spring constants
k1 = 3.0
k2 = 6.0

# Natural lengths
L1 = 0.5
L2 = 1.0

# Friction coefficients
b1 = 0.8
b2 = 0.5

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
x1 = 0.5
y1 = 0.0
x2 = 2.25
y2 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 250
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [m1, m2, k1, k2, L1, L2, b1, b2]
w0 = [x1, y1, x2, y2]

#ODE solver.
wsol = odeint(vectorfield, w0, t, args=(p,),
atol=abserr, rtol=relerr)
with open('two_springs.dat', 'w') as f:
    for t1, w1 in zip(t, wsol):
        print(f, t1, w1[0], w1[1], w1[2], w1[3])
        
# Plot the solution that was generated
from numpy import loadtxt
from pylab import figure, plot, xlabel, grid, legend, title, savefig
from matplotlib.font_manager import FontProperties
t, x1, xy, x2, y2 = loadtxt('two_springs.dat', unpack=True)
figure(1, figsize=(6, 4.5))
xlabel('t')
grid(True)
lw=1

plot(t, x1, 'b', linewidth=lw)
plot(t, x2, 'g', linewidth=lw)

legend((r'$x_1$', r'$x_2$'), prop=FontProperties(size=10))
title('Mass Displacements for the\nCoupled Spring-Mass System')
savefig('two_springs.png', dpi=50)
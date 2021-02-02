# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:37:23 2021

@author: Sabrina
"""
from numpy import arange, pi, sin, cos, sqrt, array, copy, array
from numpy.linalg import norm
from pylab import plot, show, axes

#this is just to plot the sun
rs = 0.00465
thetarange = arange(0,2*pi,0.001)
x = []
y = []

for theta in thetarange:
    x.append(rs*cos(theta))
    y.append(rs*sin(theta))
    
#define initial conditions, from astronomical data (thanks Wikipedia)
dt = 0.01
xm = 1.382
ym = 0
xj = 4.95
yj = 0 
em = 0.0934
ej = 0.0489 
vmx = 0
vmy = sqrt(4*pi**2*(1+em)/(xm*(1-em)))
vjx = 0
vjy = sqrt(4*pi**2*(1+ej)/(xj*(1-ej)))
ms = 1.99e30
mj = 1.90e27
mm = 6.42e23

#planet class for efficiency
class planet: 
    
    def __init__(self, mass, v0, e, r0):
        self.mass = mass
        self.v = copy(v0)
        self.r = copy(r0)
        self.e = e   
        self.x = [self.r[0]]
        self.y = [self.r[1]]
        

    def derivs(self, varis, Other):
        r = varis[:2]
        v = varis[2:]
        rmj = sqrt((r[0]-Other.r[0])**2 + (r[1]-Other.r[1])**2)
        xDeriv = v[0]
        yDeriv = v[1]
        vxDeriv = (-4*pi**2*r[0]/norm(r)**3-
                   4*pi**2*Other.mass/ms*(r[0]-Other.r[0])/norm(rmj)**3)
        vyDeriv = (-4*pi**2*r[1]/norm(r)**3-
                   4*pi**2*Other.mass/ms*(r[1]-Other.r[1])/norm(rmj)**3)
    
        return array([xDeriv,yDeriv,vxDeriv,vyDeriv])
    
    
    def RK4Trajectory(self, Other, dt):
        t = 0
        varis = array([self.r[0],self.r[1],self.v[0],self.v[1]])
        while t < 13:
            k1 = dt * self.derivs(varis, Other)
            k2 = dt * self.derivs(varis + 1/2 * k1, Other)
            k3 = dt * self.derivs(varis + 1/2 * k2, Other)
            k4 = dt * self.derivs(varis + k3, Other)
        
            varis += 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            self.x.append(varis[0])
            self.y.append(varis[1])

#start mars and jupiter on the x axis
dt = 0.001
mars = planet(mm, [vmx, vmy], em, [xm, ym])
jupiter = planet(mj, [vjx, vjy], ej, [xj, yj])
mars.RK4Trajectory(jupiter, dt)
jupiter.RK4Trajectory(mars, dt)

plot(mars.x, mars.y, 'tomato')
plot(jupiter.x, jupiter.y, 'springgreen')
plot(x,y,"y")
axes().set_aspect('equal')
axes().set_facecolor('k')
show()

#without jupiter, same setup
jupiterIsDead = planet(0,[vjx, vjy], ej, [xj, yj])
mars2 = planet(mm, [vmx, vmy], em, [xm, ym])
jupiterIsDead.RK4Trajectory(mars2, dt)
mars2.RK4Trajectory(jupiterIsDead, dt)
plot(mars2.x, mars2.y, "tomato")
plot(x,y,'y')
axes().set_aspect('equal')
axes().set_facecolor('k')
show()

#show the change in x and y that jupiter causes
plot((array(mars.x)-array(mars2.x)), (array(mars.y)-array(mars2.y)), 'indigo')
show()


mars3 = planet(mm, [vmy, vmx], em, [ym, xm])
jupiter3 = planet(mj, [vjx, vjy], ej, [xj, yj])
mars3.RK4Trajectory(jupiter3, dt)
jupiter3.RK4Trajectory(mars3, dt)

plot(mars3.x, mars3.y, 'tomato')
plot(jupiter3.x, jupiter3.y, 'springgreen')
plot(x,y,"y")
axes().set_aspect('equal')
axes().set_facecolor('k')
show()

#without jupiter
jupiterIsDead2 = planet(0,[vjx, vjy], ej, [xj, yj])
mars4 = planet(mm, [vmy, vmx], em, [ym, xm])
jupiterIsDead2.RK4Trajectory(mars4, dt)
mars4.RK4Trajectory(jupiterIsDead2, dt)
plot(mars4.x, mars4.y, "tomato")
plot(x,y,'y')
axes().set_aspect('equal')
axes().set_facecolor('k')
show()

plot((array(mars3.x)-array(mars4.x)), (array(mars3.y)-array(mars4.y)), 'indigo')
show()
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 15:16:47 2020

@author: Chai Wah Wu
"""
import numpy as np

class dynamical_system:
    def system_equations_ext(self,t,x):
        """ dynamical system equations extended with arc length
        """
        y = self.system_equations(t,x)
        return np.append(y,np.linalg.norm(y))

class chua_oscillator(dynamical_system):
    """ 3 dimensioanl dynamical system with a chaotic attractor.
        L. O. Chua, C. W. Wu, A. Huang, and G.-Q. Zhong, 
        "A universal circuit for studying and generating chaos part I: Routes to chaos," 
        IEEE Transactions on Circuits and Systems-I: Fundamental Theory and Applications, 
        vol. 40, pp. 732-744, Oct. 1993.
    """    
    def __init__(self):
        self.a = -1.14
        self.b = -0.714
        self.k = 1
        self.alpha = 9
        self.beta = 14
        self.gamma = 0.01
        self.startT = 0
        self.endT = 200
        self.init = [0.1,0.2,0.03]
       
    def f(self,x):
        return self.b*x+0.5*(self.a-self.b)*(abs(x+1)-abs(x-1))
    
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0] = self.k*self.alpha*(x[1]-x[0]-self.f(x[0]))
        y[1] = self.k*(x[0]-x[1]+x[2])
        y[2] = self.k*(-self.beta*x[1]-self.gamma*x[2])
        return y


class rucklidge(dynamical_system):
    """ Dynamical system model of convection 
        Rucklidge, A. M., 1992, “Chaos in Models of Double Convection,” 
        J. Fluid Mech., 237(1), pp. 209–229.
    """    
    def __init__(self):
        self.kappa = -2
        self.lambdaparam = -6.7
        self.startT = 0
        self.endT = 400
        self.init = [1,0,4.5]
    
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0]=self.kappa*x[0]-self.lambdaparam*x[1]-x[1]*x[2]
        y[1]=x[0]
        y[2]= -x[2]+x[1]**2
        return y
    
class nonauto_chaotic_system(dynamical_system):
    """ nonautonomous dynamical system with a chaotic attractor driving by a 
        sinusoidal signal. 
        C. W. Wu, G.-Q. Zhong and L. O. Chua, 
        "Synchronizing nonautonomous chaotic systems without phase-locking," 
        Journal of Circuits, Systems, and Computers, vol. 6, no. 3, pp. 227-241, 1996.  
        projecting (x,y,t) onto (x,y*cos(omega*t),y*sin(omega*t))
    """   
    def __init__(self):
        self.a = -1.37
        self.b = -0.84
        self.k = 1
        self.omega = 0.4
        self.A = 0.5
        self.beta = 0.895
        self.startT = 0
        self.endT=100*np.pi/self.omega # T is the total time
        self.init=[-9.6, 5.9, 0.38]

    def f(self,x):
        return self.b*x+0.5*(self.a-self.b)*(abs(x+1)-abs(x-1))

    def system_equations(self,t,x): 
        y = np.zeros((3,1)) 
        s = self.A*np.sin(self.omega*t)
        if np.abs(np.cos(self.omega*t)) > 0.1:
            xp = x[1]/np.cos(self.omega*t)
        else:
            xp = x[2]/np.sin(self.omega*t)
        yp = self.k*self.beta*(-x[0]-xp)
        y[0] = self.k*(xp-self.f(x[0]+s))
        y[1] = yp*np.cos(self.omega*t)-self.omega*xp*np.sin(self.omega*t)
        y[2] = yp*np.sin(self.omega*t)+self.omega*xp*np.cos(self.omega*t)
        return y
    
class rossler_hyperchaos(dynamical_system):
    """ 4 dimensional hyperchaotic dynamical system 
        O. Rossler, "An equation for hyperchaos," Physics Letters A,
        vol. 71, no. 2–3, pp. 155–157, 1979.
    """    
    def __init__(self):
        self.a = 0.25
        self.b = 3.0
        self.c = 0.5
        self.d = 0.05
        self.startT = 0
        self.endT = 200
        self.init = [-10,-6,0,10]
    
    def system_equations(self,t,x):
        y = np.zeros((4,1)) 
        y[0]=-x[1]-x[2]
        y[1]=x[0]+self.a*x[1]+x[3]
        y[2]=self.b+x[0]*x[2]
        y[3]=-self.c*x[2]+self.d*x[3]
        return y    
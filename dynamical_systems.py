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
        return np.append(y,np.linalg.norm(y[:3]))
    def get_init(self):
        """ get initial conditions """
        return self.init
    def get_Tinterval(self):
        """ get starting and ending time """
        return self.startT, self.endT

class chua_oscillator(dynamical_system):
    """ 3 dimensional dynamical system with a chaotic attractor.
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
    """ nonautonomous dynamical system driven by a 
        sinusoidal signal with a chaotic attractor. 
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
        self.endT=100*np.pi/self.omega
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

class hyperchaotic_circuit(dynamical_system):
    """ 4 dimensional hyperchaotic circuit with a chaotic attractor.
        T. Matsumoto, L. O. Chua, and K. Kobayashi, 
        "Hyperchaos: Laboratory experiment and numerical confirmation," 
        IEEE Transactions on Circuits and Systems, 
        vol. 33, pp. 1143-1147, Nov. 1986.
    """    
    def __init__(self):
        self.C1 = 0.5
        self.C2 = 0.05
        self.L1 = 1
        self.L2 = 2/3
        self.R = 1
        self.a = -0.2
        self.b = 3
        self.startT = 0
        self.endT = 200
        self.init = [0.1,0.2,0.03,0.4]
       
    def f(self,x):
        return self.b*x+0.5*(self.a-self.b)*(abs(x+1)-abs(x-1))
    
    def system_equations(self,t,x):
        y = np.zeros((4,1)) 
        y[0] = 1/self.C1*(self.f(x[1]-x[0])-x[2])
        y[1] = 1/self.C2*(-self.f(x[1]-x[0])-x[3])
        y[2] = 1/self.L1*(x[0]+self.R*x[2])
        y[3] = 1/self.L2*x[1]
        return y     

class lorenz(dynamical_system):
    """ 3 dimensional dynamical system with a chaotic attractor.
        E. N. Lorenz, "Deterministic nonperiodic flow,"
        Journal of the Atmospheric Sciences, vol. 20, no. 2, pp. 130–141, 1963.
    """    
    def __init__(self):
        self.sigma = 10.0
        self.rho = 28.0
        self.beta = 8.0/3.0
        self.startT = 0
        self.endT = 30
        self.init = [0.1,0.1,0.1]
       
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0] = self.sigma*(x[1] - x[0])
        y[1] = x[0]*(self.rho - x[2]) - x[1]
        y[2] = x[0]*x[1]-self.beta*x[2]
        return y
    
class chen(dynamical_system):
    """ 3 dimensional dynamical system with a chaotic attractor.
        G. Chen and T. Ueta, "Yet another chaotic attractor,"
        International Journal of Bifurcation and Chaos, vol. 9, no. 7, pp. 1465–1466, 1999.
    """    
    def __init__(self):
        self.a = 40.0
        self.b = 3.0
        self.c = 28.0
        self.startT = 0
        self.endT = 30
        self.init = [-0.1,0.5,0.6]
       
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0] = self.a*(x[1] - x[0])
        y[1] = (self.c-self.a)*x[0] - x[0]*x[2] +self.c*x[1]
        y[2] = x[0]*x[1]-self.b*x[2]
        return y

class arneodo(dynamical_system):
    """ 3 dimensional dynamical system with a chaotic attractor.
        A. Arneodo, P. Coullet, and E. A. Spiegel, "Chaos in a finite macroscopicsystem," 
        Physics Letters, vol. 92A, no. 8, pp. 369-373, 1982.
    """    
    def __init__(self):
        self.mu = -1.0
        self.mu0 = -5.5
        self.mu1 = 3.5
        self.mu2 = 1.0
        self.startT = 0
        self.endT = 200
        self.init = [0.2,0.2,-0.75]
       
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0] = x[1]
        y[1] = x[2]
        y[2] = self.mu*x[0]**3-self.mu0*x[0]-self.mu1*x[1]-self.mu2*x[2]
        return y
       
class brockett(dynamical_system):
    """ 3 dimensional dynamical system with a chaotic attractor.
        R. W. Brockett, "On conditions leading to chaos in feedback systems,"
        Proceedings of the IEEE Conference on decision and control, pp. 932-936, 1982.
    """    
    def __init__(self):
        self.k = 1.8
        self.a = 1.25
        self.startT = 0
        self.endT = 300
        self.init = [1.31, 0.57, 0.32]
    
    def f(self,x):
        return -self.k*x if abs(x) <= 1 else 2*self.k*x-3*self.k*np.sign(x)
       
    def system_equations(self,t,x):
        y = np.zeros((3,1)) 
        y[0] = x[1]
        y[1] = x[2]
        y[2] = -self.f(x[0])-self.a*x[1]-x[2]
        return y


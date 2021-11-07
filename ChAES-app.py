# -*- coding: utf-8 -*-
"""
Created on Sun Nov  7 10:33:46 2021

@author: Chai Wah Wu

Chaotic Attractor Explorer with Streamlit (ChAES)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from dynamical_systems import *
import streamlit as st
import pandas as pd
import plotly.express as px

def solve_for_length(sol,s,p0):
    """ find the time when the length of a curve is a given value
        Input: sol = ode solution, values can be found at any time
               s = required curve length
               p0 = initial guess 
    """
    return minimize(lambda x:(sol(x)[-1]-s)**2,p0).x[0]

@st.cache
def generate_attractor(dyn_system):
    """  based on the MATLAB file runrucklidge.m from "Modeling Dynamical Systems for 3D Printing” 
     by Stephen K. Lucas, Evelyn Sander, and Laura Taalman and available at
     http://math.gmu.edu/~sander/EvelynSite/supplementary-materials-for.html
     
     Python code for plotting chaotic
     attractors. Output is data that is
     initially equally spaced points along the
     curve, but subdivide if the angle between
     line segments is not large enough
     
     ported to Python by Chai Wah Wu, 12/14/2020
     features added in the port: 
         1. dynamical systems are encapsulated in dynamical_system class.
         2. arclength extension to the system equation are in dynamical_system base class
         3. use scipy minimize module to find the time when a certain length of the trajectory is reached.
    """
    
    NUM=1000 # NUM is the initial number of
             # pieces of equal arc length
    
    T0, T = dyn_system.get_Tinterval() # T0 is the starting time, T is the ending time
    init = dyn_system.get_init() # initial conditions
    method = dyn_system.get_integration_method() # integration method
    
    # run system from T0 to T and use ending state as initial conditions
    init = list(solve_ivp(dyn_system.system_equations_ext,(T0,T),init+[0],method=method).y[:-1,-1]) # solve the system, initial length=0
    # y is the solution at a given time
    sol = solve_ivp(dyn_system.system_equations_ext,(T0,T),init+[0],dense_output=True,method=method).sol
    y=sol(T)
    arclen=y[-1] # arclen is the length of the curve
    # Initialize data arrays, pos is position
    # data, seg is segment length data
    segnum=NUM # current number of segments
    # First point
    pos=np.zeros((segnum+1,3))
    pos[0,:]=init[:3]
    # Last point
    pos[segnum,:]=y[:3]
    # Location along arc of points
    segment=np.arange(segnum+1)*arclen/segnum
    time=np.zeros(segnum+1)
    time[0]=T0
    time[segnum]=T
    # time is time at each position
    
    for i in range(segnum-1):
       time[i+1]=solve_for_length(sol,segment[i+1], time[i])
       y=sol(time[i+1])
       pos[i+1,:]=y[:3]
    
    done=False;
    # Keep going until no segments need to be
    # split
    while not done:
       done=True
       # Identify line segment vectors and
       # lengths
       vec=np.zeros((segnum,3))
       lenvec=np.zeros(segnum)
       for i in range(segnum):
          vec[i,:]=pos[i+1,:]-pos[i,:]
          lenvec[i]=np.linalg.norm(vec[i,:])
    
       # Identify angles between segments and
       # which ones to split
       split=np.zeros(segnum).astype(int)
       for i in range(segnum-1):
          l2 = (lenvec[i]*lenvec[i+1])
          if l2 > 1e-9:
              ang=np.degrees(np.arccos(np.clip(np.dot(vec[i,:],vec[i+1,:])/l2,-1,1)))
              if ang>10:
                 split[i]=1
                 split[i+1]=1
                 done=False
       
       # Split segments at half the distance
       # between endpoints
       newn=segnum+np.sum(split)
       newpos=np.zeros((newn+1,3))
       newpos[0,:]=init[:3]
       newseg=np.zeros(newn+1)
       newseg[0]=0
       newtime=np.zeros(newn+1)
       newtime[0]=T0
       newtime[newn]=T
       j=0
       for i in range(segnum):
          if split[i]:  # Add a new midpoint
             j+=1
             news=(segment[i]+segment[i+1])/2
             newtime[j]=solve_for_length(sol,news,time[i])
             y=sol(newtime[j])
             newpos[j,:]=y[:3]
             newseg[j]=news
          
          # Add endpoint
          j +=1
          newpos[j,:]=pos[i+1,:] 
          newseg[j]=segment[i+1]
          newtime[j]=time[i+1]
       
       # Update variable list
       segment=newseg
       pos=newpos 
       time=newtime
       segnum=newn
    return time, pos

systemsdict = {
    'Chua': chua_oscillator,
    'Rucklidge': rucklidge,
    'nonauto chaos': nonauto_chaotic_system,
    'Rossler hyperchaos': rossler_hyperchaos,
    'hyperchaos': hyperchaotic_circuit,
    'Lorenz': lorenz,
    'Chen': chen,
    'Arneodo': arneodo,
    'Brockett': brockett,
    'Sparrow': sparrow,
    'Shinriki': shinriki,
    'Dmitriev': dmitriev
    }

systemslabels = sorted(systemsdict.keys())

st.title('Chaotic Attractor Explorer with Streamlit (ChAES)')

dyn_system_label = st.sidebar.radio('Select dynamical system',systemslabels,systemslabels.index('nonauto chaos')) # select dynamical system
dyn_system = systemsdict[dyn_system_label]()

time, pos = generate_attractor(dyn_system)

figwidth=1000
figheight=1000
df = pd.DataFrame(pos)
df.columns = ['x1','x2','x3']
fig = px.line_3d(df,x='x1',y='x2',z='x3',width=figwidth,height=figheight)
st.write(dyn_system.__doc__.replace('\n','  \n'))
st.plotly_chart(fig,width=figwidth,height=figheight)

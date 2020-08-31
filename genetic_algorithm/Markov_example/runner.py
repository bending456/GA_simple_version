#
# Template for simulation to be run on each process 
#
from scipy.integrate import odeint
import numpy as np 
import math

class Runner():
  def __init__(self):
    self.params={
        'k1': 0.3,
        'k2': 40000,
        'k3': 2.4,
        'k4': 50000,
        'k5': 1.58,
        'k6': 7000,
        'L1': 0.0001,
        'L2': 0.004,
        'L3': 0.1,
        'L4': 0.00000000000001}         

  def odefunc(self,ys,ts,params,ATP):
    ## Assigned your parameters applied in ODEs
    k1 = params['k1']
    k2 = params['k2']
    k3 = params['k3']
    k4 = params['k4']
    k5 = params['k5']
    k6 = params['k6']
    L1 = params['L1']
    L2 = params['L2']
    L3 = params['L3']
    L4 = params['L4']
    ## Assigned dependent variables from ys 
    C1, C2, C3, C4, Q1, Q2, Q3, Q4 = ys
    
    ### ATP switch
    rep = math.floor(ts/40)
    if (rep % 2) == 0:
        A = ATP
    elif (rep % 2) == 1:
        A = 0

    ### ODEs: There should be 8 ordinary differential equations for each state starting from C1, C2, ...
    dC1dt = L1*C4 + k1*C2       - (3*k2*A + L4)*C1
    dC2dt = 3*k2*A*C1 + 2*k3*Q1 - (k1 + 2*k4*A)*C2
    dC3dt = 2*k1*Q4 + 3*k2*A*C4 - (k1 + 2*k2*A)*C3
    dC4dt = L4*C1 + k1*C3       - (L1 + 3*k2*A)*C4
    dQ1dt = 2*k4*A*C2 + 3*k5*Q2 - (2*k3 + k6*A)*Q1
    dQ2dt = k6*A*Q1 + L2*Q3     - (L3 + 3*k5)*Q2
    dQ3dt = L3*Q2 + k2*A*Q4     - (L2 + 3*k1)*Q3
    dQ4dt = 3*k1*Q3 + 2*k2*A*C3 - (k2*A + 2*k1)*Q4
    
    ## Write out dydt
    dydt = [dC1dt,dC2dt,dC3dt,dC4dt,dQ1dt,dQ2dt,dQ3dt,dQ4dt]
    
    return dydt
    
  def simulate(self, 
               varDict     = None,      # dictionary of parameters to be used for simulation
               returnDict  = dict(),    # dictionary output to be returned by simulation
               jobDuration = 1,          # [second]
               ATP         = 32
               ):
    
    if varDict is None:
      varDict = self.params; 
  
    # range is in [s]
    ts = np.linspace(0,jobDuration,jobDuration*1000)
    # initial value of dependent variable
    y0s = [1,0,0,0,0,0,0,0]
    ys = odeint(self.odefunc,y0s,ts,args=(varDict,ATP*1e-6,))
    
    g12 = 1.5e-8
    g34 = 4.5e-8
    V   = -60e-3
    E   = 0

    I = g12*(ys[:,4] + ys[:,5])*(V-E) + g34*(ys[:,6] + ys[:,7])*(V-E)
    
    data        = dict() # This is a nuisance, but keeping backward compatible w Gotran stuff 
    data['t']   = ts          
    data['Cai'] = I
  
    returnDict['data'] = data       
    return returnDict
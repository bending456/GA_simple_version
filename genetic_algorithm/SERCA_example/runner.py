
#
# Template for simulation to be run on each process 
#
from scipy.integrate import odeint
import numpy as np 

class Runner():
  def __init__(self):
    self.params={
        'VmaxSC'  : 5,
        'kdSC'    : 0.1,
        'nSC'     : 2,
        'VmaxRyR' : 0.1,
        'kdRyR'   : 0.8,
        'nRyR'    : 4,
        'leakRate': 0.1 }       

  def odefunc(self,ys,ts,params):
    ## Assigned your parameters applied in ODEs
    VmaxSC   = params['VmaxSC']
    kdSC     = params['kdSC']
    nSC      = params['nSC']
    VmaxRyR  = params['VmaxRyR']
    kdRyR    = params['kdRyR']
    nRyR     = params['nRyR']
    leakRate = params['leakRate']
    ## Assigned dependent variables from ys 
    Cai = ys
    
    ## Write out ODEs
    dCadt = -VmaxSC/(1+(kdSC/Cai)**nSC) + VmaxRyR/(1+(Cai/kdRyR)**nRyR) + leakRate*(0.7-Cai)
    
    ## Write out dydt
    dydt = dCadt
    
    return dydt
    
  def simulate(self, 
               varDict     = None,      # dictionary of parameters to be used for simulation
               returnDict  = dict(),    # dictionary output to be returned by simulation
               jobDuration = 1          # [second]
               ):
    
    if varDict is None:
      varDict = self.params; 
  
    # range is in [s]
    ts = np.linspace(0,jobDuration,jobDuration*1000)
    # initial value of dependent variable
    y0s = [0.62]
    ys = odeint(self.odefunc,y0s,ts,args=(varDict,))
  
    data        = dict() # This is a nuisance, but keeping backward compatible w Gotran stuff 
    data['t']   = ts          
    data['Cai'] = ys[:,0]
  
    returnDict['data'] = data       
    return returnDict
    
  

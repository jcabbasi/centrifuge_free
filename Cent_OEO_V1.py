# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 09:25:36 2021

@author: 2924746
"""
import numpy as np
from CentClasses import Fluid as fl
from CentClasses import Rock as rock
from CentClasses import RockFluid as RockFluid
from CentClasses import NumericalSetting 
from CentClasses import PropertyMatrix as PM
import CentSimCore as simcore
from CentSimCore import Simulate
from CentSimCore import simProps
from DataStorageClass import *
import matplotlib.pyplot as plt


# In this file we are going to convert our matlab code to a python code -- 13 Oct 2021:

#Fluids: Density kg/m3, Viscosity p
water=fl(1000,0.001)
oil=fl(700,0.001)
#Rock:
rk=rock(0.2,1e-14,0.05,0.0009)

#Numerical Setting:
# rock, Grid Number, Initial Timestep, ErrorLimit, Max Iterartion, Maximum Number of Timestep, Minimum Saturation Change
NS=NumericalSetting(rk,200,0.1,1e-6,15,100000,1e-14)

# Rock/Fluid is needed

rf=RockFluid(rk,water,oil,0.25,0.365)
Pmax=0.34e5; Pmin=-0.1e5
rf.initializePc(1.7,0.54,Pmax,Pmin,0,4.1,500)
rf.initializeRelperm(2,2,0.5,0.5)

swt=np.zeros([5,1])
swt[0]=1
pclist=rf.pc(swt)
#Initialization using property matrix9
Init=PM()
Init.initializeuAuto(rf,NS,0)

simprops=simProps(0.25, [60],900000) # inputs: inner radius (m), vector of omega (rad/s) , maximum simulation time (sec)
satexp=rf.SwatR2Calc(10,0.25)
satZeroPc=rf.PcZeroCalc()
installation='TEO'

OEOsim=Simulate(simprops,rf,Init,NS,installation,printLog=True) 
res=OEOsim.simulateNow()

Predicted=res.ExtractSimulationResults(['Time','RF'])
fig, ax2 = plt.subplots(ncols=2,nrows=1, sharey=False)
ax2[0].scatter(Predicted[:,0], Predicted[:,1])
hh=0



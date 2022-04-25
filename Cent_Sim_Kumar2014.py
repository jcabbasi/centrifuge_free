# -*- coding: utf-8 -*-

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
import pandas as pd

# 17 Nov 2021:
#Fluids: Density kg/m3, Viscosity p
water=fl(1000,0.00071)
oil=fl(700,0.00217)
#Rock:
rk=rock(0.20,2.44E-13,0.05,0.0009)

#Numerical Setting:
# rock, Grid Number, Initial Timestep, ErrorLimit, Max Iterartion, Maximum Number of Timestep, Minimum Saturation Change
NS=NumericalSetting(rk,40,0.001,1e-7,10,25000,1e-10,simpleGridding=False)
DXv=NS.DXv
# Rock/Fluid is needed

rf=RockFluid(rk,water,oil,0.38,0.30)

#Pc WW-Base
# Pmax=0.16e5; Pmin=0; Pmid=0.005e5
# rf.initializePc(2.0,10,Pmax,Pmin,Pmid,6,1)


#Pc MW
Pmax=0.16e5; Pmin=-0.03e5; Pmid=0.005e5
rf.initializePc(2.0,10,Pmax,Pmin,Pmid,7,0.4)

#Pc SWW
# Pmax=0.16e5; Pmin=0.02; Pmid=0.02e5
# rf.initializePc(2.0,20,Pmax,Pmin,Pmid,7,1)

rf.initializeRelperm(0.7,1.8,0.15,0.35)
# rf.initializeRelperm(2.028592563,	1.432719279	,0.176004231	,0.729625955)



swt=np.zeros([5,1])
swt[0]=1
pclist=rf.pc(swt)
#Initialization using property matrix9
Init=PM()
Init.initializeuAuto(rf,NS,0)

simprops=simProps(0.25, [60],20000000) # inputs: inner radius (m), vector of omega (rad/s) , maximum simulation time (sec)
satexp=rf.SwatR2Calc(1,0.25)
satZeroPc=rf.PcZeroCalc()
installation='TEO'

OEOsim=Simulate(simprops,rf,Init,NS,installation,printLog=True) 
pd.to_pickle(OEOsim,'saveTEO_MW40.pickle')
res=OEOsim.simulateNow(FileName='OEO-MW0612_40',InteractivePlot=True,reportPeriod=5)

Predicted=res.ExtractSimulationResults(['Time','RF'])
fig, ax2 = plt.subplots(ncols=2,nrows=1, sharey=False)
ax2[0].scatter(Predicted[:,0], Predicted[:,1])
chh=0



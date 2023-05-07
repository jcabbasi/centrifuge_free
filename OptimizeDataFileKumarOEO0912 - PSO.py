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
from CentOptimizerClass import *


# 12 Dec 2021:
#Fluids: Density kg/m3, Viscosity p
water=fl(1000,0.00071)
oil=fl(700,0.00217*0.1)
#Rock:
rk=rock(0.20,2.44E-13,0.05,0.0009)

#Numerical Setting:
# rock, Grid Number, Initial Timestep, ErrorLimit, Max Iterartion, Maximum Number of Timestep, Minimum Saturation Change
NS=NumericalSetting(rk,10,0.05,1e-7,10,300,1e-9,simpleGridding=False)
DXv=NS.DXv
# Rock/Fluid is needed

rf=RockFluid(rk,water,oil,0.38,0.30)


#Pc MW
Pmax=0.16e5; Pmin=-0.03e5; Pmid=0.005e5
rf.initializePc(2.0,10,Pmax,Pmin,Pmid,7,0.4)


rf.initializeRelperm(0.7,1.8,0.15,0.35)



swt=np.zeros([5,1])
swt[0]=1
pclist=rf.pc(swt)
#Initialization using property matrix9
Init=PM()
Init.initializeuAuto(rf,NS,0)

simprops=simProps(0.25, [60],6500) # inputs: inner radius (m), vector of omega (rad/s) , maximum simulation time (sec)
satexp=rf.SwatR2Calc(1,0.25)
satZeroPc=rf.PcZeroCalc()
installation='OEO'
OEOsim=Simulate(simprops,rf,Init,NS,installation,printLog=False) 
installation='TEO'
TEOsim=Simulate(simprops,rf,Init,NS,installation,printLog=False) 

#Observed Data
filename="\\ObsData.xlsx"
observed=LoadandSave.LoadDataFromExcel(filename,'B5','C14')
observedDual=LoadandSave.LoadDataFromExcel(filename,'F5','G14')

opt=PSOoptimizer(OEOsim,rf,observed,Dual=True,obsDual=observedDual,simClassDual=TEOsim)

opt.OptimizeNowPSO(simTime=[1600,1600])




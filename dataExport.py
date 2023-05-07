import numpy as np
from CentClasses import RockFluid as RockFluid
from DataStorageClass import *
import matplotlib.pyplot as plt
import pandas as pd


simData1 = pd.read_pickle("./SeAna__muo__0.00217.pkl")
simData2 = pd.read_pickle("./SeAna__muo__0.00021700000000000002.pkl")

RFlist=[0.03,0.1,0.15,0.2,0.9]

XvD1=simData1.loc[2,'Sw'][:,0]/np.max(simData1.loc[2,'Sw'][:,0])
XvD2=simData2.loc[2,'Sw'][:,0]/np.max(simData2.loc[2,'Sw'][:,0])

npSw1=np.reshape(simData1.loc[2,'Sw'][:,1], [len(XvD1),1]) 
npP1=np.reshape(simData1.loc[2,'P'][:,1], [len(XvD1),1])  
npSw2=np.reshape(simData2.loc[2,'Sw'][:,1], [len(XvD1),1]) 
npP2=np.reshape(simData2.loc[2,'P'][:,1], [len(XvD1),1]) 
rfval=0
for ii in simData1.index:
    if simData1.loc[ii,'RF']>RFlist[rfval]:
        print(simData1.loc[ii,'RF'])
        npSw1=np.append(npSw1,  np.reshape(simData1.loc[ii,'Sw'][:,1],[len(XvD1),1])    ,axis=1)
        npP1=np.append(npP1,  np.reshape(simData1.loc[ii,'P'][:,1], [len(XvD1),1]) ,axis=1)
        rfval+=1
npSw1=np.append(npSw1,np.reshape(simData1.loc[ii,'Sw'][:,1],[len(XvD1),1]),axis=1)
npP1=np.append(npP1,np.reshape(simData1.loc[ii,'P'][:,1],[len(XvD1),1]),axis=1)

rfval=0
for ii in simData2.index:
    if simData2.loc[ii,'RF']>RFlist[rfval]:
        print(simData2.loc[ii,'RF'])
        npSw2=np.append(npSw2,np.reshape(simData2.loc[ii,'Sw'][:,1],[len(XvD1),1]),axis=1)
        npP2=np.append(npP2,np.reshape(simData2.loc[ii,'P'][:,1],[len(XvD1),1]),axis=1)
        rfval+=1
npSw2=np.append(npSw2,np.reshape(simData2.loc[ii,'Sw'][:,1],[len(XvD2),1]),axis=1)
npP2 =np.append(npP2,np.reshape(simData2.loc[ii,'P'][:,1],[len(XvD2),1]),axis=1)
   
fig= plt.figure(figsize=(10, 4), dpi=150)
ax = fig.subplots(1,2)
ax[0].set_xlabel('Dimensionless Length');ax[0].set_ylabel('Water Saturation') 
ax[0].set_xlim([0, 1])
pass
xx=11
for jj in range( len(RFlist) ):
    pass
    if jj==(len(RFlist)-1):
        ax[0].plot(np.reshape(XvD1,[len(XvD1),1]),npSw1[:,jj+1] ,'b-')
        ax[0].plot(np.reshape(XvD2,[len(XvD2),1]),npSw2[:,jj+1],'r--' )
        # ax[1].plot(np.linspace(0,1,30),npP1[:,jj+1] ,'b-')
        # ax[1].plot(np.linspace(0,1,30),npP2[:,jj+1],'r--' )
        pass
    else:
        ax[0].plot(np.reshape(XvD1,[len(XvD1),1]),npSw1[:,jj+1] ,'b-')
        ax[0].plot(np.reshape(XvD2,[len(XvD2),1]),npSw2[:,jj+1],'r--' )
        # ax[1].plot(np.linspace(0,1,30),npP1[:,jj+1] ,'b-')
        # ax[1].plot(np.linspace(0,1,30),npP2[:,jj+1],'r--' ) 
        pass           


ax[1].plot(simData1.loc[:,'Time'],simData1.loc[:,'RF'],'b-',label= 'muo=0.002 cp')
ax[1].plot(simData2.loc[:,'Time'],simData2.loc[:,'RF'],'r--',label= 'muo=0.0002 cp' )
ax[1].set_xlabel('Time (s)');ax[1].set_ylabel('Recovery Factor') 
ax[1].set_xscale('log');ax[1].set_yscale('log') 
ax[1].legend(loc="lower right")


print(npSw1)
     
     
   
        
    
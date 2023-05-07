import numpy as np
from DataStorageClass import *
import pandas as pd
from matplotlib import pyplot as plt
import tensorflow as tf

swc=0.38
snwr=0.3
muw=0.00071
munw=0.00217
filename="\\ObsData.xlsx"
observedOEO=LoadandSave.LoadDataFromExcel(filename,'B5','C113')
observedTEO=LoadandSave.LoadDataFromExcel(filename,'F5','G210')

def krw_f(Sw,coreyPara,swc,snwr):      
         krw=coreyPara[0]* np.power(((Sw-swc)/(1-swc-snwr)),coreyPara[2])
         return krw

  
def krnw_f(Sw,coreyPara,swc,snwr):  
        S_star=(1-Sw-snwr)/(1-swc-snwr)
        krnw=coreyPara[1]* np.power( S_star,coreyPara[3])
        return krnw

def fractionalFlow(krw,krnw,muw,munw):
        f=(1/(1+krnw*muw/krw/munw))
        return f
def lossf(y_true, y_pred):

      loss = tf.math.squared_difference(y_pred, y_true)**1
      loss_f=tf.reduce_mean(loss)**1
      return loss_f, loss

optdata = pd.read_pickle("optData_muochanged_v2.pkl")
# optdata = pd.read_pickle("optData_muochanged.pkl")
# optdata = pd.read_pickle("optData_muwchanged_v2.pkl")
# optdata = pd.read_pickle("optData_base_v2.pkl")

optdataSize=optdata.shape
selected=np.zeros([])
originCorey=[0.15,0.35,0.7,1.8]

sw=np.linspace(swc,1-snwr,100)

krw_base=krw_f(sw,originCorey,swc,snwr)
krnw_base=krnw_f(sw,originCorey,swc,snwr)
f_base=fractionalFlow(krw_base,krnw_base,muw,munw)

fig, ax = plt.subplots(ncols=2,nrows=2, sharey=False)
ax[0,0].plot(sw,krw_base,'--',color='black',label='Actual')
ax[0,0].plot(sw,krnw_base,'--',color='black')
ax[1,1].plot(observedTEO[:,0],observedTEO[:,1],'o',color='black',label='Actual')
# ax[0,1].plot(sw,f_base,'--',color='black',label='Actual')
ax[1,0].plot(observedOEO[:,0],observedOEO[:,1],'o',color='black',label='Actual')

ax[0,0].set_xlabel('sw');ax[0,0].set_ylabel('kr')
ax[1,0].set_xlabel('time (s)');ax[1,0].set_ylabel('RF')
ax[1,1].set_xlabel('time (s)');ax[1,1].set_ylabel('RF')
ax[0,1].set_xlabel('sw');ax[0,1].set_ylabel('MSE')
labels=['Both','OEO','TEO']
case=['ErrorVal','ErrorVal1','ErrorVal2']
col=['red','blue','green']
bg=0
end=bg+1
for jj in [0,1,2]:
        krw=0;krnw=0
        for pp in optdata.index:
                nnn=optdata.loc[pp,case[jj]]          
                if np.isnan(nnn):
                     optdata.loc[pp,case[jj]] =99999   
                
                
        optdata.sort_values(by=case[jj],inplace=True,ignore_index=True)
        # print(optdata[[case[jj]]][10:30])
        for ii in range(bg,end):
                
                coreyP=np.zeros([4])
                
                selected=np.append(selected,ii)
                coreyP[0]= optdata.loc[ii,'KrWmax']
                coreyP[1]= optdata.loc[ii,'KrNWmax']
                coreyP[2]= optdata.loc[ii,'nW']
                coreyP[3]= optdata.loc[ii,'nNW']
                errTotal=optdata.loc[ii,'ErrorVal']
                errOEO=optdata.loc[ii,'ErrorVal1']
                errTEO=optdata.loc[ii,'ErrorVal2']
                if ii==bg:
                        krw=(krw_f(sw,coreyP,swc,snwr))/1
                        krnw=(krnw_f(sw,coreyP,swc,snwr))/1                        
                else:
                        krw=(krw_f(sw,coreyP,swc,snwr)+krw)/2
                        krnw=(krnw_f(sw,coreyP,swc,snwr)+krnw)/2
                
        ax[0,0].plot(sw,krw,color=col[jj],label='Sim-'+labels[jj])
        ax[0,0].plot(sw,krnw,color=col[jj]) 
        ax[0,0].legend(loc="upper right")
        lossw,losskrw=lossf(krw, krw_base); lossnw,losskrnw=lossf(krnw, krnw_base)
        loss = 0.5*lossw+0.5*lossnw
        print(coreyP)
        print('++---  Loss Summary  ---++, errorT: {0:.3e}  OEO: {1:.3e}; TEO:{2:.3e}  '.format(errTotal, errOEO,errTEO  )   )
        print('++---  Loss Summary  ---++, krw: {0:.3e}  krnw: {1:.3e}; total:{2:.3e}  '.format(lossw, lossnw,loss  )   )
        
        rfc=optdata.loc[ii,'ProdCurve']
        ax[1,0].plot(rfc[:,0],rfc[:,1],label='Sim-'+labels[jj], color=col[jj])
        # ax[0,1].plot(sw,fractionalFlow(krw,krnw,muw,munw),'-', color='red',label='Simulation')
        ax[0,1].plot(sw,losskrw,'-', color=col[jj],label='krw-'+labels[jj])
        ax[0,1].plot(sw,losskrnw,'--', color=col[jj],label='krnw-'+labels[jj])
        ax[0,1].legend(loc="upper right")
        ax[1,0].legend(loc="lower right")
        rfcD=optdata.loc[ii,'ProdCurveDual']
        ax[1,1].plot(rfcD[:,0],rfcD[:,1],label='Sim-'+labels[jj], color=col[jj])
        ax[1,1].legend(loc="lower right")
                

  
  


  
                
fin=1                             


# for jj in [0,1,2]:

#         optdata.sort_values(by=case[jj],inplace=True,ignore_index=True)
#         for ii in range(1,5):
                
#                 coreyP=np.zeros([4])
                
#                 selected=np.append(selected,ii)
#                 coreyP[0]= optdata.loc[ii,'KrWmax']
#                 coreyP[1]= optdata.loc[ii,'KrNWmax']
#                 coreyP[2]= optdata.loc[ii,'nW']
#                 coreyP[3]= optdata.loc[ii,'nNW']
#                 errTotal=optdata.loc[ii,'ErrorVal']
#                 errOEO=optdata.loc[ii,'ErrorVal1']
#                 errTEO=optdata.loc[ii,'ErrorVal2']
#                 krw=krw_f(sw,coreyP,swc,snwr)
#                 krnw=krnw_f(sw,coreyP,swc,snwr)
#                 ax[0,0].plot(sw,krw,color=col[jj],label='Sim-'+labels[jj])
#                 ax[0,0].plot(sw,krnw,color=col[jj]) 
#                 ax[0,0].legend(loc="upper right")
#                 lossw,losskrw=lossf(krw, krw_base); lossnw,losskrnw=lossf(krnw, krnw_base)
#                 loss = 0.5*lossw+0.5*lossnw
#                 print(coreyP)
#                 print('++---  Loss Summary  ---++, errorT: {0:.3e}  OEO: {1:.3e}; TEO:{2:.3e}  \n'.format(errTotal, errOEO,errTEO  )   )
#                 print('++---  Loss Summary  ---++, krw: {0:.3e}  krnw: {1:.3e}; total:{2:.3e}  \n'.format(lossw, lossnw,loss  )   )
                
#                 rfc=optdata.loc[ii,'ProdCurve']
#                 ax[1,0].plot(rfc[:,0],rfc[:,1],label='Sim-'+labels[jj], color=col[jj])
#                 # ax[0,1].plot(sw,fractionalFlow(krw,krnw,muw,munw),'-', color='red',label='Simulation')
#                 ax[0,1].plot(sw,losskrw,'-', color=col[jj],label='krw-'+labels[jj])
#                 ax[0,1].plot(sw,losskrnw,'--', color=col[jj],label='krnw-'+labels[jj])
#                 ax[0,1].legend(loc="upper right")
#                 ax[1,0].legend(loc="lower right")
#                 rfcD=optdata.loc[ii,'ProdCurveDual']
#                 ax[1,1].plot(rfcD[:,0],rfcD[:,1],label='Sim-'+labels[jj], color=col[jj])
#                 ax[1,1].legend(loc="lower right")
                        







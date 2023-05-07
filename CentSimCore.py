# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 14:28:03 2021

@author: 2924746
"""

import numpy as np
import CentClasses as cc
from DataStorageClass import *
from FlowFncs import FlowFnc
import time as tm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from visualizationClass import *
global cutTS
import sys

class simProps():
    def __init__(self, r1, omegaVect, maxSimTime ):
        # inputs: inner radius (m), vector of omega (rad/s) , maximum simulation time (sec)
        self.r1=r1 #in meter
        self.omegaVect=omegaVect # Vector of simulation stages (rotational speeds) - rad/s
        self.maxSimTime=maxSimTime # Maximum simulation time at each simulation stage (rotation speed) - in s
        self.currentOmega=None
        

class Simulate():
    
    def __init__(self, simProps, rf, init, NSet,installation,printLog=True):
        self.rf=rf
        self.init=init
        self.NSet=NSet
        self.installation=installation #OEO or TEO#
        self.simprops=simProps
        self.Data=DataStorage()
        self.RepeatNoLowError=0
        self.__LastIterNo=0
        self.__DTnow=None
        self.__SatChangeRate=None
        self.printLog=printLog
        self.__stepNo=None
        self.cutTS=False
        self.intercativePlot=intercativeClass(NSet,simProps,rf)
    
    def refreshModel(self):
        self.Data=DataStorage()          
        self.RepeatNoLowError=0
        self.__LastIterNo=0
        self.__DTnow=None
        self.__SatChangeRate=None 

    
    def simulateNow(self,FileName='simResults',InteractivePlot=False,reportPeriod=None,ExportbySingle=True, simTime=None):
        start_time=tm.time()
        Pold=self.init.P.copy()
        Sold=self.init.Sw.copy()
        time=0
        gn=self.NSet.gn #Grid Number
        self.Data.assign(self.simprops.currentOmega,time,np.ones([gn,1]),np.ones([gn,1]),np.ones([gn,1]),0,0)
        self.__stepNo=0
        self.simTime=simTime
        self.intercativePlot.InteractivePlotInititialize(InteractivePlot)
        self.EquilData=EquilDataStorage(gn,self.simprops.omegaVect)
        stageNo=0; objTime=0
        stageSimTime=0
        
        if self.printLog==True:
            print('                                                                                         ') 
            print('                                                                                         ') 
            print('                                                                                         ') 
            print('-----------------------------------------------------------------------------------------') 
            print('----------------------------<   Simulation Is Started   >--------------------------------') 
            print('-----------------------------------------------------------------------------------------') 
            print('                                                                                         ') 

        for omega in self.simprops.omegaVect:
            stageNo+=1
            self.simprops.currentOmega=omega
            self.__SatChangeRate=np.Inf
            StageTime=0
            self.__DTnow=None
            self.__stepNo=0
            DS_guess=-0.02*np.ones([gn,1]) # saturation change
            DP_guess=-100*np.ones([gn,1])    # pressure change in Pa           
            self.__stepNo=0
            AnalyticalRF,AnalyticalSat=self.AnalyticalRFCalc(omega)
            self.intercativePlot.InteractivePlotInit(AnalyticalSat,Sold)
            thSw=self.SwDuetoRotationATx(omega,(self.simprops.r1+self.NSet.rockLength-self.NSet.DXv[-1]/2))
            if self.printLog==True:    
                #TIMESTEP CALCULATOR
                print('-----------------------------<  Changing Rotational Speed  >-----------------------------') 
            
            equil=False
            if simTime!=None:
                objTime += simTime[stageNo-1]
            DT=self.NSet.DTinit
            while equil==False:
                #Calculate DT
                DTold=DT
                DT=self.timestepCalc() 
                if ((simTime!=None) and ((time+DT)>objTime)):
                    DT=objTime-time
                Conv=False
                iter=0
                
                DS_guess*=(DT/DTold)
                # DS_guess=-0.02*np.ones([gn,1]) # saturation change
                # DP_guess=-100*np.ones([gn,1])    # pressure change in Pa
                
                
                Errold=0
                self.__stepNo=self.__stepNo+1
                cutTSno=0
                # ITERATION CALCULATOR
                while (iter < self.NSet.iterMax) and (Conv==False) or iter < 2:  #set both error and number of iters
                    iter+=1
                    
                    if (self.rf.cutTS==True and cutTSno<2):
                        iter-=1 
                        # DT=DT/2
                        # self.__DTnow=DT
                        # DS_guess=-0.02*np.ones([gn,1]) # saturation change
                        # DP_guess=-100*np.ones([gn,1])    # pressure change in Pa                        
                        print('-----------------------------<  Cutting TimeStep  >-----------------------------') 
                        self.rf.cutTS=False
                        cutTSno+=1
                        
                         
                                        #Check for error
                    for ii in range(gn):
                        if (Sold[ii]+DS_guess[ii])>(1-self.rf.snwr):
                            DS_guess[ii]=1-self.rf.snwr-Sold[ii]
                        if (Sold[ii]+DS_guess[ii])<(self.rf.swc):
                            DS_guess[ii]=self.rf.swc -Sold[ii]

                    for ii in range(gn):
                        if (Sold[ii])>(1-self.rf.snwr):
                            Sold[ii]=1-self.rf.snwr
                        if (Sold[ii])<(self.rf.swc):
                            Sold[ii]=self.rf.swc

                    for ii in range(0,gn):
                        if abs(DS_guess[ii])<1e-10:
                            DS_guess[ii]=-1e-10
                        if abs(DP_guess[ii])<1e-3:
                            DP_guess[ii]=-1    
                                                   
                     # create matrix
                    f1SN = np.zeros([gn,1])
                    f1SPO= np.zeros([gn,1])
                    f1PN = np.zeros([gn,1])
                    f2SN = np.zeros([gn,1])
                    f2SPO= np.zeros([gn,1])
                    f2PN = np.zeros([gn,1])  
                    f1PnSn = np.zeros([gn,1])  
                    f2PnSn = np.zeros([gn,1])  
                    
                    # Fi when DS/DP i-1 or i+1
                    f1Sn_neg1=   np.zeros([gn,1]) 
                    f1Pn_neg1=   np.zeros([gn,1]) 
                    f1SoPo_neg1= np.zeros([gn,1])                     
                    f2Sn_neg1=   np.zeros([gn,1]) 
                    f2Pn_neg1=   np.zeros([gn,1]) 
                    f2SoPo_neg1= np.zeros([gn,1]) 
                    f1Sn_pos1=   np.zeros([gn,1]) 
                    f1Pn_pos1=   np.zeros([gn,1])
                    f1SoPo_pos1= np.zeros([gn,1])                     
                    f2Sn_pos1=   np.zeros([gn,1]) 
                    f2Pn_pos1=   np.zeros([gn,1])                     
                    f2SoPo_pos1= np.zeros([gn,1])                     

                    ff=FlowFnc()
                    ff.initialize(self.NSet.DX,self.NSet.DXv,gn)
                                                                        
                    for ii in range(0,gn):
                        Splus=Sold.copy()                                           
                        Splus[ii]=Sold[ii]+DS_guess[ii]
                        Pplus=Pold.copy()
                        Pplus[ii]=Pold[ii]+DP_guess[ii]   
                          
                        #No Change
                        ff.buildIt(self.installation,self.rf,self.simprops,Pold,Sold,Sold,DT)
                        f1SPO[ii]=ff.F1([ii]) ; f2SPO[ii]= ff.F2([ii])                        
                        #DSw
                        ff.buildIt(self.installation,self.rf,self.simprops,Pold,Splus,Sold,DT)
                        f1SN[ii]=ff.F1([ii]) ; f2SN[ii]= ff.F2([ii])
                        #DPw
                        ff.buildIt(self.installation,self.rf,self.simprops,Pplus,Sold,Sold,DT)
                        f1PN[ii]=ff.F1([ii]) ; f2PN[ii]= ff.F2([ii])                         
                        #DPw & DSw
                        ff.buildIt(self.installation,self.rf,self.simprops,Pplus,Splus,Splus,DT)
                        f1PnSn[ii]=ff.F1([ii]) ; f2PnSn[ii]= ff.F2([ii])                          
                                                
                        
                                               
                        if ii != 0:
                            Splus=Sold.copy()                                           
                            Splus[ii-1]=Sold[ii-1]+DS_guess[ii-1]
                            Pplus=Pold.copy()
                            Pplus[ii-1]=Pold[ii-1]+DP_guess[ii-1]
                            
                            #No Change
                            ff.buildIt(self.installation,self.rf,self.simprops,Pold,Sold,Sold,DT)
                            f1SoPo_neg1[ii]=ff.F1([ii]) ; f2SoPo_neg1[ii]= ff.F2([ii])                        
                            #DSw
                            ff.buildIt(self.installation,self.rf,self.simprops,Pold,Splus,Sold,DT)
                            f1Sn_neg1[ii]=ff.F1([ii]) ; f2Sn_neg1[ii]= ff.F2([ii])
                            #DPw
                            ff.buildIt(self.installation,self.rf,self.simprops,Pplus,Sold,Sold,DT)
                            f1Pn_neg1[ii]=ff.F1([ii]) ; f2Pn_neg1[ii]= ff.F2([ii])                         
                            #DPw & DSw                        

                        if ii != (gn-1):
                            Splus=Sold.copy()                                           
                            Splus[ii+1]=Sold[ii+1]+DS_guess[ii+1]
                            Pplus=Pold.copy()
                            Pplus[ii+1]=Pold[ii+1]+DP_guess[ii+1] 
                            #
                            #No Change 
                            ff.buildIt(self.installation,self.rf,self.simprops,Pold,Sold,Sold,DT)
                            f1SoPo_pos1[ii]=ff.F1([ii]) ; f2SoPo_pos1[ii]= ff.F2([ii])                        
                            #DSw
                            ff.buildIt(self.installation,self.rf,self.simprops,Pold,Splus,Sold,DT)
                            f1Sn_pos1[ii]=ff.F1([ii]) ; f2Sn_pos1[ii]= ff.F2([ii])
                            #DPw
                            ff.buildIt(self.installation,self.rf,self.simprops,Pplus,Sold,Sold,DT)
                            f1Pn_pos1[ii]=ff.F1([ii]) ; f2Pn_pos1[ii]= ff.F2([ii])                             
                            
              
                                
                    AA=np.zeros([2*gn,2*gn])
                    BB=np.zeros([2*gn,1])
                    

                    
                    for ii in range(0,gn):
                        
                        if ii != 0:
                            AA[ii,ii-1]      =  ( f1Sn_neg1[ii] - f1SoPo_neg1[ii] ) /  ( DS_guess[ii-1] ) 
                            AA[ii,gn+ii-1]   =  ( f1Pn_neg1[ii] - f1SoPo_neg1[ii] ) /  ( DP_guess[ii-1] ) 
                            AA[gn+ii,ii-1]   =  ( f2Sn_neg1[ii] - f2SoPo_neg1[ii] ) /  ( DS_guess[ii-1] )
                            AA[gn+ii,gn+ii-1]=  ( f2Pn_neg1[ii] - f2SoPo_neg1[ii] ) /  ( DP_guess[ii-1] )

                        
                        if ii != (gn-1):
                            AA[ii,ii+1]      = ( f1Sn_pos1[ii] - f1SoPo_pos1[ii] ) /  ( DS_guess[ii+1] )
                            AA[ii,gn+ii+1]   = ( f1Pn_pos1[ii] - f1SoPo_pos1[ii] ) /  ( DP_guess[ii+1] )                            
                            AA[gn+ii,ii+1]   = ( f2Sn_pos1[ii] - f2SoPo_pos1[ii] ) /  ( DS_guess[ii+1] )
                            AA[gn+ii,gn+ii+1]= ( f2Pn_pos1[ii] - f2SoPo_pos1[ii] ) /  ( DP_guess[ii+1] )                            
                                 
                        
                        
                        AA[ii,ii]      =  ( f1SN[ii] - f1SPO[ii] )   /  ( DS_guess[ii] ) 
                        AA[ii,gn+ii]   =  ( f1PN[ii] - f1SPO[ii] )   /  ( DP_guess[ii] )                           
                        AA[gn+ii,ii]   =  ( f2SN[ii] - f2SPO[ii] )   /  ( DS_guess[ii] )
                        AA[gn+ii,gn+ii]=  ( f2PN[ii] - f2SPO[ii] )   /  ( DP_guess[ii] )                           
                         
                        BB[   ii,0]  =  -  (  f1SPO[ii]  )                      
                        BB[gn+ii,0]  =  -  (  f2SPO[ii]  )                      
                        # BB[   ii,0]  =  -  (  f1PnSn[ii]  )                      
                        # BB[gn+ii,0]  =  -  (  f2PnSn[ii]  )    
                
                    SingArrCol=self.CheckForSingularityCol(AA,BB)
                    
                    # Remove singular rows 
                    AAclc= np.delete(AA.copy(), SingArrCol.astype(int),0)
                    AAclc= np.delete(AAclc, SingArrCol.astype(int),1)
                    BBclc= np.delete(BB.copy(), SingArrCol.astype(int),0)                                                
                    
                    # calculate
                    invAAclc=   np.linalg.inv(AAclc)
                    XXmatclc=   np.matmul(invAAclc, BBclc)
                    #invAA=   np.linalg.inv(AA)
                    # XXmat=   np.matmul(invAA, BB)                    
                    SnewIt=  Sold.copy()
                    PnewIt=  Pold.copy()
                    PMprev=  self.init.copy()
                    PMprev.updateProps(Pold,Sold)
                    DS_guessN=DS_guess.copy()
                    DP_guessN=DP_guess.copy()
                    jj=0
                    for ii in np.delete(range(0,2*gn), SingArrCol):
                            
                        if ii<gn:
                            SnewIt[ii]=Sold[ii] + XXmatclc[jj]
                            DS_guessN[ii]  =  XXmatclc[jj]            
                        else:
                            PnewIt[ii-gn]=Pold[ii-gn] + XXmatclc[jj]
                            DP_guessN[ii-gn]  =  XXmatclc[jj]                                                     
                        jj=jj+1
                    # print(PnewIt)
                    # print(SnewIt)
                    PMcurrIt=self.init.copy()
                    PMcurrIt.updateProps(PnewIt,SnewIt)
                    SatChngnew=self.SatChangeCalculator(PMcurrIt,PMprev,DT)
                    Conv,ErrDiff=self.CheckForConv(SatChngnew,Errold,iter)
                    Errold=SatChngnew
                    if (self.rf.cutTS==True and cutTSno<2):
                        Conv=False
                    
                    if np.isnan(SnewIt[1]):
                        wait=0
                        break
                    
                    DS_guess=DS_guessN.copy()*0.92
                    DP_guess=DP_guessN.copy()*0.92
                    if self.installation=='TEO':
                        PnewIt[0]=self.PnwEqDuetoRotationMAX(omega)-self.rf.pc(Sold[0])

                    #print(SnewIt)
                ###############################################################

                self.__LastIterNo=iter
                self.__SatChangeRate=abs(SatChngnew)
                
                StageTime+=DT; time+=DT 
                
                Sold=SnewIt; Pold=PnewIt
                # Calculate RF:
                RF=self.RFcalc(self.rf,self.init.Sw,SnewIt)
                Pc=self.rf.pc(SnewIt)
                
                
                equil=self.CheckEquil(StageTime,time,objTime)                
                prstr=('--- Step: {7}; Om: {0} rad/s, Time: {1:.2f} s; TS: {2:.2f} s; RF: {5:.4f}, SatChange: {6:.2e}; Errs: {3:.2e} ; Iter: {4} '.format(omega, time, DT, ErrDiff,iter, RF,self.__SatChangeRate,self.__stepNo) )   

                if self.printLog==True:
                    print(prstr)   
                else:
                    print( prstr+ "\r", end='', flush=True)
                    sys.stdout.flush()
                
                if reportPeriod!=None:
                    if np.remainder(self.__stepNo,reportPeriod)==0 or Conv==True:
                        
                        self.Data.assign(self.simprops.currentOmega,time,np.concatenate((self.NSet.Xv,PnewIt),axis=1),np.concatenate((self.NSet.Xv,SnewIt),axis=1),Pc,RF,iter)


                strPlt=  (' Time: {0:.2f} s; \n Timestep: {1:.2f} s; \n RF: {2:.2f} %;  \n IterNo: {3}; \n SatChange: {4:.3e}'.format(time,DT, RF*100,iter, self.__SatChangeRate )   )

                self.intercativePlot.updateInteractivePlot(InteractivePlot,self.__stepNo,Sold,time,RF,ErrDiff, strPlt)
            
            if self.printLog==True:    
                print('++++ Rotation Stage Finished ++++, Omega: {0} rad/s, Time: {1:.2f} s; RF: {2:.2f}; PcMax: {3:.2f} psi; SwMin: {4:.2f} \n'.format(omega, time, RF, self.rf.pc(SnewIt).max(),SnewIt.min() )   )
                print('+++++++: Step: {7}; Time: {1:.2f} s; TS: {2:.2f} s; RF: {5:.4f}, SatChange: {6:.2e}; Errs: {3:.2e} ; Iter: {4} '.format(omega, time, DT, ErrDiff,iter, RF,self.__SatChangeRate,self.__stepNo) )   

            
            self.EquilData.assign(self.simprops.currentOmega,time,StageTime,AnalyticalRF,PnewIt,SnewIt,Pc,RF,AnalyticalSat)
            

            
            ## End of RotationStage

        if self.printLog==True:    
            print('--------  Simulation is Finished Successfully  -----------')
        endtime=tm.time()
        print('++++ Simulation is Finished ++++, Last Omega: {0} rad/s, Last Time: {1:.2f} s; RF: {2:.2f}; RunTime: {3:.2f} s; Time Vector: {4:.2f} \n'.format(omega, time, RF, abs(start_time-endtime),StageTime )   )
    
        if ExportbySingle==True:
            # Export Results
            self.Data.SavetoPickle(FileName)
            self.Data.SaveRFtoCSV(FileName)
            self.EquilData.SavetoPickle(FileName+"Equil")
            self.EquilData.SavetoCSV(FileName+"Equil")

        
        return   self.Data  
                        
                
                          
    
    
    def RFcalc(self,rockfluid,Swi,Swnow):
 
        RFact=np.mean( abs(Swi-Swnow) )/(abs(1-rockfluid.snwr))
        return RFact
 
    def CheckForConv(self,ErrNew,ErrOld,iter):
        ErrDiff=1000000
        Conv=False
        if iter>1:  
            
            if (abs(ErrNew)+abs(ErrOld))==0:
                Conv=True
                ErrDiff=0

            if (abs(ErrNew)+abs(ErrOld))>0:
                ErrDiff=abs(abs(ErrNew)-abs(ErrOld))/(abs(ErrNew)+abs(ErrOld))
                if ErrDiff<self.NSet.errorLimit:
                    Conv=True
                    
        return Conv, ErrDiff
           

    def PcEqDuetoRotationMAX(self,omega):
        r2=self.simprops.r1+self.rf.rock.length
        pceqMax= 0.5 * self.rf.dDen() * omega * omega * ( r2 * r2 - self.simprops.r1 * self.simprops.r1 )
        return pceqMax

    def PcEqDuetoRotation(self,omega,x):
        r2=self.simprops.r1+self.rf.rock.length        
        pceqMax= 0.5 * self.rf.dDen() * omega * omega * ( r2 * r2 - x * x )
        return pceqMax       

    def SwDuetoRotationATx(self,omega,x):
        r2=self.simprops.r1+self.rf.rock.length        
        pceqMax= 0.5 * self.rf.dDen() * omega * omega * ( r2 * r2 - x * x )
        satPc=self.rf.PcNonZeroCalc(pceqMax)
        
        return satPc      

    def PnwEqDuetoRotationMAX(self,omega):
        r2=self.simprops.r1+self.rf.rock.length
        pnwMax= - 0.5 * self.rf.fnw.den * omega * omega * ( r2 * r2 - self.simprops.r1 * self.simprops.r1 )
        return pnwMax      
    
    # This function calculate the expected reacovery factor using analytical relationship
    def AnalyticalRFCalc(self,omega):
        #Sat Profile
       
        xx=self.NSet.Xv+self.simprops.r1
        pceq=self.PcEqDuetoRotation(omega,xx)
        
        #Pc Curve
        satrange = np.linspace(self.rf.swc,1-self.rf.snwr, num=500, endpoint=True)
        pccurve=self.rf.pc(satrange) 
        f = interp1d(pccurve.squeeze(), satrange)
        
        #RF Calc
        satExp=f(pceq)
        RF=self.RFcalc(self.rf,self.init.Sw,satExp)
    
        return RF,satExp
    
    
    
    def timestepCalc(self):
        
        if self.__DTnow==None:
            self.__DTnow= self.NSet.DTcalc()    
        elif self.__LastIterNo<4:
            if self.__DTnow<self.NSet.DTmax:
                if self.__LastIterNo==2:
                    self.__DTnow=self.__DTnow*1.2
                else:
                    self.__DTnow=self.__DTnow*1.15

            else:
                pass
        elif self.__LastIterNo<5:
            pass
        else:
            if self.__DTnow>self.NSet.DTmin:
                self.__DTnow=self.__DTnow*0.9
            else:
                pass
        

        
        return self.__DTnow
            
    def CheckForSingularity(self,AA,BB):
        SingularRows=np.array([])
        for ii in range(0,len(BB)):
            
            if AA[ii,:].any()!=0:
                pass
            else:
                SingularRows=np.append(SingularRows, ii)
        SingularRows=SingularRows.astype(int)
        return SingularRows
    
    def CheckForSingularityCol(self,AA,BB):
        SingularRows=np.array([])
        for ii in range(0,len(BB)):
            
            if AA[:,ii].any()!=0:
                pass
            else:
                SingularRows=np.append(SingularRows, ii)
        SingularRows=SingularRows.astype(int)
        return SingularRows
    
        
    def SatChangeCalculator(self, StatusCurrIter,StatusPrevIter,Dt):
        
        #ErrorP= (StatusCurrIter.P-StatusPrevIter.P)/2
        ChangeS= (StatusCurrIter.Sw-StatusPrevIter.Sw)/1
        TotalChange=ChangeS
        ChangeMean=TotalChange.mean()/Dt
        return ChangeMean
    
    
    def updateRockFluid(self,RockFluid):
        self.rf=RockFluid
    
    def ExportResults(self):
        return self.Data
    
    def copy(self):
        newClass=Simulate(self.simprops,self.rf,self.init,self.NSet,self.installation,printLog=True)
        newClass.self=self
        return newClass
    
    def CheckEquil(self,StageTime,time,objTime):
        equil=False

                     
        if self.simTime!=None: 
            
            if (time+2.2e-6)>objTime:
                equil=True
            else:
                equil=False
        
        else:
            if StageTime > self.simprops.maxSimTime:
                equil=True
            
            elif abs(self.__SatChangeRate)<self.NSet.ProdRateMin: 
                equil=True
                
        if self.__stepNo>self.NSet.tsMax:
            equil=True            
            
            
        return equil

    
    
        
        
    
    
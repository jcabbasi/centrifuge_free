# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

global cutTS

import numpy as np
from scipy.interpolate import interp1d

class Fluid:
    
    def __init__(self,den,mu):
        self.den=den
        self.mu=mu
    
    
class Rock:
      def __init__(self,phi,k,length,area):         
        self.phi=phi
        self.k=k
        self.length=length
        self.area=area
        

class RockFluid:
    
      def __init__(self,rock,fW,fNW,swc,snwr):
        self.swc=swc
        self.snwr=snwr
        self.rock=rock
        self.fw=fW
        self.fnw=fNW
        self.cutTS=False
        
      def initializePc(self,n1,n2,Pcmax,Pcmin,Pcmid,k1,k2):
          self.n1=n1
          self.n2=n2
          self.a1=Pcmax
          self.a2=Pcmin
          self.a3=Pcmid
          self.k1=k1
          self.k2=k2
          
      def initializeRelperm(self,nW,nNW,krWMAX,krNWMAX):
          self.nW=nW
          self.nNW=nNW
          self.krWMAX=krWMAX
          self.krNWMAX=krNWMAX  
        
      def  pc(self,Sw):
     # From: Andersen et al. 2017, SPE-188291-MS
     
          Swstar=(Sw-self.swc)/(1-self.swc-self.snwr)
          aa1= ((self.a1-self.a3)*np.power((1+self.k2),self.n2)  +  (self.a3-self.a2))  *  np.power((1+self.k1),self.n1) / ( np.power(1+self.k1,self.n1)* np.power(1+self.k2,self.n2)-1 )
          aa2= ((self.a1-self.a3)+(self.a3-self.a2)* np.power(1+self.k1,self.n1))* np.power(1+self.k2,self.n2) / ( np.power(1+self.k1,self.n1) * np.power(1+self.k2,self.n2) - 1)       
          capP=np.zeros([len(Sw),1])
          for jjjj in range(0,len(Sw)):
                capP[jjjj,0]= aa1  /  np.power(  1+self.k1*Swstar[jjjj]  ,  self.n1) - aa2/np.power(1+self.k2*(1-Swstar[jjjj]),self.n2)  +  self.a3              
          return capP
      
      # Maximum water saturation:
      def Swmax(self):
        return 1-self.snwr

      def Swmin(self):
        return self.swc   

      def dDen(self):
        return np.abs(self.fw.den-self.fnw.den)    
    
      def krw_f(self,Sw):      
         krw=self.krWMAX* np.power(((Sw-self.swc)/(1-self.swc-self.snwr)),self.nW)
         return krw

  
      def krnw_f(self,Sw):  
        S_star=(1-Sw-self.snwr)/(1-self.swc-self.snwr)
        if S_star<0:
          if abs(S_star)>1e-4:
            #print('WARNING: Negative Sat Values@ W Phase')
            S_star=0            
            self.cutTS=True

          if abs(S_star)<1e-4:
            S_star=0
         
        krnw=self.krNWMAX* np.power( S_star,self.nNW)
        return krnw
     
      def lambdaW_f(self,Sw):
        lambdaW=np.zeros([len(Sw),1])
        for jjjj in range(0,len(Sw)):
            lambdaW[jjjj,0]=self.rock.k*self.krw_f(Sw[jjjj])/self.fw.mu
        return lambdaW
    
      def lambdaNW_f(self,Sw):
        if type(Sw) is int or type(Sw) is float:
          Sw=np.array([Sw])
        lambdaNW=np.zeros([len(Sw),1])
        for jjjj in range(0,len(Sw)):
            lambdaNW[jjjj,0]=self.rock.k*self.krnw_f(Sw[jjjj])/self.fnw.mu
        return lambdaNW   
      
      
      def SwatR2Calc(self,omega,r1):
        r2=r1+self.rock.length        
        pcatboundary= 0.5 * self.dDen() * omega * omega * ( r2 * r2 - r1 * r1 ) 
        # satrange= range(self.swc,1-self.snwr,0.001) 
        satrange = np.linspace(self.swc,1-self.snwr, num=100, endpoint=True)
        pccurve=self.pc(satrange) 
        f = interp1d(pccurve.squeeze(), satrange)
        return f(pcatboundary)
        
      def PcZeroCalc(self):
        pcwanted=0 
        satrange = np.linspace(self.swc,1-self.snwr, num=100, endpoint=True)
        pccurve=self.pc(satrange) 
        f = interp1d(pccurve.squeeze(), satrange)
        
        try:
          zeropc= f(pcwanted)
        except:
          zeropc=np.max(satrange)
          #print('AN EXCEPTION OCCURED')
        return zeropc

      def PcNonZeroCalc(self,pcwanted):
        
        satrange = np.linspace(self.swc,1-self.snwr, num=100, endpoint=True)
        pccurve=self.pc(satrange) 
        f = interp1d(pccurve.squeeze(), satrange)
        
        try:
          zeropc= f(pcwanted)
        except:
          zeropc=np.max(satrange)
          print('AN EXCEPTION OCCURED')
        return zeropc
      
      # Quality control and modification of Corey rel perm variables: used for GA optimization
      def QCRelPermCoefs(self,nW,nNW,KrWmax,KrNWmax):
          
          if nW<0:
            nW=0  
          if nNW<0:
            nNW=0  
           
          if KrWmax<0:
            KrWmax=0
          if KrNWmax<0:
            KrNWmax=0
          
          if KrWmax>1:
            KrWmax=1
          if KrNWmax>1:
            KrNWmax=1           
          return   nW,nNW,KrWmax,KrNWmax 



class NumericalSetting:
       # rock, Grid Number, Initial Timestep, ErrorLimit, Max Iterartion, Maximum Number of Timestep, Minimum Saturation Change
       def __init__(self,rock,gn,DTinit,errorLimit,iterMax,tsMax,ProdRateMin,simpleGridding=True):
        self.gn=gn  #Grid Number
        self.DTinit=DTinit # Initial Time Step
        self.errorLimit=errorLimit
        self.ProdRateMin=ProdRateMin
        self.iterMax=iterMax
        self.tsMax=tsMax
        self.DX=rock.length/self.gn  #Grid Size in m
        self.rockLength=rock.length
        self.DT=None
        self.DTmax=10000
        self.DTmin=0.001
        self.DXv=None
        self.ExpFactor=0
        self.minCellLength=rock.length/self.gn/5
        self.simpleGridding=simpleGridding
        self.DXvSet()
        
       def DTcalc(self):
           DTret=self.DTinit
           return DTret
         
       def DXvSet(self):
        self.DXv=np.zeros([self.gn,1])  
        
        if self.simpleGridding==True:
          self.DXv=(self.rockLength/self.gn)*np.ones([self.gn,1])
  
        else:
          dist=np.exp(np.linspace(1,0,num=(self.gn-1))*self.ExpFactor)
          self.DXv[0:self.gn-1]=np.reshape( dist/np.sum(dist)*(self.rockLength-self.minCellLength) , [self.gn-1,1] )        
          self.DXv[-1]=self.minCellLength          
      
      
        self.Xv=np.zeros([self.gn,1]) 
        xx=0
        for ii in range(self.gn):
            xx = np.sum(self.DXv[0:ii+1] )- self.DXv[ii]/2 
            self.Xv[ii]= xx         
        
         
        
        
   
        
class PropertyMatrix:
    
      def __init__(self):
        self.P=None
        self.Sw=None

        
      def initialize(self, NumericalSetting , Pressure , WaterSat):
              
          self.P= Pressure  *  np.ones([NumericalSetting.gn,1])
          self.Sw=WaterSat  *  np.ones([NumericalSetting.gn,1])
          
      def initializeuAuto(self,rf, NumericalSetting , Pressure):
          self.P= Pressure  *  np.ones([NumericalSetting.gn,1])
          self.Sw=(1-rf.snwr)  *  np.ones([NumericalSetting.gn,1])
                        


      def updateProps(self,Pressure,WaterSat):
          self.P=  Pressure #Matrix
          self.Sw= WaterSat #Matrix
      
      def copyProps(self, PM):
          self.P=  PM.P  #Matrix
          self.Sw= PM.Sw #Matrix  
      
      def copy(self):
          P=np.zeros(np.shape(self.P))
          Sw=np.zeros(np.shape(self.P))
          PMnew=PropertyMatrix()
          PMnew.updateProps(P,Sw)
          
          return PMnew
          
 
          
        

  
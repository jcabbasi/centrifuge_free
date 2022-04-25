
import numpy as np
import CentClasses as cc
from DataStorageClass import *
from CentSimCore import *

class FlowFnc():
    
    
    def __init__(self):
        self.lambdaW =None
        self.lambdaNW=None
        self.lambdaT=None
        self.lambdaWH =None
        self.lambdaNWH=None
        self.lambdaTH=None 
        self.PcX=None
        self.Pw=None
        self.SwNew=None
        self.SwOld=None 
        self.rf=None
        self.simprops=None
        
        self.__DX=None
        self.__DXv=None

        self.__gn=None
        self.__DT=None
        self.__phi=None
        self.__CT=None
        
        
        
                    
    
    
    
    def initialize(self,DX,DXv,gn):
        self.__gn=gn
        self.__DXv=DXv
        
        
        DXX=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        DXX[1:self.__gn+1]=self.__DXv     #allocation of saturation to the real cells 
        DXX[0]=0.001 #self.__DXv[0] 
         
        DXX[self.__gn+1]=0.001 #self.__DXv[-1]/250 
        self.DXX=DXX
        self.__DXvH= np.zeros([self.__gn+1,1])
        self.__DXvH=(DXX[0:self.__gn+1]+DXX[1:self.__gn+2])/2

        self.__XvH= np.zeros([self.__gn+1,1])
        
        xx=0
        for ii in range(gn):
            xx += self.__DXv[ii] 
            if ii==0:
                self.__XvH[ii+1]= 0
            elif ii==gn-2:
                pass
                self.__XvH[ii+1]=xx # np.sum(DXv)  
            else:
                self.__XvH[ii+1]= xx  
        xx=3
        
        
    
    def buildIt(self,installation,rf,simprops, Pw,SwNew,SwOld,DT):
        self.SwNew = SwNew
        self.SwOld=SwOld
        self.Pw=Pw
        self.installation=installation
        self.SwNew = SwNew
        self.SwOld = SwOld          
        self.__DT=DT
        self.rf=rf
        self.simprops=simprops
        
        # Constant values, extracted for code simplicity
        self.__CT=(1/rf.rock.phi)*self.__DT/self.__DXv
        omega=simprops.currentOmega

        
        
        SwX=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        SwX[1:self.__gn+1]=SwNew     #allocation of saturation to the real cells 
        SwXold=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        SwXold[1:self.__gn+1]=SwOld     #allocation of saturation to the real cells         
        
         
        #determining lambda
        lambdaWX=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        lambdaWX[1:self.__gn+1]=rf.lambdaW_f(SwNew) #allocation of saturation to the real cells 
        #
        lambdaNWX=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        lambdaNWX[1:self.__gn+1]=rf.lambdaNW_f(SwNew) #allocation of saturation to the real cells 
    
        #
        PX=np.zeros([self.__gn+2,1]) #including real and imaginary cells
        PX[1:self.__gn+1]=Pw #allocation of saturation to the real cells

        PcX=np.zeros([self.__gn+2,1])  
        PcX[1:self.__gn+1]=rf.pc(SwNew)
        #
 
        
        if installation=='OEO':
            PX[0]=0 ; PX[self.__gn+1]=0
            SwX[0]=0; SwX[self.__gn+1]=rf.Swmin()
            SwXold[0]=0;  SwXold[self.__gn+1]=rf.Swmin()
            
            lambdaWX[0]=0 #closed
            lambdaWX[self.__gn+1]=lambdaWX[self.__gn] #
            lambdaNWX[0]=0 #closed
            lambdaNWX[self.__gn+1]=rf.lambdaNW_f(rf.Swmin())    # Outside the core
            # lambdaNWX[self.__gn  ]=rf.lambdaNW_f(rf.Swmin())    # Outside the core

            PcX[0]=0
            PcX[self.__gn+1]=0        
      

                   
        elif installation=='TEO':
            SwX[0]=rf.Swmin()
            SwX[self.__gn+1]=rf.Swmax()    
            SwXold[0]=rf.Swmin()
            SwXold[self.__gn+1]=rf.Swmax()  
            PX[0]=self.PnwEqDuetoRotationMAX(omega) 
            PX[1]=self.PnwEqDuetoRotationMAX(omega)-rf.pc(SwX[1])
            PX[self.__gn+1]=0
            lambdaWX[0]=0 
            lambdaWX[self.__gn+1]=rf.lambdaW_f([rf.Swmax()])  # Inside the core
            lambdaNWX[0]=rf.lambdaNW_f([rf.Swmin()])
            lambdaNWX[self.__gn+1]=0     # Outside the core - 
            PcX[0]=0
            PcX[self.__gn+1]=0               

        # Mobility Properties on the grid faces (boundaries)
        lambdaWH=np.zeros([self.__gn+1,1])
        lambdaNWH=np.zeros([self.__gn+1,1])
        
        self.DXT=self.DXX.copy()
        # self.DXT[-1]=self.DXT[-1]*200
        
        for ii in range(0,self.__gn+1):        
            lambdaWH[ii] = ( lambdaWX[ii]  * self.DXT[ii]  + lambdaWX[ii+1]  * self.DXT[ii+1]  ) / (self.DXT[ii]+self.DXT[ii+1] )    #Half properties: The vector is defined for the edge of cells (inner edge)
            lambdaNWH[ii]= ( lambdaNWX[ii] * self.DXT[ii]  + lambdaNWX[ii+1] * self.DXT[ii+1]  ) / (self.DXT[ii]+self.DXT[ii+1] )   #Half properties: The vector is defined for the edge of cells (inner edge)
        
        
        lambdaTH=lambdaWH+lambdaNWH

        if installation   == 'OEO':
            lambdaWH[0] =0
            lambdaNWH[0]=0 
            lambdaTH[0]=0
            lambdaNWH[self.__gn]=lambdaNWX[self.__gn+1]
            lambdaWH[self.__gn]= rf.lambdaW_f([self.rf.PcZeroCalc()])
            lambdaTH[self.__gn]= lambdaWH[self.__gn]+lambdaNWH[self.__gn]                      
            
        elif installation == 'TEO':
            lambdaWH[0]= 0
            #lambdaWH[self.__gn]= lambdaWX[self.__gn]
            lambdaNWH[self.__gn]= 0            
            lambdaNWH[0]=lambdaNWX[0]
            lambdaTH[0]=lambdaNWH[0]
            lambdaTH[self.__gn]=lambdaWH[self.__gn]

               
            
            
        self.lambdaW =lambdaWX
        self.lambdaNW=lambdaNWX
        self.lambdaT=None
        self.lambdaWH =lambdaWH
        self.lambdaNWH=lambdaNWH
        self.lambdaTH=lambdaTH 
        self.PcX=PcX
        self.Pw=PX
        
        

        
    def F1(self,cellList):
        
        omega=self.simprops.currentOmega
        rhoW=self.rf.fw.den
        CT=self.__CT
        lambdaWH=self.lambdaWH
        PX=self.Pw  
        DXvH=self.__DXvH
        f1=np.zeros([len(cellList),1])
        jj=0
        for ii in cellList:
 
            # if ii==self.__gn-1:
                
            #     f1[jj]= self.SwNew[ii]-self.SwOld[ii] - CT[ii]*( (  lambdaWH[ii+1]*( (PX[ii+2]-PX[ii+1])/DXvH[ii+1] -rhoW*omega*omega*( self.__XvH[ii+1] + self.simprops.r1)) 
            #             -(lambdaWH[ii])*( (PX[ii+1]-PX[ii])/DXvH[ii] -rhoW*omega*omega*( self.__XvH[ii+1] + self.simprops.r1 )) ))
            
                        
            
            # else:


            f1[jj]= self.SwNew[ii]-self.SwOld[ii] - CT[ii]*( (  lambdaWH[ii+1]*( (PX[ii+2]-PX[ii+1])/DXvH[ii+1] -rhoW*omega*omega*( self.__XvH[ii+1] + self.simprops.r1)) 
                        -(lambdaWH[ii])*( (PX[ii+1]-PX[ii])/DXvH[ii] -rhoW*omega*omega*( self.__XvH[ii] + self.simprops.r1 )) ))
            
            
            
            jj=jj+1
        
        return f1
     
     
        
    def F2(self,cellList):
        omega=self.simprops.currentOmega

        lambdaWH=self.lambdaWH
        lambdaNWH=self.lambdaNWH
        lambdaTH=self.lambdaTH
        PX=self.Pw  
        PcX=self.PcX 
        rhoW=self.rf.fw.den
        rhoNW=self.rf.fnw.den
        
        DXvH=self.__DXvH

        f2=np.zeros([len(cellList),1])
        jj=0
        for ii in cellList:
     
            # if ii==self.__gn-1:
 
            #     f2[jj]=( (  lambdaNWH[ii+1]  )*( (PcX[ii+2]-PcX[ii+1])/DXvH[ii+1] ) + (  lambdaTH [ii+1]    )*( (PX[ii+2]- PX[ii+1]) /DXvH[ii+1] )
            # - omega * omega*  (self.__XvH[ii+1] + self.simprops.r1)*(   lambdaWH[ii+1] * rhoW +  lambdaNWH[ii+1] * rhoNW)
            # -(lambdaNWH[ii]  )*( (PcX[ii+1]-PcX[ii])/DXvH[ii] ) - (lambdaTH [ii]  )*( (PX[ii+1]- PX[ii]) /DXvH[ii] )
            # + omega * omega * (self.__XvH[ii+1] + self.simprops.r1) * (  lambdaWH[ii] * rhoW + lambdaNWH[ii]  * rhoNW ) )
                           
            
            # else:

            f2[jj]=( (  lambdaNWH[ii+1]  )*( (PcX[ii+2]-PcX[ii+1])/DXvH[ii+1] ) + (  lambdaTH [ii+1]    )*( (PX[ii+2]- PX[ii+1]) /DXvH[ii+1] )
            - omega * omega*  (self.__XvH[ii+1] + self.simprops.r1)*(   lambdaWH[ii+1] * rhoW +  lambdaNWH[ii+1] * rhoNW)
            -(lambdaNWH[ii]  )*( (PcX[ii+1]-PcX[ii])/DXvH[ii] ) - (lambdaTH [ii]  )*( (PX[ii+1]- PX[ii]) /DXvH[ii] )
            + omega * omega * (self.__XvH[ii] + self.simprops.r1) * (  lambdaWH[ii] * rhoW + lambdaNWH[ii]  * rhoNW ) )
                
            jj=jj+1
        
        return f2

    
    
        
    def PcEqDuetoRotation(self,omega,x):
        r2=self.simProps.r1+self.rf.rock.length        
        pceqMax= 0.5 * self.rf.dDen() * omega * omega * ( r2 * r2 - x * x )
        return pceqMax       
     
    def PnwEqDuetoRotationMAX(self,omega): 
        r2=self.simprops.r1+self.rf.rock.length
        pnwMax= - 0.5 * self.rf.fnw.den * omega * omega * ( r2 * r2 - self.simprops.r1 * self.simprops.r1 )
        return pnwMax     
         
            

    
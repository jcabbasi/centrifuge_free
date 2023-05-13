import time as tm
import matplotlib.pyplot as plt
import numpy as np


class intercativeClass():
    
    def __init__(self,NSet,simprops,rf):
        

        self.Line1=None
        self.Line2=None
        self.Line3=None 
        self.Line4=None 
        
        self.fig=None
        self.gn=NSet.gn
        self.simprops=simprops
        self.ErrorsVec=np.zeros([20,1])
        self.rf=rf 
        self.updateSeq=10
        self.__unpdatingRun=False
        self.ax=None
        self.t=[0]
        self.RF=[0]
        self.XvD=NSet.Xv/NSet.rockLength
        self.Xv=NSet.Xv
        self.axoptimize=None
        self.figopt=None
        self.optObserved=None
        
        
        
    def InteractivePlotInititialize(self,InteractivePlot):
       self.InteractivePlot=InteractivePlot
       if self.InteractivePlot==True:
         # enable interactive mode
            plt.ion()
            # creating subplot and figure
            self.fig = plt.figure(figsize=(13, 12), dpi=80)
            #ax = fig.add_subplot(111)
            self.ax = self.fig.subplots(2,2)
    
    
    
    def InteractivePlotInit(self,AnalyticalSat,Sold):
        if self.InteractivePlot==True:
         # enable interactive mode
            plt.ion()
            self.line1, = self.ax[0,0].plot(self.XvD, Sold)
            self.ax[0,0].scatter(self.XvD,   AnalyticalSat,color='black')
            self.ax[0,0].set_xlim([0, 1])
            self.ax[0,0].set_ylim([np.min(AnalyticalSat)*0.99, np.max(AnalyticalSat)*1.01])
            # setting labels
            self.ax[0,0].set_xlabel("Dimensionless Length")
            self.ax[0,0].set_ylabel("Water Saturation")
            self.ax[0,0].set_title("Saturation Plot: Updating ...")  

            r2 = self.simprops.r1 + self.rf.rock.length

            x = self.simprops.r1+self.Xv
            
            r2=self.simprops.r1+self.rf.rock.length        
            pceq= 0.5 * self.rf.dDen() * self.simprops.currentOmega * self.simprops.currentOmega * ( r2 * r2 - x * x )
            
            self.line2,= self.ax[0,1].plot(Sold, pceq , 'ro')            
            satrange = np.linspace(self.rf.swc,1-self.rf.snwr, num=500, endpoint=True)
            pccurve=self.rf.pc(satrange)    
            self.ax[0,1].plot(satrange, pccurve)
                        
            self.ax[0,1].set_xlabel("Water Saturation")
            self.ax[0,1].set_ylabel("Capillary Pressure")
        
            self.ax[0,1].set_title("Capillary Pressure: Updating ...")              
            
            
            self.line3, = self.ax[1,0].plot(0.1, 0.0000001 , 'ro')
            self.ax[1,0].set_xlabel("Time (s)")
            self.ax[1,0].set_ylabel("Recovery Factor")
            self.ax[1,0].set_title("Recovery Factor: Updating ...")             
            self.ax[1,0].set_xscale('symlog'); self.ax[1,0].set_yscale('symlog')
            
            
            self.line4, = self.ax[1,1].plot(np.linspace(1,len(self.ErrorsVec),len(self.ErrorsVec)), self.ErrorsVec , 'o')
            self.ax[1,1].set_xlabel("last TS")
            self.ax[1,1].set_ylabel("Converged Error Values")
            self.ax[1,1].set_title("Convergance Analysis: Updating ...")             
            self.ax[1,1].set_yscale('symlog')            
            self.ax[1,1].set_xlim(0,len(self.ErrorsVec))
            self.txt=self.ax[1,1].text(1, 1, '--', style='italic',
                            bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 10})
            # re-drawing the figure
            self.fig.canvas.draw()
            # to flush the GUI events
            self.fig.canvas.flush_events()            
        
    def updateInteractivePlot(self,InteractivePlot,stepNo,Sold, t, RF,Errs,str):
        
        if np.remainder(stepNo,self.updateSeq)==0:
            if InteractivePlot==True:
                self.line1.set_ydata(Sold)
                self.line2.set_xdata(Sold)

                
                self.t.append(t);self.RF.append(RF)
                self.line3.set_xdata(self.t); self.line3.set_ydata(self.RF)
                self.ax[1,0].set_xlim([0, t*1.2]);  self.ax[1,0].set_ylim([0, RF*2])
                self.ax[1,0].set_xscale('log'); self.ax[1,0].set_yscale('log')


                self.ErrorsVec[0:-1]=self.ErrorsVec[1:]
                self.ErrorsVec[-1]=Errs
                self.line4.set_ydata(self.ErrorsVec)
                self.ax[1,1].set_ylim([np.min(self.ErrorsVec)*0.8, np.max(self.ErrorsVec)*1.2])
                self.txt.set_text(str)
                self.txt.set_position([1, np.max(self.ErrorsVec)*0.85])

                if self.__unpdatingRun==False:
                    self.ax[0,0].set_title("Saturation Plot: Updating ..",color='black')
                    self.ax[0,1].set_title("Saturation Plot: Updating ..",color='black')
                    self.ax[1,0].set_title("Recovery Factor: Updating ..",color='black')                 
                    self.__unpdatingRun=True 
                else:
                    self.ax[0,0].set_title("Saturation Plot: Updating ....",color='green')  
                    self.ax[0,1].set_title("Saturation Plot: Updating ....",color='green')  
                    self.ax[1,0].set_title("Recovery Factor: Updating ....",color='green')                 
                    
                    self.__unpdatingRun=False 
                    
                    
                # re-drawing the figure
                self.fig.canvas.draw()
                # to flush the GUI events
                self.fig.canvas.flush_events()
                
                
                
                tm.sleep(0.000)
                
             
           
    def InteractivePlotInitOptimizer(self,obs,optPlt):
       self.optPlt=optPlt
       if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            self.optObserved=obs
            # creating subplot and figure
            self.figopt = plt.figure(figsize=(7, 12), dpi=80)
            #ax = fig.add_subplot(111)
            self.axoptimize = self.figopt.subplots(2,2)
            self.figopt.canvas.draw()

    
    def visualizeOptimization(self,predict,optErr):
        self.optErr=optErr
        if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            #self.ax[0,0].cla() 
            self.line1, = self.axoptimize[0,0].plot(predict[:,0], predict[:,1])
            self.axoptimize[0,0].scatter(self.optObserved[:,0],   self.optObserved[:,1],color='black')
            self.axoptimize[0,0].set_xlim([0, 1])
            self.axoptimize[0,0].set_ylim([0, np.max(self.optObserved[:,1])*1.01])
            # setting labels
            self.axoptimize[0,0].set_xlabel("Time (s)")
            self.axoptimize[0,0].set_ylabel("Recovery Factor")
            self.axoptimize[0,0].set_title("Updating ...")  


            self.line2,= self.axoptimize[0,1].plot([], [] , 'ro')            
                        
            self.axoptimize[0,1].set_xlabel("Water Saturation")
            self.axoptimize[0,1].set_ylabel("Capillary Pressure")
            self.axoptimize[0,1].set_title("Capillary Pressure: Updating ...")              
            
            
            self.line3, = self.axoptimize[1,0].plot(0.1, 0.0000001 , 'ro')
            self.axoptimize[1,0].set_xlabel("Time (s)")
            self.axoptimize[1,0].set_ylabel("Recovery Factor")
            self.axoptimize[1,0].set_title("Recovery Factor: Updating ...")             
            self.axoptimize[1,0].set_xscale('symlog'); self.axoptimize[1,0].set_yscale('symlog')
            
            
            self.line4, = self.axoptimize[1,1].plot(self.optErr[:,0], self.optErr[:,1] , 'o')
            self.axoptimize[1,1].set_xlabel("Iteration")
            self.axoptimize[1,1].set_ylabel("Error Value")
            self.axoptimize[1,1].set_title("Updating ...")             
            self.axoptimize[1,1].set_yscale('symlog')            
            # re-drawing the figure
            # self.figopt.canvas.draw()
            # to flush the GUI events
            # self.figopt.canvas.flush_events()     
                
        
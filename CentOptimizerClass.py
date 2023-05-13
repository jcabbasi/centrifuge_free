
import numpy as np
import CentClasses as cc
from DataStorageClass import *
import time as tm
from  CentSimCore import *
from scipy.interpolate import interp1d
import os
from visualizationClass import *
from geneticalgorithm import geneticalgorithm as ga
import matplotlib.pyplot as plt
import matplotlib
import pyswarms as ps
from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)
import tensorflow as tf

#matplotlib.use('Agg')
# plt.switch_backend('Agg')

class GAoptimizer():
    
    def __init__(self,SimClass,RockFluid,observed):
        self.SimClass=SimClass
        self.rf=RockFluid
        self.observed=observed
        self.runNo=0
        self.optData=  pd.DataFrame(index=[], columns=['nW','nNW','KrWmax','KrNWmax','ErrorVal','FinalTime','FinalRF','Generation','ProdCurve'])
        self.PopSize=20
        self.optErr=np.zeros([])
        self.figOpt=None
        self.intercativePlot=intercativeClass(SimClass.NSet,SimClass.simprops,RockFluid)

    #
    def assignResults(self,nW,nNW,KrWmax,KrNWmax,ErrorVal,FinalTime,FinalRF,gen,Predicted):
        xx=len(self.optData)
        self.optData.loc[xx+1,'nW']=nW
        self.optData.loc[xx+1,'nNW']=nNW
        self.optData.loc[xx+1,'KrWmax']=KrWmax
        self.optData.loc[xx+1,'KrNWmax']=KrNWmax
        self.optData.loc[xx+1,'ErrorVal']=ErrorVal
        self.optData.loc[xx+1,'Generation']=gen
        self.optData.loc[xx+1,'FinalTime']=FinalTime
        self.optData.loc[xx+1,'FinalRF']=FinalRF
        self.optData.loc[xx+1,'ProdCurve']=Predicted

    
    
    def OptimizeNowGA(self,optPlt=True,simTime=None):
        self.simTime=simTime
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        self.intercativePlot.InteractivePlotInitOptimizer(self.observed,optPlt=False)

        def fitness_func(x):
            
            gen=np.floor(self.runNo/self.PopSize)
            self.runNo=self.runNo+1

            nW=x[0]
            nNW=x[1]
            KrWmax=x[2]
            KrNWmax=x[3]        
            ErrorVal=np.zeros([np.size(nW)])   
            nW,nNW,KrWmax,KrNWmax=cc.RockFluid.QCRelPermCoefs(self,nW,nNW,KrWmax,KrNWmax)
            print('NEXT:')
            print('++++ Optimization Variables ++++, RunNo: {4}  nW: {0:.2f} , nNW: {1:.2f} ; KrWmax: {2:.2f}; KrNWmax: {3:.2f} ; \n'.format(nW, nNW, KrWmax, KrNWmax,self.runNo )   )
            self.rf.initializeRelperm(nW,nNW,KrWmax,KrNWmax)
            self.SimClass.refreshModel()
            self.SimClass.updateRockFluid(self.rf)
            results=self.SimClass.simulateNow(InteractivePlot=False,reportPeriod=2,ExportbySingle=False,simTime=self.simTime)
            Predicted=results.ExtractSimulationResults(['Time','RF'])
            FinalTime=Predicted[-1,0]
            FinalRF  =Predicted[-1,1]
            ErrorVal, ErrorCurve=self.findErrorBetweenTwoVector(self.observed, Predicted)
            self.optErr=np.append(self.optErr,ErrorVal)
            print( '++++++++ Error Value:::::  {0}  \n'.format(ErrorVal )   )
            self.assignResults(nW,nNW,KrWmax,KrNWmax,ErrorVal,FinalTime,FinalRF,gen,Predicted)
            self.intercativePlot.visualizeOptimization(Predicted,ErrorCurve)
            self.optData.to_pickle("./optData.pkl")
            self.optData.loc[:,['nW','nNW','KrWmax','KrNWmax','FinalTime','FinalRF','Generation','ErrorVal']].to_csv("./optData.csv")
            if FinalRF<0.38:
                ErrorVal=100
            return ErrorVal  


        algorithm_param = {'max_num_iteration': None,\
                   'population_size':50,\
                   'mutation_probability':0.15,\
                   'elit_ratio': 0.01,\
                   'crossover_probability': 0.3,\
                   'parents_portion': 0.2,\
                   'crossover_type':'uniform', 
                   'max_iteration_without_improv':None}

        varbound=np.array([[0.2,4],[0.2,4],[0.1,1],[0.1,1]])
        # varbound=np.array([[0.7,0.701],[1.8,1.8001],[0.15,0.1501],[0.35,0.3501]])
        # varbound=np.array([[2.028,2.028592],[1.432,1.43271],[0.17,0.176],[0.729,0.729625]])

        model=ga(function=fitness_func,dimension=4,variable_type='real',variable_boundaries=varbound,function_timeout=5e5,algorithm_parameters=algorithm_param)

        model.run()
        
        model

    
 
    
    def findErrorBetweenTwoVector(self,Obsrved,Predicted):
        

        f = interp1d(Predicted[:,0], Predicted[:,1])
        maxTobs=np.max(Obsrved[:,0]);maxTpred=np.max(Predicted[:,0])
        if maxTobs>maxTpred:
            delt = np.where(Obsrved[:,0] > maxTpred)
            ObsrvedN=np.delete(Obsrved,delt,axis=0)
        else:
            ObsrvedN=Obsrved.copy()
        ErrorCurve=np.zeros(np.shape(ObsrvedN))
        ErrorCurve[:,0]=ObsrvedN[:,0]
        InterPolatedObservedValues=f(ObsrvedN[:,0])
        ErrorVal=np.power(sum(abs(ObsrvedN[:,1]*ObsrvedN[:,1]-InterPolatedObservedValues*InterPolatedObservedValues)),0.5)/np.power(sum(abs(ObsrvedN[:,1]*ObsrvedN[:,1])),0.5)
        ErrorVal=ErrorVal**4
        ErrorCurve[:,1]=np.power((abs(ObsrvedN[:,1]*ObsrvedN[:,1]-InterPolatedObservedValues*InterPolatedObservedValues)),0.5)
        
        return ErrorVal, ErrorCurve
    
    
    
           
    def InteractivePlotInititialize(self,optPlt):
       self.optPlt=optPlt
       if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            # creating subplot and figure
            self.fig = plt.figure(figsize=(7, 12), dpi=80)
            #ax = fig.add_subplot(111)
            self.ax = self.fig.subplots(2,2)
            plt.switch_backend('Agg')

    
    def visualizeOptimization(self,predict):
        if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            #self.ax[0,0].cla() 
            self.line1, = self.ax[0,0].plot(predict[:,0], predict[:,1])
            self.ax[0,0].scatter(self.observed[:,0],   self.observed[:,1],color='black')
            self.ax[0,0].set_xlim([0, 1])
            self.ax[0,0].set_ylim([0, np.max(self.observed[:,1])*1.01])
            # setting labels
            self.ax[0,0].set_xlabel("Time (s)")
            self.ax[0,0].set_ylabel("Recovery Factor")
            self.ax[0,0].set_title("Updating ...")  


            self.line2,= self.ax[0,1].plot([], [] , 'ro')            
                        
            self.ax[0,1].set_xlabel("Water Saturation")
            self.ax[0,1].set_ylabel("Capillary Pressure")
            self.ax[0,1].set_title("Capillary Pressure: Updating ...")              
            
            
            self.line3, = self.ax[1,0].plot(0.1, 0.0000001 , 'ro')
            self.ax[1,0].set_xlabel("Time (s)")
            self.ax[1,0].set_ylabel("Recovery Factor")
            self.ax[1,0].set_title("Recovery Factor: Updating ...")             
            self.ax[1,0].set_xscale('symlog'); self.ax[1,0].set_yscale('symlog')
            
            plt.switch_backend('Agg')
            
            self.line4, = self.ax[1,1].plot(np.linspace(1,len(self.optErr),len(self.optErr)), self.optErr , 'o')
            self.ax[1,1].set_xlabel("Iteration")
            self.ax[1,1].set_ylabel("Error Value")
            self.ax[1,1].set_title("Updating ...")             
            self.ax[1,1].set_yscale('symlog')            
            # re-drawing the figure
            self.fig.canvas.draw()
            # to flush the GUI events
            self.fig.canvas.flush_events()     
        
        
        
        

class PSOoptimizer():
    
    def __init__(self,SimClass,RockFluid,observed,Dual=False,obsDual=None,simClassDual=None):
        self.SimClass=SimClass
        self.rf=RockFluid
        self.observed=observed
        self.runNo=0
        self.optData=  pd.DataFrame(index=[], columns=['nW','nNW','KrWmax','KrNWmax','ErrorVal','FinalTime','FinalRF','Generation','ProdCurve','ErrorVal1','ErrorVal2','ProdCurveDual'])
        self.PopSize=20
        self.optErr=np.zeros([])
        self.figOpt=None
        self.intercativePlot=intercativeClass(SimClass.NSet,SimClass.simprops,RockFluid)
        self.dual=Dual
        self.obsDual=obsDual
        self.simClassDual=simClassDual
    #
    def assignResults(self,nW,nNW,KrWmax,KrNWmax,ErrorVal,FinalTime,FinalRF,gen,Predicted,ErrorVal1=0,ErrorVal2=0,ProdCurveDual=0):
        xx=len(self.optData)
        self.optData.loc[xx+1,'nW']=nW
        self.optData.loc[xx+1,'nNW']=nNW
        self.optData.loc[xx+1,'KrWmax']=KrWmax
        self.optData.loc[xx+1,'KrNWmax']=KrNWmax
        self.optData.loc[xx+1,'ErrorVal']=ErrorVal
        self.optData.loc[xx+1,'Generation']=gen
        self.optData.loc[xx+1,'FinalTime']=FinalTime
        self.optData.loc[xx+1,'FinalRF']=FinalRF
        self.optData.loc[xx+1,'ProdCurve']=Predicted
        self.optData.loc[xx+1,'ErrorVal1']=ErrorVal1
        self.optData.loc[xx+1,'ErrorVal2']=ErrorVal2
        self.optData.loc[xx+1,'ProdCurveDual']=ProdCurveDual

    
    
    def OptimizeNowPSO(self,optPlt=True,simTime=None):
        self.simTime=simTime
        os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
        self.intercativePlot.InteractivePlotInitOptimizer(self.observed,optPlt=False)
        

            
        
        
        # PSO Algorithm:
        maxbound=[1.5,2.9,0.6,0.6]
        minbound=[0.3,1.0,0.05,0.1]
        
        # maxbound=[0.7,1.8,0.15,0.35]
        # minbound=maxbound
        
        bounds=(minbound,maxbound)      
        # Initialize swarm
        options = {'c1': 0.6, 'c2': 0.5, 'w':0.61}
        n_particles=40
        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles, dimensions=4, options=options, bounds=bounds)
        self.round=0
        ErrorValDual=0

        def fitfunc(xx):
            self.round+=1
            scoreReg=np.zeros(n_particles)
            for ii in range(len(scoreReg)):
                ErrorVal=fitness_func(xx[ii,:])
                scoreReg[ii]=ErrorVal
            return scoreReg
            

        def fitness_func(x):
            
            self.runNo=self.runNo+1
            
            nW=x[0]
            nNW=x[1]
            KrWmax=x[2]
            KrNWmax=x[3]        
            ErrorVal=np.zeros([np.size(nW)])   
            nW,nNW,KrWmax,KrNWmax=cc.RockFluid.QCRelPermCoefs(self,nW,nNW,KrWmax,KrNWmax)
            print('NEXT: Optimization Variables; RunNo: {4}  nW: {0:.2f} , nNW: {1:.2f} ; KrWmax: {2:.2f}; KrNWmax: {3:.2f} ; \n'.format(nW, nNW, KrWmax, KrNWmax,self.runNo )   )
            
            self.rf.initializeRelperm(nW,nNW,KrWmax,KrNWmax)
            self.SimClass.refreshModel()
            self.SimClass.updateRockFluid(self.rf)
            results=self.SimClass.simulateNow(InteractivePlot=False,reportPeriod=2,ExportbySingle=True,simTime=self.simTime)
            Predicted=results.ExtractSimulationResults(['Time','RF'])
            FinalTime=Predicted[-1,0]
            FinalRF  =Predicted[-1,1]
            ErrorValsingle, ErrorCurve=self.findErrorBetweenTwoVector(self.observed, Predicted)
            ErrorValsingle=ErrorValsingle*1000

            if self.dual:
                self.rf.initializeRelperm(nW,nNW,KrWmax,KrNWmax)
                self.simClassDual.refreshModel()
                self.simClassDual.updateRockFluid(self.rf)
                results=self.simClassDual.simulateNow(InteractivePlot=False,reportPeriod=2,ExportbySingle=True,simTime=self.simTime)
                PredictedDual=results.ExtractSimulationResults(['Time','RF'])
                FinalTime=PredictedDual[-1,0]
                FinalRF  =PredictedDual[-1,1]
                ErrorValDual, ErrorCurve=self.findErrorBetweenTwoVector(self.obsDual, PredictedDual)
                ErrorValDual=ErrorValDual*1000
                
            ErrorVal=ErrorValsingle+ErrorValDual
            
            self.optErr=np.append(self.optErr,ErrorVal)
            print( '++++++++ Error Value::::: Total: {0:.3f}; First: {2:.3f};  Second: {1:.3f} \n'.format(ErrorVal,ErrorValDual,ErrorValsingle )   )

            if self.dual:            
                self.assignResults(nW,nNW,KrWmax,KrNWmax,ErrorVal,FinalTime,FinalRF,self.round,Predicted,ErrorVal1=ErrorValsingle,ErrorVal2=ErrorValDual,ProdCurveDual=PredictedDual)
            else:
                self.assignResults(nW,nNW,KrWmax,KrNWmax,ErrorVal,FinalTime,FinalRF,self.round,Predicted)
            
            self.intercativePlot.visualizeOptimization(Predicted,ErrorCurve)
            self.optData.to_pickle("./optData.pkl")
            self.optData.loc[:,['nW','nNW','KrWmax','KrNWmax','FinalTime','FinalRF','Generation','ErrorVal']].to_csv("./optData.csv")
            # if FinalRF<0.38:
            #     ErrorVal=100
            return ErrorVal 
        
        #Optimization:
        x0=[40,5.0,0.06]
        kwargs={"ga":x0[1],"ep":x0[2]}
        cost, pos = optimizer.optimize(fitfunc, 40)
        plot_cost_history(cost_history=optimizer.cost_history)
        plt.show()

        
        

    
 
    
    def findErrorBetweenTwoVector(self,Obsrved,Predicted):
        

        
        f = interp1d(Predicted[:,0], Predicted[:,1])
        maxTobs=np.max(Obsrved[:,0]);maxTpred=np.max(Predicted[:,0])
        if maxTobs>maxTpred:
            delt = np.where(Obsrved[:,0] > maxTpred)
            ObsrvedN=np.delete(Obsrved,delt,axis=0)
        else:
            ObsrvedN=Obsrved.copy()
        
        # dtv=np.zeros(len(ObsrvedN[:,0]))
        # for ll in range(len(ObsrvedN[:,0])-1):
        #    dtv[ll+1]= ObsrvedN[ll,0]-ObsrvedN[ll+1,0]
        # meandt=np.mean(dtv)
        
        ErrorCurve=np.zeros(np.shape(ObsrvedN))
        ErrorCurve[:,0]=ObsrvedN[:,0]
        predInterp=f(ObsrvedN[:,0]) #InterPolatedObservedValues
            
        loss,lossc=self.lossf(ObsrvedN[:,1], predInterp); 
        

        return loss, lossc
    
  
    def lossf(self,y_true, y_pred):
        loss = tf.math.squared_difference(y_pred, y_true)**1
        loss_f=tf.reduce_mean(loss)**1
        return loss_f, loss  
    
           
    def InteractivePlotInititialize(self,optPlt):
       self.optPlt=optPlt
       if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            # creating subplot and figure
            self.fig = plt.figure(figsize=(7, 12), dpi=80)
            #ax = fig.add_subplot(111)
            self.ax = self.fig.subplots(2,2)
            plt.switch_backend('Agg')

    
    def visualizeOptimization(self,predict):
        if self.optPlt==True:
         # enable interactive mode
            plt.ion()
            #self.ax[0,0].cla() 
            self.line1, = self.ax[0,0].plot(predict[:,0], predict[:,1])
            self.ax[0,0].scatter(self.observed[:,0],   self.observed[:,1],color='black')
            self.ax[0,0].set_xlim([0, 1])
            self.ax[0,0].set_ylim([0, np.max(self.observed[:,1])*1.01])
            # setting labels
            self.ax[0,0].set_xlabel("Time (s)")
            self.ax[0,0].set_ylabel("Recovery Factor")
            self.ax[0,0].set_title("Updating ...")  


            self.line2,= self.ax[0,1].plot([], [] , 'ro')            
                        
            self.ax[0,1].set_xlabel("Water Saturation")
            self.ax[0,1].set_ylabel("Capillary Pressure")
            self.ax[0,1].set_title("Capillary Pressure: Updating ...")              
            
            
            self.line3, = self.ax[1,0].plot(0.1, 0.0000001 , 'ro')
            self.ax[1,0].set_xlabel("Time (s)")
            self.ax[1,0].set_ylabel("Recovery Factor")
            self.ax[1,0].set_title("Recovery Factor: Updating ...")             
            self.ax[1,0].set_xscale('symlog'); self.ax[1,0].set_yscale('symlog')
            
            plt.switch_backend('Agg')
            
            self.line4, = self.ax[1,1].plot(np.linspace(1,len(self.optErr),len(self.optErr)), self.optErr , 'o')
            self.ax[1,1].set_xlabel("Iteration")
            self.ax[1,1].set_ylabel("Error Value")
            self.ax[1,1].set_title("Updating ...")             
            self.ax[1,1].set_yscale('symlog')            
            # re-drawing the figure
            self.fig.canvas.draw()
            # to flush the GUI events
            self.fig.canvas.flush_events()     
        
  

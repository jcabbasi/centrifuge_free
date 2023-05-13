import pandas as pd
from CentSimCore import Simulate
import matplotlib.pyplot as plt
import time as tm



class SentivityAnalysis():
    

    
    def DoSensAnalysis(fileName,parameter,Values,ExportbySingle=False):
        para={'muo','muw','omega','gridno','installation'}
        fileName=fileName
    
        base=pd.read_pickle(fileName)
        out1=pd.DataFrame(index=range(5000), columns=[])
        
        #Output:
        fig, ax=plt.subplots(1)
         # enable interactive mode
        plt.ion()
        ax.set_xscale('symlog')            
        ax.set_yscale('symlog')  


        sensClass=base.copy()
        for val in Values:
            col1='Time_'+parameter+'='+str(val)
            col2='RF___'+parameter+'='+str(val)
            out1[col1]=None
            out1[col2]=None
                        
            if parameter=='muo':
                sensClass.rf.fnw.mu=val
                
            if parameter=='muw':
                sensClass.rf.fw.mu=val            
            
            if parameter=='krWMAX':
                sensClass.rf.krWMAX=val
                         
            if parameter=='omega':
                sensClass.simprops.omegaVect=[val]
                
            if parameter=='gridno':
                sensClass.NSet.gn=val
                sensClass.NSet.DXvSet()
            if parameter=='installation':
                sensClass.installation=val
                
            sensClass.init.initializeuAuto(sensClass.rf,sensClass.NSet,0)

            name='SeAna__' + parameter + '__' + str(val)
            print('------------------------------  Sensitivity Analysis  ----------------------------')
            print('++++ Now the value of   {0}  is: {1}, .......Lets Go....... \n'.format(parameter, val )   )
            senClassTry=Simulate(sensClass.simprops,sensClass.rf,sensClass.init,sensClass.NSet,sensClass.installation,printLog=True) 
            res=senClassTry.simulateNow(FileName=name,InteractivePlot=True,reportPeriod=1,ExportbySingle=ExportbySingle)
            Predicted=res.ExtractSimulationResults(['Time','RF'])
            out1.loc[range(len(Predicted[:,0])),col1]=Predicted[:,0]
            out1.loc[range(len(Predicted[:,0])),col2]=Predicted[:,1]
            ax.scatter(Predicted[2:,0], Predicted[2::,1])
            ax.set_xscale('symlog')            
            ax.set_yscale('symlog')          
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("RF - Sensitivity")
            # re-drawing the figure
            fig.canvas.draw()
            # to flush the GUI events
            fig.canvas.flush_events()    
            tm.sleep(0.000)
        
        fn= "./"+fileName+"__Sens__.csv"
        # out1.to_pickle(fn)   
        out1.to_csv(fn)
        print('Finish')
       
        
        
        
        
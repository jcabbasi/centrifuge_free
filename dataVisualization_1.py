import numpy as np
from CentClasses import RockFluid as RockFluid
from DataStorageClass import *
from CentOptimizerClass import *
import matplotlib.pyplot as plt
import pandas as pd


optData = pd.read_pickle("./optData.pkl")

nWvector=np.array(optData.loc[:,'nW'])
nNWvector=np.array(optData.loc[:,'nNW'])
ErrorVal=np.array(optData.loc[:,'ErrorVal'])
gen=np.array(optData.loc[:,'Generation'],dtype=int)

colorList=np.array(['blue','red','orange','green','black','purple','grey'])
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(nWvector, nNWvector, ErrorVal,color=colorList[gen])
ax.set_xlabel('nW');ax.set_ylabel('nNW');ax.set_zlabel('Error')
vv=8
fig, ax2 = plt.subplots(ncols=2,nrows=1, sharey=False)
filename="\\ObsData.xlsx"
observed=LoadandSave.LoadDataFromExcel(filename,'B5','C110')
ax2[0].scatter(observed[:,0], observed[:,1],color="red",s=150,alpha=0.7)
ax2[0].set_xlabel('Time (s)');ax2[0].set_ylabel('RF')
ax2[1].set_xlabel('Run Number');ax2[1].set_ylabel('Error')

for ii in optData.index:
    if optData.loc[ii,'ErrorVal']<0.84:
        prodCurve=optData.loc[ii,'ProdCurve']
        ax2[0].scatter(prodCurve[:,0], prodCurve[:,1],color="blue")
        ax2[1].scatter(ii,optData.loc[ii,'ErrorVal'],color="blue")

cc=3




from CentSensitivityAnalysisClass import SentivityAnalysis as sa



# sa.DoSensAnalysis('saveOEO_MW40.pickle', 'gridno',[10,20,30,40,50])

sa.DoSensAnalysis('saveOEO_MW40.pickle', 'muo',[0.002170,0.002170*0.1],ExportbySingle=True)

# sa.DoSensAnalysis('saveOEO_MW40.pickle', 'muw',[0.00071*0.1],ExportbySingle=True)

# sa.DoSensAnalysis('saveOEO_MW40.pickle', 'krwMAX',[0.1,0.5,0.9],ExportbySingle=True)

# sa.DoSensAnalysis('saveOEO.pickle', 'installation',['OEO'],ExportbySingle=True)



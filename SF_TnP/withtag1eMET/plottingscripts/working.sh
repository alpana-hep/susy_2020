
# ## MC-Data comparisons 1 electron CR
## For SF                                                                                                          
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,5)'


root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,5)'

## Calc FR in Data and MC as a function of different kinematics and processes
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/data",1)'
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/",0)'
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/SF_checks",3)'

## purity calculation
root -b -q 'v1_purityCalc.C("Results/PurityFact",0)'
# calculate SF
root -b -q 'CorrSF_calc.C("Results/SF_plots",0,0)'


## systematic checks
root -b -q 'v1_CorrSF_calc.C("Results/SF_plots/Checks",0,1)'

root -b -q 'CorrSF_calc.C("Results/SF_plots/Compare_FR",0,2)' ###compare FR
root -b -q 'CorrSF_calc.C("Results/SF_plots/Compare_SF",0,3)' ## compare SF btw Zee with W like selections and Tnp method
root -b -q 'CorrSF_calc.C("Results/SF_plots/withTTJets",0,4)'
root -b -q 'CorrSF_calc.C("Results/SF_plots/withWJets",0,5)'

# # ## MC-Data comparisons 1 electron CR - Stacked plots
                                                                                                          
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,0)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,1)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,2)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,3)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,4)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,5)'


# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,0)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,1)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,2)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,3)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,4)'
# # root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,5)'

root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,5)'

# ## attempts at calculating TF from SR (Znunu) to CR (Zll) - not used in this method
# root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",2)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",3)' ## for both muon and electron

# root -b -q 'SRBins_compare.C("Results/SRBins")' ## comparison of CR and SR in SR bins - not to be used yet
root -b  -q 'purityFactor_calc.C("Results/purityFactor_Calc/",0)' ## Calculation of purity factor - MC
root -b -q  'purityFactor_calc.C("Results/purityFactor_Calc/data",1)'  ## same but using Data
# root -b -q 'purityFactor_Plots.C("Results/purityFactor_Calc/overlayPlots",0)' ## Purity plots -


# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar",0)' ## SR vs CR shape comparisons
# root -b -q 'DataMC_SRvsCR_kinem.C("Results/SF_CR_McvsData",0)' ## SF calculation

root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,0)' ## Data MC comparisons in CR after weighting data with purity
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,1)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,2)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,3)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,4)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,5)'


#root -b -q 'CorrSF_calc.C("Results/SF_plots",0)'

## systematic studies
root -b -q 'DataMC_SRvsCR_kinem.C("puSysUp",1)'
root -b -q 'DataMC_SRvsCR_kinem.C("puSysDown",2)'
root -b -q 'DataMC_SRvsCR_kinem.C("JetSys_JECup",3)'
root -b -q 'DataMC_SRvsCR_kinem.C("JetSys_JECdown",4)'
root -b -q 'DataMC_SRvsCR_kinem.C("JetSys_JERup",5)'
root -b -q 'DataMC_SRvsCR_kinem.C("JetSys_JERdown",6)'
root -b -q 'DataMC_SRvsCR_kinem.C("btagSFdown",7)'
root -b -q 'DataMC_SRvsCR_kinem.C("btagSFup",8)'

root -b -q 'CompareNominalSF_system.C("puSysUp",1,"puSysUp")'
root -b -q 'CompareNominalSF_system.C("puSysDown",2,"puSysDown")'
root -b -q 'CompareNominalSF_system.C("JetSys_JECup",3,"JetSys_JECup")'
root -b -q 'CompareNominalSF_system.C("JetSys_JECdown",4,"JetSys_JECdown")'
root -b -q 'CompareNominalSF_system.C("JetSys_JERup",5,"JetSys_JERup")'
root -b -q 'CompareNominalSF_system.C("JetSys_JERdown",6,"JetSys_JERdown")'
root -b -q 'CompareNominalSF_system.C("btagSFdown",7,"btagSFdown")'
root -b -q 'CompareNominalSF_system.C("btagSFup",8,"btagSFup")'


## predictions
root -b -q 'oneDPred_compareSys.C("Results/",1,"puSysUp")'
root -b -q 'oneDPred_compareSys.C("Results/",2,"puSysDown")'
root -b -q 'oneDPred_compareSys.C("Results/",3,"JetSys_JECup")'
root -b -q 'oneDPred_compareSys.C("Results/",4,"JetSys_JECdown")'
root -b -q 'oneDPred_compareSys.C("Results/",5,"JetSys_JERup")'
root -b -q 'oneDPred_compareSys.C("Results/",6,"JetSys_JERdown")'
root -b -q 'oneDPred_compareSys.C("Results/",7,"btagSFdown")'
root -b -q 'oneDPred_compareSys.C("Results/", 8,"btagSFup")'



root -b -q 'oneDPred_compareSys.C("Results/",1,9,"SF_uncert_up")'
root -b -q 'oneDPred_compareSys.C("Results/",1, 10,"SF_uncert_down")'


root -l 'finaloneD_Pred_LL.C("Results",3,1,1,"final")'

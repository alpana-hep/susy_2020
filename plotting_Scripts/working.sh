# ## MC-Data comparisons 1 electron CR - Stacked plots
                                                                                                          
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,0)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,1)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,2)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,3)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,4)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Electron/",0,5)'


# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,0)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,1)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,2)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,3)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,4)'
# root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Muon/",1,5)'

root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/Lepton/",2,5)'

## attempts at calculating TF from SR (Znunu) to CR (Zll) - not used in this method
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",1)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",2)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TFbins_v1",3)' ## for both muon and electron

root -b -q 'SRBins_compare.C("Results/SRBins")' ## comparison of CR and SR in SR bins - not to be used yet
root -b  -q 'purityFactor_calc.C("Results/purityFactor_Calc/",0)' ## Calculation of purity factor - MC
root -b -q  'purityFactor_calc.C("Results/purityFactor_Calc/data",1)'  ## same but using Data
root -b -q 'purityFactor_Plots.C("Results/purityFactor_Calc/overlayPlots",0)' ## Purity plots -


root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar",0)' ## SR vs CR shape comparisons
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SF_CR_McvsData",0)' ## SF calculation

root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,0)' ## Data MC comparisons in CR after weighting data with purity
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,1)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,2)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,3)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,4)'
root -b -q 'SF_CRcompare_DatavsMC.C("Results/data_mcComparisons/SF_withPurity",2,5)'


#root -b -q 'CorrSF_calc.C("Results/SF_plots",0)'


## For fakerate
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots",0)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/HEm_veto",1)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/L1TrigProb",2)'


# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic",0)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/HEM_veto",1)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/L1TrigProb",2)'

# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar",0)'
# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/HEM_veto",1)'
root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/",0)'

# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots",0)'
# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/HEM_veto",1)'
root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/",0)'

## Kinematics comparisons SR vs CR
#root -b 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar")'

## Transfer factors
#root -b 'SRvsCRStacked_LL_varRatio.C("Results/TransferFactors/TFbins_v1_nJets_BJets",1)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TransferFactors/TFbins_v3_phopt_qmulti",3)'

## MC-Data comparisons 1 electron CR
                                                                                                          
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,0)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,1)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,2)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,3)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,4)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,5)'


# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,5)'


# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,5)'

## validation on MC
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/",1,1,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti/v1",1,1,1)'
# root -b 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/Bins_inPhopT_Qmulti/v1",1,1,1)'
root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/",1,1,1)'


## new SR bins
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti",1,4,1)'
# extra plots
root -b 'extraPlots_CRvsSR.C("Results/ExtraPlots_validationChecks/SRvsCR",0)'
root -b 'extraPlots_valid.C("Results/ExtraPlots_validationChecks/Validation_MC",0)'
root -b 'extraPlots_2D.C("Results/ExtraPlots_validationChecks/2Dplots",0)'


                                                                            
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,0,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,1,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,2,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,3,1)'                                                                                      \

root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,4,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,5,1)'


### systematic studies

root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/pileup",3,1,"pileup_sys_up")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/pileup",3,2,"pileup_sys_down")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/btagSF",3,3,"btagSFdown")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/btagSF",3,4,"btagSFup")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/Jetsys_JER",3,6,"JetSys_JERup")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/Jetsys_JER",3,5,"JetSys_JERdown")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/Jetsys_JEC",3,8,"JetSys_JECup")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/Jetsys_JEC",3,7,"JetSys_JECdown")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/CrossSection",3,9,"CrossSecUp")'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/CrossSection",3,10,"CrossSecDown")'



###
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,1,1,"btag_SF_up")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,2,1,"btag_SF_down")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,3,1,"JEC_Sys_up")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,4,1,"JEC_Sys_down")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,5,1,"JER_Sys_up")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,6,1,"JER_Sys_down")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,7,1,"Pileup_Sys_up")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,8,1,"Pileup_Sys_down")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,9,1,"CrossSection_Sys_up")'
root -b -q 'oneDTF_compareSys.C("Results/TFcompare/",3,10,1,"CrossSection_Sys_down")'


#### For predictions
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,1,1,"btag_SF_up")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,2,1,"btag_SF_down")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,3,1,"JEC_Sys_up")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,4,1,"JEC_Sys_down")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,5,1,"JER_Sys_up")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,6,1,"JER_Sys_down")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,7,1,"Pileup_Sys_up")'
root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,8,1,"Pileup_Sys_down")'
# root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,9,1,"CrossSection_Sys_up")'
# root -b -q 'oneDPred_compareSys.C("Results/Pred_compare/",3,10,1,"CrossSection_Sys_down")'


root -b -q 'finaloneD_Pred_LL.C("Results/",3,1,1,"btag_SF_up")'

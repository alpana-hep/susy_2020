## to fix overflow in ratio plot -- set last bin as overflow before you take a ratio of SR/CR

## data MC comparisons
# root -l 'plotAlps_RatioPlots.C("path to save plots",which lepton (1 for electron, 0 muon, 2 for e+mu),which year (0=2016 (pre+post VFP), 1=2017, 2=2018, 3=2016pre,4=2016post,5=fullrun2))'
## total lepton
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,0)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,1)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,2)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,3)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,4)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,5)'


## Electron                                                                                                                 
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,0)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,1)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,2)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,3)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,4)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,5)'

## Muon                                                                                                                 
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,0)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,1)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,2)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,3)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,4)'
root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,5)'

## showing TF as 1 SR and 1 CR (as it is in Bhumika's case)
                                                          
root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",1,1,1)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",2,1,1)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",3,1,1)'

# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",1,2,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",2,2,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",3,2,1)'

# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",1,3,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",2,3,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",3,3,1)'

# ## without qapplying trigger effciency
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/woTrigEff",3,1,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/woTrigEff",2,1,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/woTrigEff",1,1,1)'

## calculating TF for lost leptons v1 - 8 bins

root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",3,1,1)'
root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",2,1,1)'
root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",1,1,1)'
                                                                                                     
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",3,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",2,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",1,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",3,3,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",2,3,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",1,3,1)'

## TF vs different categories for lost lepton (failing acceptance, isolation and identification)
root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v1",1,1,1)'
root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v1",2,1,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v2",1,2,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v2",2,2,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v3",1,3,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v3",2,3,1)'


## TF in varbins for different kinematics
root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",1,3,1)'
root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",2,3,1)'
root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",3,3,1)'
## validation of TF on MC in search bins
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",1,1,1)' #electron for TFBins v1
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",2,1,1)' # muon
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",3,1,1)' # all leptons

                                                                                              
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",1,2,1)' #electron for TFBins v2                                                                          
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",2,2,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",3,2,1)'

                                                                                              
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",1,3,1)' #electron for TFBins v2                                                                          
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",2,3,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",3,3,1)'


## validation of TF in MC - different kinematic variables
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",3,3,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",2,3,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",1,3,1)'

                                                                      
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",3,2,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",2,2,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",1,2,1)'

                                                                      
root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",3,1,1)'
root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",2,1,1)'
root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",1,1,1)'



## prediction in data and MC expected comparisons in search bins
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,0,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,1,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,2,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,3,1)'                                                                                           
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,4,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,5,1)'                                                                              
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,0,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,1,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,2,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,3,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,4,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,5,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,0,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,1,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,2,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,3,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,4,1)'
root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,5,1)'

                                                                                      
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,5,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,5,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,5,2)'


                                                                                      
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,5,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,5,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,5,3)'


#####
                                                                        
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,0,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,1,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,2,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,3,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,4,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",3,5,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,0,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,1,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,2,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,3,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,4,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",2,5,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,0,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,1,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,2,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,3,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,4,1)'
root -b -q 'plot_DataValida_SRBins_Stacked.C("Validation_data/Kinematics/TFbins_v1",1,5,1)'


## validation of TF in Data and MC - different kinematic variables
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,0,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,1,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,2,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,3,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,4,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,5,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,0,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,1,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,2,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,3,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,4,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,5,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,0,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,1,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,2,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,3,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,4,1)'
root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,5,1)'

                                                                              
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,5,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,5,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,5,2)'

# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,5,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,5,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,5,3)'


## adding copying steps to CERN cluster

## copying data-MC comparisons

scp data_mcComparisons/Lepton_LL_* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/data_mcComparisons/
scp data_mcComparisons/Electron_LL_* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/data_mcComparisons/Electron
scp data_mcComparisons/Muon_LL_* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/data_mcComparisons/Muon

## transfer factors

scp TransferFactors/SRvsCR/TFbins_v1/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/SRvsCR/TFbins_v1
scp TransferFactors/SRvsCR/TFbins_v2/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/SRvsCR/TFbins_v2
scp TransferFactors/SRvsCR/TFbins_v3/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/SRvsCR/TFbins_v3


scp TransferFactors/SRvsCR/TFbins_v1/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/SRvsCR/TFbins_v1/Electron
scp TransferFactors/SRvsCR/TFbins_v2/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/SRvsCR/TFbins_v2/Electron
scp TransferFactors/SRvsCR/TFbins_v3/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/SRvsCR/TFbins_v3/Electron

scp TransferFactors/SRvsCR/TFbins_v1/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/SRvsCR/TFbins_v1/Muon
scp TransferFactors/SRvsCR/TFbins_v2/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/SRvsCR/TFbins_v2/Muon
scp TransferFactors/SRvsCR/TFbins_v3/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/SRvsCR/TFbins_v3/Muon


## event cateogiry wise
scp TransferFactors/EventCateg_TF/TFbins_v1/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/EventCateg_TF/TFbins_v1/Muon
scp TransferFactors/EventCateg_TF/TFbins_v2/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/EventCateg_TF/TFbins_v2/Muon
scp TransferFactors/EventCateg_TF/TFbins_v3/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/EventCateg_TF/TFbins_v3/Muon

scp TransferFactors/EventCateg_TF/TFbins_v1/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/EventCateg_TF/TFbins_v1/Electron
scp TransferFactors/EventCateg_TF/TFbins_v2/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/EventCateg_TF/TFbins_v2/Electron
scp TransferFactors/EventCateg_TF/TFbins_v3/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/EventCateg_TF/TFbins_v3/Electron

scp TransferFactors/EventCateg_TF/TFbins_v1/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/EventCateg_TF/TFbins_v1
scp TransferFactors/EventCateg_TF/TFbins_v2/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/EventCateg_TF/TFbins_v2
scp TransferFactors/EventCateg_TF/TFbins_v3/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/EventCateg_TF/TFbins_v3


## transfer factors - electron muon tau
scp TransferFactors/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/TFbins_v1/Muon
scp TransferFactors/TFbins_v2/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/TFbins_v2/Muon
scp TransferFactors/TFbins_v3/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TransferFactors/TFbins_v3/Muon

scp TransferFactors/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/TFbins_v1/Electron
scp TransferFactors/TFbins_v2/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/TFbins_v2/Electron
scp TransferFactors/TFbins_v3/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostElectron_bkg/TransferFactors/TFbins_v3/Electron

scp TransferFactors/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/TFbins_v1/Lepton
scp TransferFactors/TFbins_v2/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/TFbins_v2/Lepton
scp TransferFactors/TFbins_v3/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostLepton_bkg/TransferFactors/TFbins_v3/Lepton


## TF in variables
scp TF_inKinematics/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TF_inKinematics/Muon
scp TF_inKinematics/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TF_inKinematics/Electron
scp TF_inKinematics/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/TF_inKinematics/


## validation data
scp Validation_data/Kinematics/TFbins_v1/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v1/
scp Validation_data/Kinematics/TFbins_v2/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v2
scp Validation_data/Kinematics/TFbins_v3/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v3

scp Validation_data/Kinematics/TFbins_v1/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v1/Electron
scp Validation_data/Kinematics/TFbins_v2/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v2/Electron
scp Validation_data/Kinematics/TFbins_v3/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v3/Electron
scp Validation_data/Kinematics/TFbins_v1/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v1/Muon
scp Validation_data/Kinematics/TFbins_v2/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v2/Muon
scp Validation_data/Kinematics/TFbins_v3/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/Kinematics/TFbins_v3/Muon

## SR bins
scp Validation_data/SRbins/TFbins_v1/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v1/
scp Validation_data/SRbins/TFbins_v2/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v2
scp Validation_data/SRbins/TFbins_v3/Lepton_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v3

scp Validation_data/SRbins/TFbins_v1/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v1/Electron
scp Validation_data/SRbins/TFbins_v2/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v2/Electron
scp Validation_data/SRbins/TFbins_v3/Electron_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v3/Electron
scp Validation_data/SRbins/TFbins_v1/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v1/Muon
scp Validation_data/SRbins/TFbins_v2/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v2/Muon
scp Validation_data/SRbins/TFbins_v3/Muon_LL* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_data/SRbins/TFbins_v3/Muon

## validation in MC
scp Validation_MC/Sbins/TFBIns_v1/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v1/
scp Validation_MC/Sbins/TFBIns_v2/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v2
scp Validation_MC/Sbins/TFBIns_v3/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v3

scp Validation_MC/Sbins/TFBIns_v1/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v1/Electron
scp Validation_MC/Sbins/TFBIns_v2/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v2/Electron
scp Validation_MC/Sbins/TFBIns_v3/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v3/Electron
scp Validation_MC/Sbins/TFBIns_v1/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v1/Muon
scp Validation_MC/Sbins/TFBIns_v2/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v2/Muon
scp Validation_MC/Sbins/TFBIns_v3/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Sbins/TFBIns_v3/Muon



## kinematics
scp Validation_MC/Kinematics_Valid/TFBins_v1/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v1/
scp Validation_MC/Kinematics_Valid/TFBins_v2/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v2
scp Validation_MC/Kinematics_Valid/TFBins_v3/Lepton_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v3

scp Validation_MC/Kinematics_Valid/TFBins_v1/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v1/Electron
scp Validation_MC/Kinematics_Valid/TFBins_v2/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v2/Electron
scp Validation_MC/Kinematics_Valid/TFBins_v3/Electron_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v3/Electron
scp Validation_MC/Kinematics_Valid/TFBins_v1/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v1/Muon
scp Validation_MC/Kinematics_Valid/TFBins_v2/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v2/Muon
scp Validation_MC/Kinematics_Valid/TFBins_v3/Muon_LL*W+TTBar* kalpana@lxplus9.cern.ch:/eos/user/k/kalpana/www/folder/HGCAL_TDAQ/Plots/Susy_Analysis/lostMuon_bkg/Validation_MC/Kinematics_Valid/TFBins_v3/Muon




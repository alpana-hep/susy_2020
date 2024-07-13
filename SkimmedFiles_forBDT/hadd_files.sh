#!/bin/bash                                                                                                                                      
#path=/store/user/kalpana/ul_rootop_Analys_May2/v1 phoID_loose_runList_Summer16v3_signal_T5bbbbZg_pt100_MET200_job0.root

# hadd -f Signal_T5bbbbZg_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/signal | grep 'phoID_loose_runList_Summer16v3_signal_T5bbbbZg'`
# hadd -f Signal_T5qqqqHg_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/signal | grep 'phoID_loose_runList_Summer16v3_signal_T5qqqqHg'`

# hadd -f Signal_T5ttttZg_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/signal | grep 'phoID_loose_runList_Summer16v3_signal_T5ttttZg'`

# hadd -f Signal_T6ttZg_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/signal | grep 'phoID_loose_runList_Summer16v3_signal_T6ttZg'`


hadd -f Signal_TChiWG_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24\
/signal | grep 'phoID_loose_runList_Summer16v3_signal_TChiWG'`
hadd -f Signal_TChiNG_v18_PhoIdloose_phopt40_BL_skimsForBDT.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/signal | grep 'phoID_loose_runList_Summer16v3_signal_TChiNG'`

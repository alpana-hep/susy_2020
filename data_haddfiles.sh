#!/bin/bash                                                                                                                                      
path=/FR_UL_Oct24/Skimmed//store/user/kalpana/FR_UL_Oct24/Skimmed
hadd -f out_Data_UL2018_Run2018D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2018_Run2018A_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
hadd -f out_Data_UL2018_Run2018B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
hadd -f out_Data_UL2018_Run2018C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

hadd -f out_Data_UL2017_Run2017B_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2017_Run2017C_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2017_Run2017D_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2017_Run2017E_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2017_Run2017F_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


hadd -f out_Data_UL2016_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016_Run2016G_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
hadd -f out_Data_UL2016_Run2016H_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loose_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`

# ### sys studies

# hadd -f out_Data_UL2018_Run2018D_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_looseJetSys_JECup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_looseJetSys_JECup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_looseJetSys_JECup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_looseJetSys_JECup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_looseJetSys_JECup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_looseJetSys_JECup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECup_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2018_Run2018D_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_looseJetSys_JECdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_looseJetSys_JECdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_looseJetSys_JECdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_looseJetSys_JECdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_looseJetSys_JECdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_looseJetSys_JECdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JECdown_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`


# ## JER up
# hadd -f out_Data_UL2018_Run2018D_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_looseJetSys_JERdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_looseJetSys_JERdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_looseJetSys_JERdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_looseJetSys_JERdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_looseJetSys_JERdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_looseJetSys_JERdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERdown_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2018_Run2018D_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_looseJetSys_JERup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_looseJetSys_JERup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_looseJetSys_JERup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_looseJetSys_JERup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_looseJetSys_JERup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_looseJetSys_JERup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_looseJetSys_JERup_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`



# ### btagSF
# hadd -f out_Data_UL2018_Run2018D_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_loosebtagSFdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_loosebtagSFdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_loosebtagSFdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_loosebtagSFdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_loosebtagSFdown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET2
# 00'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_loosebtagSFdown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFdown_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`



# hadd -f out_Data_UL2018_Run2018D_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_loosebtagSFup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_loosebtagSFup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_loosebtagSFup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_loosebtagSFup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_loosebtagSFup_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_loosebtagSFup_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosebtagSFup_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`


# ### pileup sys
# hadd -f out_Data_UL2018_Run2018D_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_loosepuSysDown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_loosepuSysDown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_loosepuSysDown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_loosepuSysDown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_loosepuSysDown_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_loosepuSysDown_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysDown_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`



# hadd -f out_Data_UL2018_Run2018D_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018A_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018B_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2018_Run2018C_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_loosepuSysUp_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_loosepuSysUp_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_loosepuSysUp_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_loosepuSysUp_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_loosepuSysUp_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`


# hadd -f out_Data_UL2016_Run2016F_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_loosepuSysUp_pt40_MET200.root `xrdfsls -u /store/user/kalpana/FR_UL_Oct24/Skimmed | grep 'phoID_loosepuSysUp_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`




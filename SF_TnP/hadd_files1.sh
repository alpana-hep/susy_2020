#!/bin/bash                                                                                                                                      
path=/store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp



hadd -f out_Data_UL2018_Run2018A_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2018_Run2018A'`

hadd -f out_Data_UL2018_Run2018B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2018_Run2018B'`

hadd -f out_Data_UL2018_Run2018C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2018_Run2018C'`

hadd -f out_Data_UL2018_Run2018D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2018_Run2018D'`

hadd -f out_Data_UL2017_Run2017B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2017_Run2017B'`
hadd -f out_Data_UL2017_Run2017C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2017_Run2017C'`
hadd -f out_Data_UL2017_Run2017D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2017_Run2017D'`
hadd -f out_Data_UL2017_Run2017E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2017_Run2017E'`
hadd -f out_Data_UL2017_Run2017F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2017_Run2017F'`

hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016APV_Run2016B'`
hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016APV_Run2016C'`
hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016APV_Run2016D'`
hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016APV_Run2016E'`
hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016APV_Run2016F'`

hadd -f out_Data_UL2016_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016_Run2016F'`
hadd -f out_Data_UL2016_Run2016G_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016_Run2016G'`
hadd -f out_Data_UL2016_Run2016H_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_UL2016_Run2016H'`






















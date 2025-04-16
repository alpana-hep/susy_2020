#!/bin/bash
#/store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET
# hadd -f Summer20UL18_ZNuNuGJets_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL18_ZNuNuGJets_v1_pt100_MET200'`
# hadd -f Summer20UL17_ZNuNuGJets_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL17_ZNuNuGJets_v1_pt100_MET200'`
# hadd -f Summer20UL16_ZNuNuGJets_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL16_ZNuNuGJets_v1_pt100_MET200'`
# hadd -f Summer20UL16APV_ZNuNuGJets_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_ZNuNuGJets_v1_pt100_MET200'`

# hadd -f Summer20UL18_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL18_ZJetsToNuNu_HT_v1_pt100_M
# ET200'`
# hadd -f Summer20UL17_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL17_ZJetsToNuNu_HT_v1_pt100_MET200'`
# hadd -f Summer20UL16_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL16_ZJetsToNuNu_HT_v1_pt100_M
# ET200'`
# hadd -f Summer20UL16APV_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp/Skims | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_ZJetsToNuNu_HT_v1_p
# t100_MET200'`


hadd -f Summer20UL18_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL18_ZLLGJets_MonoPhoton'`
hadd -f Summer20UL17_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL17_ZLLGJets_MonoPhoton'`
hadd -f Summer20UL16APV_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL16APV_ZLLGJets_MonoPhoton'`
hadd -f Summer20UL16_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL16_ZLLGJets_MonoPhoton'`


hadd -f Summer20UL18_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL18_DYJets_Mt50'`
hadd -f Summer20UL17_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL17_DYJets_Mt50'`
hadd -f Summer20UL16_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL16_DYJets_Mt50'`
hadd -f Summer20UL16APV_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_loose_skimmed_runList_Summer20UL16APV_DYJets_Mt50'`



# hadd -f Summer20UL18_ZLLGJets_MonoPhoton_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL18_ZLLGJets_MonoPhoton'`
# hadd -f Summer20UL17_ZLLGJets_MonoPhoton_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL17_ZLLGJets_MonoPhoton'`
# hadd -f Summer20UL16APV_ZLLGJets_MonoPhoton_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL16APV_ZLLGJets_MonoPhoton'`
# hadd -f Summer20UL16_ZLLGJets_MonoPhoton_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL16_ZLLGJets_MonoPhoton'`


# hadd -f Summer20UL18_DYJets_Mt50_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL18_DYJets_Mt50'`
# hadd -f Summer20UL17_DYJets_Mt50_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL17_DYJets_Mt50'`
# hadd -f Summer20UL16_DYJets_Mt50_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL16_DYJets_Mt50'`
# hadd -f Summer20UL16APV_DYJets_Mt50_PhoIdlooseSF_check_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_1tageAsMET | grep 'phoID_looseSF_check_skimmed_runList_Summer20UL16APV_DYJets_Mt50'`



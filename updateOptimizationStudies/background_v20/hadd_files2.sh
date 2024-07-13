#!/bin/bash                                                                                                                                      
path=/store/user/kalpana/ul_rootop_Analys_May2/v1


hadd -f Summer20UL18_ZNuNuGJets_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_ZNuNuGJets'`
hadd -f Summer20UL17_ZNuNuGJets_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_ZNuNuGJets'`
hadd -f Summer20UL16_ZNuNuGJets_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_ZNuNuGJets'`
hadd -f Summer20UL16APV_ZNuNuGJets_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_ZNuNuGJets'`

hadd -f Summer20UL18_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_ZJetsToNuNu_HT'`
hadd -f Summer20UL17_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_ZJetsToNuNu_HT'`
hadd -f Summer20UL16_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_ZJetsToNuNu_HT'`
hadd -f Summer20UL16APV_ZJetsToNuNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_ZJetsToNuNu_HT'`

hadd -f Summer20UL18_QCD_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_QCD_HT'`
hadd -f Summer20UL17_QCD_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_QCD_HT'`
hadd -f Summer20UL16_QCD_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_QCD_HT'`
hadd -f Summer20UL16APV_QCD_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_QCD_HT'`

hadd -f Summer20UL18_GJets_DR-0p4_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_GJets_DR-0p4'`
hadd -f Summer20UL17_GJets_DR-0p4_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_GJets_DR-0p4'`
hadd -f Summer20UL16_GJets_DR-0p4_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_GJets_DR-0p4'`
hadd -f Summer20UL16APV_GJets_DR-0p4_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_GJets_DR-0p4'`

hadd -f Summer20UL17_TTJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_TTJets_inc'`
hadd -f Summer20UL16_TTJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_TTJets_inc'`
hadd -f Summer20UL18_TTJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_TTJets_inc'`
hadd -f Summer20UL16APV_TTJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_TTJets_inc'`


hadd -f Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_WJetsToLNu_HT'`
hadd -f Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_WJetsToLNu_HT'`
hadd -f Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_WJetsToLNu_HT'`
hadd -f Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_WJetsToLNu_HT'`


hadd -f Summer20UL17_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_WGJets_MonoPhoton_PtG-130'`
hadd -f Summer20UL16_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_WGJets_MonoPhoton_PtG-130'`
hadd -f Summer20UL18_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-130'`
hadd -f Summer20UL16APV_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_WGJets_MonoPhoton_PtG-130'`

hadd -f Summer20UL17_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_WGJets_MonoPhoton_PtG-40to130'`
hadd -f Summer20UL16_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_WGJets_MonoPhoton_PtG-40to130'`
hadd -f Summer20UL18_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-40to130'`
hadd -f Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130'`



hadd -f Summer20UL17_TTJets_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_TTJets_HT'`
hadd -f Summer20UL16_TTJets_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_TTJets_HT'`
hadd -f Summer20UL18_TTJets_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_TTJets_HT'`
hadd -f Summer20UL16APV_TTJets_HT_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_TTJets_HT'`

hadd -f Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16APV_TTGJets_inc'`
hadd -f Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL16_TTGJets_inc'`
hadd -f Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL17_TTGJets_inc'`
hadd -f Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_BL_BDTwithT5bbbbZg_7variables.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24 | grep 'phoID_loose_runList_skimmed_Summer20UL18_TTGJets_inc'`

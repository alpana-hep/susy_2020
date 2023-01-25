#!/bin/bash                                                                                                                                    

path=/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/updatedbkgNtuples_skimmed_wphotPt20_MET100_Njets2_Jan2023/updatedbkgNtuples_skimmed_wphotPt20_MET100_Njets2_Jan2023
#/store/user/kalpana/updatedbkgNtuples_skimmed_wphotPt20_MET100_Njets2_Jan2023/ #/store/user/kalpana/myProduction/ 
#/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/GMSB_lowPtPhot_2021/Samples/Wjets/skims/


for year in  Autumn18
do
    ls ${path} | grep SR_${year}.WGJets > runList_WGJets_${year}_v18.txt
    ls ${path} | grep SR_${year}.WJets > runList_WJets_${year}_v18.txt

    #ls ${path}/SR_${year}.WGJets_MonoPhoton_PtG-*_TuneCP5_13TeV-madgraph_v18.root > runList_WGJets_${year}_v18.txt

    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep SR_Autumn18.WGJets_MonoPhoton_PtG >runList_skimmed_Autumn18_lowPhopT_WGJetsPtG_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep SR_Autumn18.WJetsToLNu_HT >runList_skimmed_Autumn18_lowPhopT_WJetsToLNu_HT_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WGJets_MonoPhoton_PtG-130 >runList_Autumn18_lowPhopT_WGJetsPtG-130_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WGJets_MonoPhoton_PtG-40to130 >runList_Autumn18_lowPhopT_WGJetsPtG-40to130_v18.txt
   #  xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-100To200 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-200To400 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-400To600 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-600To800 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-800To1200 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-1200To2500 >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_v18.txt
#     xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-2500ToInf >runList_Autumn18_lowPhopT_Autumn18.WJetsToLNu_HT-2500ToInf_TuneCP5\
# _13TeV-madgraphMLM-pythia8_v18.txt
    #xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-600To800 >temp_runList_Autumn18_lowPhopT_WJetsToLNu_HT_600to800_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-800To1200 >temp_runList_Autumn18_lowPhopT_WJetsToLNu_HT_800to1200_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-1200To2500 >temp_runList_Autumn18_lowPhopT_WJetsToLNu_HT_1200to2500_v18.txt
    # xrdfs root://cmseos.fnal.gov/ ls ${path} | grep WJetsToLNu_HT-2500ToInf >temp_runList_Autumn18_lowPhopT_WJetsToLNu_HT_2500toInf_v18.txt


    # ls ${path}/${year}.TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_v18.root > runList_TTGJets_${year}_v18.txt
    # ls ${path}/SR_${year}.TTJets_HT*.root > runList_TTJets-HT_${year}_v18.txt
    # ls ${path}/SR_${year}.TTJets_*Lept*.root > runList_TTJets-inc_${year}_v18.txt
    # eosls ${path}/${year}.WGJets_MonoPhoton_PtG-*.root > runList_WGJets_${year}_v18.txt
    # eosls ${path}/${year}.WJetsToLNu_HT-*.root > runList_WJets_${year}_v18.txt
    # ls ${path}/SR_${year}.ZJetsToNuNu_HT*.root > runList_ZJets_${year}_v18.txt
    # ls ${path}/SR_${year}.ZNuNuGJets_MonoPhoton_PtG-*.root > runList_ZGJets_${year}_v18.txt
    # ls ${path}/SR_${year}.QCD_HT*.root > runList_QCD_${year}_v18.txt
    # ls ${path}/SR_${year}.GJets_DR*.root > runList_GJets_DR_${year}_v18.txt
done

# year=Summer16v3
   
# ls ${path}/SR_${year}.TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v18.root > runList_TTGJets_${year}_v18.txt
# ls ${path}/SR_${year}.TTJets_HT*.root > runList_TTJets-HT_${year}_v18.txt
# ls ${path}/SR_${year}.TTJets_*Lept*.root > runList_TTJets-inc_${year}_v18.txt
# ls ${path}/SR_${year}.WGJets_MonoPhoton_PtG-*_TuneCUETP8M1_13TeV-madgraph_v18.root > runList_WGJets_${year}_v18.txt
# ls ${path}/SR_${year}.WJetsToLNu_HT-*.root > runList_WJets_${year}_v18.txt
# ls ${path}/SR_${year}.ZJetsToNuNu_HT*.root > runList_ZJets_${year}_v18.txt
# ls ${path}/SR_${year}.ZNuNuGJets_MonoPhoton_PtG-*.root > runList_ZGJets_${year}_v18.txt
# ls ${path}/SR_${year}.QCD_HT*.root > runList_QCD_${year}_v18.txt
# ls ${path}/SR_${year}.GJets_DR*.root > runList_GJets_DR_${year}_v18.txt

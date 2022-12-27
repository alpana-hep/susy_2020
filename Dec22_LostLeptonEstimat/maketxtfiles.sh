#!/bin/bash                                                                                                                                    

path=/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/GMSB_lowPtPhot_2021/Samples/Wjets/skims/


for year in Fall17 Autumn18
do
    ls ${path}/SR_${year}.TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_v18.root > runList_TTGJets_${year}_v18.txt
    ls ${path}/SR_${year}.TTJets_HT*.root > runList_TTJets-HT_${year}_v18.txt
    ls ${path}/SR_${year}.TTJets_*Lept*.root > runList_TTJets-inc_${year}_v18.txt
    ls ${path}/SR_${year}.WGJets_MonoPhoton_PtG-*_TuneCP5_13TeV-madgraph_v18.root > runList_WGJets_${year}_v18.txt
    ls ${path}/SR_${year}.WJetsToLNu_HT-*.root > runList_WJets_${year}_v18.txt
    ls ${path}/SR_${year}.ZJetsToNuNu_HT*.root > runList_ZJets_${year}_v18.txt
    ls ${path}/SR_${year}.ZNuNuGJets_MonoPhoton_PtG-*.root > runList_ZGJets_${year}_v18.txt
    ls ${path}/SR_${year}.QCD_HT*.root > runList_QCD_${year}_v18.txt
    ls ${path}/SR_${year}.GJets_DR*.root > runList_GJets_DR_${year}_v18.txt
done

year=Summer16v3
   
ls ${path}/SR_${year}.TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v18.root > runList_TTGJets_${year}_v18.txt
ls ${path}/SR_${year}.TTJets_HT*.root > runList_TTJets-HT_${year}_v18.txt
ls ${path}/SR_${year}.TTJets_*Lept*.root > runList_TTJets-inc_${year}_v18.txt
ls ${path}/SR_${year}.WGJets_MonoPhoton_PtG-*_TuneCUETP8M1_13TeV-madgraph_v18.root > runList_WGJets_${year}_v18.txt
ls ${path}/SR_${year}.WJetsToLNu_HT-*.root > runList_WJets_${year}_v18.txt
ls ${path}/SR_${year}.ZJetsToNuNu_HT*.root > runList_ZJets_${year}_v18.txt
ls ${path}/SR_${year}.ZNuNuGJets_MonoPhoton_PtG-*.root > runList_ZGJets_${year}_v18.txt
ls ${path}/SR_${year}.QCD_HT*.root > runList_QCD_${year}_v18.txt
ls ${path}/SR_${year}.GJets_DR*.root > runList_GJets_DR_${year}_v18.txt

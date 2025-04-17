#!/bin/sh
#analyzeLightBSM 800 127 FastSim T5qqqqHg
executable=$1
gluinoMass=$2
nlspMass=$3
anaArg=$4
outName=$5
#ext=$6
hist=$6 
#h_Sbins_v6_withOnlyBL_Selec_HT1TeV_Met250
hist1=$6
currDir=$(pwd)

echo "root://cmseos.fnal.gov//store/group/lpcsusyhad/SusyPhotonMET/Run2ProductionV18/signals/unskimmed_SortedSignalscans/scan/${anaArg}/${outName}_${anaArg}_${gluinoMass}_${nlspMass}_unSkimmed.root" 

echo "root://cmseos.fnal.gov//store/group/lpcsusyhad/SusyPhotonMET/Run2ProductionV18/signals/unskimmed_SortedSignalscans/scan/${anaArg}/${outName}_${anaArg}_${gluinoMass}_${nlspMass}_unSkimmed.root" > runFileList.txt

######################################
# SETUP CMSSW STUFF...
######################################


echo "======== ${PWD} ==============="

source /cvmfs/cms.cern.ch/cmsset_default.sh
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#  export SCRAM_ARCH=slc7_amd64_gcc820
# scram p CMSSW CMSSW_11_1_0_pre3
# cd CMSSW_11_1_0_pre3/src
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src


# export SCRAM_ARCH=slc7_amd64_gcc700
# scram p CMSSW CMSSW_10_2_21
# cd CMSSW_10_2_21/src
eval `scramv1 runtime -sh`
echo $CMSSW_RELEASE_BASE
cd $currDir
echo $currDir

outRootFile="${anaArg}_${outName}_${gluinoMass}_${nlspMass}_v18"
./$executable runFileList.txt ${outRootFile}.root Run2 ${anaArg}signal loose
#exit
#_${outName}_2016
# export SCRAM_ARCH=el9_amd64_gcc12
# scram p CMSSW CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src

rm runFileList.txt
#rm -rf CMSSW_11_1_0_pre3/
# mv runFileList.txt tmp/runFileList.txt

echo "========== HiggsAnalysis/CombinedLimit ==============="
cd $currDir
echo "${PWD}"

# export SCRAM_ARCH=slc7_amd64_gcc700
# cmsrel CMSSW_10_2_13
# cd CMSSW_10_2_13/src
# eval `scramv1 runtime -sh`
# cmsenv
# echo "${PWD}"
# echo "ls curr dir"
                                           
cd ${currDir}
echo "curr dir"
ls ${currDir}
rm -rf dataCards/

echo "making datacards"
mkdir dataCards
mkdir dataCards/${outRootFile}_${hist1}
echo "mkdir dataCards/${outRootFile}_${hist1}"
root -l -q -b 'PredmakeDatacard_SBins.C('${gluinoMass}','${nlspMass}',"'${outRootFile}'.root","'${hist1}'","'${hist}'","'${anaArg}'")'
#exit
#rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin1.txt
#exit
#!/bin/bash
#a=$"Hello i am pass";
# if [ $(echo $hist1 | grep -c "MET_200") -gt 0 ]
# then
#     echo "Success"
#     for i in {35..52}
#     do
# 	rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin${i}.txt
#     done
# elif [ $(echo $hist1 | grep -c "MET_100") -gt 0 ]
# then
#     echo "Success"
#     for i in {40..52}
#     do
#         rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin${i}.txt
#     done

# elif [ $(echo $hist1 | grep -c "MET_300") -gt 0 ]
# then
#     echo "Success"
#     for i in {28..52}
#     do
#         rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin${i}.txt
#     done

# else
#     echo "Fail";
# fi
#exit
# for i in {35..52}
# do
#     rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin${i}.txt
# done

# for i in {0..2}
# do
#     rm dataCards/${outRootFile}_${hist1}/${outRootFile}_${hist}_bin${i}.txt
# done

echo "ls dataCards/${outRootFile}_${hist1}/*"
ls dataCards/${outRootFile}_${hist1}/*

rm ${outRootFile}.root



# ls ${currDir}
pwd
# cd /uscms/home/kalpana/nobackup/public/work/Susy_lowPho_analysis/CMSSW_14_0_0_pre0/src/CAT_workingDirectory/SignalRegionGetLimits/datacards/
#cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
#cd /uscms/home/kalpana/nobackup/public/work/Susy_lowPho_analysis/CMSSW_14_0_0_pre0/src/CAT_workingDirectory/SignalRegionGetLimits/datacards/
tar -xvf ${currDir}/higgsAnalysis.tar
cd HiggsAnalysis/CombinedLimit
scramv1 b clean
scramv1 b
echo "cmssw"
                                           
#cd HiggsAnalysis/CombinedLimit #CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit
eval `scramv1 runtime -sh`
echo $CMSSW_RELEASE_BASE
pwd
#combine --help
echo "combineCards.py ${currDir}/dataCards/${outRootFile}_${hist1}/*.txt > dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt"
combineCards.py ${currDir}/dataCards/${outRootFile}_${hist1}/*.txt > dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt

echo "calculating limit"
mH="$(echo "${gluinoMass}+${nlspMass}*0.0001" | bc)"
echo $mH
#mH=2200.020
# echo "combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -t -1 -n ${outName}_${gluinoMass}_${nlspMass}_${hist} -m ${mH}"
# combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -t -1 -n ${outName}_${gluinoMass}_${nlspMass}_${hist} -m ${mH}

                  
echo "combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -n ${outName}_${gluinoMass}_${nlspMass}_${hist} -m ${mH}"
combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -t -1 -n ${outName}_${gluinoMass}_${nlspMass}_${hist} -m ${mH}

#xrdcp -f higgsCombine${outName}_${gluinoMass}_${nlspMass}_${hist}*.root root://cmseos.fnal.gov//store/user/bkansal/myProduction/limits_rootout/T5qqqqHg/${hist1}/
#xrdcp -f higgsCombine${outName}_${gluinoMass}_${nlspMass}_${hist}*.root root://cmseos.fnal.gov//store/user/bkansal/myProduction/limits_rootout/T6ttZg/${hist1}/
# cd /uscms/home/kalpana/nobackup/public/work/Susy_lowPho_analysis/CMSSW_14_0_0_pre0/src/CAT_workingDirectory/v10_combinetool/CMSSW_14_1_0_pre4/src/SignalRegion_GetLimits
# mkdir /eos/uscms/store/user/kalpana/Susy_phoMet/limit_rootout/v17_June2024_optimization/${anaArg}/${hist}/
xrdcp -f higgsCombine${outName}_${gluinoMass}_${nlspMass}_${hist}*.root root://cmseos.fnal.gov//store/user/kalpana/Susy_phoMet/limit_rootout/v17_June2024_optimization/${anaArg}/
#higgsCombine${outRootFile}_${gluinoMass}_${nlspMass}_${hist}*.root root://cmseos.fnal.gov//store/user/kalpana/Susy_phoMet/limit_rootout/T5bbbbZg/${hist1}/
rm higgsCombine${outName}_${gluinoMass}_${nlspMass}_${hist}*.root
#cmssw-el7 -- exit
#exit

#rm *.txt
cd ${currDir}
rm -rf dataCards/
#rm datacards/*txt
# cd ${currDir}
# rm signal_T5bbbbZg_*.root

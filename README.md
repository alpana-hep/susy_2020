# susy_2020

##To make exclusion plots for strong productions :
Links : 
https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/

*copy files from 
Git clone -b  GetSignalLimits git@github.com:alpana-hep/susy_2020.git .

mkdir SignalRegion_GetLimits
Cd SignalRegion_GetLimits
Copy AnalyzeLightBSM.cc, ​​AnalyzeLightBSM.h, NtupleVariables.cc, NtupleVariables.h, Makefile, bkg samples root files after analyzer code to current directory.
Make
Setup higgs combined tool (one time setup):
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv

git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

cd HiggsAnalysis/CombinedLimit

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
scramv1 b clean; scramv1 b # always make a clean build

cd $CMSSW_BASE/src/
tar cvf higgsAnalysis.tar HiggsAnalysis
cp higgsAnalysis.tar ../../.
rm -rf  CMSSW_10_2_13/

mkdir datacards
Script for making datacards :
makeDatacard_SBins.C  .
//Check script and make sure file name and directories are correctly mentioned
root -l 'makeDatacard_SBins.C(2350,10,"Out_T5bbbbZG_2350_10_v18.root", "h_Sbins_v6_withOnlyBL_Selec","h_Sbins_v6_withOnlyBL_Selec")'
Arguements : Mg, mX, output file containing all the histograms after basleine selections, histogram names


Commands :

./worker_SP.sh analyzeLightBSM 2200 200 FastSim T5bbbbZg Summer16v3 <folder name>
./worker_SP.sh analyzeLightBSM 2350 10 FastSim T5bbbbZG Summer16v3 h_Sbins_v6_withOnlyBL_Selec

Condor job submission for making data-cards and calculating r value for various mass points :

./calcLimit.sh T5bbbbZg_MassScan.txt

Make limit plots : 
Check :
https://github.com/alpana-hep/susy_2020/tree/Exc_plotsLimit







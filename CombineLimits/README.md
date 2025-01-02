# This repository contains the instructions to run and setup the combine tool
## Setup Combine tool 
```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.0
scramv1 b clean; scramv1 b
tar cvf higgsAnalysis.tar HiggsAnalysis

```
This should setup the working directory for combined
```
mkdir SignalRegion_GetLimits
cd SignalRegion_GetLimits
git clone -b  GetSignalLimits git@github.com:alpana-hep/susy_2020.git .
cp AnalyzeLightBSM.* SignalRegion_GetLimits && cp NupleVariables.* SignalRegion_GetLimits && cp Makefile SignalRegion_GetLimits && bkg samples root files after analyzer code to current directory.
make
cp  ../HiggsAnalysis/CombinedLimit/higgsAnalysis.tar .
```
Create data cards -
```
mkdir datacards
makeDatacard_SBins.C
```
Check script and make sure file name and directories are correctly mentioned
```
root -l -q -b 'makeDatacard_SBins.C(${gluinoMass},${nlspMass},'${outRootFile}.rootq,'${hist1}','${hist}')'
```
Example to run the script - 
```
root -l -q -b 'makeDatacard_SBins.C(2200,200,”out_T5bbbbZg_2200_200.root”,”h_Sbins_LL_MET_200”)’
```
Where hist and hist1 are the name of the histograms with different search bin
To run the full chain - 
```
./worker_SP.sh analyzeLightBSM 2200 200 T5bbbbZg Summer16v3   h_Sbins_LL_MET_200
```

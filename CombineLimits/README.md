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
To submit for full mass scan
```
./calcLimit.sh list_T6ttZg.txt T6ttZg Summer16v3 h_Sbins_LL_v4_MET_200_withMvaCut 1000 BL_BDTwith13variables_T6ttZg TMVAClassification_T6ttZg_Phopt40_MET200_13va\
riables_200trees_2maxdepth.weights.xml
```
In case you have different file names or use files for background with different BDT training
```
source temp.sh T5bbbbZg
```

Once you get the limits calculated, to get the plots, follow these instructions
```
cd plotLimits
```
First hadd files for all mass points
```
hadd -f higgsCombineSummer16v3_TChiWG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut_v1.root `xrdfsls -u  /store/user/kalpana/Susy_phoMet/limit_rootout/v17_June2024_optimization/TChiWG/ | grep 'h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.AsymptoticLimits.'`
```
Then run the analyzer script for strong production
```
ls ${Path}/higgsCombineSummer16v3_T5ttttZg_h_Sbins_LL_newSbins_v3_MET_200_v1.root >input_combine_T5ttttZg_h_Sbins_LL_newSbins_v3_MET_200.txt
./plotlimit input_combine_T5ttttZg_h_Sbins_LL_newSbins_v3_MET_200.txt out_T5ttttZg_combine_h_Sbins_LL_newSbins_v3_MET_200_v1.root T5ttttZg
root -b -q 'getExclusion.C("out_T5ttttZg_combine_h_Sbins_LL_newSbins_v3_MET_200_v1.root")'


```
For electroweakino samples
```
python3 limitplotter_TChiWG_combine.py higgsCombineSummer16v3_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root
```




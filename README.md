To copy samples from fermilab -
```
xrdcp -r root://cmseos.fnal.gov//store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/InputFiles/ .
```

# Machinery to analyze the input files
NtupleVariables.h - header file where tree branches are defined
AnalyzerLightBSM.cc - where event level analysis happens

## Instructions to run the code -
```
git clone -b Skeleton  https://github.com/alpana-hep/susy_2020.git .
make 
```
Arguments for ./analyzeLightBSM :
```
./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
filelist: see under the 'inputFiles' directory
outfile: as you want to name your file
year : which year to process
process: MC smaples for which job is running or data for data files

### Example to run the script

```
./analyzeLightBSM infile.txt out_wjets.root 2018 WJets

```

# Machinery to train BDTs
use script 'BDT_strong_TTJets.py'
Description of arguments is provided in the script.

Example to run the script -
```
python3 BDT_strong_TTJets.py -y all -s1 T5gg_2200_deltaM10 -b1 FullRun2  -m RandS -nt 200 -md 2 -n Equalweight_T5gg_2200_DeltaM10_v1phopt40_MET200_13variables -cuts 'MET>200'
```
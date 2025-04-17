## Description - this repository contains the codes/instructions to run optimization, bdt trainings, and combined tools

## Instructions to run the code -
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_0_0_pre0
cd CMSSW_14_0_0_pre0/src
eval `scramv1 runtime -sh`
git clone -b UL_BDTOptmizationStudies  https://github.com/alpana-hep/susy_2020.git .
```
## Prepare skims for BDT

You can use the already skimmed files stored on fermilab - all the baseline selections are applied
```
/store/user/kalpana/ul_rootop_Analys_May2/OptmizationStudies_June24/BDTTraining_July24/InputFiles/
```
Or create for yourself following these steps -

1. GMSB signal (v17) skims for training the BDT
```
cd SkimmedFiles_forBDT
cp Filelist/*.txt .
make
./analyzeLightBSM  runList_unSkimmed_Autumn18_T5bbbbZg_2200_10.txt out_T5bbbbZg_2200_10.root Run2 T5bbZgsignal loose
```
Or submit for all the input mass points -
```
source runInt.sh
```
combine the output of different mass points in different models
```
source hadd_files.sh
```

2. Skims for T5gg files

```
cd SkimmedFiles_forBDT/v20_T5ggSkims
make
source submany.sh
source haddFiles.sh
```
3. Skims for v20 background samples
```
cd SkimmedFiles_forBDT/background_v20/
cp Filelist*.txt .
source submitMany1.sh && source submitMany2.sh
source combine.sh
```

## Training a boosted decision tree
```
cd BDT_training/
python3 BDT_strong_TTJets.py -y all -s1 T5gg_2200_deltaM10 -b1 FullRun2  -m RandS -nt 200 -md 2 -n Equalweight_T5gg_2200_DeltaM10_v1phopt40_MET200_13variables -cuts 'MET>200'
```
Description of the arguments to BDT_strong_TTJets.py
```
-y = year; "all" for Run2 luminosity
-s1 = signal string - indicate which signal files to pick; "T5bbbbZg" or "T5gg" pr "TChiWg" or etc.
-b1 = bakcground string - indicating which background files to pick; "FullRun2"
-m = RandS; (not used explicitly in the code)- can be ignored
-nt = number of trees in BDT;
-md = maximum depth in each tree BDT;
-n = name string; can be anything you want
-cuts = string containing logical expression of preselection to implement before training a BDT on the files; it can be "MET>200 && photon_pt>40" or something like this
```

See "runme.sh" for more such combinations submitted

All the xml files saved in this directory are the output of training done on individual signal models

## optimization studies - baseline selection plus BDT response
In this, we create histogram based on different preselections used and compared background and signal

These are also the files and script used as input while running combine tool (see later instructions)

```
cd updateOptimizationStudies
make
```
This script reads xml file of BDT weights based on the input training features used.

And depending on your signal process string provided in the job, it will read a specific XML file - more details please check in AnalyzeLightBSM.cc

For v17 GMSB signal models - follow the instructions listed below
```
./analyzeLightBSM  runList_unSkimmed_Autumn18_T5bbbbZg_2200_10.txt out_T5bbbbZg_2200_10.root Run2 T5bbbbZgsignal loose
```
OR submit jobs for few mass points in different signal models at once
```
source runInt.sh
```

For v20 T5gg model - follow the instructions listed below :
```
cd  v20_T5ggModel/
make
source submany.sh
```

For v20 Background - follow the instructions listed below :
```
cd background_v20/
cp ../*.xml .
make
source submitMany2.sh
source hadd_files2.sh
source combine.sh
```

The output from the above step can be used to create stacked plots from thesis - overlay plots etc..

Follow these instructions to produce such plots
```
cd PlottingScripts
source working.sh
```


## Running combine tool and making limit plots -
```
cd CombineLimits
```
Description of different files in the directory -
```
makeDatacard_SBins.C - create the data cards
make sure you don't give empty datacards , combine won't run.
worker_SP.sh - shell scripts run your analyzer script over signal models files one by one, create datacards, and run combine tool and save the output tree in a root file

combine.sh - submit the jobs for all the grid points and take input the mass scan
```

Setup the combine tool using the instruction given here-  https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/

```                                                                                                                                                                  cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.1
scramv1 b clean; scramv1 b
cd CMSSW_14_1_0_pre4/src
tar -xvf  higgsAnalysis.tar HiggsAnalysis
```

Commands to execute-
```
make
./analyzeLightBSM  runList_unSkimmed_Autumn18_T5bbbbZg_2200_10.txt out_T5bbbbZg_2200_10.root Run2 T5bbZgsignal loose

```

Create data cards-
```
root -l -q -b 'makeDatacard_SBins.C(2700,1600,"T5qqqqHg_Summer16v3_2700_1600_v18.root","h_Sbins_LL_newSbins_v7_MET_200","h_Sbins_LL_newSbins_v7_MET_200","T5qqqqHg")'
```

Run a job interactively which runs over a signal mass point, create data cards and run the combine tool

```
./worker_SP.sh >executable> <mg> <mnlsp> <signal model> <Extension to read the file > <hist name>
```

Example -
```
./worker_SP.sh analyzeLightBSM 2200 200 T5bbbbZg Summer16v3   h_Sbins_LL_MET_200

```
Submit jobs for a signal model - taking list of different mass point as input 
```
./calcLimit.sh <mas scan txt file> <signal model> <extension> <hist>
```

Example -
```
./calcLimit.sh T5bbbbZg_MassScan.txt T5bbbbZg Summer16v3 h_Sbins_LL_MET_200

```

Combine all the files after running the combine tool successfully -

```
cd plotLimits
source hadd_files.sh
./plotlimit in_file.txt out.root T5bbbbZg
root -b 'getExclusion.C("out.root")'

```
To make all the plots -
```
cd PlottingScripts
source runme.sh
```
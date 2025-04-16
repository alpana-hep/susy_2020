### Z invisible background estimation

## Instructions to run the code -
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_0_0_pre0
cd CMSSW_14_0_0_pre0/src
eval `scramv1 runtime -sh`
git clone -b UL_ZinvisibleEstimation  https://github.com/alpana-hep/susy_2020.git .
make
cp inputfiles/*.txt .
./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>

```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.

photon ID: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)

Also, along with photon ID, this string also indicate which systematic you are studying, and corresponding which TF to use, For example -
```
loose --> looseJetSys_JECup for JET sys studies - JEC and up
loose --> looseJetSys_JECdown for	JET sys	studies	- JEC and down
```
See sys_submit.sh for more details

filelist: see under the 'inputfiles' directory

outfile: as you want to name your file

process: MC smaples for which job is running or data for data files


There are some flags in Analyzer code which one needs to be careful as they are switching on and off some of the corrections. Brief description is given below :
```
applyTrgEff=true - to apply trigger efficiency
applyHEMveto=true - to apply HEM veto to 2-18 and 2017
applyL1TrigFire_prob=true - to apply L1 trigger prefire correction to 2016 and 2017
applyPUwt = true - apply pileup weights
applybTagSFs=true - apply btag SF
applysys=false - apply or not do systematic studies - should be false in default case - true when you are calculating TF for that systematic
```

Example to run the machinery -

```
./analyzeLightBSM inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt out_DYJetsToLL_M-50_HT-100to200.root 2017 DYJetsToLL_M-50UL loose
```

To submit the condor jobs:

executable is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted.

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

```
root -l -q 'splitRunList.C("runList_Summer20UL17_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","WJetsToLNu_HT-100To200UL","loose")'

```

sys_submit.sh - submit multiple jobs at a time for MC - dilepton
```
source sys_submit.sh
```
Hadding the output of above jobs (please change the path)
```
source hadd_files_sys.sh
```

For data - submission of jobs in bulk
```
source submitMany2.sh
```
combining the output
```
source  hadd_files_dataele.sh && source hadd_files_datamu.sh
```
For ZNuNu process - submission of jobs in bulk
```
source pred_sys_submit.sh
```
combining the output -
```
source Predsys_hadd.sh
```

clean*.sh - to clean the log files of the condor jobs

run 

```
source cleanupBatchfiles.sh

```

To make the plots the plotsse combine.sh to hadd files and get the overall MC and data contribution, and use  working.sh to get the plots (it also is explaining the role of each script)

```
cd plottingScripts
source combine.sh
source working.sh

```

The input skims files for MC -
```
/store/user/lpcsusyphotons/kalpana/SkimsUL_June2023/
```
Unskimmed data and MC files
```
/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/
```

The output files for this studies are on lxplus
```
/eos/user/k/kalpana/Susy_outputFiles/Zinv_bkg
```
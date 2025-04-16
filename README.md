# Fake rate background estimation

## Instructions to run the code -
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_0_0_pre0
cd CMSSW_14_0_0_pre0/src
eval `scramv1 runtime -sh`
git clone -b UL_FakeRate_bkgEstimation  https://github.com/alpana-hep/susy_2020.git .
make
./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.

<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'

(Note-  first three are cutbased ID recommended by Egamma group.)

Also, along with photon ID, this string also indicate which systematic you are studying, and corresponding which TF to use, For example

```
loose --> looseJetSys_JECup for JET sys studies - JEC and up
loose --> looseJetSys_JECdown for	JET sys	studies	- JEC and down
```

filelist: see under the 'inputFiles' directory

outfile: as you want to name your file

process: MC smaples for which job is running or data for data files - "UL" should always be added in the end of this string

To run a single job interactively -
```
./analyzeLightBSM inputfiles/runList_Summer20UL18_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8.txtp out_Summer20UL18_TTJets_HT_1200to2500.root 2018  TTJets_HTUL  loose

```
There are some flags in Analyzer code which one needs to be careful as they are switching on and off some of the corrections. Brief description is given below :
```
applyTrgEff=true - to apply trigger efficiency
applyHEMveto=true - to apply HEM veto to 2-18 and 2017
applyL1TrigFire_prob=true - to apply L1 trigger prefire correction to 2016 and 2017
bool apply_pixelveto=true; - apply pixel veto
applyPUwt = true - apply pileup weights
applybTagSFs=true - apply btag SF
applysys=false - apply or not do systematic studies - should be false in default case - true when you are calculating TF for that systematic
```

Default TF of FR file to be read from -
```
Electron_FR_TFbins_v3_phopt_qmulti_phoID_loose_09Jan24.root
```
Default - SF
```
out_SF_FR_Data_MC_Default.root
```

To submit the condor jobs:

executable is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

```
root -l -q 'splitRunList.C("runList_Summer20UL18_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",25,"2018","WJetsToLNu_HT-100To200UL","loose")'

```
OR add this string for all MC data samples in a shell script and submit all at once submit multiple jobs at a time.

skimmed_submit.sh - submit multiple jobs at a time for MC
```
source skimmed_submit.sh
```
to hadd the output files
```
source skim_hadd.sh
```
* please update the path for files as to what one is using

To submit jobs for data (full sample - data files are not skimmed) (for systematics studies - comment out the rest of the jobs in it)
```
source submitMany_data.sh
```
to hadd the output files
```
source data_haddfiles.sh
```


clean*.sh - to clean the log files of the condor jobs

run 

```
source cleanupBatchfiles.sh

```
To make the plots: use combine.sh to hadd files and get the overall MC and data contribution, and use  working.sh to get the plots (it also is explaining the role of each script)
```
cd plottingScripts
source combine.sh
source wroking.sh

```
Please make sure the directory path exists in your area.

The input skims files for MC -
```
/store/user/lpcsusyphotons/kalpana/SkimsUL_June2023/FR/
```
Unskimmed data and MC files
```
/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/
```

The output files for this studies are on lxplus
```
/eos/user/k/kalpana/Susy_outputFiles/FR_bkg
```

## Scale factor calculation
We will be using Single electron and Egamms dataset instead of MET and DY Jets (ZLL Gamma )samples
```
cd SF_TnP/
make
./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.

photon ID: 'loose', 'medium', 'tight','mva_wp90','mva_p80'

(Note-  first three are cutbased ID recommended by Egamma group.)

filelist: see under the 'inputFiles' directory

outfile: as you want to name your file

process: MC smaples for which job is running or data for data files

Example to run the machinery -
```
./analyzeLightBSM inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt out_DYJetsToLL_M-50_HT-100to200.root 2017 DYJetsToLL_M-50UL loose
```

To submit the condor jobs:

executable is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted.

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

```
root -l -q 'splitRunList.C("runList_Summer20UL18_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",25,"2018","WJetsToLNu_HT-100To200UL","loose")'

```

skim_submit.sh - submit multiple jobs at a time for MC 
```
source skim_submit.sh
```
for hadding the files
```
source hadd_files_skim.sh
```

For data
```
source submitMany1.sh
```

To hadd files
```
source hadd_files1.sh
```
clean*.sh - to clean the log files of the condor jobs

run

```
source cleanupBatchfiles.sh

```

To make the plotsse combine.sh to hadd files and get the overall MC and data contribution, and use  working.sh to get the plots (it also is explaining the role of each script)
```
cd SF_TnP/plottingScripts
source combine.sh
source working.sh

```
working.sh contains the detasils of the scripts

Please make sure the directory path exists in your area.

The input skims files for MC -
```
/store/user/lpcsusyphotons/kalpana/SkimsUL_June2023/FR/
```
Unskimmed data and MC files
```
/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/
```

The output files for this studies are on lxplus
```
/eos/user/k/kalpana/Susy_outputFiles/FR_bkg/SF_TnP
```

### SF studies with considering 1 tag electron as MET
```
cd SF_TnP/withtag1eMET
make
```

skim_submit.sh - submit multiple jobs at a time for MC
```
source skim_submit.sh
```
for hadding the files
```
source hadd_files_skim.sh
```
For WJets & TTJets
```
source submitMany2.sh
```
To hadd files
```
source hadd_files1.sh
```
For data
```
source submitMany1.sh
```

To hadd files
```
source hadd_files1.sh
```


To make the plotsse combine.sh to hadd files and get the overall MC and data contribution, and use  working.sh to get the plots (it also is explaining the role of each script)
```
cd SF_TnP/withtag1eMET/plottingScripts
source combine.sh
source working.sh

```
working.sh contains the details of the scripts

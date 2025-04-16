#Lost lepton background estimation

## Instructions to run the code -
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_0_0_pre0
cd CMSSW_14_0_0_pre0/src
eval `scramv1 runtime -sh`
git clone -b LostLeptonbkg_studies https://github.com/alpana-hep/susy_2020.git .
make (run make everytime you change anything in any of the source/header file)
./analyzeLightBSM <filelist> <outfile> <year> <process> <which_lepton> <photon ID>
```

Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples. - not needed anymore.

<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group and last two are MVA based IDs)
Also, along with photon ID, this string also indicate which systematic you are studying, and corresponding which TF to use, For example
```
loose --> looseJetSys_JECup for JET sys studies - JEC and up
loose --> looseJetSys_JECdown for	JET sys	studies	- JEC and down
```
See full list of examples in submitMany1.sh

<file_list>: containing the path for the samples to run over - see under inputFiles directory
<year>: which year dataset you are running on : "2016postVFP", "2016preVFP", "2017","2018"
<process>: MC sample name string or data for all year data files - should always contain a string "UL"
<which_lepton>: "Electron" or "Muon" 

Example to run the job interactively for a case
```
./analyzeLightBSM inputFiles/runList_skimmed_Summer20UL16_TTGJets_inc.txt out_Summer20UL16_TTGJets_inc_v20_lostElectron.root 2016postVFP TTGJetsUL Electron loose
```
OR for muon background
```
./analyzeLightBSM inputFiles/runList_skimmed_Summer20UL16_TTGJets_inc.txt out_Summer20UL16_TTGJets_inc_v20_lostMuon.root 2016postVFP TTGJetsUL Muon loose

```

There are some flags in Analyzer code which one needs to be careful as they are switching on and off some of the corrections. Brief description is given below :
```
applyTrgEff=true - to apply trigger efficiency
applyHEMveto=true - to apply HEM veto to 2-18 and 2017
applyL1TrigFire_prob=true - to apply L1 trigger prefire correction to 2016 and 2017
applyPUwt = true - apply pileup weights
applybTagSFs=true - apply btag SF
applysys=false - apply or not do systematic studies - should be false in default case - true when you are calculating TF for that systematic
```
Default TF file to be read from -
```
Lepton_LL_TFv7_HT_bjets_phopT_PhoIdloose_phoID_loose_09Jan24.root
```

To submit the condor jobs:
<executable> is 'worker2.sh' (change or add destination path for output files in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the input files which you want to transfer which your code will be using interactively)
To submit multiple jobs for a given samples at a time:
```
root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2018","WGJets_MonoPhoton_PtG-40to130UL","Electron","loose")'

OR add this string for all MC data samples in a shell script and submit all at once  -- submitMany1.sh - submit multiple jobs at a time.

```
source submitMany1.sh
```
the above files also has jobs for systematic studies (which are commented out)
To hadd the output
```
source hadd_files_final.sh
```
* please update the path for files as to what one is using

To submit jobs for data (full sample - data files are not skimmed) (for systematics studies - comment out the rest of the jobs in it)
```
source submitMany_data.sh
```
to hadd the output files
```
source hadd_files_data.sh
```

clean*.sh - to clean the log files of the condor jobs
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
/store/user/lpcsusyphotons/kalpana/SkimsUL_June2023/
```
Unskimmed data and MC files
```
/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/
```

The output files for this studies are on lxplus
```
/eos/user/k/kalpana/Susy_outputFiles/LL_bkg
```


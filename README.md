## Script to skim UL samples -

# Instructions to run the code -
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
scram p CMSSW CMSSW_14_0_0_pre0
cd CMSSW_14_0_0_pre0/src
eval `scramv1 runtime -sh`
git clone -b ULSkims_FinalV3 https://github.com/alpana-hep/susy_2020.git .
cd LL_bkgskims_withbtagEff or cd FR_bkg_skims_withbtagEff or cd DY_skim/ 
make 
./analyzeLightBSM <filelist> <outfile> <year <process> Electron <photon ID>
```
String 'Electron' is not used inside the code anywhere, so keep it constant.

Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.

photon ID: 'loose', 'medium', 'tight','mva_wp90','mva_p80'

(Note-  first three are cutbased ID recommended by Egamma group.)

Example to run the script
```
./analyzeLightBSM runList_Summer20UL18_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt out.root 2018 DYJetsUL Electron loose
```

To submit the condor jobs:

<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

submitMany1.sh - submit multiple jobs at a time. 

clean*.sh - to clean the log files of the condor jobs

Submit multiple condor jobs
```
source submitMany2.sh (or other submit shell scripts as per the directory)
```
Note: with this method you can only submit 1 skimming job per file. Can not merge multiple files.

to combine the files -
```
source hadd_files1.sh
```
# Script to skim UL samples -

## Instructions to run the code -
1. git clone -b ULSkims_FinalV3 https://github.com/alpana-hep/susy_2020.git .
2. cd LL_bkgskims_withbtagEff or cd FR_bkg_skims_withbtagEff (whichever skims one wants to create)
3. make (run make after you change anything in any of the source/header file)
4. ./analyzeLightBSM <filelist> <outfile> <year <process> Electron <photon ID>

String 'Electron' is not used inside the code anywhere, so keep it constant.
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.
<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)

Example to run the script
```
./analyzeLightBSM temp_QCD.txt out_Summer20UL18_WJetsToLNu_HT-100To200_v20.root 2018 WJetsUL Electron loose
```

To submit the condor jobs:
<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)
submitMany2.sh - submit multiple jobs at a time. 
clean*.sh - to clean the log files of the condor jobs

Submit multiple condor jobs
```
source submitMany2.sh (or other submit shell scripts)
```
Note: with this method you can only submit 1 skimming job per file. Can not merge multiple files.

Example to run the script interactively
```
./analyzeLightBSM temp_QCD.txt out_Summer20UL18_WJetsToLNu_HT-100To200_v20.root 2018 WJetsUL Electron loose
```
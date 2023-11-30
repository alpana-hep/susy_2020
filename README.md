## This code is for implementing different photon IDs (reading from the ntupls) & study its impact on prompt vs non prompt photons.

## Instructions to run the code -
1. git clone -b UL_photonIDStudy https://github.com/alpana-hep/susy_2020.git .
2. make (run make after you change anything in any of the source/header file)
3. ./analyzeLightBSM <filelist> <outfile> <year <process> <Lepton estimation you are interested in> <photon ID>

Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.

<Lepton> : 'Electron' or 'Muon'
<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)

To submit the condor jobs:
<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)
submitMany2.sh - submit multiple jobs at a time.
clean*.sh - to clean the log files of the condor jobs

Example to run the script
```
./analyzeLightBSM temp_QCD.txt out_Summer20UL18_WJetsToLNu_HT-100To200_v20.root 2018 WJetsUL Electron loose
```


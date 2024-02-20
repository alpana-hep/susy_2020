#Lost lepton background estimation

## Instructions to run the code -
```
1. git clone -b LostLeptonbkg_studies https://github.com/alpana-hep/susy_2020.git .
2. make (run make everytime you change anything in any of the source/header file)
3. ./analyzeLightBSM <filelist> <outfile> <year> <process> <which_lepton> <photon ID>
```

Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.
<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)

<file_list>: containing the path for the samples to run over - see under inputFiles directory
<year>: which year dataset you are running on : 2016postVFP, 2016preVFP, 2017,2018
<process>: MC sample string or data for all year data files
<which_lepton>: Electron or Muon

To submit the condor jobs:
<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)
To submit multiple jobs for a given samples at a time:
```
root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2018","WGJets_MonoPhoton_PtG-40to130UL","Electron","loose")'

```
submitMany1.sh - submit multiple jobs at a time.
OR submit all at a time:
```
source submitMany1.sh
```

clean*.sh - to clean the log files of the condor jobs

```
source cleanupBatchfiles.sh
```

Example to run the script interactively (I suggest not to do until or unless only one file you are looping over)

```
./analyzeLightBSM inputFiles/runList_skimmed_Summer20UL16_TTGJets_inc.txt out_Summer20UL16_TTGJets_inc_v20_lostElectron.root 2016postVFP TTGJetsUL Electron loose
```
OR for muon background
```
./analyzeLightBSM inputFiles/runList_skimmed_Summer20UL16_TTGJets_inc.txt out_Summer20UL16_TTGJets_inc_v20_lostMuon.root 2016postVFP TTGJetsUL Muon loose

```

To make the plots:
```
cd plottingScripts
source wroking.sh

```
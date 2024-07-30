# Machinery to analyze the input files
## Instructions to run the code -
```
1. git clone -b UL_optmizationStudies  https://github.com/alpana-hep/susy_2020.git .
2. make (run make after you change anything in any of the source/header file)
3. ./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in
 'map_crosssection_SMprocess_v1.txt' and should contain UL in it if you are analyzing UL ntuples.
<filelist>: see under the 'inputFiles' directory
<outfile>: as you want to name your file
<year> : which year to process
<process>: MC smaples for which job is running or data for data files
Example to run the script

```
./analyzeLightBSM infile.txt out_wjets.root 2018 WJets

```

# Machinery to train BDTs

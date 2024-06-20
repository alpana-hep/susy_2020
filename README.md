# Running combine tool and making limit plots
makeDatacard_SBins.C - create the data cards
make sure you don't give empty datacards , combine won't run.
worker_SP.sh - shell scripts run your analyzer script over signal models files one by one, create datacards, and run combine tool and save the output tree in a root file

combine.sh - submit the jobs for all the grid points and take input the mass scan

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
git clone -b UL_CombineTool_Limits  https://github.com/alpana-hep/susy_2020.git .
make
```
Create data cards-
```
root -l -q -b 'makeDatacard_SBins.C(2200,200,'${outRootFile}.root','${hist1}','${hist}')'
```

Where hist and hist1 are the search bins  histogram to be read

Run a job interactively -
```
./worker_SP.sh >executable> <mg> <mnlsp> <signal model> <Extension to read the file > <hist name>
```

Example -
```
./worker_SP.sh analyzeLightBSM 2200 200 T5bbbbZg Summer16v3   h_Sbins_LL_MET_200

```
Submit jobs for a signal model
```
./calcLimit.sh <mas scan txt file> <signal model> <extension> <hist>
```
Example -
```
./calcLimit.sh T5bbbbZg_MassScan.txt T5bbbbZg Summer16v3 h_Sbins_LL_MET_200

```

Combine all the files after running the combine tool successfully
```
cd plotLimits
source hadd_files.sh
./plotlimit in_file.txt out.root T5bbbbZg
root -b 'getExclusion.C("out.root")'

```

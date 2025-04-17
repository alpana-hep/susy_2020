# This code is for implementing all baseline selections and checking their impact on different MC samples.

## Instructions to run the code -
```
git clone -b UL_BDTOptimizationStuides https://github.com/alpana-hep/susy_2020.git .
make 
./analyzeLightBSM input_gl2200_X1600_deltaM10.txt out_T5gg_2200_1600_delM10.root FullRun2 T5bbbbZg_pprovasignalUL loose
```

To submit the condor jobs:
<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 

spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)
submitMany2.sh - submit multiple jobs at a time.
clean*.sh - to clean the log files of the condor jobs

Example to run the script
For signals  
```
source submany.sh
```

For backgrounds
```
cd background_v20
source submitmany2.sh
```
Please check the content of these .sh file and comment in and out as per the need.

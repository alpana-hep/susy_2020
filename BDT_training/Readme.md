# Machinery to train BDTs
use script 'BDT_strong_TTJets.py'
Description of arguments is provided in the script.
Example to run the script -
```
python3 BDT_strong_TTJets.py -y all -s1 T5gg_2200_deltaM10 -b1 FullRun2  -m RandS -nt 200 -md 2 -n Equalweight_T5gg_2200_DeltaM10_v1phopt40_MET200_13variables -cuts 'MET>200'
```

To open GUI tab after training -
```
root -l
TMVA::TMVAGUI("./BDT_output_all_Equalweight_T5gg_2200_DeltaM100_phopt40_MET200_13variables_200trees_2maxdepth.root")
```

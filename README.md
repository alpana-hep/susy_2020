## Instructions

*you may need to change the path of the directory in maketxtfiles.sh.


To run the Analyzer code, follow these Commands :
```
make

./analyzeLightBSM runList/<txt file name> <output rootfile name> <dataset name>

```

for example : 
```
./analyzeLightBSM runList_TTGJets_2017_v17.txt TTGJets_2017_v17.root TTG_v17_2016
```

After getting output root files, you can plot the distributions using plotKinStack.C code. 

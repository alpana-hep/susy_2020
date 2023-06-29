# Script to read UL ntuples & make plots to compare with pre legacy 
## (code for pre-legacy is under the directory "pre-legacy")
to run the script -
./analyzeLightBSM <runlist.txt> <outputfile> <year> <dataset_name> <Lepton for_which_you_are_analysing > <phoID>

<dataset> should contain UL in it if you are analyzing UL ntuples and should not contain UL while analyzing pre-legacy ntuples.

Example to run the script -
./analyzeLightBSM temp_QCD.txt out_Summer20UL18_WJetsToLNu_HT-100To200_v20.root 2018 WJetsUL Electron loose


The folder 'plotting_scripts' contains the script to make the plots to compare pre-legacy and UL

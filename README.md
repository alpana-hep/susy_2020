# BDT Training

To run the script
python BDT_strong_v3.py -y 2016 -s1 T5bbbbZg  -s2 T5bbbbZg -b TTJets -m RandS -nt 200 -md 2

See BDT_strong_v3.py for more details on the arguements

And also change the path for input file txt accordingly in BDT_strong_v3.py

After running this script you get the root files containing all the information for training. To run the example macros of TMVA,
root -l 
root[0] TMVA::TMVAMultiClassGui("path/to/tmva/output/file.root")
You will get the GUI corresponding to this command and which can generate all the plots by one click

# susy_2020
Code can be run after invoking 'make'
# command line
./analyzeLightBSM filelist outputfile dataset sample_tag
  
filelist - .txt file containing the path of file for the sample

dataset - 2016,2017,2018 (for respective background samples)
              - signal (for signal samples)
      
Sample_tag - options are "QCD_Jets" & "GJets_DR" or "random"
          to combine QCD and Gamma+jet samples.


The path for the files should be replaced depending upon the cluster you are working on. The current txt files contains the path for TIFR2 (india-cms) cluster

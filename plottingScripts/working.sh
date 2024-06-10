## plots for different modesl -
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5bbbbZg","T5bbbbZg","2500")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5bbbbZg","T5bbbbZg","2400")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5bbbbZg","T5bbbbZg","2300")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5bbbbZg","T5bbbbZg","2200")'

root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5qqqqHg","T5qqqqHg","2500")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5qqqqHg","T5qqqqHg","2400")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5qqqqHg","T5qqqqHg","2300")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5qqqqHg","T5qqqqHg","2200")'

root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2500")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2600")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2400")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2300")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2200")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T5ttttZg","T5ttttZg","2100")'

root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T6ttZg","T6ttZg","1300")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T6ttZg","T6ttZg","1400")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T6ttZg","T6ttZg","1600")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/T6ttZg","T6ttZg","1800")'

root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/TChiWG","TChiWG","0")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Signals/TChiNG","TChiNG","0")'

# background
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/WGJets","WGJets","0")'
root -b	-q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/WJets","WlnuJets","0")'
root -b	-q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/TTGJets","ttbarG","0")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/TTJets","ttbarJets","0")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/GJets","GJets","0")'
root -b -q 'overlay_1file.C("Results/withBaselineSelection/Backgrounds/ZNuNuGJets","ZnunuGJets","0")'
## overlay Plots
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5bbbbZg","T5bbbbZg","2500")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5bbbbZg","T5bbbbZg","2400")'

root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5bbbbZg","T5bbbbZg","2300")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5bbbbZg","T5bbbbZg","2200")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5qqqqHg","T5qqqqHg","2500")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5qqqqHg","T5qqqqHg","2400")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5qqqqHg","T5qqqqHg","2300")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5qqqqHg","T5qqqqHg","2200")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5ttttZg","T5ttttZg","2500")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5ttttZg","T5ttttZg","2400")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5ttttZg","T5ttttZg","2300")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T5ttttZg","T5ttttZg","2200")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T6ttZg","T6ttZg","1300")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T6ttZg","T6ttZg","1400")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T6ttZg","T6ttZg","1600")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/T6ttZg","T6ttZg","1800")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/TChiWG","TChiWG","0")'
root -b -q 'overlay_multifile.C("Results/withBaselineSelection/Signals/OverlayPlots/TChiNG","TChiNG","0")'
## stacked plots
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5bbbbZg","T5bbbbZg","2300")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5bbbbZg","T5bbbbZg","2400")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5bbbbZg","T5bbbbZg","2500")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5bbbbZg","T5bbbbZg","2200")'

root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5qqqqHg","T5qqqqHg","2500")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5qqqqHg","T5qqqqHg","2400")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5qqqqHg","T5qqqqHg","2300")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5qqqqHg","T5qqqqHg","2200")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5ttttZg","T5ttttZg","2500")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5ttttZg","T5ttttZg","2400")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5ttttZg","T5ttttZg","2300")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T5ttttZg","T5ttttZg","2200")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T6ttZg","T6ttZg","1300")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T6ttZg","T6ttZg","1400")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T6ttZg","T6ttZg","1600")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/T6ttZg","T6ttZg","1800")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/TChiWG","TChiWG","0")'
root -b -q 'StackPlots_multifile.C("Results/withBaselineSelection/Signals/StackedPlots/TChiNG","TChiNG","0")'





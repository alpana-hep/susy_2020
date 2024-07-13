#Do relevent imports
from ROOT import TChain, gROOT, TFile, TCut, TMVA, TLorentzVector
import numpy as np
import sys
import argparse

#  An example of running this would be
#  python BDT_strong_v3.py -y 2016 -s1 SMS-T5Wg_m19' -s2 SMS-T6Wg_m17 -b GJets_DR -m RandS -nt 200 -md 10 -n sometext
#  This would train a BDT with filenames listed in filelist_mvaprepped.txt that contain the string "SMS-T5Wg_m19" or "SMS-T6Wg_m17" used as the signal and filenames that contain the string "GJets_DR" as the background.
#  The resulting BDT would have 200 trees and a maxdepth of 10.
#  The output file would have the name "BDT_output_2016_sometext_200trees_10maxdepth.root". 
#  You can make the "filelist_mvaprepped.txt" file by doing 'xrdfsls /store/group/target/directory > filelist_mvaprepped.txt'.

gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description="Training BDT to separate real and fake MET")



parser.add_argument('-s1', '--signal1', dest='sig_name1', default='SMS-T5Wg_m19', help='This is the sample used for signal, ie ZGGtonunuGG ')
#parser.add_argument('-s2', '--signal2', dest='sig_name2', default='SMS-T6Wg_m17', help='This is the sample used for signal, ie ZGGtonunuGG ')
parser.add_argument('-y', dest='year', default='2016', help='This is the year or you can put "all" for full Run2 or "phase1" for only Phase1')
parser.add_argument('-b1', '--background1', dest='back_name1', default='GJets_DR', help='This is the sample used for the backbaround, ie GJets or QCD')
#parser.add_argument('-b2', '--background2', dest='back_name2', default='QCD', help='This is the sample used for the backbaround, ie GJets or QCD')
#parser.add_argument('-b3', '--background3', dest='back_name3', default='TTJets', help='This is the sample used for the backbaround, ie GJets or QCD')
parser.add_argument('-m', '--mode', dest='mode', default='RandS', help='If training on RandS samples, you should put RandS here.')
parser.add_argument('-n', '--name', dest='name', default='test', help='This is a string appended to the output file name to serve as an identifier')
parser.add_argument('-nt', '--ntrees', dest='Ntrees', type = int, help = 'This is the number of trees in the bdt')
parser.add_argument('-md', '--maxdepth', dest='maxdepth', type = int, help = 'This is the maximum depth for each tree')
parser.add_argument('-cuts', '--cuts_selection', dest='Cuts_string', default='MET>200',help = 'Sets of preselection to be applied')

args = parser.parse_args()

year = args.year
sig_name1 = args.sig_name1
#sig_name2 = args.sig_name2
back_name1 = args.back_name1
#back_name2 = args.back_name2
#back_name3 = args.back_name3

Cuts = args.Cuts_string
#Read file lists to make TChains.  There should be one TChain for signal and one for background
sigchainlist = []
sigcountlist = []

bkgchainlist = []
bkgcountlist = []

if args.mode!='RandS':
    print( "Rebalance and smear has been selected")
    sig = TChain("TreeMaker2/PreSelection")
    bkg = TChain("TreeMaker2/PreSelection")
else:
    sig = TChain("PreSelection")
    bkg = TChain('PreSelection')

Isfastsim1 = False
Isfastsim2 = False

if 'SMS' in sig_name1:
     Isfastsim1 = True
# if 'SMS' in sig_name2:
#      Isfastsim2 = True


if '2016' in year:
     sig_name1 = 'Summer16v3Fast_{0}'.format(sig_name1)
     sig_name2 = 'Summer16v3Fast_{0}'.format(sig_name2)
     back_name = 'Summer16v3.{0}'.format(back_name)
     lumiInfb=35.922
if '2017' in year:
     sig_name1 = 'Fall17Fast.{0}'.format(sig_name1)
     sig_name2 = 'Fall17Fast.{0}'.format(sig_name2)
     back_name = 'Fall17.{0}'.format(back_name)
     lumiInfb=41.529
if '2018' in year:
     sig_name1 = 'Autumn18Fast.{0}'.format(sig_name1)
     sig_name2 = 'Autumn18Fast.{0}'.format(sig_name2)
     back_name = 'Autumn18.{0}'.format(back_name)
     lumiInfb= 59.74
if 'all' in year:
     sig_name1 = '{0}'.format(sig_name1)
     #sig_name2 = '{0}'.format(sig_name2)
     back_name1 = '{0}'.format(back_name1)
     #back_name2 = '{0}'.format(back_name2)
     #back_name3 = '{0}'.format(back_name3)
     lumiInfb= 137.19
if 'phase1' in year:
     sig_name1 = '{0}'.format(sig_name1)
     sig_name2 = '{0}'.format(sig_name2)
     back_name = '{0}'.format(back_name)
     lumiInfb= 137.19
if not Isfastsim1:
     sig_name1 = sig_name1.replace('Fast','')
# if not Isfastsim2:
#      sig_name2 = sig_name2.replace('Fast','')



flist = open('./filelist_sig_bkg.txt')
for line in flist:
    if 'phase1' in year and 'Summer16v3' in line: continue
    if sig_name1 in line: # or sig_name2 in line:
        chtemp = TChain("Pre_Selection")
        chtemp.Add(format(line.strip()))#'root://cmseos.fnal.gov/{0}'.format(line.strip()))                                                                    
        sigchainlist.append(chtemp)
        counttemp = TChain("Pre_Selection")
        counttemp.Add(format(line.strip()))#'root://cmseos.fnal.gov/{0}'.format(line.strip()))                                                                 
        sigcountlist.append(counttemp)
        print ('Loading signal file:{0}'.format(line.strip()))
    if back_name1 in line: # or back_name2 in line: # or back_name3 in line:
        chtemp = TChain("Pre_Selection")
        chtemp.Add(format(line.strip()))
        bkgchainlist.append(chtemp)
        counttemp = TChain("Pre_Selection")
        counttemp.Add(format(line.strip()))
        bkgcountlist.append(counttemp)
        print ('Loading bkg file:{0}'.format(line.strip()))

        # bkg.Add(format(line.strip()))#'root://cmseos.fnal.gov/{0}'.format(line.strip()))format(line.strip())                                                  
        # print ('Loading background file: {0}'.format(line.strip()))
        #     if 'Summer16v3.GJet_Pt' in line:                             
        #          bkg.Add('root://cmseos.fnal.gov/{0}'.format(line.strip())                                                    
flist.close()     

fout = TFile.Open("./Results/BDT_output_{3}_{0}_{1}trees_{2}maxdepth.root".format(args.name, args.Ntrees, args.maxdepth, year),"RECREATE")

TMVA.Tools.Instance()

factory = TMVA.Factory("TMVAClassification", fout, "V:!Silent:Color:Transformations=I:DrawProgressBar:AnalysisType=Classification")


loader = TMVA.DataLoader("dataset_bdt_{3}_{0}_{1}trees_{2}maxdepth".format(args.name, args.Ntrees, args.maxdepth, year))
#evtwt=0.0
#Add TChains to factory
for i, sig in enumerate(sigchainlist):
    tree= sigcountlist[i]#.Get("")
    # print(tree)    
    # for entry in tree:
    #     if(entry>1):
    #         break
    #     else:
    #         evtwt=tree.evtwt#sigcountlist[i].GetBranch("evtwt")
    #tree.SetBranchAddress("evtwt", &evtwt)
    #tree.GetEntry(i);
    # w=evtwt#.MET;
    # print('alpas',w)
    w = sigcountlist[i].GetEntries()
    print(w,'alps')
     
    # #     print 'weight for this mass point is {0}'.format(weight)
    weight = 1.0/w
    loader.AddSignalTree(sig, weight)
for i, bkg in enumerate(bkgchainlist):
    # tree= bkgcountlist[i]#.Get("")                                                                                                              
    # print(tree)
    # for entry in tree:
    #     if(entry>1):
    #         break
    #     else:
    #         evtwt=tree.evtwt
    w=bkgcountlist[i].GetEntries()
    # print(evtwt)
    weight =1.0/w
    loader.AddBackgroundTree(bkg, weight)

if args.mode != 'RandS':
#     loader.AddVariable('mva_Ngoodjets','I')
#     loader.AddVariable('mva_ST','F')
     loader.AddVariable('mva_Pt_jets/mva_ST','F')
     loader.AddVariable('abs(mva_dPhi_GG)','F')
     loader.AddVariable('mva_Photons0Et/mva_ST','F')
     loader.AddVariable('mva_Photons1Et/mva_ST','F')
     loader.AddVariable('mva_HardMET/mva_ST','F')
     loader.AddVariable('mva_Pt_GG/mva_ST','F')
     loader.AddVariable('mva_ST_jets/mva_ST','F')
     loader.AddVariable('mva_min_dPhi','F')
     loader.AddVariable('mva_dPhi1','F')
     loader.AddVariable('mva_dPhi2','F')
     loader.AddVariable('abs(mva_dPhi_GGHardMET)','F')
     loader.AddVariable('mva_dRjet1photon1','F')
     loader.AddVariable('mva_dRjet1photon2','F')
     loader.AddVariable('mva_dRjet2photon1','F')
     loader.AddVariable('mva_dRjet2photon2','F')
#     loader.AddVariable('JetsClean[0].Pt()/PhotonsLoose[0].Pt()','F')
#     loader.AddVariable('mva_Pt_jets/mva_Pt_GG')
     loader.AddSpectator('mva_Ngoodjets', 'I')
#     loader.AddSpectator('EventWeight', 'F')

#     loader.SetWeightExpression("EventWeight")
#     loader.AddVariable('mass_GG','F')




else:
#Add variables and spectators to factory
     loader.AddVariable('MET','F')
     loader.AddVariable('NhadJets','F')
     loader.AddVariable('NbJets','F')
     loader.AddVariable('mTPhoMET_','F')
     loader.AddVariable('dPhi_PhoMET_','F')
     loader.AddVariable('dPhi_Met_Jet1','F')
     loader.AddVariable('ST','F')
     loader.AddVariable('PhoPt','F')
     loader.AddVariable('dPhi_Met_Jet2','F')
     loader.AddVariable('dPhi_Pho_Jet2','F')
     loader.AddVariable('dPhi_Pho_Jet1','F')
     loader.AddVariable('LeadJets1Pt','F')
     loader.AddVariable('LeadJets2Pt','F')


#loader.AddVariable('HardMET','F')
#      loader.AddVariable('Pt_GG','F')
#      loader.AddVariable('ST_jets','F')

# #--Some deltaPhi variables:
#      loader.AddVariable('min_dPhi','F')
#      loader.AddVariable('dPhi_GGHardMET','F')

#Define cuts to be applied
# sigcuts = TCut('min_dPhi>-1 && min_dPhi<10 && Photons_isEB[0]==1 && Photons_isEB[1]==1 && PhotonsLoose[0].Et()>80 && PhotonsLoose[1].Et()>80 && mva_Ngoodjets>1') 
# bkgcuts = TCut('min_dPhi>-1 && min_dPhi<10 && Photons_isEB[0]==1 && Photons_isEB[1]==1 && PhotonsLoose[0].Et()>80 && PhotonsLoose[1].Et()>80 && mva_Ngoodjets>1')
#Cuts = 'PhoPt>200'
sigcuts = TCut(Cuts)#'MET>200')
bkgcuts =  TCut(Cuts)#'PhoPt>200')
#sigcuts_RandS = TCut('Pho1_hasPixelSeed==0 && Pho2_hasPixelSeed==0 && IsRandS==0 && HardMETPt>100 && mva_Ngoodjets>1 && mva_Photons0Et>80 && mva_Photons1Et>80 && abs(HardMetMinusMet)<90 && NVtx>0 && NPhotons==2 && (abs(analysisPhotons[0].Eta())<1.442 || abs(analysisPhotons[1].Eta())<1.442)')
#bkgcuts_RandS = TCut('Pho1_hasPixelSeed==0 && Pho2_hasPixelSeed==0 && IsRandS==1 && HardMETPt>100 && mva_Ngoodjets>1 && mva_Photons0Et>80 && mva_Photons1Et>80 && abs(HardMetMinusMet)<90 && NVtx>0 && NPhotons==2 && (abs(analysisPhotons[0].Eta())<1.442 || abs(analysisPhotons[1].Eta())<1.442)')

# sigcuts_RandS = TCut('Pho1_hasPixelSeed==0 && Pho2_hasPixelSeed==0 && IsRandS==0 && HardMETPt>130 && mva_Ngoodjets>1 && mva_Photons0Et>80 && mva_Photons1Et>80 && abs(HardMetMinusMet)<90 && NVtx>0 && NPhotons==2 && (abs(analysisPhotons[0].Eta())<1.44 || abs(analysisPhotons[1].Eta())<1.44) && abs(analysisPhotons[0].Eta())<2.4 && abs(analysisPhotons[1].Eta())<2.4')
# bkgcuts_RandS = TCut('Pho1_hasPixelSeed==0 && Pho2_hasPixelSeed==0 && IsRandS==1 && HardMETPt>130 && mva_Ngoodjets>1 && mva_Photons0Et>80 && mva_Photons1Et>80 && abs(HardMetMinusMet)<90 && NVtx>0 && NPhotons==2 && (abs(analysisPhotons[0].Eta())<1.44 || abs(analysisPhotons[1].Eta())<1.44) && abs(analysisPhotons[0].Eta())<2.4 && abs(analysisPhotons[1].Eta())<2.4')

#Prepare trees for training and testing (I think this just applies the cuts)
if args.mode != 'RandS':
     loader.PrepareTrainingAndTestTree(sigcuts_RandS, bkgcuts_RandS, "SplitMode=random:!V")

else:
     loader.PrepareTrainingAndTestTree(sigcuts, bkgcuts, "SplitMode=random:!V")

#Use BookMethod to specify BDT

#factory.BookMethod(loader, TMVA.Types.kBDT, 'BDT_otherthing', '!H:!V:NTrees=200:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex')

factory.BookMethod(loader, TMVA.Types.kBDT, 'BDT_{0}trees_{1}maxdepth'.format(args.Ntrees, args.maxdepth), "NTrees={0}:MaxDepth={1}:nCuts={2}".format(args.Ntrees, args.maxdepth, -1))

#factory.BookMethod(loader, TMVA.Types.kBDT, 'BDTB_{0}trees_{1}maxdepth'.format(args.Ntrees, args.maxdepth), "NTrees={0}:MaxDepth={1}:BoostType=Bagging".format(args.Ntrees, args.maxdepth))

#factory.BookMethod(loader, TMVA.Types.kBDT, 'BDTRA_{0}trees_{1}maxdepth'.format(args.Ntrees, args.maxdepth), "NTrees={0}:MaxDepth={1}:BoostType=RealAdaBoost".format(args.Ntrees, args.maxdepth))

#factory.BookMethod(loader, TMVA.Types.kBDT, 'BDTF{0}trees_{1}maxdepth'.format(args.Ntrees, args.maxdepth), "NTrees={0}:MaxDepth={1}:UseFisherCuts".format(args.Ntrees, args.maxdepth))

#factory.BookMethod(loader, TMVA.Types.kBDT, 'BDT_15D', "NTrees=200:MaxDepth=15")

#Train, test, and evaluate
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
#close output file
fout.Close()

print ('TMVA Classification is done!')

c = factory.GetROCCurve(loader)

name = "ROC_%s_%s_%s.png"%(args.sig_name1,args.back_name1,args.name)#,args.back_name3)
c.Draw()
c.SaveAs(name) #"example_%s.png",name)
#sys.stdin.readline()


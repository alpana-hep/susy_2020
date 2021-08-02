#ifndef AnalyzeLightBSM_H
#define AnalyzeLightBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"

class AnalyzeLightBSM : public NtupleVariables{

 public:
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *N2_mass="mass");
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *, const char *);
  TLorentzVector getBestPhoton();
  
  void print(Long64_t);
  //  void findObjMatchedtoG(TLorentzVector);
  double NumEvents;

  //  double   Weight;
  double wt,lumiInfb=35.9;//35.86;//36.814;//35.862824;//36.814;
  int N_0emt=0,N_all=0,N_1e=0,N_2e=0,N_1m=0,N_2m=0,N_1t=0,N_2t=0;
  int n_electrons,n_muon,n_tau;
  Int_t zeroleptons =0, Nzeroleptons =0, Nelectrons =0,Nmuons = 0,Ntaus =0;
  Int_t lept_zeroleptons =0, lept_Nzeroleptons =0, lept_Nelectrons =0,lept_Nmuons = 0,lept_Ntaus =0;
  Int_t Iso_zeroleptons =0, Iso_Nzeroleptons =0, Iso_Nelectrons =0,Iso_Nmuons = 0,Iso_Ntaus =0;
  bool isSignal=false;
  /* TH1D *h_nEvts; */
  /* TH1I *h_RunNum; */
  /* TH1D *h_intLumi; */
  /* TH1F *h_njets; */
  /* TH1F *h_bjets; */
  /* TH1F *h_MET; */
  /* TH1F *h_photonPt */;
  TH1F *h_weights;
  TH1F *h_events;
  //  TH1F *h_events;
  TH1F *h_isotrack;
  TH1F *h_zerolepton;
  TFile *oFile;
  //histogram exclusive of all the cuts;
  TH1F *h_NJets_exclusive;
  TH1F *h_NbJets_exclusive;
  TH1F *h_Met_exclusive;
  TH1F *h_PhoPt_exclusive;
  TH1F *h_PhoPhi_exclusive;

  // histogram after applying baseline selections :: lepton-veto, isotrac veto, 
  TH1F *h_NJets_LepIsoVeto;
  TH1F *h_NbJets_LepIsoVeto;
  TH1F *h_Met_LepIsoVeto;
  TH1F *h_PhoPt_LepIsoVeto;
  TH1F *h_PhoPhi_LepIsoVeto;
  TH1F *h_HT_LepIsoVeto;
  // histograms after applying MET cut,
  TH1F *h_NJets_aMetcut;
  TH1F *h_NbJets_aMetcut;
  TH1F *h_Met_aMetcut;
  TH1F *h_PhoPt_aMetcut;
  TH1F *h_PhoPhi_aMetcut;
  TH1F *h_HT_aMetcut;
  //histogram after applying pt cut
  TH1F *h_NJets_PhoPt;
  TH1F *h_NbJets_PhoPt;
  TH1F *h_Met_PhoPt;
  TH1F *h_PhoPt_PhoPt;
  TH1F *h_PhoPhi_PhoPt;
  TH1F *h_HT_PhoPt;
  // histograms after applying Njet>=2 cut,                                                                                                                                                                   
  TH1F *h_NJets_aJetcut;
  TH1F *h_NbJets_aJetcut;
  TH1F *h_Met_aJetcut;
  TH1F *h_PhoPt_aJetcut;
  TH1F *h_PhoPhi_aJetcut;
  TH1F *h_HT_aJetcut;

  // histograms after applying Njet>=2 cut,                                                                                                                                                
  TH1F *h_NJets_ExMET;
  TH1F *h_NbJets_ExMET;
  TH1F *h_Met_ExMET;
  TH1F *h_PhoPt_ExMET;
  TH1F *h_PhoPhi_ExMET;
  TH1F *h_HT_ExMET;
  //histogram inclusive of all the cuts;
  TH1F *h_NJets_inclusive;
  TH1F *h_NbJets_inclusive;
  TH1F *h_Met_inclusive;
  TH1F *h_PhoPt_inclusive;
  TH1F *h_PhoPhi_inclusive;
  TH1F *h_HT_inclusive;

};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
  int chi2_mass= atoi(N2_mass);
  //  char hname[200], htit[200];
  double xlow = 0.0,  xhigh = 2*(2350-chi2_mass), xhighh = 4.0*(2350-chi2_mass);
  //  int nbins = 2000;
  char name[100],title[100];
  
  oFile = new TFile(outFileName, "recreate");
  ////  TH1::SetDefaultSumw2(1);
  h_NJets_inclusive= new TH1F("h_NJets_inclusive","Total number of Jets (inclusive of Baseline cuts)",50,0,50);
  h_NbJets_inclusive= new TH1F("h_NbJets_inclusive","Number of B-tagged jets (inclusive of Baseline cuts)",50,0,50);
  h_Met_inclusive=new TH1F("h_Met_inclusive","Total MET (inclusive of Baseline cuts)",150,0,xhigh);
  h_PhoPt_inclusive= new TH1F("h_PhoPt_inclusive","Photon pt (inclusive of Baseline cuts)",150,0,xhigh);
  h_PhoPhi_inclusive= new TH1F("h_PhoPhi_inclusive","Photon pt (inclusive of Baseline cuts)",150,-15,15);
  h_HT_inclusive = new TH1F("h_HT_inclusive","",150,0,xhighh);
  h_NJets_exclusive= new TH1F("h_NJets_exclusive","Total number of Jets (exclusive of Baseline cuts)",50,0,50);
  h_NbJets_exclusive= new TH1F("h_NbJets_exclusive","Number of B-tagged jets (exclusive of Baseline cuts)",50,0,50);
  h_Met_exclusive= new TH1F("h_Met_exclusive","Total MET (exclusive of Baseline cuts)",150,0,xhigh);
  h_PhoPt_exclusive= new TH1F("h_PhoPt_exclusive","Photon pt (exclusive of Baseline cuts)",150,0,xhigh);
  h_PhoPhi_exclusive= new TH1F("h_PhoPhi_exclusive","Photon pt (exclusive of Baseline cuts)",150,-15,15);


  h_NJets_LepIsoVeto= new TH1F("h_NJets_LepIsoVeto","Total number of Jets after LepIsoVeto cut",50,0,50);
  h_NbJets_LepIsoVeto= new TH1F("h_NbJets_LepIsoVeto","Number of B-tagged jets after LepIsoVeto cut ",50,0,50);
  h_Met_LepIsoVeto= new TH1F("h_Met_LepIsoVeto","Total MET after LepIsoVeto cut",150,0,xhigh);
  h_PhoPt_LepIsoVeto= new TH1F("h_PhoPt_LepIsoVeto","Photon pt after LepIsoVeto cut",150,0,xhigh);
  h_PhoPhi_LepIsoVeto= new TH1F("h_PhoPhi_LepIsoVeto","Photon pt after LepIsoVeto cut",150,-15,15);
  h_HT_LepIsoVeto = new TH1F("h_HT_LepIsoVeto","",150,0,xhighh);

  h_NJets_aMetcut= new TH1F("h_NJets_aMetcut","Total number of Jets after Metcut",50,0,50);
  h_NbJets_aMetcut= new TH1F("h_NbJets_aMetcut","Number of B-tagged jets after Metcut ",50,0,50);
  h_Met_aMetcut= new TH1F("h_Met_aMetcut","Total MET after Metcut ",150,0,xhigh);
  h_PhoPt_aMetcut= new TH1F("h_PhoPt_aMetcut","Photon pt after Metcut ",150,0,xhigh);
  h_PhoPhi_aMetcut= new TH1F("h_PhoPhi_aMetcut","Photon pt after Metcut ",150,-15,15);
  h_HT_aMetcut = new TH1F("h_HT_aMetcut","",150,0,xhighh);
  h_NJets_PhoPt= new TH1F("h_NJets_PhoPt","Total number of Jets after pho pt cut",50,0,50);
  h_NbJets_PhoPt= new TH1F("h_NbJets_PhoPt","Number of B-tagged jets after pho pt cut ",50,0,50);
  h_Met_PhoPt= new TH1F("h_Met_PhoPt","Total MET after pho pt cut ",150,0,xhigh);
  h_PhoPt_PhoPt= new TH1F("h_PhoPt_PhoPt","Photon pt after pho pt cut ",150,0,xhigh);
  h_PhoPhi_PhoPt= new TH1F("h_PhoPhi_PhoPt","Photon pt after pho pt cut ",150,-15,15);
  h_HT_PhoPt = new TH1F("h_HT_PhoPt","",150,0,xhighh);

  h_NJets_aJetcut= new TH1F("h_NJets_aJetcut","Total number of Jets after Jet cut",50,0,50);
  h_NbJets_aJetcut= new TH1F("h_NbJets_aJetcut","Number of B-tagged jets after Jet cut ",50,0,50);
  h_Met_aJetcut= new TH1F("h_Met_aJetcut","Total MET after Jet cut",150,0,xhigh);
  h_PhoPt_aJetcut= new TH1F("h_PhoPt_aJetcut","Photon pt after Jet cut",150,0,xhigh);
  h_PhoPhi_aJetcut= new TH1F("h_PhoPhi_aJetcut","Photon pt after Jet cut",150,-15,15);
  h_HT_aJetcut = new TH1F("h_HT_aJetcut","",150,0,xhighh);
  h_NJets_ExMET= new TH1F("h_NJets_ExMET","Total number of Jets after Jet cut",50,0,50);
  h_NbJets_ExMET= new TH1F("h_NbJets_ExMET","Number of B-tagged jets after Jet cut ",50,0,50);
  h_Met_ExMET= new TH1F("h_Met_ExMET","Total MET after MET cut",150,0,xhigh);
  h_PhoPt_ExMET= new TH1F("h_PhoPt_ExMET","Photon pt after Jet cut",150,0,xhigh);
  h_PhoPhi_ExMET= new TH1F("h_PhoPhi_ExMET","Photon pt after Jet cut",150,-15,15);
  h_HT_ExMET = new TH1F("h_HT_ExMET","",150,0,xhighh);
  /* h_nEvts=new TH1D("nEvents","no. of events in this tree",4,0,4); */
  /* h_RunNum=new TH1I("runs","Run nos.",300000,0,300000); */

  /* h_intLumi=new TH1D("intLumi","integrated luminosity in /fb",10000,25,200);  */
  /* //declare histogram here                                                                                                                      */
  /* h_njets = new TH1F("h_njets","Njets distribution",30,0,30); */
  /* h_njets->GetXaxis()->SetTitle("number of jets in an event"); */
  /* h_njets->SetLineColor(kBlue); */
  /* h_photonPt = new TH1F("h_photon_Pt","pt distribution of photon at gen level",100,0,2400); */
  /* h_photonPt->GetXaxis()->SetTitle("pt of photon"); */
  /* h_photonPt->SetLineColor(kBlue); */

  /* h_MET = new TH1F("h_MET","Met distribution",100,0,2400); */
  /* h_MET->GetXaxis()->SetTitle("Total MET"); */
  /* h_MET->SetLineColor(kBlue); */

  /* h_bjets = new TH1F("h_bjets","NbJets distribution",30,0,30); */
  /* h_bjets->GetXaxis()->SetTitle("number of b-tagged jets in an event"); */
  /* h_bjets->SetLineColor(kBlue); */

  h_events = new TH1F("h_events","binwise filling of events",10,0,10);
  
  const char *ee[10] = {"h-h","e-h","t-h","m-h","e-e","t-t","m-m","e-m","e-t","t-m"};
  for(Int_t j =1; j<11;j++)
    {h_events->GetXaxis()->SetBinLabel(j,ee[j-1]);
    }
  h_isotrack = new TH1F("h_isotrack","binwise filling of events",10,0,10);

  //  const char *ee[10] = {"h-h","e-h","t-h","m-h","e-e","t-t","m-m","e-m","e-t","t-m"};
  for(Int_t j =1; j<11;j++)
    {h_isotrack->GetXaxis()->SetBinLabel(j,ee[j-1]);
    }
  h_zerolepton = new TH1F("h_zerolepton","binwise filling of events",10,0,10);

  //  const char *ee[10] = {"h-h","e-h","t-h","m-h","e-e","t-t","m-m","e-m","e-t","t-m"};
  for(Int_t j =1; j<11;j++)
    {h_zerolepton->GetXaxis()->SetBinLabel(j,ee[j-1]);
    }

}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* dataset, const char* N2_mass) {
  string nameData=dataset;//vvv

  TChain *tree = new TChain("PreSelection");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName, N2_mass);
  
}

Bool_t AnalyzeLightBSM::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeLightBSM::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

AnalyzeLightBSM::~AnalyzeLightBSM() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif


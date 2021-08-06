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
  //Long64_t transMass(float , float, float, float);
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
  //
  TH1F *h_dR_Electron;
  TH1F *h_dR_Ele_zoomed;
  TH1F *h_dR_Ele_lastBin;

  TH1F *h_dEta_Electron;
  TH1F *h_dPhi_Electron;
  TH1F *h_dR_DiffMethod;
  TH1F *h_mindR_Electron;
  TH2F *h_dPhiVsdEta_Electron;
  TH2F *h_dEtaVsdR;
  TH2F *h_dPhiVsdR;
  TH2F *h_dEtaVsdR_zoomed;
  TH2F *h_dPhiVsdR_zoomed;
  /* TH2F *h_dEtaVsdR_last; */
  /* TH2F *h_dPhiVsdR; */

  // transverse mass distribution
  TH1F *h_Mt;
  TH2F *h_MtvsEta;
  TH1F *h_Mt_zoomed;
  TH1F *h_weights;
  TH1F *h_events;
  //  TH1F *h_events;
  TH1F *h_isotrack;
  TH1F *h_zerolepton;
  TFile *oFile;
  // distribution for W-jets
  TH1F *h_GenW_Pt;
  TH1F *h_GenW_P;
  TH1F *h_GenW_Eta;
  TH1F *h_GenW_Phi;
  TH2F *h_GenW_PtvsP;
  TH2F *h_GenW_PtvsEta;
  TH2F *h_GenW_PvsEta;
  TH2F *h_GenW_PtvsPhi;
  //histogram exclusive of all the cuts;
  TH1F *h_GenEle_Pt;
  TH1F *h_GenEle_Pt_zoomed;
  TH1F *h_GenEle_Pt_lastbin;

  TH1F *h_GenEle_Mt;
  TH1F *h_GenEle_Mt_zoomed;
  TH1F *h_GenEle_Mt_lastbin;
  
  TH1F *h_GenEle_Phi;
  TH1F *h_GenEle_Eta;
  TH1F *h_RecEle_Pt;
  TH1F *h_GenEle_Pt_reBin;
  TH1F *h_RecEle_Pt_zoomed;
  TH1F *h_RecEle_Pt_lastbin;
  TH1F *h_RecEle_Phi;
  TH1F *h_RecEle_Eta;
  //2d histogram;
  TH2F *h_NGenVsReco_Elec;
  TH2F *h_PtVsEta_GenElec;
  TH2F *h_PtVsEta_RecElec;
  TH2F *h_PtVsEta_RecElec_zoomed;
  TH2F *h_PtVsEta_GenElec_zoomed;

  TH2F *h_PtVsPhi_GenElec;
  TH2F *h_PtVsPhi_RecElec;
  TH2F *h_EtavsPhi_RecElec;
  TH2F *h_EtavsPhi_GenElec;
  TH2F *h_temp;
  TH2F *h_Pt_eVsnu;
  TH2F *h_Eta_eVsnu;
  TH2F *h_Phi_eVsnu;
  //Gen-reco electron matching plots
  TH1F *h_GenRecEle_Pt;
  TH1F *h_GenRecEle_Pt_rebin;
  TH1F *h_GenRecEle_Pt_reBin;

  TH1F *h_GenRecEle_Eta;
  //Gen e- matching with photon
  TH1F *h_GenElec_RecPho_Pt;
  TH1F *h_GenElec_RecPho_Pt_rebin;
  TH1F *h_GenElec_RecPho_Pt_reBin;
  TH1F *h_GenElec_RecPho_Eta;
  TH1F *h_dR_GenEleVsPhoton;
  TH1F *h_mindR_GenEleVsPhoton;
  TH2F *h_Pt_WVselec;
  TH2F *h_Pt_WVsMet;
  TH2F *h_ElecPtVs_MET;
  TH2F *h_Pt_WVsGen_nu;
  TH2F *h_Pt_WVs_GenMet;
  TH2F *h_Pt_W_vsMt;
  TH1F *h_nu_pt;
  
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
  int chi2_mass= atoi(N2_mass);
  //  char hname[200], htit[200];
  double xlow = 0.0,  xhigh = 3200, xhigh1 = 3500,xhigh2=300;//4.0*(2350-chi2_mass);
  //  int nbins = 2000;
  char name[100],title[100];
  
  oFile = new TFile(outFileName, "recreate");
  //  TH1::SetDefaultSumw2(1);
  h_Pt_W_vsMt = new TH2F("h_Pt_W_vsMt","W-pt vs Mt",700,0,xhigh1,200,0,1000);
  h_Pt_W_vsMt->SetMinimum(0);
  h_nu_pt= new TH1F("h_nu_pt","nuetrino Pt",700,0,xhigh1);
  h_Pt_WVs_GenMet= new TH2F("h_Pt_WVs_GenMet","W pt vs GenMET",700,0,xhigh1, 700,0,xhigh1);
  h_Pt_WVs_GenMet->SetMinimum(0);
  h_Pt_WVsGen_nu=new TH2F("h_Pt_WVsGen_nu"," w-pt vs nuetrino pt",700,0,xhigh1, 700,0,xhigh1);
  h_Pt_WVsGen_nu->SetMinimum(0);
  h_Pt_WVselec= new TH2F("h_Pt_WVselec","Pt for W vs that of e-",700,0,xhigh1,700,0,xhigh1);
  h_Pt_WVselec->SetMinimum(0);
  h_Pt_WVsMet = new TH2F("h_Pt_WVsMet","W-Pt vs MET",700,0,xhigh1, 700,0,xhigh1);
  h_Pt_WVsMet->SetMinimum(0); 
  h_ElecPtVs_MET = new TH2F("h_ElecPtVs_MET","Elec Pt vs MET",700,0,xhigh1, 700,0,xhigh1);
  h_ElecPtVs_MET->SetMinimum(0);
  h_mindR_GenEleVsPhoton=new TH1F("h_mindR_GenEleVsPhoton","min dR between gen e- & photon",200,0,8);
  h_dR_GenEleVsPhoton= new TH1F("h_dR_GenEleVsPhoton","dR between gen e- & photon",200,0,8);
  
  h_GenRecEle_Pt= new TH1F("h_GenRecEle_Pt","Pt of gen e- matching with reco e-",700,0,xhigh1);
  h_GenElec_RecPho_Pt= new TH1F("h_GenElec_RecPho_Pt","Pt of gen e- faking as photon", 700,0, xhigh1);
  h_GenRecEle_Pt_rebin= new TH1F("h_GenRecEle_Pt_rebin","Pt of gen e- matching with reco e-",320,0,xhigh);
  h_GenRecEle_Pt_reBin= new TH1F("h_GenRecEle_Pt_reBin","Pt of gen e- matching with reco e-",110,0,xhigh);

  h_GenElec_RecPho_Pt_rebin= new TH1F("h_GenElec_RecPho_Pt_rebin","Pt of gen e- faking as photon",320,0,xhigh);
  h_GenElec_RecPho_Pt_reBin= new TH1F("h_GenElec_RecPho_Pt_reBin","Pt of gen e- faking as photon",110,0,xhigh);
  h_GenRecEle_Eta= new TH1F("h_GenRecEle_Eta","Eta of gen e- matching with reco e-",150,-8,8);
  h_GenElec_RecPho_Eta = new TH1F("h_GenElec_RecPho_Eta","Eta of gen e- faking as photon",150,-8,8);

  h_Pt_eVsnu = new TH2F ("h_Pt_eVsnu","Tranverse mom of E- vs nu",300,0,1500,300,0,1500);
  h_Pt_eVsnu->SetMinimum(0);
  h_Eta_eVsnu = new TH2F("h_Eta_eVsnu","Eta for e- vs that for nu",150,-8,8,150,-8,8);
  h_Eta_eVsnu->SetMinimum(0);
  h_Phi_eVsnu = new TH2F("h_Phi_eVsnu","Phi for e- vs that for nu",150,-6,6,150,-6,6);
  h_Phi_eVsnu->SetMinimum(0);
  h_MtvsEta= new TH2F("h_MtvsEta","transverse mass vs eta",150,-8,8,300,0,1500);
  h_MtvsEta->SetMinimum(0);
  h_Mt = new TH1F("h_Mt","Transverse mass",300,0,1500);
  h_Mt_zoomed=new TH1F("h_Mt_zoomed","Transverse mass",250,0,500);
  h_GenW_Pt=new TH1F("h_GenW_Pt","Pt distribution for W-bosons at gen level",700,0,xhigh1);
  h_GenW_P=new TH1F("h_GenW_P","Momentum distribution for W-bosons at gen level",700,0,xhigh1);
  h_GenW_Eta=new TH1F("h_GenW_Eta","Eta distribution for W-bosons at gen level",150,-8,8);
  h_GenW_Phi=new TH1F("h_GenW_Phi","Phi distribution for W-bosons at gen level",150,-6,6);
  h_GenW_PtvsP= new TH2F("h_GenW_PtvsP","Pt vs total momentum for W at gen level",700,0,xhigh1,700,0,xhigh1);
  h_GenW_PtvsP->SetMinimum(0);
  h_GenW_PtvsEta= new TH2F("h_GenW_PtvsEta","Pt vs ETa  for W at gen level",150,-8,8,700,0,xhigh1);
  h_GenW_PtvsEta->SetMinimum(0);
  h_GenW_PvsEta= new TH2F("h_GenW_PvsEta","total momentum vs ETa  for W at gen level",150,-8,8,700,0,xhigh1);
  h_GenW_PvsEta->SetMinimum(0);
  h_GenW_PtvsPhi= new TH2F("h_GenW_PtvsPhi","pt vs phi for w bosons",150,-8,8,700,0,xhigh1);
  h_GenW_PtvsPhi->SetMinimum(0);
  h_GenEle_Pt= new TH1F("h_GenEle_Pt","Pt for electrons at Gen Level",320,0,xhigh);
  h_GenEle_Pt_reBin= new TH1F("h_GenEle_Pt_reBin","Pt for electrons at Gen Level",110,0,xhigh);
  h_GenEle_Pt_zoomed= new TH1F("h_GenEle_Pt_zoomed","Pt for electrons at Gen Level",60,0,xhigh2);
  h_GenEle_Pt_lastbin= new TH1F("h_GenEle_Pt_lastbin","Pt for electrons at Gen Level",700,0,xhigh1);
  h_GenEle_Phi= new TH1F("h_GenEle_Phi","Phi distribution for eletcron at Gen level",150,-6,6);
  h_GenEle_Eta=new TH1F("h_GenEle_Eta","Eta distribution for eletcron at Gen level",150,-8,8);

  h_RecEle_Pt= new TH1F("h_RecEle_Pt","Pt for electrons at Gen Level",200,0,xhigh);
  h_RecEle_Pt_zoomed= new TH1F("h_RecEle_Pt_zoomed","Pt for electrons at Rec Level",60,0,xhigh2);
  h_RecEle_Pt_lastbin= new TH1F("h_RecEle_Pt_lastbin","Pt for electrons at Rec Level",700,0,xhigh1);

  h_RecEle_Phi= new TH1F("h_RecEle_Phi","phi distribution for eletcron at Gen level",150,-6,6);
  h_RecEle_Eta=new TH1F("h_RecEle_Eta","Eta distribution for eletcron at Gen level",100,-8,8);

  h_NGenVsReco_Elec= new TH2F("h_NGenVsReco_Elec"," Gen electrons vs Reco electrons",8,0,8,4,0,4);
  h_PtVsEta_GenElec = new TH2F("h_PtVsEta_GenElec"," Pt vs Eta for gen electrons",100,-8,8,200,0,xhigh);
  h_PtVsEta_GenElec->SetMinimum(0);
  h_PtVsEta_GenElec_zoomed = new TH2F("h_PtVsEta_GenElec_zoomed"," Pt vs Eta for Gen electrons",100,-8,8,700,0,xhigh1);
  h_PtVsEta_GenElec_zoomed->SetMinimum(0);
  h_PtVsEta_RecElec_zoomed = new TH2F("h_PtVsEta_RecElec_zoomed"," Pt vs Eta for Rec electrons",100,-8,8,700,0,xhigh1);
  h_PtVsEta_RecElec_zoomed->SetMinimum(0);

  h_PtVsPhi_GenElec = new TH2F("h_PtVsPhi_GenElec"," Pt vs phi for gen electrons",100,-6,6,200,0,xhigh);
  h_PtVsPhi_GenElec->SetMinimum(0);

  h_EtavsPhi_GenElec = new TH2F("h_EtaVsPhi_GenElec","Eta vs Phi for gen electron",200,-8,8,200,-6,6);
  h_EtavsPhi_GenElec->SetMinimum(0);
  
  h_PtVsEta_RecElec = new TH2F("h_PtVsEta_RecElec"," Pt vs Eta for rec electrons",100,-8,8,200,0,xhigh);
  h_PtVsEta_RecElec->SetMinimum(0);
  h_PtVsPhi_RecElec = new TH2F("h_PtVsPhi_RecElec"," Pt vs phi for rec electrons",150,-6,6,200,0,xhigh);
  h_PtVsPhi_RecElec->SetMinimum(0);
  h_EtavsPhi_RecElec = new TH2F("h_EtaVsPhi_RecElec","Eta vs Phi for rec electron",200,-8,8,200,-6,6);
  h_EtavsPhi_RecElec->SetMinimum(0);
  h_temp = new TH2F("h_temp","tem",100,-8,8,200,0,xhigh);
  h_dPhiVsdEta_Electron= new TH2F("h_dPhiVsdEta_Electron","d-phi vs d-eta",100,-8,8,100,-8,8);
  h_dPhiVsdEta_Electron->SetMinimum(0);
  h_dEtaVsdR= new TH2F("h_dEtaVsdR","d-eta vs dR :d_eta:dR",100,-8,8,200,0,8);
  h_dEtaVsdR->SetMinimum(0);
  h_dPhiVsdR= new TH2F("h_dPhiVsdR","d-phi vs dR :d_phi:dR",150,-8,8,200,0,8);
  h_dPhiVsdR->SetMinimum(0);
  h_dEtaVsdR_zoomed= new TH2F("h_dEtaVsdR_zoomed","d-eta vs dR :d_eta:dR",100,-8,8,200,0,2);
  h_dEtaVsdR_zoomed->SetMinimum(0);

  h_dPhiVsdR_zoomed= new TH2F("h_dPhiVsdR_zoomed","d-phi vs dR :d_phi:dR",150,-8,8,200,0,2);
  h_dPhiVsdR_zoomed->SetMinimum(0);
  h_dR_Electron = new TH1F("h_dR_Electron","dR (Gen e- && Rec e-)",200,0,8);
  h_dR_Ele_zoomed = new TH1F("h_dR_Ele_zoomed","dR (Gen e- && Rec e-)",200,0,2);
  h_dR_Ele_lastBin = new TH1F("h_dR_Ele_lastBin","dR (Gen e- && Rec e-)",200,0,2);

  h_dEta_Electron=new TH1F("h_dEta_Electron" ,"dEta(Gen e- & Rec e-)",100,-8,8);
  h_dPhi_Electron=new TH1F("h_dPhi_Electron" ,"dPhi(Gen e- & Rec e-)",150,-8,8);
  h_dR_DiffMethod= new TH1F("h_dR_DiffMethod","dR (Gen e- && Rec e-)",200,0,8);  
  h_mindR_Electron= new TH1F("h_mindR_Electron","min_dR (Gen e- && Rec e-)",200,0,8);

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


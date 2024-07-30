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

#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;

class AnalyzeLightBSM : public NtupleVariables{

 public:
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *sample="sample");
  //  std::cout<<"alpana"<<std::endl;
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *,const char *);
  void     BookHistogram(const char *, const char *);
  int getBinNoV4(int);
  int getBinNoV7(int,int);
  int getBinNoV6(int);
  int getBinNoV6_WithOnlyBLSelec(int,int);
  int getBinNoV6_WithOnlyBLSelec_v1(int,int);
  int getBinNoV6_WithOnlyBLSelec_v2(int,int);
  int getBinNo_v1FR(double , int );
  int getBinNo_v0FR(double , double, double );
  int getBinNo_v2FR(double, int,int);
  TLorentzVector getBestPhoton(int);
  TLorentzVector getPhoton_withoutfullID();
  vector <TLorentzVector> getLorentzVector(int, Float_t[],Float_t[],Float_t[],Float_t[]);
  void FillHistogram_Kinematics(int ,int , int, double , double , double, double , double, int , double , double , double );//int ,TLorentzVector,int, int, double , double, double, double,float,float,double,float,double,double,double,double,double,TLorentzVector,vector<TLorentzVector>,TLorentzVector, int, double, double);
  void FillHistogram_Kinematics_varBin(int ,double, int , int , double , double,double, double);//float,float,double,float,double,double,double,double );
  void FillTFBins_Valid(int, TLorentzVector,int, int, double , double, double, double,float,float,double,float,double,double,double,double,double,TLorentzVector,vector<TLorentzVector>,TLorentzVector, int, double,double,double, double);
  
  int Photons_OriginType();
  //  <vector>
  double getGendRLepPho(int);
  double getGendRElecPho(int);
  double getdR_GenPho_RecoPho(TLorentzVector);
  void CrossSection_Map_Init();
  double getGenLep1(TLorentzVector, bool);
  double getGenLep(TLorentzVector);
  double getGenRecodRLep(TLorentzVector);
  int getBinNoV7_le(int , int);
  int getBinNoV1_le(int , int);
  int getBinNoV2_st(int,int, int, int);
  int   getBinNoV7_highMET(int, int);
  int getBinNoV16_le(int, int, double);
  std::vector<int>dR_recoPho_GenParticle(TLorentzVector);
    
  //Long64_t transMass(float , float, float, float);
  void print(Long64_t);
  //  void findObjMatchedtoG(TLorentzVector);
  double NumEvents;
  //  double   Weight;
  double wt,lumiInfb=35.9;//35.86;//36.814;//35.862824;//36.814;
  int bestPhotonIndxAmongPhotons=-100;
  int N_0emt=0,N_all=0,N_1e=0,N_2e=0,N_1m=0,N_2m=0,N_1t=0,N_2t=0;
  int n_electrons,n_muon,n_tau;
  Int_t zeroleptons =0, Nzeroleptons =0, Nelectrons =0,Nmuons = 0,Ntaus =0;
  Int_t lept_zeroleptons =0, lept_Nzeroleptons =0, lept_Nelectrons =0,lept_Nmuons = 0,lept_Ntaus =0;
  Int_t Iso_zeroleptons =0, Iso_Nzeroleptons =0, Iso_Nelectrons =0,Iso_Nmuons = 0,Iso_Ntaus =0;

  vector<TLorentzVector> GenElectrons_v1;
  vector<TLorentzVector> GenParticles_v1;
  vector<TLorentzVector> GenMuons_v1;
  vector<TLorentzVector> GenTaus_v1;
  vector<TLorentzVector> GenJets_v1;
  vector<TLorentzVector> Electrons_v1;
  vector<TLorentzVector> Photons_v1;
  vector<TLorentzVector> Muons_v1;
  vector<TLorentzVector> Taus_v1;
  vector<TLorentzVector> Jets_v1;
  vector<TLorentzVector>HLTElectronObjects_v1;
  vector<TLorentzVector> TAPElectronTracks_v1;
  int HLTElectronObjects_ =0;
  int  TAPElectronTracks_=0;
  int Jets_ =0;
  int Taus_=0;
  int Muons_=0;
  int Photons_=0;
  int Electrons_=0;
  int GenJets_=0;
  int GenTaus_ =0;
  int  GenMuons_=0;
  int GenParticles_=0;
  int GenElectrons_=0;

  int BTags;
  bool isSignal=false;
  /* vector<double> METLowEdge={200,270,350,450,750,2000}; */
  /* vector<double> METLowEdge1={100,200,270,350,450,750,2000}; */
  /* vector<double> METLowEdge2={100,200,270,350,450,2000}; */
  /* vector<double> METLowEdge_v3={200,300,370,450,600,750,900,2000}; */
  /* vector<double> METLowEdge_v3_1={200,300,370,450,600,900,2000}; */
  vector<double> METLowEdge_lowMET={100,370,450,600};
  vector<double> METLowEdge_highMET={300,370,450,600};


  vector<double> METLowEdge_v2={100,200,300,370,450,600,750,900};//{100,200,,270,350,450,600,750,900,2000};
  vector<double> METLowEdge_v2_1={100,200,300,370,450,600,750};
  vector<double> METLowEdge_v2_2={100,200,300,370,450,600};
  vector<double> METLowEdge_v3={200,300,370,450,600,750,900};
  vector<double> METLowEdge_v3_1={200,300,370,450,600,750};
  vector<double> METLowEdge_v3_2={200,300,370,450,600};
  vector<double> METLowEdge_v1={300,370,450,600,750,900};
  vector<double> METLowEdge_v1_1={300,370,450,600,750};
  vector<double> METLowEdge_v1_2={300,370,450,600};

  vector<double> BestPhotonPtBinLowEdge={40,70,100,120,140,160,200,240,300,450,600,1000};
  vector<double> QMultLowedge={0,2,4,7,100};
  vector<double>  nJetsLowedge={2,5,10,20};
  vector<double>  nbtagsLowedge={0,1,10};
  //define histograms here
  TH1F *h_selectBaselineYields_;
  TH1F *h_selectBaselineYields_v2;
   TH1F *h_selectBaselineYields_v1;

  TFile *oFile;
  /* TH1F *h_events; */
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
  int chi2_mass= atoi(N2_mass);
  //  char hname[200], htit[200];
  double xlow = 0.0,  xhigh = 3200, xhigh1 = 3500,xhigh2=300;//4.0*(2350-chi2_mass);
  //  int nbins = 2000;
  char name[100],title[100];
  char hname[1000],hname1[1000], hname1_2d[1000],hname_2d[10000],hname_njets[10000],hname_nBjets[10000], hname_Met[10000],hname_PhoPt[10000],hname_Mt_phopt[10000],hname_dPhi[10000],hname_st[1000],hname_ht[1000],hname_njet_vs_ST[1000],hname_njet_vs_HT[1000],hname_ST_vs_ptPho[1000];
  vector<string> baseline = {"No_selection","PreSelection","MET_200"};//,"lepton_veto","veto_chargTracks","photon_pT40","MET_100","nJets_2","ST_300","Trig_eff","Even  cout<<"size of baseline vector"<<"\t"<<baseline.size()<<endl;
  
  vector<string> checks={"NoSelection","PreSelection","Elec_CR","Pho_SR"};//,"HEM_veto_Elec_CR","HEM_veto_Pho_SR","L1Trig_Elec_CR","L1Trig_Pho_SR","ProbL1Trig_Elec_C
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //Initialize histogram here
  h_selectBaselineYields_ = new TH1F("cutflows","cutflows",60,-0.5,60.5);
  h_selectBaselineYields_v2 = new TH1F("cutflows_BL","cutflows_LL",60,-0.5,60.5);
  h_selectBaselineYields_v1 = new TH1F("cutflows_LL","cutflows_LL",60,-0.5,60.5);

}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* dataset, const char* N2_mass) {
  string nameData=dataset;//vvv
  //TDirectory * dir = new TDirectory("TreeMaker2");
  TChain *tree = new TChain("Pre_Selection");
    //  TChain *tree = new TChain("PreSelection");
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
  CrossSection_Map_Init();

  //Jets = 0;
}
void AnalyzeLightBSM::CrossSection_Map_Init()
{
  char *f_name_EH = new char[2000];
  sprintf(f_name_EH,"./map_crosssection_SMprocess_v1.txt");//,chi2_method);
  std::ifstream in_EH(f_name_EH);
  if(!in_EH) {
    cout<<"ERROR => "<<f_name_EH<<" Not found"<<endl;
    //return;                                                                                                                                
    exit(0);
  }
  string process_name;
  float value, entries;
  cout<<"File name = "<<f_name_EH<<endl;
  while(in_EH>>process_name>>value>>entries){
    std::pair<std::string, float> temp_pair;
    /* std::vector<float> temp_vector; */
    /* temp_vector.push_back(w1); */
    /* temp_vector.push_back(w2); */
    /* temp_vector.push_back(w3); */
    float weight =value/entries;

    temp_pair = std::make_pair(process_name,weight);
    cross_sectionValues.insert(temp_pair);
  }
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


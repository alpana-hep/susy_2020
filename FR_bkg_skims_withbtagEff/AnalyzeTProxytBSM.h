#ifndef ANALYZETPROXYTBSM_H
#define ANALYZETPROXYTBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVarsTProxy.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"

//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;
//#pragma link C++ class NtupleVarsTProxy+;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myLV;

class AnalyzeTProxytBSM : public NtupleVarsTProxy{

 public:
  AnalyzeTProxytBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *sample="sample", const char* LostlepFlag ="Flag", const char* phoID="phoID");
  //  std::cout<<"alpana"<<std::endl;
  ~AnalyzeTProxytBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  //void     EventLoop(const char *,const char *,const char *,const char *, const char*, const char*);
  //  void     EventLoop();
  void     BookHistogram(const char *);//, const char *);
  Bool_t  Process(Long64_t entry);
myLV getBestPhoton(int);
 void storeBTagEff(const char *);
 int bestPhotonIndxAmongPhotons=-100;
 void     EventLoop(const char *,const char *);

 TH1F *h_selectBaselineYields_;
 vector<float> xbins = {-1,0, 30, 50, 70, 100, 140, 200, 300, 600, 1000,99999};//20,30,40,50,60,70,80,100,120,160,210,260,320,400,500,600,800,99999}, 
 vector<float> ybins = {0.0,0.8,1.6,2.4,3.0};
 TH2F *d_eff_b, *n_eff_b;
 TH2F *d_eff_c, *n_eff_c;
 TH2F *d_eff_udsg, *n_eff_udsg;
 
 TFile *oFile;
 /* TTree* outtree; */
};
#endif


#ifdef ANALYZETPROXYTBSM_cxx

//void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
void AnalyzeTProxytBSM::BookHistogram(const char *outFileName) {
  std::cout << "AnalyzeLightBSM::BookHistogram " << std::endl;
  TH1::SetDefaultSumw2(1);

  oFile = new TFile(outFileName, "recreate");
  //outtree = new TTree("PreSelection");
  //outtree = new TTree();
  std::cout << "AnalyzeTProxytBSM::BookHistogram " << fChain->GetEntries() << std::endl;
  char name[100],title[100];
  h_selectBaselineYields_ = new TH1F("cutflows","cutflows",10,-0.5,9.5);
  d_eff_b = new TH2F("d_eff_b_","d_eff_b",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));
  n_eff_b = new TH2F("n_eff_b","n_eff_b",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));
  d_eff_c = new TH2F("d_eff_c","d_eff_c",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));
  n_eff_c = new TH2F("n_eff_c","n_eff_c",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));
  d_eff_udsg = new TH2F("d_eff_udsg","d_eff_udsg",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));
  n_eff_udsg = new TH2F("n_eff_udsg","n_eff_udsg",xbins.size()-1,&(xbins[0]), ybins.size()-1,&(ybins[0]));

  //  outtree = fChain->CloneTree(0);
  //outtree->Print();
}

AnalyzeTProxytBSM::AnalyzeTProxytBSM(const TString &inputFileList, const char *outFileName,const char *dataset, const char *sample, const char* LostlepFlag, const char* phoID) {
  
  std::cout << outFileName << std::endl;

  string nameData=dataset;//vvv

  //TDirectory * dir = new TDirectory("TreeMaker2");
  TChain *tree = new TChain("TreeMaker2/PreSelection");

  // for skimmed tree
  //TChain *tree = new TChain("PreSelection"); 

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }
  
  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  //NtupleVarsTProxy::Init(tree,nameData);
  NtupleVarsTProxy::Init(tree);
  
  BookHistogram(outFileName); //, N2_mass);
  //CrossSection_Map_Init();

  //Jets = 0;

}

Bool_t AnalyzeTProxytBSM::FillChain(TChain *chain, const TString &inputFileList) {

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

Long64_t AnalyzeTProxytBSM::LoadTree(Long64_t entry) {
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


AnalyzeTProxytBSM::~AnalyzeTProxytBSM() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  /* std::cout << "AnalyzeTProxytBSM::~AnalyzeTProxytBSM " << outtree->GetEntries() << std::endl;  */
  oFile->Write();
  oFile->Close();

}

#endif // AnalyzeTProxytBSM_cxx

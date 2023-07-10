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
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *sample="sample", const char* LostlepFlag ="Flag", const char* phoID="phoID");
  //  std::cout<<"alpana"<<std::endl;
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *,const char *, const char*, const char*);
  void     BookHistogram(const char *, const char *);
  int getBinNoV4(int);
  int getBinNoV7(int,int);
  int getBinNoV6(int);
  int getBinNoV6_WithOnlyBLSelec(int,int);
  TLorentzVector getBestPhoton(int);
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
  vector<TLorentzVector> GenJetsAK8_v1;
  vector<TLorentzVector> GenJetsAK15_v1;
  vector<TLorentzVector> JetsAK8_v1;
  vector<TLorentzVector> JetsAK15_v1;
  vector<TLorentzVector> JetsAK15_subjets_v1;
  vector<TLorentzVector> JetsAK8_subjets_v1;
  vector<TLorentzVector>JetsConstituents_v1;
  vector<TVector3> PrimaryVertices_v1;
  vector<TVector3> Tracks_v1;
  vector<TVector3> Tracks_referencePoint_v1;
  vector<TVector3> GenVertices_v1;
  vector<TLorentzVector>TAPElectronTracks_v1;
  vector<TLorentzVector>TAPMuonTracks_v1;
  vector<TLorentzVector>TAPPionTracks_v1;

  

  int BTags;
  bool isSignal=false;
  vector<double> METLowEdge={200,270,350,450,750,2000};
  vector<double> METLowEdge1={100,200,270,350,450,750,2000};
  vector<double> METLowEdge2={100,200,270,350,450,2000};
  vector<double> METLowEdge_v3={100,300,370,450,600,750,900};
  vector<double> METLowEdge_v3_1={100,300,370,450,600,750};
  TH1F *h_selectBaselineYields_;
  TH1F *h_selectBaselineYields_v1;
  /* TH1F *h_selectBaselineYields_SR; */
  /* TH1F *h_selectBaselineYields_CR;  */
  /* TH1F *h_selectBaselineYields_MuSR; */
  /* TH1F *h_selectBaselineYields_MuCR; */
  /* TH2F *h_2d_genPromptPho_VsGenPho; */
  /* TH1F *h_isotrack; */
  /* TH1F *h_zerolepton; */
  TFile *oFile;
  /* TH1F *h_events; */
  /* TH1D *h_nEvts; */

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
  //  const char *baseline[25]= {"Nocut","skim","bkg_comp","Met-filter","veto-lep","veto-iso","dPhi_Met","pt_jet","ST_300","Met_250","pT_100","nocut","HT_1TeV_Met100","HT_1TeV_Met250","HT_1TeV_Met250_pt_100","HT_15TeV_Met100","HT_15TeV_Met250","HT_15TeV_Met250_pt_100","HT_175TeV_Met100","HT_175TeV_Met250","HT_175TeV_Met250_pt_100","HT_2TeV_Met100","HT_2TeV_Met250","HT_2TeV_Met250_pt_100","nocut_sam"};//"st_300_Met100","pt_st_Met_250","st_300_Met250","nocut"};
  //Double_t xbins_PhotPt[110]={};//{20,25,30,35,40,,7,10,20,30,40,50,80,90,100,150};
  //const char *baseline[24]= {"Nocut","skim","bkg_comp","Met-filter","veto-lep","veto-iso","dPhi_Met","pt_jet","ST_300","Met_250","pT_100","nocut","HT_1TeV_Met100","HT_1TeV_Met250","HT_1TeV_Met250_pt_100","HT_15TeV_Met100","HT_15TeV_Met250","HT_15TeV_Met250_pt_100","HT_175TeV_Met100","HT_175TeV_Met250","HT_175TeV_Met250_pt_100","HT_2TeV_Met100","HT_2TeV_Met250","HT_2TeV_Met250_pt_100"};
  //  const char *baseline[25]={"Nocut","photon_selec","Phot_pT_20","nHadJets_2","MET_100","ST_300","bkg_comp","Met_cleaning","lept_veto","veto_chargedTracks","dPhi_MET","jet_pT_Pho_pT","MET_250","pho_pt_100","Final","Pho_pT_30","HT_1TeV_Met250","HT_1TeV_Met250_pt_100","HT_15TeV_Met100","HT_15TeV_Met250","HT_15TeV_Met250_pt_100","HT_175TeV_Met100","nocut_sam","basic_sam"};//"st_300_Met100","pt_st_Met_250","st_300_Met250","nocut"
  //const char *baseline[48]={"Nocut","SignalRegion","ControlRegion","lostElec_SR","lostMu_SR","lostTau_SR","else_pho_SR","lostElec_CR","lostMu_CR","lostTau_CR","else_pho_CR","else_SR","else_CR","lostElecTau_SR","lostMuTau_SR","Tau_hadronic_SR","lostElec_SR_Accept","lostElec_SR_ident","SignalRegion_BDTcut1","ControlRegion_BDTcut1","lostElec_SR_BDTcut1","lostMu_SR_BDTcut1","lostTau_SR_BDTcut1","else_pho_SR_BDTcut1","else_SR_BDTcut1","lostElecTau_SR_BDTcut1","lostMuTau_SR_BDTcut1","Tau_hadronic_SR_BDTcut1","lostElec_CR_BDTcut1","lostMu_CR_BDTcut1","lostTau_CR_BDTcut1","else_pho_CR_BDTcut1","else_CR_BDTcut1","SignalRegion_BDTcut2","ControlRegion_BDTcut2","lostElec_SR_BDTcut2","lostMu_SR_BDTcut2","lostTau_SR_BDTcut2","else_pho_SR_BDTcut2","else_SR_BDTcut2","lostElecTau_SR_BDTcut2","lostMuTau_SR_BDTcut2","Tau_hadronic_SR_BDTcut2","lostElec_CR_BDTcut2","lostMu_CR_BDTcut2","lostTau_CR_BDTcut2","else_pho_CR_BDTcut2","else_CR_BDTcut2"};

  vector<string> baseline = {"NoSelection","PreSelection","woStitch", "PreSelection","Elec_CR","Elec_SR","FailAcep_ElecSR","FailId_ElecSR","FailIso_ElecSR","Mu_CR","Mu_SR","FailAcep_MuSR","FailId_MuSR","FailIso_MuSR","Elec_CR_BDTcut1","Elec_SR_BDTcut1","FailAcep_ElecSR_BDTcut1","FailId_ElecSR_BDTcut1","FailIso_ElecSR_BDTcut1","Mu_CR_BDTcut1","Mu_SR_BDTcut1","FailAcep_MuSR_BDTcut1","FailId_MuSR_BDTcut1","FailIso_MuSR_BDTcut1","Elec_failEtacut_SR","Elec_failpTcut_SR","Mu_failEtacut_SR","Mu_failpTcut_SR","Elec_failEtacut_SR_BDTcut1","Elec_failpTcut_SR_BDTcut1","Mu_failEtacut_SR_BDTcut1","Mu_failpTcut_SR_BDTcut1","Elec_CR_BDTcut2","Elec_SR_BDTcut2","FailAcep_ElecSR_BDTcut2","FailId_ElecSR_BDTcut2","FailIso_ElecSR_BDTcut2","Mu_CR_BDTcut2","Mu_SR_BDTcut2","FailAcep_MuSR_BDTcut2","FailId_MuSR_BDTcut2","FailIso_MuSR_BDTcut2","Elec_SR_valid","Elec_SR_valid_BDTcut1","Mu_SR_valid","Mu_SR_valid_BDTcut1","FailAcep_TauHadronic_MuSR","FailAcep_TauHadronic_MuSR_BDTcut1","FailAcep_Not_TauHadronic_MuSR","FailAcep_Not_TauHadronic_MuSR_BDTcut1","Mu_failAccep_tauHadronic_SR","Mu_failAccep_tauHadronic_SR_BDTcut1"};
  cout<<"size of baseline vector"<<"\t"<<baseline.size()<<endl;
  vector <string> baseline1={"Elect_Inc","Mu_Inc","Tau_Inc","Elect_Inc_v1","Mu_Inc_v1","Tau_Inc_v1","Elect_SR_bin0","Elect_SR_v1_bin1","Elect_SR_bin1","Mu_SR_bin0","Mu_SR_v1_bin0","Mu_SR_bin1","Mu_SR_v1_bin1","Mu_CR","Tau_Inc","Tau_SR","Tau_CR","Tau_ElecSR","Tau_ElecCR","Tau_MuSR","Tau_MuCR","ElecNu_Inc","ElecNu_SR","MuNu_Inc","MuNu_SR","TauNu_Inc","TauNu_SR","ElecNu_CR","MuNu_CR","TauNu_CR",""};

  /* vector<string> baseline_v1 = {"PreSelection","Elec_CR_FR","Elec_CR_promptPho","Elec_CR_NonpromptPho","Elec_CR","Elec_SR_promptPho","Elec_SR_promptPho_GenMatch","Elec_SR_promptPho_NoGenMatch","Elec_SR_promptPho_NoGenMatch_noFR","Elec_SR_promptPho_NoGenMatch_FR","Elec_SR_promptPho_NoGenMatch_Else","Elec_SR_noPromptPho","Elec_SR_noPromptPho_noFR","Elec_SR_noPromptPho_FR","Elec_SR_noPromptPho_FR_else"}; */

  vector<string> baseline_v2 = {"NoSelection","PreSelection","Elec_CR_FR","Elec_CR_promptPho","Elec_CR_NonpromptPho","Elec_CR","Elec_SR_promptPho","Elec_SR_promptPho_GenMatch","Elec_SR_promptPho_noGenPhoInfo","Elec_SR_promptPho_NoGenMatch","Elec_SR_promptPho_NoGenMatch_noFR","Elec_SR_promptPho_NoGenMatch_noFR_tauMatch","Elec_SR_promptPho_NoGenMatch_noFR_qgmatch","Elec_SR_promptPho_NoGenMatch_Else","Elec_SR_promptPho_NoGenMatch_FR","2e_SR_noFR","2e_SR_FR","2e_SR_noFR_else","Elec_SR_noPromptPho","Elec_SR_noPromptPho_noFR","Elec_SR_noPromptPho_noFR_tauMatch","Elec_SR_noPromptPho_noFR_qgmatch","Elec_SR_noPromptPho_Else","Elec_SR_noPromptPho_FR","2e_SR_noPrompt_noFR","2e_SR_noPrompt_FR","2e_SR_noPrompt_noFR_else","Mu_CR_FR","Mu_CR_promptPho","Mu_CR_NonpromptPho","Mu_CR","Mu_SR_promptPho","Mu_SR_promptPho_GenMatch","Mu_SR_promptPho_noGenPhoInfo","Mu_SR_promptPho_NoGenMatch","Mu_SR_promptPho_NoGenMatch_noFR","Mu_SR_promptPho_NoGenMatch_noFR_qgmatch","Mu_SR_promptPho_NoGenMatch_Else","Mu_SR_promptPho_NoGenMatch_FR","Mu_SR_promptPho_1tauHad","Mu_SR_promptPho_2tauHad","Mu_SR_noPromptPho","Mu_SR_noPromptPho_noFR","Mu_SR_noPromptPho_NoGenMatch_noFR_qgmatch","Mu_SR_noPromptPho_NoGenMatch_Else","Mu_SR_noPromptPho_NoGenMatch_FR","Mu_SR_noPromptPho_1tauHad","Mu_SR_noPromptPho_2tauHad","Elec_CR_0Pho","Elec_CR_Phol40","Elec_SR_0Pho","Elec_SR_Phol40","Elec_CR_promptPho_v1","Elec_CR_NonpromptPho_v1","Elec_CR_ElecFake","Elec_CR_Else","Elec_CR_Else_v1","Elec_SR_promptPho_v1","Elec_SR_NonpromptPho_v1","Elec_SR_ElecFake","Elec_SR_Else","Elec_SR_Else_v1","Elec_CR_promptPho_hasTrue","Elec_CR_NonpromptPho_hasTrue","Elec_CR_ElecFake_hasTrue","Elec_CR_Else_hasTrue","Elec_CR_promptPho_hasNot_True","Elec_CR_NonpromptPho_hasNot_True","Elec_CR_ElecFake_hasNot_True","Elec_CR_Else_hasNot_True","Elec_SR_promptPho_hasTrue","Elec_SR_NonpromptPho_hasTrue","Elec_SR_ElecFake_hasTrue","Elec_SR_Else_hasTrue","Elec_SR_promptPho_hasNot_True","Elec_SR_NonpromptPho_hasNot_True","Elec_SR_ElecFake_hasNot_True","Elec_SR_Else_hasNotTrue","Mu_CR_promptPho_v1","Mu_CR_NonpromptPho_v1","Mu_CR_ElecFake","Mu_CR_Else","Mu_CR_Else_v1","Mu_SR_promptPho_v1","Mu_SR_NonpromptPho_v1","Mu_SR_ElecFake","Mu_SR_Else","Mu_SR_Else_v1","TauHad_SR_promptPho_v1","TauHad_SR_NonpromptPho_v1","TauHad_SR_ElecFake","TauHad_SR_Else","TauHad_SR_Else_v1","preSelection_promptPho_v1","preSelection_NonpromptPho_v1","preSelection_ElecFake","preSelection_Else","preSelection_Else_v1","preSelection"};
  cout<<"length of baseline_v2 vector"<<"\t"<<baseline_v2.size()<<endl;
  const char *baseline3[18]={"Wjets_Inclusive","Wjets_SR","Wjets_CR","Wjets_lostElec","Wjets_lostMu","Wjets_lostTau","Wjets_photon","Wjets_unidentified","Elect_Inc","Elect_SR","Elect_CR","Mu_Inc","Mu_SR","Mu_CR","Tau_Inc","Tau_SR","Tau_CR"};//

  const char *baseline2[5]={"WjetsVsElec","WjetsVsMu","WjetsVsTau","TauVsElec","TauVsMu"};
  
  vector<string> checks={"NoSelection","PreSelection","Elec_CR_promptPho_v1","Elec_CR_NonpromptPho_v1","Elec_CR_ElecFake","Elec_CR_Else","Elec_CR_Else_v1","Elec_SR_promptPho_v1","Elec_SR_NonpromptPho_v1","Elec_SR_ElecFake","Elec_SR_Else","Elec_SR_Else_v1","Mu_CR_promptPho_v1","Mu_CR_NonpromptPho_v1","Mu_CR_ElecFake","Mu_CR_Else","Mu_CR_Else_v1","Mu_SR_promptPho_v1","Mu_SR_NonpromptPho_v1","Mu_SR_ElecFake","Mu_SR_Else","Mu_SR_Else_v1","TauHad_SR_promptPho_v1","TauHad_SR_NonpromptPho_v1","TauHad_SR_ElecFake","TauHad_SR_Else","TauHad_SR_Else_v1"};

  /* char check1[1000]; */
  /* for(int c =0; c<checks.size();c++) */
  /*   { */
  /*     sprintf(check1,"h_mindR_elec_pho_%s",checks[c].str()); */
  /*     h_mindR_elec_pho[c] = new TH1F(check1,"mindR(#gamma, e), where e is gen e if SR or reco e if CR",1000,0,10); */
  /*     sprintf(check1,"h_pTratio_elec_pho_%s",checks[c].str()); */
  /*     h_pTratio_elec_pho[c] = new TH1F(check1,"#frac{p_{T}^{#gamma}}{p_{T}^{e}} where e is gen e if SR or reco e if CR",1000,0,10); */
  /*     sprintf(check1,"h_mindR_Vs_pTratio_elec_pho_%s",checks[c].str()); */
  /*     h_mindR_Vs_pTratio_elec_pho[c] = new TH2F(check1,"mindR(#gamma, e) vs #frac{p_{T}^{#gamma}}{p_{T}^{e}} where e is gen e if SR or reco e if CR",1000,0,10,1000,0,10); */
  /*   } */

  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  vector<double> xbins_PhotPt;
  vector<double> xbins_Eta;
  vector<double> xbins_En;
  vector<double> xbins_Phi;
  
  


  /* TH1F *h_GenpT[25]; */
  /* TH1F *h_GenEta[25]; */
  /* TH1F *h_GenPhi[25]; */
  /* TH1F *h_GenEn[25]; */
  for(int i_bin=0;i_bin<200;i_bin++)
    {
      if(i_bin<=20) xbins_PhotPt.push_back(0+(1*(i_bin)));
      if(i_bin>20 && i_bin<=40) xbins_PhotPt.push_back(20+(5*(i_bin-20)));
      if(i_bin>40 && i_bin<=100) xbins_PhotPt.push_back(120+(10*(i_bin-40)));
      if(i_bin>100) xbins_PhotPt.push_back(720+(20*(i_bin-70)));     
      
    }

  //  vector<double> METLowEdge={200,270,350,450,750,2000};
  vector<double> METLowEdge_LL={100,200,300,370,450,600,750,900,2000};
  vector<double> Pho_pt_Edge_LL={30,40,50,70,90,100,150,200,300,500,700,1000};
  vector<double> NHadjets_Edge_LL={2,4,5,6,7,8,9,10,11,12,13,14};
  vector<double> NBjets_Edge_LL={0,1,2,4,5,6,7,8,9,10,11,12,13,14};
  vector<double> ST_Edge_LL={300,400,500,700,900,1000,1500};

  /* for (int i =0; i<14;i++) */
  /*   { */
  /*     sprintf(hname,"h_GenpT_%s",baseline1[i].c_str()); */
  /*     h_GenpT[i] = new TH1F(hname,hname,3000,0,6000);///xbins_PhotPt.size()-1,&(xbins_PhotPt[0])); */
  /*     sprintf(hname,"h_GenEta_%s",baseline1[i].c_str()); */
  /*     h_GenEta[i] = new TH1F(hname,hname,500,-8,8); */
  /*     /\* sprintf(hname,"h_GenPhi_%s",baseline1[i].c_str()); *\/ */
  /*     /\* h_GenPhi[i] = new TH1F(hname,hname,500,-10,10); *\/ */
  /*     /\* sprintf(hname,"h_GenEn_%s",baseline1[i].c_str()); *\/ */
  /*     /\* h_GenEn[i] = new TH1F(hname,hname,3000,0,6000);//xbins_PhotPt.size()-1,&(xbins_PhotPt[0])); *\/ */
  /*     /\* sprintf(hname,"h_GenMET_%s",baseline1[i].c_str()); *\/ */
  /*     /\* h_GenMET[i] = new TH1F(hname,hname,3000,0,6000);//xbins_PhotPt.size()-1,&(xbins_PhotPt[0])); *\/ */
  /*     /\* sprintf(hname,"h_RecoMET_%s",baseline1[i].c_str()); *\/ */
  /*     /\* h_RecoMET[i] = new TH1F(hname,hname,3000,0,6000);//xbins_PhotPt.size()-1,&(xbins_PhotPt[0])); *\/ */

  /*   } */

  /* sprintf(hname,"h_mindR_genVsreco_Elec"); */
  /* h_mindR_genVsreco_Elec  = new TH1F(hname,"mindR(gen e-, reco e-)",1000,0,10); */
  /* sprintf(hname,"h_mindRVspT_recoElec"); */
  /* h_mindRVspT_recoElec = new TH2F(hname,"mindR(gen e-, reco e-)[X-axis] vs p^{reco elec}_{T}",1000,0,10,1000,0,2000); */
  /* sprintf(hname,"h_mindRVspT_genElec"); */
  /* h_mindRVspT_genElec = new TH2F(hname,"mindR(gen e-, reco e-)[X-axis] vs p^{gen elec}_{T}",1000,0,10,1000,0,2000); */
  /* sprintf(hname,"h_mindRVsEta_recoElec"); */
  /* h_mindRVsEta_recoElec = new TH2F(hname,"mindR(gen e-, reco e-)[X-axis] vs Eta^{reco elec}",1000,0,10,1000,-4,4); */
  /* sprintf(hname,"h_mindRVsEta_genElec"); */
  /* h_mindRVsEta_genElec = new TH2F(hname,"mindR(gen e-, reco e-)[X-axis] vs Eta^{gen elec}",1000,0,10,1000,-4,4); */
  /* sprintf(hname,"h_nele"); */
  /* h_nele = new TH1F(hname,"nelectrons_reco",2,0,2); */
  /* sprintf(hname,"h_pT_recovsGen"); */
  /* h_pT_recovsGen = new TH2F(hname,"p_{T}: Gen electron vs reco electron",1000,0,2000, 1000,0,2000); */
  /* sprintf(hname,"h_eta_recoVsGen"); */
  /* h_eta_recoVsGen = new TH2F(hname,"Eta: Gen electron vs reco electron",500,-4,4,500,-4,4); */
  
      
  h_selectBaselineYields_ = new TH1F("cutflows","cutflows",60,-0.5,60.5);
  h_selectBaselineYields_v1 = new TH1F("cutflows_LL","cutflows_LL",60,-0.5,60.5);

}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* dataset, const char* N2_mass, const char* LostlepFlag, const char* phoID) {
  string nameData=dataset;//vvv
  //TDirectory * dir = new TDirectory("TreeMaker2");
    TChain *tree = new TChain("TreeMaker2/PreSelection");
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


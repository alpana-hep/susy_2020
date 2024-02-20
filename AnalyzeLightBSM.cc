#define AnalyzeLightBSM_cxx
#include "AnalyzeLightBSM.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;
using namespace TMVA;
int main(int argc, char* argv[])
{

  if (argc < 6) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"which Lostlep bkg"<< " "<<"Which pho_ID"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  const char *elec = argv[5];
  const char *phoID = argv[6];
  //TString pho_ID = phoID;

  AnalyzeLightBSM ana(inputFileList, outFileName, data,sample, elec,phoID);
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  ana.EventLoop(data,inputFileList,sample,outFileName,elec,phoID);
  Tools::Instance();
  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
  cout<<"inside event loop"<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  TString Elec_flag = elec;
  TString pho_ID_str = phoID;
  
  int pho_ID=-1;
  if(pho_ID_str.Contains("loose"))
    pho_ID=0;
  else if (pho_ID_str.Contains("medium"))
    pho_ID=1;
  else if (pho_ID_str.Contains("tight"))    
    pho_ID=2;
  else if (pho_ID_str.Contains("mva_wp90"))
    pho_ID=3;
  else if (pho_ID_str.Contains("mva_wp80"))
    pho_ID=4;
  cout<<"which photon ID "<<"\t"<<"case"<<"\t"<<pho_ID<<endl;
  // TString pho_D = phoID;
  bool lost_elec_flag=false;
  if(Elec_flag.Contains("Electron"))
    lost_elec_flag = true;
  else
    lost_elec_flag = false;
  
  float counter=0.0;
  TString s_data=data;
 char* s_cross = new char[100];
  sprintf(s_cross,"%s.%s",data,sample);
  TString s_sample= sample;
  
  Long64_t nbytes = 0, nb = 0; 
  int decade = 0;
  int count=0;
  double lumiInfb=137.19;
  int count_QCD=0;
  Long64_t nSurvived = 0,bkg_comp=0,MET_rej=0,nocut=0;
  int genphomatch_after=0,genphomatch_before=0;
  bool genphocheck=false;
  bool v17=true, v12=false;
  bool EWselection=true;
  bool Debug=false;
  bool higMET=true, highdphi=false;
  double deepCSVvalue=0,p0=0,p1=0,p2=0;
  bool applyTrgEff=false;
  if(s_data.Contains("2016preVFP")){ lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV
  if(s_data.Contains("2016postVFP")) { lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016

  if(s_data.Contains("2017")) {lumiInfb=41.48; deepCSVvalue = 0.4506;}
  if(s_data.Contains("2018")){ lumiInfb=59.83;deepCSVvalue = 0.4168;}
  if(s_data.Contains("signal"))lumiInfb= 137.19; //since we only have 2018 dataset

  if(!s_sample.Contains("data"))
    {
      cout<<" analyzing MC bkgs" <<endl;
      if(s_data.Contains("2016preVFP") && highdphi) {lumiInfb=19.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.6001;}
      if(s_data.Contains("2016postVFP") && highdphi) {lumiInfb=16.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.5847;}
      
      if(s_data.Contains("2017") && highdphi) {lumiInfb=41.48; p0=1.82e+02; p1=6.336e+01; p2=9.171e-01; deepCSVvalue = 0.4506;}
      if(s_data.Contains("2018") && highdphi) {lumiInfb=59.83; p0=1.787e+02; p1=6.657e+01; p2=9.47e-01; deepCSVvalue = 0.4168;}

      if(s_data.Contains("2016preVFP") && !highdphi) {lumiInfb=19.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.6001;}
      if(s_data.Contains("2016postVFP") && highdphi) {lumiInfb=16.5;p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.5847;}

      if(s_data.Contains("2017") && !highdphi) {lumiInfb=41.48; p0=1.82e+02; p1=6.336e+01; p2=9.171e-01; deepCSVvalue = 0.4506;}
      if(s_data.Contains("2018") && !highdphi) {lumiInfb=59.83; p0=1.787e+02; p1=6.657e+01; p2=9.47e-01; deepCSVvalue = 0.4168;}

      if(s_data.Contains("FastSim")) lumiInfb=137.19;
      cout<<"Trigger efficiency flag "<<applyTrgEff<<" p0 "<<p0<<"  p1  "<<p1<<"  p2 "<< p2<<endl;

    }
  // if(!s_data.Contains("data"))
  //   {

  //     if(s_data.Contains("2016") && highdphi) {lumiInfb=35.922; p0=1.66539e+02; p1=8.13254e+01; p2=9.71152e-01; };//deepCSVvalue = 0.6321;}
  //     if(s_data.Contains("2017") && highdphi) {lumiInfb=41.529; p0=1.86744e+02; p1=6.74978e+01; p2=9.65333e-01;};// deepCSVvalue = 0.4941;}
  //     if(s_data.Contains("2018") && highdphi) {lumiInfb=59.74; p0=1.89868e+02; p1=6.60434e+01; p2=9.79618e-01;};// deepCSVvalue = 0.4184;}

  //     if(s_data.Contains("2016") && !highdphi) {lumiInfb=35.922; p0=1.67229e+02; p1=8.52729e+01; p2=8.29784e-01; deepCSVvalue = 0.6321;}
  //     if(s_data.Contains("2017") && !highdphi) {lumiInfb=41.529; p0=1.67641e+02; p1=1.21487e+02; p2=9.23864e-01; deepCSVvalue = 0.4941;}
  //     if(s_data.Contains("2018") && !highdphi) {lumiInfb=59.74; p0=1.45529e+02; p1=1.08431e+02; p2=9.27220e-01; deepCSVvalue = 0.4184;}

  //     if(s_data.Contains("FastSim") && s_data.Contains("2016")) lumiInfb=137.19;                                                                                    
  //    }
    if(s_sample.Contains("data"))
    {
      cout<<" analyzing data" <<endl;

      if(s_data.Contains("2016preVFP")) {deepCSVvalue = 0.6001;}
      if(s_data.Contains("2016postVFP")) {deepCSVvalue = 0.5847;}

      if(s_data.Contains("2017")) {deepCSVvalue = 0.4506;}
      if(s_data.Contains("2018")) {deepCSVvalue = 0.4168;}
    }

  cout<<"Calculating the weight for process "<<s_cross<<" deepCSVvalue  "<<deepCSVvalue<<"\t"<<"sdata  "<<s_data<<"\t"<<"s_sample "<<s_sample<<endl;
 cout<<"Trigger efficiency flag "<<applyTrgEff<<" p0 "<<p0<<"  p1  "<<p1<<"  p2 "<< p2<<endl;
  std::string s_process = s_cross;
  TString s_Process = s_process;
  double cross_section = getCrossSection(s_process);
  cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_cross<<"\t"<<s_process<<endl;
 //  cout<<"Event"<<"\t"<<"par-ID "<<"\t"<<"parentID"<<"\t"<<"GenMET"<<"\t"<<"MET"<<"\t"<<"pT"<<"\t"<<"Eta"<<"\t"<<"Phi"<<"\t"<<"E"<<endl;
 //  float met=0.0,st=0.0, njets=0, btags=0,mTPhoMET=0.0,dPhi_PhoMET=0.0,dPhi_MetJet=0.0;
  
 //  // TMVA::Reader *reader1 = new TMVA::Reader();
 //  // reader1->AddVariable( "MET", &met );
 //  // reader1->AddVariable( "NhadJets", &njets );
 //  // reader1->AddVariable( "BTags", &btags );
 //  // reader1->AddVariable("mTPhoMET_",&mTPhoMET);
 //  // reader1->AddVariable("dPhi_PhoMET_",&dPhi_PhoMET);
 //  // reader1->AddVariable("dPhi_Met_Jet",&dPhi_MetJet);
 //  // reader1->AddVariable( "ST", &st );
 //  // reader1->BookMVA( "BDT_100trees_2maxdepth method", "TMVAClassification_BDT_100trees_2maxdepth.weights.xml");
  float nsurVived=0.0;
  int searchBin=0, Tfbins=0;
  float nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0,wt_LL=0.0;

 //  // counters for events yields after each selection/rejection
  float nEvents_Selec[100]={};
  const char* out_nEventsTags[100] ={};
  const char* filetag[42]={"TTGJets_2018","TTGJets_2017","TTGJets_2016postVFP","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016postVFP","Run2_WGJets","TTJets_2018","TTJets_2017","TTJets_2016postVFP","Run2_TTJets","WJets_2018","WJets_2017","WJets_2016postVFP","Run2_WJets","WGJets+WJets_2018","WGJets+WJets_2017","WGJets+WJets_2016postVFP","Run2_WGJets+WJets","TTGJets+TTJets_2018","TTGJets+TTJets_2017","TTGJets+TTJets_2016postVFP","Run2_TTGJets+TTJets","W+TTBar_2018","W+TTBar_2017","W+TTBar_2016postVFP","W+TTBar_FullRun2","TTGJets_2016preVFP","WGJets_2016preVFP","TTJets_2016preVFP","WJets_2016preVFP","WGJets+WJets_2016preVFP","TTGJets+TTJets_2016preVFP","W+TTBar_2016preVFP","TTGJets_2016","WGJets_2016","TTJets_2016","WJets_2016","WGJets+WJets_2016","TTGJets+TTJets_2016","W+TTBar_2016"};

  //const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
  char* hname = new char [200];
  sprintf(hname, "Lepton_LL_TFbins_v1_nJetsBjets_phoID_loose_09Jan24.root");
  char* hname1 = new char [200];
  sprintf(hname1, "Lepton_LL_TFbins_v2_nJetsBjets_PhoPt_phoID_loose_09Jan24.root");//"TF_Electron_LLEstimation_binsV0_phoID_%s_08Aug23.root",phoID);
  char* hname2 = new char [200];
  sprintf(hname2, "Lepton_LL_TFbins_v3_nJetsBjets_MET_phoID_loose_09Jan24.root");// "TF_allin1_LLEstimation_binsV0_phoID_%s_08Aug23.root",phoID);
  TFile* f_TFbins= new TFile(hname);
  TFile* f_TFbins_v1 =new TFile(hname1);
  TFile* f_TFbins_v2= new TFile(hname2);
  char* histname = new char[2000];
  TH1F* h_TF;
  TH1F* h_TF1;
  TH1F* h_TF2;
 auto sample1="";
 if(s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-130")|| s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-40to130"))
   sample1 = "WGJets";
 else
   sample1 = sample;
 sprintf(histname,"h_TFbins_LL_W+TTBar_FullRun2");//h_TFbins_LL_W+TTBar_FullRun2");//2018",data);                                                                   
 h_TF = (TH1F*)f_TFbins->Get(histname);
 h_TF1 = (TH1F*)f_TFbins_v1->Get(histname);
 h_TF2 = (TH1F*)f_TFbins_v2->Get(histname);
 cout<<"Reading TF for lost electron estimation:  "<<"\t"<<histname<<endl;
 for(int i=0; i<h_TF->GetNbinsX();i++)
   {cout<<"TFBIns_v0"<<"\t"<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}

 for(int i=0; i<h_TF1->GetNbinsX();i++)
   {cout<<"TFBIns_v1"<<"\t"<<i<<"\t"<<h_TF1->GetBinContent(i)<<endl;}
 for(int i=0; i<h_TF2->GetNbinsX();i++)
   {cout<<"TFBIns_v2"<<"\t"<<i<<"\t"<<h_TF2->GetBinContent(i)<<endl;}

 // if(lost_elec_flag)
 //   {
 //     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_Elec__WGJets_2018");//%s_%s",sample1,data);
 //     h_TF = (TH1F*)f_elec->Get(histname);
 //     cout<<"Reading TF for lost electron estimation:  "<<"\t"<<histname<<endl;
 //     for(int i=0; i<h_TF->GetNbinsX();i++)
 //       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}
 //   }
 // else
 //   {
 //     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_FailAcep_TauHadronic_WGJets_2018");//%s_%s",sample1,data);
 //     cout<<"Reading TF for lost muon estimation:  "<<"\t"<<histname<<endl;
 //     h_TF = (TH1F*)f_muon->Get(histname);
 //     for(int i=0; i<h_TF->GetNbinsX();i++)
 //       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}
 //   }

  int coutt=0;
  //nentries=1000;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
     {
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
        cout << 10 * k << " %" << endl;
      decade = k;
      
      // ===============read this entry == == == == == == == == == == ==                                                                        
      Long64_t ientry = LoadTree(jentry);
      if(Debug)
	cout<<"===load tree entry ==="<<jentry<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(Debug)
        cout<<"===load tree entry ==="<<jentry<<endl;

      // ==========================================================================================\\
      // ====================================== Calculating weight per event ======================\\
      // ==========================================================================================\\
      //cout<<"weight setting"<<"\t"<<jentry<<endl;
      // if(s_sample.Contains("MCMC_86_7257_1750"))
      //   wt = (0.002991*lumiInfb*1000)/nentries; //gluino mass = 1750 GeV                                                                                                      
      // else if(s_sample.Contains("MCMC_86_7257_1400"))
      //   wt = (0.0284*lumiInfb*1000)/nentries; //need to check luminosity                                                                                                      
      // else if(s_sample.Contains("MCMC_86_7257_1500"))
      //   wt = (0.0157*lumiInfb*1000)/nentries;
      // else if(s_sample.Contains("MCMC_86_7257_1600"))
      // 	wt = (0.00887*lumiInfb*1000)/nentries;
      // else if(s_sample.Contains("MCMC_86_7257_1700"))
      //   wt = (0.00507*lumiInfb*1000)/nentries;
      // else  if(s_sample.Contains("MCMC_86_7257_1800"))
      //   wt = (0.00293*lumiInfb*1000)/nentries;
      // else if(s_sample.Contains("MCMC_86_7257_1900"))
      //   wt = (0.00171*lumiInfb*1000)/nentries;
      // else if(s_sample.Contains("MCMC_86_7257_2000"))
      //   wt = (0.00101*lumiInfb*1000)/nentries;
      //cout<"Alps check "<<s_Process<<endl;
      //if(MET>100)
      //h_MET->Fill(MET);
      if(s_Process.Contains("2018.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2018.WGJets_MonoPhoton_PtG-130UL")|| s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-40to130UL") ||s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-130UL") || s_Process.Contains("2017.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2017.WGJets_MonoPhoton_PtG-130UL") ||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-130UL")||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-40to130UL"))//|| s_Process.Contains("2018.TTGJets_incUL")) //(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018") || s_sample.Contains("WGJets"))
      	{	  
	  if(jentry==0)
      	    cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_process<<endl;
      	  wt = cross_section*lumiInfb*1000.0;//)/nentries; // Weight*lumiInfb*1000.0; //(cross_section*lumiInfb*1000.0)/nentries;
	  //if (jentry<10) cout<<"Alpana urgent check "<<wt<<"\t"<<cross_section<<"\t"<<lumiInfb<<endl;
      	}
      else// if (s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))// || s_sample.Contains("WGJets"))
      	{
      	  wt = Weight*lumiInfb*1000.0;
      	}
     if (jentry<10) cout<<"Printing weight befor apply trig eff "<<wt<<"\t"<<cross_section<<"\t"<<lumiInfb<<endl;
      if(s_sample.Contains("data")) wt =1;
      if(jentry<10)
	{
	  cout<<"Event "<<jentry<<" before event trig  "<<wt<<"\t"<<Weight<<" TriggerPass->size() "<<TriggerPass->size()<<endl;
	  // for(int i=0;i<TriggerPass->size();i++)
	  //   cout<<"index "<<i<<" (*TriggerPass)[i] "<<(*TriggerPass)[i]<<"\t"<<(*TriggerVersion)[i]<<endl;
	}
      bool tighte_trgpass=true;
      if(s_sample.Contains("data") && applyTrgEff)
	{
	  wt=1;
	  if(TriggerPass->size()!=213) continue;
	  if((*TriggerPass)[183]==1 || (*TriggerPass)[164]==1 || (*TriggerPass)[166]==1 || (*TriggerPass)[170]==1 || (*TriggerPass)[172]==1 || (*TriggerPass)[174]==1) tighte_trgpass=true;
	  else tighte_trgpass=false;	  
	  if(tighte_trgpass==false)  continue;
	}
      if(!s_sample.Contains("data") && !s_sample.Contains("signal") && applyTrgEff)
	{
	  wt = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
	}
      if(Debug)
	cout<<"alpa"<<"\t"<<jentry<<"\t"<<endl;
      if(jentry<10)
	{    
	  cout<<"Event "<<jentry<<" event weight"<<"\t"<<wt<<"\t"<<cross_section<<"\t"<<Weight<<"lumiInfb  "<<lumiInfb<<endl;
	  cout<<"Event "<<jentry<<"\t"<<wt<<"\t"<<"MET  "<<MET<<"\t"<<(((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2)<<endl;
	}
      h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
      nocut++;
      if(Debug)
    	cout<<"Before filling lorentz vectors"<<endl;
      //filling the lorentz vector variables 
      
      if(s_sample.Contains("UL")){
    	//set objetc pointer
    	Electrons_v1.clear();
    	GenElectrons_v1.clear();
    	GenParticles_v1.clear();
    	GenMuons_v1.clear();
    	GenTaus_v1.clear();
    	GenJets_v1.clear();//_v1.clear();
    	Electrons_v1.clear();
    	Photons_v1.clear();
    	Muons_v1.clear();
    	Taus_v1.clear();
    	Jets_v1.clear();
    	if(Debug)
    	  cout<<jentry<<"\t"<<"After setting pointer  lorentz vectors"<<endl;
	Electrons_v1 = getLorentzVector(Electrons_,Electrons_fCoordinates_fPt,Electrons_fCoordinates_fEta,Electrons_fCoordinates_fPhi,Electrons_fCoordinates_fE);
	Muons_v1 = getLorentzVector(Muons_,Muons_fCoordinates_fPt,Muons_fCoordinates_fEta,Muons_fCoordinates_fPhi,Muons_fCoordinates_fE);
	Photons_v1 = getLorentzVector(Photons_,Photons_fCoordinates_fPt,Photons_fCoordinates_fEta,Photons_fCoordinates_fPhi,Photons_fCoordinates_fE);
	Jets_v1 = getLorentzVector(Jets_,Jets_fCoordinates_fPt,Jets_fCoordinates_fEta,Jets_fCoordinates_fPhi,Jets_fCoordinates_fE);
	
	if(!s_sample.Contains("data")) {
	  GenElectrons_v1 = getLorentzVector(GenElectrons_,GenElectrons_fCoordinates_fPt,GenElectrons_fCoordinates_fEta,GenElectrons_fCoordinates_fPhi,GenElectrons_fCoordinates_fE);
	  GenMuons_v1 = getLorentzVector(GenMuons_,GenMuons_fCoordinates_fPt,GenMuons_fCoordinates_fEta,GenMuons_fCoordinates_fPhi,GenMuons_fCoordinates_fE);
	  GenTaus_v1 = getLorentzVector(GenTaus_,GenTaus_fCoordinates_fPt,GenTaus_fCoordinates_fEta,GenTaus_fCoordinates_fPhi,GenTaus_fCoordinates_fE);
	  GenJets_v1 = getLorentzVector(GenJets_,GenJets_fCoordinates_fPt,GenJets_fCoordinates_fEta,GenJets_fCoordinates_fPhi,GenJets_fCoordinates_fE);
	  GenParticles_v1 = getLorentzVector(GenParticles_,GenParticles_fCoordinates_fPt,GenParticles_fCoordinates_fEta,GenParticles_fCoordinates_fPhi,GenParticles_fCoordinates_fE);
	}

      if(Debug)
    	cout<<jentry<<"\t"<<"After filling  reco Jets"<<endl;

      }
      
      int count=0;
      int branch_size = GenParticles_v1.size();
      // vector<TLorentzVector> v1_GenParticles;
      // v1_GenParticles = getLorentzVector(GenParticles_,GenParticles_fCoordinates_fPt,GenParticles_fCoordinates_fEta,GenParticles_fCoordinates_fPhi,GenParticles_fCoordinates_fE);
      // if(Debug){
      // cout<<"Entry "<<jentry<<"\t"<<"Event  "<<EvtNum<<" Run  "<<RunNum<<" Gen particles  "<<"\t"<<branch_size<<"\t"<<endl;
      // 	//	cout<<"Entry "<<jentry<<"\t"<<"Method 2 Gen particles  "<<"\t"<<v1_GenParticles.size()<<"\t"<<endl;
      
      // for (int i=0;i<branch_size;i++)
      // 	{
      // 	  cout<<GenParticles_v1[i].Pt()<<"\t"<<GenParticles_v1[i].Eta()<<"\t"<<GenParticles_v1[i].Phi()<<"\t"<<GenParticles_v1[i].E()<<endl;
      // 	  //	  cout<<"Method-2 "<<v1_GenParticles[i].Pt()<<"\t"<<v1_GenParticles[i].Eta()<<"\t"<<v1_GenParticles[i].Phi()<<"\t"<<v1_GenParticles[i].E()<<endl;
      // 	}

      // }

      int ele_branch=Electrons_;
      if(Debug)
        cout<<"Electrons "<<"\t"<<Electrons_<<endl;

      int pho_branch = Photons_;
      if(Debug)
    	        cout<<"Photons "<<"\t"<<Photons_<<endl;
           
    //   //variables to be used in getting dR
      float iGenEle_eta =99999.0, iGenEle_phi=99999.0, iRecEle_eta =99999.0, iRecEle_phi=99999.0, iRecphot_eta =99999.0, iRecphot_phi=99999.0,iGen_Wpt=99999.0,iGen_Wp=99999.0,iGen_WEta=99999.0,iGen_WPhi=99999.0;
      float dR_Ele=0, min_dR=9999, min_dR_pho=9999;
      int pdgID=0, parentId=0, status;
      int Ids[4]={11,13,15,22};
      const char* ids[4]={"e","mu","tau","pho"};
      // loop over gen particles
      int count_genElec=0,count_genMu=0,count_genTau=0, igenPho=0, count_genEl_tau=0,count_genMu_tau=0,count_haddecay=0;
      TLorentzVector genPho_W,genElec_W,genMu_W,genTau_W,genElec_Tau,genMu_Tau,genW,genNuElec_W,genNuMu_W,genNuTau_W,genElecNu_Tau,genMuNu_Tau;
      int elec_reco=0,elec_reco0_before=0,elec_reco1_before=0,muon_reco=0,elec_gen3=0,elec_gen2=0, elec_gen=0, muon_gen=0,elec_reco0=0,elec_reco1=0,evtSurvived_preselec=0,elec_reco2=0,elec2_reco=0,survived_vetohad=0,elec_reco1_CR=0,survived_elecge1=0,events_cr=0,events_sr=0,total=0,remain=0,elec_reco0_genel=0,efakepho=0,ele=0,genphomatch_after=0,genphomatch_before=0,elec_gen4=0,gentauhad2=0,lep2=0,lep=0;
      TLorentzVector leadGenPho,genPho,genLep,v_lep1,v_lep2,v_genMu,v_genNuMu,v_genElec,v_genNuElec,v_genTau,v_genNuTau,v_genW,v_genNu;
    
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////   Making gen level collections of all particles  ////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(Debug)
    	cout<<"cp - before gen particles"<<endl;
      // for (int igen=0; igen<branch_size; igen++)
      // 	{	  
      // 	  pdgID= (*GenParticles_PdgId)[igen];
      // 	  parentId=(*GenParticles_ParentId)[igen];
      // 	  // if(Debug)
      // 	  //   {	    cout<<"printing Gen particles parent ID and pdgID"<<endl;
      // 	  //     cout<<parentId<<"\t"<<pdgID<<endl;
      // 	  //     cout<<GenParticles_fCoordinates_fPt[igen]<<"\t"<<GenParticles_fCoordinates_fEta[igen]<<endl;}
	  
      // 	  if(abs(pdgID)==24)
      // 	    {h_W_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_W_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_W_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_W_E->Fill(GenParticles_v1[igen].E(),wt); v_genW = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==11 && abs(parentId)==24){h_elec_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_elec_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_elec_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_elec_E->Fill(GenParticles_v1[igen].E(),wt);v_genElec = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==12 && abs(parentId)==24){h_Nu_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_Nu_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_Nu_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_Nu_E->Fill(GenParticles_v1[igen].E(),wt);v_genNuElec = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==13 && abs(parentId)==24){h_mu_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_mu_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_mu_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_mu_E->Fill(GenParticles_v1[igen].E(),wt); v_genMu=GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==15 && abs(parentId)==24){h_tau_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_tau_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_tau_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_tau_E->Fill(GenParticles_v1[igen].E(),wt); v_genTau=GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==14 && abs(parentId)==24) v_genNuMu = GenParticles_v1[igen];
      // 	  if(abs(pdgID)==16 && abs(parentId)==24) v_genNuTau = GenParticles_v1[igen];
      // 	  if((abs(pdgID)==14 || abs(pdgID)==16 || abs(pdgID)==12) && abs(parentId)==24) v_genNu= GenParticles_v1[igen];
      // 	}

      // h_2d_elec_nu_pT->Fill(v_genElec.Pt(),v_genNuElec.Pt(),wt);
      // h_2d_mu_nu_pT->Fill(v_genMu.Pt(),v_genNuMu.Pt(),wt);
      // h_2d_Tau_nu_pT->Fill(v_genTau.Pt(),v_genNuTau.Pt(),wt);
      // h_2d_WvsElec_pT->Fill(v_genW.Pt(),v_genElec.Pt(),wt);
      // h_2d_WvsElecNu_pT->Fill(v_genW.Pt(),v_genNuElec.Pt(),wt);
      // h_2d_WvsMu_pT->Fill(v_genW.Pt(),v_genMu.Pt(),wt);
      // h_2d_WvsMuNu_pT->Fill(v_genW.Pt(),v_genNuMu.Pt(),wt);
      // h_2d_WvsTau_pT->Fill(v_genW.Pt(),v_genTau.Pt(),wt);
      // h_2d_WvsTauNu_pT->Fill(v_genW.Pt(),v_genNuTau.Pt(),wt);
      
      // h_2d_elec_nu_Eta->Fill(v_genElec.Eta(),v_genNuElec.Eta(),wt);
      // h_2d_mu_nu_Eta->Fill(v_genMu.Eta(),v_genNuMu.Eta(),wt);
      // h_2d_Tau_nu_Eta->Fill(v_genTau.Eta(),v_genNuTau.Eta(),wt);
      // h_2d_WvsElec_Eta->Fill(v_genW.Eta(),v_genElec.Eta(),wt);
      // h_2d_WvsElecNu_Eta->Fill(v_genW.Eta(),v_genNuElec.Eta(),wt);
      // h_2d_WvsMu_Eta->Fill(v_genW.Eta(),v_genMu.Eta(),wt);
      // h_2d_WvsMuNu_Eta->Fill(v_genW.Eta(),v_genNuMu.Eta(),wt);
      // h_2d_WvsTau_Eta->Fill(v_genW.Eta(),v_genTau.Eta(),wt);
      // h_2d_WvsTauNu_Eta->Fill(v_genW.Eta(),v_genNuTau.Eta(),wt);
      // h_2d_nu_MET_pT->Fill(v_genNu.Pt(),MET,wt);
      //h_MET->Fill(MET,wt);
      // h_GenMET->Fill(GenMET,wt);
      // h_GenJets->Fill(GenJets_,wt);
      if(Debug)
      cout<<"===load tree entry check1 - after filling gen level at entry ==="<<"\t"<<jentry<<endl;           
      
    //   //selecting electron and muons
      int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
      TLorentzVector genPho1,genEle1,genMu1,genTau1,recElec, recMu;//,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
      vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_genLep1,v_recEle,v_recMu;
      v_genEle2.clear();
      v_genEle1.clear();
      v_genMu1.clear();
      v_genMu2.clear();
      v_genTau1.clear();
      v_genTau2.clear();

           int leadGenPhoIdx=-100, mu_index=-100;
      int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
      int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
      if(!s_sample.Contains("data")){
    	for(int i=0;i<GenParticles_;i++){
    	  if(GenParticles_v1[i].Pt()!=0){
    	    if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=100) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 )
    	      {
    		leadGenPhoIdx = i;
    		genPho1 = (GenParticles_v1[i]);
    		v_genPho1.push_back(genPho1);
    	      }
    	  }
    	}
    	bool hadtau = false;
    	for(int i=0 ; i < GenElectrons_; i++)
    	  {
    	    if(GenElectrons_v1[i].Pt()!=0)// && (*GenElectrons)[i].Pt()!=inf)                                                                                                                                          
    	      {
    		nGenEle1++;
    		genEle1 = (GenElectrons_v1[i]);
    		v_genEle2.push_back(genEle1);
    		v_genLep2.push_back(genEle1);
    		//h_GenpT[3]->Fill(genEle1.Pt()); h_GenEta[3]->Fill(genEle1.Eta());
    	      }
    	  }
    	for(int i=0 ; i < GenMuons_; i++)
    	  {
    	    if(GenMuons_v1[i].Pt()!=0)// && (*GenMuons)[i].Pt()!=inf)                                                                                                                                                  
    	      {
    		nGenMu1++;
    		genMu1 = (GenMuons_v1[i]);
    		v_genMu2.push_back(genMu1);
    		v_genLep2.push_back(genMu1);
    		//h_GenpT[4]->Fill(genMu1.Pt()); h_GenEta[4]->Fill(genMu1.Eta());
    	      }
    	  }
    	//cout<<"entry:"<<"\t"<<jentry<<"\t"<<nGenMu1<<"\t"<<"genMuon size"<<"\t"<<v_genMu2.size()<<endl;                                                                                                              

    	for(int i=0 ; i < GenTaus_; i++)
    	  {
    	    if(GenTaus_v1[i].Pt()!=0) // && (*GenTaus)[i].Pt()!=inf)                                                                                                                                                   
    	      {
    		nGenTau1++;
    		genTau1 = (GenTaus_v1[i]);
    		v_genTau2.push_back(genTau1);
    		//v_genLep2.push_back(genTau1);                                                                                                                                                                        
    		//h_GenpT[5]->Fill(genTau1.Pt()); h_GenEta[5]->Fill(genTau1.Eta());
    		if((*GenTaus_had)[i])
    		  nGentau_had1++;

    	      }
    	  }
      }
      sortTLorVec(&v_genTau2);
      sortTLorVec(&v_genMu2);
      sortTLorVec(&v_genEle2);
      sortTLorVec(&v_genLep2);

      if(Debug)
    	cout<<"=== load entry "<<jentry<<"\t"<<"GenTaus_  "<<GenTaus_<<endl;
      //      int nlep=0;
      int total_lost_el = 0,cr_el=0,sr_el,e_index=-1,nlep=0, NgenElec=0,leadGenPhoIdx1=0;
      bool Elec_passEtacut=false,Elec_passpTcut=false,Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false, genElec_passEtacut=false,genElec_passpTcut=false,genElec_passAccep = false, genElec_passId= false, genElec_passIso = false,genMu_passEtacut=false,genMu_passpTcut=false,genMu_passAccep = false, genMu_passId= false, genMu_passIso = false;
      bool Mu_passEtacut=false,Mu_passpTcut=false, Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
      int fail_isoElec=0;
      nlep=0; nelec_reco=0;
      sortTLorVec(&Electrons_v1);
      sortTLorVec(&Muons_v1);

      for(int i=0;i<Electrons_;i++)
    	{
    	  if(nelec_reco>0) continue;
    	  if((Electrons_v1[i].Pt()>10) && abs(Electrons_v1[i].Eta()) < 2.5){
	    Elec_passAccep = true;
	    Elec_passEtacut=true; Elec_passpTcut=true;
	    if( (*Electrons_passIso)[i]==1)
    	    {	      
    	      pass_isoElec++; Elec_passIso = true; nelec_reco++; nlep++; e_index=i; recElec=Electrons_v1[i]; v_recEle.push_back(Electrons_v1[i]);
    	      // h_recoElec_pT->Fill(recElec.Pt(),wt);  h_recoElec_Eta->Fill(recElec.Eta(),wt);	     h_recoElec_Eta->Fill(recElec.Eta(),wt);
    	      // h_recoElec_Phi->Fill(recElec.Phi(),wt);
    	    }	    
    	  else if((*Electrons_passIso)[i]!=1)
    	    {	      fail_isoElec++; Elec_passIso = false;}
	  }
	  else
	    {
	      if((Electrons_v1[i].Pt()<10) )       	 Elec_passpTcut = false;
	      if(abs(Electrons_v1[i].Eta()>2.5))		Elec_passEtacut=false;
	    }
    	}
      sortTLorVec(&v_recEle);
      if(Debug)      cout<<"entry: "<<"\t"<<jentry<<" "<<nlep<<"\t"<<"reco e size"<<" "<<v_recEle.size()<<" "<<Electrons_v1.size()<<"\t"<<NElectrons<<endl;
      for(int i=0;i<Muons_;i++)
    	{
	  if(nmu_reco>0) continue;
    	  if((Muons_v1[i].Pt()>10) && abs(Muons_v1[i].Eta()) < 2.5 && nmu_reco<1){
	    Mu_passAccep = true;
            Mu_passEtacut=true; Mu_passpTcut=true;
	    if((*Muons_passIso)[i]==1)
    	    {   	      pass_isoMu++; Mu_passIso = true;nmu_reco++; nlep++; mu_index=i; recMu=Muons_v1[i]; v_recMu.push_back(Muons_v1[i]);
    	      // h_recoMu_pT->Fill(recMu.Pt(),wt);   	      h_recoMu_Eta->Fill(recMu.Eta(),wt);
    	      // h_recoMu_Eta->Fill(recMu.Eta(),wt);    	      h_recoMu_Phi->Fill(recMu.Phi(),wt);
    	    }
    	  else if((*Muons_passIso)[i]!=1)
    	    { 	fail_isoMu++; Mu_passIso = false;}
	  }
	  else{
	    if((Muons_v1[i].Pt()<10) )             Mu_passpTcut = false;
	    if(abs(Muons_v1[i].Eta()>=2.5))                Mu_passEtacut=false;
	  }

	}
      if(Debug && v_recMu.size()==1)      cout<<"entry: "<<"\t"<<jentry<<" "<<nlep<<"\t"<<"reco e size"<<" "<<v_recEle.size()<<" "<<Muons_v1.size()<<"\t"<<NMuons<<endl;
      if(Electrons_==0 && nelec_reco!=0) //checking if the contribution to SR with no reco e- is coming up or not?                              
    	cout<<"Entry "<< jentry<<"\t"<<nelec_reco<<endl;
      float ratio =0.0, ratio1=0.0, mindr_genElecPho=-999, mindr_=-9999;
      int count_genEle=0,count_recEle=0;
      vector<TLorentzVector> goodPho_n;
      vector<int> goodPhoIndx_n;
      // if(Debug){
      // 	cout<<" Photons_cutBasedIDFall17V2 "<<endl;
      // 	cout<<" Photons_cutBasedIDFall17V2 "<<Photons_cutBasedID->size()<<endl;
      // }
      for(int iPho=0;iPho<Photons_;iPho++){
    	if((*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
    	  goodPho_n.push_back( Photons_v1[iPho] );
    	  goodPhoIndx_n.push_back(iPho);
    	}
      }
      int nPhotons = goodPho_n.size();
      TLorentzVector bestPhoton=getBestPhoton(pho_ID);      
      // *******************  Selecting Jet objects ********************************//
       int minDRindx=-100,phoMatchingJetIndx=-100,hadJetID=-999,bJet1Idx=-100,nHadJets=0;
      double minDR=99999,ST=0,Ht=0;
      vector<TLorentzVector> hadJets,bjets;
      hadJets.clear();
      bjets.clear();
      vector<int> jetMatchindx;
      bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
      if(Debug)
        cout<<"===load tree entry check2 at entry ==="<<"\t"<<jentry<<endl;
      for(int i=0;i<Jets_;i++)
    	{
    	  if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
    	    double dR=bestPhoton.DeltaR(Jets_v1[i]);
    	    if(dR<minDR){minDR=dR;minDRindx=i;}
    	  }
      }
      if(Debug)
    	cout<<"===load tree entry  ==="<<"\t"<<jentry<<"\t"<<"Jets check == "<<minDR<<endl;
      
    for(int i=0;i<Jets_;i++)
      { if(Debug)
    	cout<<"  = Jets.Pt()  ==  "<<Jets_v1[i].Pt()<<"\t"<< " = Jets.Eta() == "<<Jets_v1[i].Eta()<<endl;
    	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
    	if(Debug)
	  cout<< "==== loadjets ==="<<"\t"<<i<<"\t"<<minDR<<endl;
          if( !(minDR < 0.3 && i==minDRindx) )
    	    { 
	      hadJetID= (*Jets_ID)[i];	      
    	      //if(true){
	      if(hadJetID){
		hadJets.push_back(Jets_v1[i]);
    		  // hadJets_hadronFlavor.push_back((*Jets_hadronFlavor)[i]);
    		  // hadJets_HTMask.push_back((*Jets_HTMask)[i]);
    		  // hadJets_bJetTagDeepCSVBvsAll.push_back((*Jets_bJetTagDeepCSVBvsAll)[i]);
    		  // if(q==1) leadjet_qmulti=(*Jets_chargedMultiplicity)[q];
    		  // if(q==1) leadjet_Pt=(*Jets)[q].Pt();
	      //cout<<"  (*Jets_bJetTagDeepCSVBvsAll)[i]  "<<(*Jets_bJetTagDeepCSVBvsAll)[i]<<endl;
		if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
		  bjets.push_back(Jets_v1[i]); bJet1Idx = i;}
		// hadJets.push_back((*Jets)[i]);
		jetMatchindx.push_back(i);
	      }
    	    }
    	}
      }
    if(hadJets.size()==0) continue;
    if(Debug)
      cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;
    sortTLorVec(&hadJets);
    BTags = bjets.size();
    // if(BTags!= BTagsDeepCSV){
    //   coutt+=wt;
    //   cout<<"Entry: "<<jentry<<"\t"<<"BTags "<<BTags<<"\t"<<" BTagsDeepCSV "<<BTagsDeepCSV<<"\t deepCSVvalue "<<deepCSVvalue<<"\t"<<hadJets.size()<<endl;}
    if( minDR<0.3 ) {phoMatchingJetIndx=minDRindx; recoJetMatch_recoPho=true;}
    double genmindr=99999, recojetmindr=99999;
    if(Debug)
      cout<<"====load entry === "<<jentry<<"\t"<<"phoMatchingJetIndx ===  "<<phoMatchingJetIndx<<endl;
    //genmindr = MinDr(bestPhoton,GenJets_v1);
    if(Debug)
      cout<<"===load tree entry === "<<"\t"<<jentry<<" hadJets.size()  "<< hadJets.size()<<endl;    
    for(int i=0;i<hadJets.size();i++){
      if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());Ht=Ht+(hadJets[i].Pt());}
      if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
    }
    sortTLorVec(&hadJets);
    if(hadJets.size()<BTags)
      cout<<"Entry: "<<jentry<<"\t"<<"njets "<<nHadJets<<"\t"<<" hadJets.size() "<<hadJets.size()<<endl;
    // // ********* This is to account into all visible energy: Adding photon matching with Jet ******************//
    if( minDR<0.3 ) ST=ST+bestPhoton.Pt(); 
    // //transverse mass of best photon+MET                                                                                  
    double mTPhoMET=0.0;//sqrt(2*(bestPhoton.Pt())*MET*(1-cos(DeltaPhi(METPhi,bestPhoton.Phi()))));
    // // ******* Dphi between MET and Best Photon     ******************** //                                                                 
    double dPhi_PhoMET= 0.0;//abs(DeltaPhi(METPhi,bestPhoton.Phi()));
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" mTPhoMET "<<mTPhoMET<<" METPhi  "<<bestPhoton.Eta()<<" Photons_ "<<Photons_<<" hadJets.size() " <<"\t"<<hadJets.size()<<" NElectrons "<<NElectrons<< "  Jets  "<<Jets_<<"\t"<< "Jets size  "<<Jets_v1.size()<<endl;
    // h_leadJet_pT->Fill(hadJets[0].Pt(),wt);
    // h_subleadJet_pT->Fill(hadJets[1].Pt(),wt);    
    // h_leadJet_Eta->Fill(hadJets[0].Eta(),wt);
    // h_subleadJet_Eta->Fill(hadJets[1].Eta(),wt);
    h_Njets[0]->Fill(nHadJets,wt);
    h_Nbjets[0]->Fill(BTags,wt);
    h_MET_[0]->Fill(MET,wt);
    //h_PhotonPt[0]->Fill(bestPhoton.Phi(),wt);
    h_Mt_PhoMET[0]->Fill(mTPhoMET,wt);
    //h_dPhi_PhoMET[0]->Fill(dPhi_PhoMET,wt);
    h_St[0]->Fill(ST,wt);
    h_HT[0]->Fill(HT,wt);
    if(jentry<10)
    cout<<"madMinPhotonDeltaR  "<<madMinPhotonDeltaR<<endl;
    h_madminPhotonDeltaR_noSelection->Fill(madMinPhotonDeltaR,wt);
    
    TLorentzVector Met;
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    Double_t deta_jet_pho= 0.0,deta_jet_met=0.0,deta_met_pho=0.0;
    if(Debug)
      cout<<"METPhi "<<METPhi<<endl;
    if(    goodPho_n.size()==0 || Photons_==0) continue;
    
    // h_recoPho_pT->Fill(bestPhoton.Pt(),wt); h_recoPho_Eta->Fill(bestPhoton.Eta(),wt);
    // h_recoPho_Phi->Fill(bestPhoton.Phi(),wt);
    h_PhotonPt[0]->Fill(bestPhoton.Pt(),wt);
    //deta_jet_pho = abs(bestPhoton.Eta()-hadJets[0].Eta());
    double mT= 0.0, dPhi_METjet1=5, dPhi_METjet2=5, dPhi_phojet1=5, dPhi_phojet2=5, dPhi_phoMET=5;
    // if(hadJets.size() > 0) dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    // if(hadJets.size() > 1) dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    // dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    // dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));

    dPhi_METjet1 = abs(Met.DeltaPhi(hadJets[0]));
    dPhi_METjet2 = abs(Met.DeltaPhi(hadJets[1]));
    //    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    // dPhi_METjet1=3.8,dPhi_METjet2=3.8;                                                                                                         
    // Met.SetPtEtaPhiE(MET,0,METPhi,0);
    // // dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    // // dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    // //    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    // dPhi_METjet1 = abs(DeltaPhi(METPhi,(hadJets)[0].Phi()));
    // dPhi_METjet2 = abs(DeltaPhi(METPhi,(hadJets)[1].Phi()));
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" hadJets.size() "<<hadJets.size()<<endl;



    h_selectBaselineYields_->Fill("all",wt);
    if(bestPhotonIndxAmongPhotons<0) continue;
    bool bestPhoHasPxlSeed=true;
    if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
    if( bestPhoHasPxlSeed ) continue;
    //    cout<<"Event : " <<jentry <<"\t"<<(*Photons_mvaValuesID)[bestPhotonIndxAmongPhotons] << " (*Photons_mvaValuesID)[iPho]"<< "\t"<<bestPhotonIndxAmongPhotons<< " Photons _ "<<Photons_ <<"\t"<<Photons_mvaValuesID->size()<<endl;

    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" bestPhoHasPxlSeed "<<bestPhoHasPxlSeed<<endl;
    if(hadJets.size() == 0) continue;

    h_selectBaselineYields_->Fill("hadsize!=0",wt);
    if(bestPhoton.Pt()>40)    h_selectBaselineYields_->Fill("#gamma_{pt}>40",wt);
    else continue;

    //cout<<"Event : " <<jentry <<"\t"<<(*Photons_mvaValuesID)[bestPhotonIndxAmongPhotons] << " (*Photons_mvaValuesID)[iPho]"<< "\t"<<bestPhotonIndxAmongPhotons<< " Photons _ "<<Photons_ <<"\t Photons_mvaValuesID->size() "<<Photons_mvaValuesID->size()<<"\t"<<Photons_v1.size()<<"\t"<<Photons_fullID->size()<<"\t"<<Photons_cutBasedID->size()<<endl;
    if (nHadJets>=2)h_selectBaselineYields_->Fill("njets>=2",wt);
    else continue;
    if(MET>200) 
      h_selectBaselineYields_->Fill("MET>200",wt);
    else continue;
    if(ST>300)
      h_selectBaselineYields_->Fill("ST>300",wt);
    else continue;
    //recommended MET filtser for UL - taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_should
    //check nvtx filter - less efficiency for 2016 - https://indico.cern.ch/event/1057110/#27-met-filters-performance-stu
    if(PFCaloMETRatio >=  5) continue;
    if (s_data.Contains("2017") || s_data.Contains("2018"))
      if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 && NVtx>0)) continue;
    if (s_data.Contains("2016")){
      if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 &&EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 )) continue;
    }

    h_selectBaselineYields_->Fill("MetCleaning",wt);
    
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
        h_selectBaselineYields_->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    if(phoMatchingJetIndx>=0 && (Jets_v1[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
    if(phoMatchingJetIndx<0) continue;
    if(higMET)
      { 
    	if (bestPhoton.Pt()<=100)continue;
    	if(MET<=200) continue;                                                                                              
      }
    //    FillHistogram_Kinematics(2,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,1);
	 
	 if(!s_sample.Contains("data")){
    if(Debug)
      cout<<"===loading tree entry === "<<jentry<<"\t"<<"  ===after preselection  ==="<<endl;
    h_madminPhotonDeltaR_beforeStitching->Fill(madMinPhotonDeltaR,wt);
    //remove the overlapped events from different MC samples as discussed in these slides - https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf

    if((s_sample.Contains("TTJets_HT")||s_sample.Contains("TTJets-HT")) && madHT<600) continue;        
    if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept")) && madHT>600) continue;
    if(!genphocheck)      {        genphomatch_before++;
        double mindr_Pho_genlep=getGenLep(bestPhoton);
        if( s_sample.Contains("TTG") )
          {  if(!hasGenPromptPhoton)              {		
    		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;	    }
            else if(hasGenPromptPhoton)              {
	      h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
	      if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
		if(madMinPhotonDeltaR<0.5)		    h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt);
		  if(mindr_Pho_genlep<0.5)		    h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt);		
    		  continue;}
    		else
    		  {    		    if(madMinPhotonDeltaR >= 0.5)    		      h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
    		    if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);}
	    }
	  }                                             
        if(s_sample.Contains("WG"))
          {            if(!hasGenPromptPhoton){                  		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;              }
            else if(hasGenPromptPhoton)              {	      h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
	      if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
		  if(madMinPhotonDeltaR<0.5)                    h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt);
                  if(mindr_Pho_genlep<0.5)                    h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt);
    		  continue;}
    		else
                  {           if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
                  }}}
        if(s_sample.Contains("WJets"))
          { if(!hasGenPromptPhoton){h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
	      if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}              
            else if(hasGenPromptPhoton){  h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);   
    		  continue;}
    		else{                  
                    if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);
		}}}
                        
        if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets") || s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept"))
          {            if(!hasGenPromptPhoton){           
    		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}              
            else if(hasGenPromptPhoton){             
    		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){// h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
    		  continue;}
    		else{                  
                    if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)                      h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);}                          }}
        if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
          {            if(!(madMinPhotonDeltaR>0.4)) continue;}
          
        if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
          {            if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;}         
    	if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG")) || s_sample.Contains("ZNuNuGJets")))          {
	  if(!(madMinPhotonDeltaR>0.5)) continue;}
          
        if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
          {            if(!(madMinPhotonDeltaR<=0.5)) continue;}          
        genphomatch_after++;
      }
    h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt); 
	 }
    h_selectBaselineYields_v1->Fill("Pre-Selection",wt);
    FillHistogram_Kinematics(1,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
     	 //    h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt);
    if(Debug)
      cout<<" just before photon identification - prompt/non-prompt ========   ===="<<jentry<<endl;
    ////////////////////////////////////////////////////////////////////////////////////
    //// ======= ======= photon identification - prompt/non-prompt ======== ////////////
    ////// based on the following studies   ///////////////////////////////////////////
    // https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf
    int flag_phoOrigin = -1;
    if(!s_sample.Contains("data")){
      flag_phoOrigin =Photons_OriginType();
      if(isoMuonTracks ==0 &&  isoPionTracks==0 && isoElectronTracks==0) {
	if(NElectrons==0 && NMuons==0)
	  {
	    if(flag_phoOrigin==0)
	      FillHistogram_Kinematics(2,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	    else if(flag_phoOrigin==1)
	      FillHistogram_Kinematics(3,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	    else if(flag_phoOrigin==2)
	      FillHistogram_Kinematics(4,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	    else if(flag_phoOrigin==3)
	      FillHistogram_Kinematics(5,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);	  
	  }
      }
    }
    // if(jentry<10000)
    //   cout<<"Entry: "<<jentry<<"\t"<<bestPhotonIndxAmongPhotons<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;
    
    if(Debug)
      cout<<"just before 1 lepton CR + with 1photon ====== "<<jentry<<endl;
    ////////////////////////////////////////////////////////////////////////////////////
    //  ====  ========== ======== 1 lepton CR + with 1photon  =====================  ///
    ///////////////////////////////////////////////////////////////////////////////////
    
    bool elec_CR= false; bool mu_CR = false;
    bool elec_SR= false; bool mu_SR = false;
    int TFbins = getBinNoV1_le(nHadJets,BTags);
    int TFbins1 = getBinNoV16_le(nHadJets,BTags,bestPhoton.Pt());;
    int TFbins2= getBinNoV7_le(nHadJets,BTags); ;
    double wt_LL = wt*h_TF->GetBinContent(TFbins+1);
    double wt_LL1 = wt*h_TF1->GetBinContent(TFbins1+1);
    double wt_LL2 = wt*h_TF2->GetBinContent(TFbins2+1);
    
    if (lost_elec_flag && nelec_reco == 1 && nmu_reco == 0){
      if(isoMuonTracks !=0 || isoPionTracks!=0) continue; // veto muon/pions from 1 electron CR
      if(Debug)
	cout<<"in Electron CR + with 1photon ====== "<<jentry<<endl;

      double dr2 = bestPhoton.DeltaR(Electrons_v1[e_index]); // FR elec-CR(1e,0pho), Photon SR (0e, 1pho and dR(gen e,pho)<0.2)
      if(dr2<=0.2)continue;
      if(Debug)
        cout<<"in Electron CR + with 1photon after dR cut ====== "<<jentry<<endl;

      double mTElecMET=sqrt(2*(Electrons_v1[e_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,Electrons_v1[e_index].Phi()))));
      if(mTElecMET>100) { continue;}//h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);continue;} // remove signal contamination 
      elec_CR = true;      
      if(elec_CR){
      FillHistogram_Kinematics(6, nHadJets, BTags, bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      FillHistogram_Kinematics_varBin(2,nHadJets, BTags, bestPhoton.Pt(),ST,wt);
      //      int TFbins = getBinNoV1_le(nHadJets,BTags);
      int searchBin = getBinNoV6_WithOnlyBLSelec(nHadJets,BTags);
      if(nHadJets>=7)
	cout<<"Entry: "<<jentry<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<bestPhoton.Pt()<<"\t"<<wt<<endl;
      //double wt_LL = wt*h_TF->GetBinContent(TFbins+1);
      FillTFBins_Valid(6, nHadJets, BTags,wt_LL, wt_LL1, wt_LL2, bestPhoton.Pt(),ST);//
      // h_TFbins_ElecLL_validation[6]->Fill(TFbins+1,wt_LL);
      // h_TFbins_ElecLL_validation_v1[6]->Fill(TFbins+1,wt_LL);
      // h_Sbins_LL_Validation[6]->Fill(searchBin,wt_LL);
      //      cout<<"TFbins " <<TFbins<<"\t"<<h_TF->GetBinContent(TFbins+1)<<"\t"<<wt_LL<<endl;
      }
      elec_CR = true;
    }
    
    if (!lost_elec_flag && nelec_reco == 0 && nmu_reco == 1){
      if(Debug) cout<<"=== inside muon cr condition === "<<jentry<<endl;
      if(isoElectronTracks !=0 || isoPionTracks!=0) continue; // veto muon/pions from 1 electron CR                                           
      if(Debug)
        cout<<"in Muon CR + with 1photon ====== "<<jentry<<endl;

      double mTElecMET=sqrt(2*(Muons_v1[mu_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,Muons_v1[mu_index].Phi()))));
      if(mTElecMET>100) { continue;}//h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);continue;} // remove signal contamination               
      double dr2 = bestPhoton.DeltaR(Muons_v1[mu_index]); 
      if(dr2<=0.2)continue;
      if(Debug)
        cout<<"in Muon CR + with 1photon after dR cut ====== "<<jentry<<endl;

      mu_CR = true;
      if(mu_CR){     
      FillHistogram_Kinematics_varBin(4,nHadJets, BTags, bestPhoton.Pt(),ST,wt);
      FillHistogram_Kinematics(7, nHadJets, BTags, bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      //int TFbins = getBinNoV1_le(nHadJets,BTags);
      int searchBin = getBinNoV6_WithOnlyBLSelec(nHadJets,BTags);
      if(searchBin==35 || searchBin==34)
	cout<<"entry: "<<jentry<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<searchBin<<endl;
      if(nHadJets>=7)
	cout<<"Entry: "<<jentry<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<bestPhoton.Pt()<<"\t"<<wt<<endl;
   
      //double wt_LL = wt*h_TF->GetBinContent(TFbins+1);
      FillTFBins_Valid(7, nHadJets, BTags,wt_LL,wt_LL1, wt_LL2, bestPhoton.Pt(),ST);
      // h_TFbins_ElecLL_validation[7]->Fill(TFbins+1,wt_LL);
      // h_TFbins_ElecLL_validation_v1[7]->Fill(TFbins+1,wt_LL);
      // h_Sbins_LL_Validation[7]->Fill(searchBin,wt_LL);
      } 
    }
    
     
    // filling the histograms for various different scenarios
    if(elec_CR && MET>100){ //filling events based on reco-photon origin
      if(Debug)
        cout<<"in Electron CR + with 1photon  before mindR ====== "<<jentry<<"\t"<<" v_genEle2.size()  "<<v_genEle2.size()<<"GenElectrons_ "<<GenElectrons_<<"Electrons_ "<<Electrons_<<endl;
      h_mindR_elec_pho[2]->Fill(bestPhoton.DeltaR(Electrons_v1[e_index]),wt);
      h_pTratio_elec_pho[2]->Fill(Electrons_v1[e_index].Pt()/bestPhoton.Pt(),wt);
      h_mindR_Vs_pTratio_elec_pho[2]->Fill( bestPhoton.DeltaR(Electrons_v1[e_index]),Electrons_v1[e_index].Pt()/bestPhoton.Pt(),wt);

	//Transfer Factor bins 
	// int TFbins_v1 = getBinNoV1_le(BTags,nHadJets);
	// int TFbins_v2 = getBinNoV1_le(BTags,nHadJets);
	// int searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
      if(!s_sample.Contains("data")){
      if(flag_phoOrigin==0)
	FillHistogram_Kinematics(23,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==1)
	FillHistogram_Kinematics(24,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==2)
	FillHistogram_Kinematics(25,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==3)
	FillHistogram_Kinematics(26,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      }
    }
    
    if(mu_CR && MET>100){ //filling events based on reco-photon origin
      if(Debug)
        cout<<"in Muon CR + with 1photon  before mindR ====== "<<jentry<<"\t"<<" v_genEle2.size()  "<<v_genEle2.size()<<"GenMuons_ "<<GenMuons_<<"Muons_ "<<Muons_<<endl;

      if(v_genEle2.size()!=0){h_mindR_elec_pho[4]->Fill(bestPhoton.DeltaR(Muons_v1[mu_index]),wt);	
	h_pTratio_elec_pho[4]->Fill(Muons_v1[mu_index].Pt()/bestPhoton.Pt(),wt);
	h_mindR_Vs_pTratio_elec_pho[4]->Fill( bestPhoton.DeltaR(Muons_v1[mu_index]),Muons_v1[mu_index].Pt()/bestPhoton.Pt(),wt);}
      if(!s_sample.Contains("data")){
      if(flag_phoOrigin==0)
        FillHistogram_Kinematics(32,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==1)
        FillHistogram_Kinematics(33,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==2)
        FillHistogram_Kinematics(34,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==3)
        FillHistogram_Kinematics(35,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      }
    }

    if(Debug) cout<<"=== Just before 0 lepton SR === "<<jentry<<endl;
    
    ////////////////////////////////////////////////////////////////////////////////////                                                        
    //  ====  ========== ======== 0 lepton + 1photon SR for LL  =====================  ///                                                       
    ///////////////////////////////////////////////////////////////////////////////////     
    bool fakeElectron=false, proc=false, proc1=false;
    if (!s_sample.Contains("data") && lost_elec_flag && nelec_reco == 0 && nmu_reco == 0){
      if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
      if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue; //reject w->qq events or all hadronic decay of TTbar
      if(nGenMu1>0) continue;
      if(nGenEle1==0) continue;
      if(((s_sample.Contains("WGJets") || s_sample.Contains("WJets"))&& nGentau_had1>0) ||( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")) && nGentau_had1>1)) continue; //removing tau-hadronic from lost electron SR - these events will be considered in lost muon SR -- allowing genTaHad+1 electron SR cases.
      
      sortTLorVec(&v_genEle2);
      v_genEle1 =v_genEle2;
      if(Debug) cout<<"=== Just before fakerate condition in elec SR === "<<jentry<<endl;

      proc = v_genEle2.size()==1;//(s_sample.Contains("WGJets")|| s_sample.Contains("WJets"));                                                  
      proc1 = v_genEle2.size()>1;//( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")));                                            
      fakeElectron =(proc && (bestPhoton.DeltaR(v_genEle2[0])<0.2 && bestPhoton.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[0].Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && bestPhoton.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[0].Pt()>0.8)|| (v_genEle2.size()>1 && (bestPhoton.DeltaR(v_genEle2[1])<0.2 && bestPhoton.Pt()/v_genEle2[1].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[1].Pt()>0.8))));

      //      fakeElectron =(proc && (bestPhoton.DeltaR(v_genEle2[0])<0.2 && v_genEle2[0].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[0].Pt()/bestPhoton.Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && v_genEle2[0].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[0].Pt()/bestPhoton.Pt()>0.8)|| (v_genEle2.size()>1 && (bestPhoton.DeltaR(v_genEle2[1])<0.2 && v_genEle2[1].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[1].Pt()/bestPhoton.Pt()>0.8))));
      //if(fakeElectron) continue;// should contribute in FR estimation
      if(!fakeElectron) elec_SR = true;      
      if(!fakeElectron && MET>100) {h_mindR_elec_pho[3]->Fill(bestPhoton.DeltaR(v_genEle2[0]),wt);
	h_pTratio_elec_pho[3]->Fill(v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
	h_mindR_Vs_pTratio_elec_pho[3]->Fill( bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);}
      if(elec_SR && MET>100){ FillHistogram_Kinematics(8, nHadJets, BTags, bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	 FillHistogram_Kinematics_varBin(3,nHadJets, BTags, bestPhoton.Pt(),ST,wt);}
      // Events falling into e SR due to failing different acceptance, ID, ISo criteria
      double dR = 9999;
      if(elec_SR){
	for(int i=0;i<v_genEle2.size();i++){
	  if(v_genEle2[i].Pt()>10 && abs(v_genEle2[i].Eta())<2.5){ genElec_passAccep=true; genElec_passEtacut=true; genElec_passpTcut=true;}
	  else {genElec_passAccep = false; if(v_genEle2[i].Pt()<10) genElec_passpTcut = false; if(abs(v_genEle2[i].Eta())>2.5) genElec_passEtacut=false;}		  
	  if(genElec_passAccep){ if(Electrons_v1.size()!=0) dR = MinDr(v_genEle2[i],Electrons_v1);
	    if(dR>0.2) genElec_passId= false; else  genElec_passId= true;
	    if(dR<0.2 && Elec_passIso) genElec_passIso=true; else genElec_passIso=false;
	  }
	}    
      }
    }
     
//mu_SR
    bool tauHad_SR=false;
    if(Debug) cout<<"=== Just before muon SR === "<<jentry<<endl;
    if (!s_sample.Contains("data") && !lost_elec_flag && nelec_reco == 0 && nmu_reco == 0){
      if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue;//to reject W->qq' type of events
      if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
      if(nGenEle1>1) continue;
      if(nGenMu1==0 && v_genTau2.size()==0 && nGenEle1==1) continue;
      if(nGenMu1==0 && v_genTau2.size()==1 && nGenEle1==1) continue; //considered in lost electron SR
      if(Debug) cout<<"=== Just before FR condition muon SR === "<<jentry<<endl;

      proc = (v_genEle2.size()==1);
      proc1 = v_genEle2.size()>1;
      bool fakeElectron=false;
      if(v_genEle2.size()!=0 ){
	fakeElectron =(proc && (bestPhoton.DeltaR(v_genEle2[0])<0.2 && bestPhoton.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[0].Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && bestPhoton.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[0].Pt()>0.8)|| (v_genEle2.size()>1 && (bestPhoton.DeltaR(v_genEle2[1])<0.2 && bestPhoton.Pt()/v_genEle2[1].Pt()<1.2 && bestPhoton.Pt()/v_genEle2[1].Pt()>0.8))));

      }	//	fakeElectron =(proc && (bestPhoton.DeltaR(v_genEle2[0])<0.2 && v_genEle2[0].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[0].Pt()/bestPhoton.Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && v_genEle2[0].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[0].Pt()/bestPhoton.Pt()>0.8)|| (v_genEle2.size()>1 && (bestPhoton.DeltaR(v_genEle2[1])<0.2 && v_genEle2[1].Pt()/bestPhoton.Pt()<1.2 && v_genEle2[1].Pt()/bestPhoton.Pt()>0.8))));}
      //checking minDr condition for fake electron
      if(!fakeElectron && v_genEle2.size()!=0 && MET>100) {h_mindR_elec_pho[5]->Fill(bestPhoton.DeltaR(v_genEle2[0]),wt);
      h_pTratio_elec_pho[5]->Fill(v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
      h_mindR_Vs_pTratio_elec_pho[5]->Fill( bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);}
      //if(fakeElectron) continue;
      if(!fakeElectron) mu_SR = true;
      if(mu_SR && MET>100) FillHistogram_Kinematics_varBin(5,nHadJets, BTags, bestPhoton.Pt(),ST,wt);

      if(((s_sample.Contains("WGJets") || s_sample.Contains("WJets"))&& nGentau_had1>0) ||( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")) && (nGentau_had1>1 || (nGentau_had1>0 && nGenMu1==0)))) // nGenTauHad>0 scenario for cases where T-->hadTau and T--> qq
	tauHad_SR = true;
      if(mu_SR && !tauHad_SR && MET>100) FillHistogram_Kinematics(14, nHadJets, BTags, bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      if(mu_SR && tauHad_SR && MET>100) FillHistogram_Kinematics(20, nHadJets, BTags, bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      // Events falling into e SR due to failing different acceptance, ID, ISo criteria                                                          
      double dR = 9999;
      if(mu_SR){
      for(int i=0;i<v_genMu2.size();i++){
	if(v_genMu2[i].Pt()>10 && abs(v_genMu2[i].Eta())<2.5){ genMu_passAccep=true; genMu_passEtacut=true; genMu_passpTcut=true;}
        else {genMu_passAccep = false; if(v_genMu2[i].Pt()<10) genMu_passpTcut = false; if(abs(v_genMu2[i].Eta())>2.5) genMu_passEtacut=false;}
	if(genMu_passAccep){ if(Muons_v1.size()!=0) dR = MinDr(v_genMu2[i],Muons_v1);
	  if(dR>0.2) genMu_passId= false; else  genMu_passId= true;
	 if(dR<=0.2 && Mu_passIso) genMu_passIso=true; else genMu_passIso=false;
        }
      }
      }
      }
    
    // filling the histograms for various different scenarios based on photon origin in that event
	 if(!s_sample.Contains("data")){
    if(elec_SR && MET>100){
      if(flag_phoOrigin==0)
        FillHistogram_Kinematics(27,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==1)
        FillHistogram_Kinematics(28,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==2)
        FillHistogram_Kinematics(29,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==3)
        FillHistogram_Kinematics(30,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
    }

    if(mu_SR && !tauHad_SR && MET>100){
      if(flag_phoOrigin==0)
        FillHistogram_Kinematics(35,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==1)
        FillHistogram_Kinematics(36,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==2)
        FillHistogram_Kinematics(37,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(flag_phoOrigin==3)
        FillHistogram_Kinematics(38,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
    }
    if(tauHad_SR && mu_SR && MET>100)
      {
	if(flag_phoOrigin==0)
	  FillHistogram_Kinematics(39,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	else if(flag_phoOrigin==1)
	  FillHistogram_Kinematics(40,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	else if(flag_phoOrigin==2)
	  FillHistogram_Kinematics(41,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	else if(flag_phoOrigin==3)
	  FillHistogram_Kinematics(42,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      }

    ///   filling histogram for avrious sub-categories of different signal regions
    if(elec_SR && MET>100) { if(!genElec_passAccep) FillHistogram_Kinematics(9,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	// if(!genElec_passpTcut)  FillHistogram_Kinematics(13,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	// else if(!genElec_passEtacut) FillHistogram_Kinematics(12,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);}
	else if(!genElec_passId && genElec_passAccep) FillHistogram_Kinematics(10,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	else if(!genElec_passIso && genElec_passId && genElec_passAccep) FillHistogram_Kinematics(11,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt); }
                                
    if(mu_SR && !tauHad_SR && MET>100) { if(!genMu_passAccep) FillHistogram_Kinematics(15,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	// if(!genMu_passpTcut)  FillHistogram_Kinematics(19,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
	// else if(!genMu_passEtacut) FillHistogram_Kinematics(18,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);}
      else if(!genMu_passId && genMu_passAccep) FillHistogram_Kinematics(16,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt);
      else if(!genMu_passIso && genMu_passId && genMu_passAccep) FillHistogram_Kinematics(17,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,wt); }
	 }


      
    nsurVived+=wt;
    } //loop over entries

  if(Debug)
    cout<<"filling the branches in tree"<<endl;
   //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
  cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;

  cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<"\t"<< coutt<<endl;
  cout << "============   Electron     ===============" <<endl;

  cout<<"applied weights "<<" "<<wt<<endl;
  cout<<"CR electron :==  "<<"\t"<<nCR_elec<<"\t"<<endl; 
  cout<<"SR electron :==  "<<"\t"<<nSR_elec<<"\t"<<endl;
  cout<<"SR e- : Fail acceptance "<<"\t"<<FailAccept_Elec<<"\t"<<endl;
  cout<<"SR e- : Fail Id " <<"\t"<<FailId_Elec<<"\t"<<endl;
  cout<< "SR e- : Fail Iso" <<"\t"<<FailIso_Elec<<"\t"<<endl;

  cout << "============   Muon     ===============" <<endl;

  // cout<<"applied weights "<<" "<<wt<<endl;
  // cout<<"CR muon :==  "<<"\t"<<nCR_mu<<"\t"<<endl;
  // cout<<"SR muon :==  "<<"\t"<<nSR_mu<<"\t"<<endl;
  // cout<<"SR mu : Fail acceptance "<<"\t"<<FailAccept_Mu<<"\t"<<endl;
  // cout<<"SR mu : Fail Id " <<"\t"<<FailId_Mu<<"\t"<<endl;
  // cout<< "SR mu : Fail Iso" <<"\t"<<FailIso_Mu<<"\t"<<wt*FailIso_Mu<<endl;

  // for (int i =0;i<49;i++)
  //   {
  //     cout<<out_nEventsTags[i]<<" :"<<"\t"<<"=====: "<<"\t"<<nEvents_Selec[i]<<endl;
      
  //     cout<<" "<<"\t"<<endl;
      
  //   }

  // for (int i =0;i<49;i++)
  //   {
  //     cout<<nEvents_Selec[i]<<endl;
  //   }  
  // cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;
  // outfile->cd();
  // skim_tree->Write();
  // outfile->Close();
  // cout<<"outFile: "<<outfileName<<" written!!"<<endl;

}


int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nHadJets,int nbjets)
{
  
  int sBin=-100,m_i=0;
  if(nbjets==0 ){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=7;}
    else if(nHadJets>=7)               { sBin=13;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=18;}
    else if(nHadJets==5 || nHadJets==6){ sBin=23;}
    else if(nHadJets>=7)               { sBin=28;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<199.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
	break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7         ;
        break; }
    }
  }
  else if(sBin==7 || sBin==33 || sBin==39){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
      if(METLowEdge_v3_1[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;break;}
      else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6; break; }
    }
  }

  else 
    {
      for(int i=0;i<METLowEdge_v3_2.size()-1;i++){
	if(METLowEdge_v3_2[i]<199.99) continue;
	m_i++;
	if(MET >= METLowEdge_v3_2[i] && MET < METLowEdge_v3_2[i+1]){ sBin = sBin+m_i;break; }
	else if(MET >= METLowEdge_v3_2[METLowEdge_v3_2.size()-1])  { sBin = sBin+5; break; }
      }
    }
  // if(sBin==0){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<199.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  //       break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()])  { sBin = 6;
  //       break; }
  //   }
  //   else if (sBin==7){
  //     for(int i=0;i<METLowEdge_v3.size()-1;i++){
  // 	if(METLowEdge_v3[i]<199.99) continue;
  // 	int sBin1=sBin;
  // 	m_i++;
  // 	if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  // 	  break; }
  // 	else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 6;
  // 	  break; }

  //   }
    
  // }
  // else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3_1[i]<199.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
  //       break;}
  //     else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
  //       break; }
  //   }
  // }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }                                                                                                                    

  //   }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }                                                                                                                    

//     }
//   }
// -
  // int sBin=-100,m_i=0;
  // if(nbjets==0 ){
  //   if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
  //   else if(nHadJets==5 || nHadJets==6){ sBin=8;}
  //   else if(nHadJets>=7)               { sBin=15;}
  // }
  // else{
  //   if(nHadJets>=2 && nHadJets<=4)     { sBin=22;}
  //   else if(nHadJets==5 || nHadJets==6){ sBin=29;}
  //   else if(nHadJets>=7)               { sBin=36;}
  // }
  // if(sBin==0){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  // 	break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
  // 	break; }
  //   }
  // }
  // else if(sBin==8 || sBin==15 || sBin==22 || sBin==29 || sBin==36){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
  //     if(METLowEdge_v3_1[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
  // 	break;}
  //     else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
  // 	break; }
  //   }
  // }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }

  // }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }

  //   }
  //}

  return sBin;
}
int AnalyzeLightBSM::getBinNoV7_le(int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(MET<300)sBin=1; else if (MET>300) sBin=2;}
    else if(nHadJets==3)     { if(MET<300)sBin=3; else if (MET>300) sBin=4;}
    else if(nHadJets==4)     { if(MET<300)sBin=5; else if (MET>300) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<300)sBin=7; else if (MET>300) sBin=8;}
    else if(nHadJets>=7)   { if(MET<300)sBin=9; else if (MET>300) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(MET<300)sBin=11; else if (MET>300) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<300)sBin=13; else if (MET>300) sBin=14;}
    else if(nHadJets>=7)   { if(MET<300)sBin=15; else if (MET>300) sBin=16;}
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV16_le(int nHadJets, int nbjets, double photon_pT){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(photon_pT<100)sBin=1; else if (photon_pT>100) sBin=2;}
    else if(nHadJets==3)     { if(photon_pT<100)sBin=3; else if (photon_pT>100) sBin=4;}
    else if(nHadJets==4)     { if(photon_pT<100)sBin=5; else if (photon_pT>100) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(photon_pT<100)sBin=7; else if (photon_pT>100) sBin=8;}
    else if(nHadJets>=7)   { if(photon_pT<100)sBin=9; else if (photon_pT>100) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(photon_pT<100)sBin=11; else if (photon_pT>100) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(photon_pT<100)sBin=13; else if (photon_pT>100) sBin=14;}
    else if(nHadJets>=7)   { if(photon_pT<100)sBin=15; else if (photon_pT>100) sBin=16;}
  }
  return sBin;
}


int AnalyzeLightBSM::getBinNoV1_le( int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2)     { sBin=1;}
    else if(nHadJets==3)     { sBin=2;}
    else if(nHadJets==4)     { sBin=3;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=4;}
    else if(nHadJets>=7)   { sBin=5;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      { sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=7;}
    else if(nHadJets>=7)   { sBin=8;}
  }
  return sBin;
}
int AnalyzeLightBSM::getBinNoV7(int nHadJets, int btags){
  int sBin=-100,m_i=0;
  if(btags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=4;}
    else if(nHadJets>=7)               { sBin=8;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=12;}
    else if(nHadJets>=5){ sBin=16;}
    //    else if(nHadJets>=7)               { sBin=26;}
  }
  //  if(sBin==-99)//
  //  cout<<sBin<<"\t"<<nHadJets<<"\t"<<btags<<endl;
  if(sBin==0){
    for(int i=0;i<METLowEdge_lowMET.size()-1;i++){
      if(METLowEdge_lowMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_lowMET[i] && MET < METLowEdge_lowMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_lowMET[METLowEdge_lowMET.size()-1])  { sBin = 4         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge_lowMET.size()-1;i++){
      if(METLowEdge_lowMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_lowMET[i] && MET < METLowEdge_lowMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_lowMET[METLowEdge_lowMET.size()-1])  { sBin = sBin+4   ;break; }
    }
  }
  return sBin;
}
int AnalyzeLightBSM::getBinNoV7_highMET(int nHadJets, int btags){
  int sBin=-100,m_i=0;
  if(btags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=4;}
    else if(nHadJets>=7)               { sBin=8;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=12;}
    else if(nHadJets>=5){ sBin=16;}
    //    else if(nHadJets>=7)               { sBin=26;}                                                                                                                     
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge_highMET.size()-1;i++){
      if(METLowEdge_highMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_highMET[i] && MET < METLowEdge_highMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_highMET[METLowEdge_highMET.size()-1])  { sBin = 4         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge_highMET.size()-1;i++){
      if(METLowEdge_highMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_highMET[i] && MET < METLowEdge_highMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_highMET[METLowEdge_highMET.size()-1])  { sBin = sBin+4   ;break; }
    }
  }
  return sBin;
}

TLorentzVector AnalyzeLightBSM::getBestPhoton(int pho_ID){
  // TLorentzVector v1;
  // vector<TLorentzVector> goodPho;
  // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
  //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
  // }

  // if(goodPho.size()==0) return v1;
  // sortTLorVec(&goodPho);
  // return goodPho[0];
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;
  //  cout<<"\t"<<(*Photons_mvaValuesID).size()<<endl;
  for(int iPho=0;iPho<Photons_;iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    //    cout<<(*Photons_hasPixelSeed)[iPho] <<"\t (*Photons_hasPixelSeed)[iPho]  "<<" (*Photons_fullID)[iPho]"<<"\t"<<(*Photons_fullID)[iPho]<<"\t"<<Photons_v1[iPho].Eta()<<endl;
    if(abs(Photons_v1[iPho].Eta())<2.4 && (*Photons_hasPixelSeed)[iPho]<0.001 && ( ((*Photons_fullID)[iPho] && pho_ID==0) || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42)))
      {
      goodPho.push_back(Photons_v1[iPho] );
      goodPhoIndx.push_back(iPho);
      //      cout<<(*Photons_mvaValuesID)[iPho] << " (*Photons_mvaValuesID)[iPho]"<<endl;
    }
  }
  //  cout<<(*Photons_hasPixelSeed)[iPho] <<"\t(*Photons_hasPixelSeed)[iPho]  "<<endl;
  //cout<<" goodPho.size()" <<"\t"<<  goodPho.size()<<endl;
  int highPtIndx=-100;
  for(int i=0;i<goodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
  }

  if(highPtIndx>=0){
    bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
  }
  else bestPhotonIndxAmongPhotons = -100;
  if(highPtIndx==-100){TLorentzVector v0;return v0;}
  else return goodPho[highPtIndx];

  
}

// TLorentzVector AnalyzeLightBSM::getBestPhoton_tightID(){
//   vector<TLorentzVector> goodPho;
//   vector<int> goodPhoIndx;
//   for(int iPho=0;iPho<Photons->size();iPho++){
//     if((*Photons_hadTowOverEM)[iPho]<0.02148 && (*Photons_sigmaIetaIeta)[iPho]<0.00996 &&  (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
//       goodPho.push_back( (*Photons)[iPho] );
//       goodPhoIndx.push_back(iPho);
//     }
//   }

//   int highPtIndx=-100;
//   for(int i=0;i<goodPho.size();i++){
//     if(i==0) highPtIndx=0;
//     else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
//   }

//   if(highPtIndx>=0){
//     bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
//   }
//   else bestPhotonIndxAmongPhotons = -100;
//   if(highPtIndx==-100){TLorentzVector v0;return v0;}
//   else return goodPho[highPtIndx];
// }



double AnalyzeLightBSM::getGendRLepPho(int pho_ID){//MC only
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  // vector<TLorentzVector> goodPho;
  // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
  //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back((*Photons)[iPhoton]);
  // }
  // if(goodPho.size()!=0) 
  genPho1 =getBestPhoton(pho_ID);
  
  for(int i=0;i<GenParticles_;i++){
     if(GenParticles_v1[i].Pt()!=0){
       if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	 if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
       }
     }
  }//for
  if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
  else return 1000.0;
}

 double AnalyzeLightBSM::getdR_GenPho_RecoPho(TLorentzVector bestPhoton)
 {
   
   TLorentzVector genPho;
   int leadGenPhoIdx=-100;
   int minDR = 9999;
   vector<TLorentzVector> v_genPho;
   for (int igen=0; igen<GenParticles_; igen++)
     {
       if(GenParticles_v1[igen].Pt()!=0){
	 if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
	   genPho = (GenParticles_v1[igen]);
	   v_genPho.push_back(genPho);
	 }
       }
     }
   return MinDr(bestPhoton,v_genPho);   
 }
double AnalyzeLightBSM::getGendRElecPho(int pho_ID){//MC only                                                                                             
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  genPho1 =getBestPhoton(pho_ID);
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      if( abs((*GenParticles_PdgId)[i])==11 && (abs((*GenParticles_ParentId)[i])<=25 )&& (abs((*GenParticles_ParentId)[i])!=15)){
	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
      }
    }
  }//for                                                                                                                                       
  if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
  else return 1000.0;

}
double AnalyzeLightBSM::getGenLep1(TLorentzVector bestPhoton, bool flag){                                     
  vector<TLorentzVector> v_genLep2;
  TLorentzVector genMu1, genEle1;
  if(flag)
    {
      for(int i=0 ; i < GenElectrons_; i++)
	{
	  if(GenElectrons_v1[i].Pt()!=0)
	    {
	      genEle1 = (GenElectrons_v1[i]);
	      v_genLep2.push_back(genEle1);
	    }

	}
    }
  else
    {
      for(int i=0 ; i < GenMuons_; i++)
	{
	  if(GenMuons_v1[i].Pt()!=0)
	    {
	      genMu1 = (GenMuons_v1[i]);
	      v_genLep2.push_back(genMu1);
	    }
	}
    }
  return MinDr(bestPhoton,v_genLep2);
}


double AnalyzeLightBSM::getGenLep(TLorentzVector bestPhoton){
  vector<TLorentzVector> v_genLep2;
  TLorentzVector genMu1, genEle1;
  // if(flag)
  //   {
      for(int i=0 ; i < GenElectrons_; i++)
        {
          if(GenElectrons_v1[i].Pt()!=0)
            {
              genEle1 = (GenElectrons_v1[i]);
              v_genLep2.push_back(genEle1);
            }

        }
  //   }
  // else
  //   {
      for(int i=0 ; i < GenMuons_; i++)
        {
          if(GenMuons_v1[i].Pt()!=0)
            {
              genMu1 = (GenMuons_v1[i]);
              v_genLep2.push_back(genMu1);
            }
        }
      //  }
  return MinDr(bestPhoton,v_genLep2);
}

double AnalyzeLightBSM::getGenRecodRLep(TLorentzVector recoLep1){//MC only                                                                                                
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
      }
    }
  }//for                                                                                                                                          
  if(recoLep1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(recoLep1);
  else return 1000.0;
}

std::vector<int> AnalyzeLightBSM::dR_recoPho_GenParticle(TLorentzVector recoPho){//MC only                                                                     
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100, pdgID =-99999, parentID=-99999, count=0;
  int flag = 0;
  std::vector<int> array;//={};
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      //if(abs((*GenParticles_PdgId)[i])==22 || abs((*GenParticles_ParentId)[i])==22) continue;
      if(recoPho.DeltaR(GenParticles_v1[i])<0.2 )
	{
	  //cout<<"You got to be kidding me"<<endl;	  
	  //h_parID_matchRecoPho->Fill((*GenParticles_PdgId)[i],wt);
	  //h_parentID_vsPartId_matchPho->Fill((*GenParticles_PdgId)[i],(*GenParticles_ParentId)[i],wt);
	  if(abs((*GenParticles_PdgId)[i])==22)
	    {
	      //h_phopT_BeamRemenant->Fill(recoPho.Pt(),wt);
	      continue;
	    }
	  if(abs((*GenParticles_ParentId)[i])==22) continue;	  
	  pdgID =(*GenParticles_PdgId)[i];
	  count++;
	  flag = 1;	  
	  //cout<<"You got to be kidding me"<<"\t"<<flag<<"\t"<<pdgID<<"\t"<<(*GenParticles_ParentId)[i]<<"\t"<<(*GenParticles_Status)[i]<<endl;
	}
    }
  }
  array.push_back(pdgID);
  array.push_back(flag);
  //cout<<"after kidding me"<<"\t"<<array[0]<<"\t"<<array[1]<<endl;
  return array;//pdgID, flag;
 
}
vector <TLorentzVector> AnalyzeLightBSM::getLorentzVector(int size, Float_t Pt_size[],Float_t Eta_size[],Float_t Phi_size[],Float_t E_size[])
{
  //  cout<<"filling the vector of lorentz vector"<< "\t"<< "vector <TLorentzVector> AnalyzeLightBSM::getLorentzVector"<<endl;
  vector<TLorentzVector> Temp;
  for(int i=0; i<size;i++)
    {
      TLorentzVector P1;
      P1.SetPtEtaPhiE(Pt_size[i],Eta_size[i],Phi_size[i],E_size[i]);
      //cout<< "Pt "<<Pt_size[i]<<" Eta "<<Eta_size[i]<< " Phi "<<Phi_size[i]<< "E "<<E_size[i]<<endl;
      Temp.push_back(P1);
    }
  //cout<< "Out lorentz vectors size "<<Temp.size() <<endl;
  return Temp;
}
void AnalyzeLightBSM::FillHistogram_Kinematics(int i, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;
  h_Njets[i]->Fill(Njets,wt);
  h_Nbjets[i]->Fill(btags,wt);
  h_MET_[i]->Fill(MET,wt);
  h_PhotonPt[i]->Fill(pho_Pt,wt);
  h_Mt_PhoMET[i]->Fill(mt_phoMET,wt);
  h_dPhi_PhoMET[i]->Fill(dPhi,wt);
  h_St[i]->Fill(ST,wt);
  h_HT[i]->Fill(HT,wt);
  int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  h_TFbins_LL_v4[i]->Fill(TFbins_v4,wt);
  h_TFbins_LL_v2[i]->Fill(TFbins_v2,wt);
  
  h_Sbins_LL[i]->Fill(searchBin,wt);
  h_TFbins_LL_v3[i]->Fill(TFbins_v3,wt);

  //cout<<"Alps "<<i<<"\t"<<"TFbins_v1  "<<TFbins_v1<<" TFbins_v2 "<< TFbins_v2<<"\t MET "<<MET<<"\t NJets  "<<Njets<<"  btags "<<btags<<"\t"<<h_TFbins_LL_v1[i]->GetBinContent(TFbins_v1)<<"\t"<<wt<<"\t"<<searchBin<<endl;
  //cout<<h_Njets[i]->GetMean()<<endl;

}
void AnalyzeLightBSM::FillTFBins_Valid(int i, int Njets, int btags, double wt,double wt1, double wt2, double pho_Pt, double ST){
  int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
  h_TFbins_ElecLL_validation[i]->Fill(TFbins_v2,wt);
  h_TFbins_ElecLL_validation_v1[i]->Fill(TFbins_v2,wt);
  h_Sbins_LL_Validation[i]->Fill(searchBin,wt);
  h_Njets_validation[i]->Fill(Njets,wt);
  h_Nbjets_validation[i]->Fill(btags,wt);
  h_MET_validation[i]->Fill(MET,wt);
  h_PhotonPt_validation[i]->Fill(pho_Pt,wt);
  h_St_validation[i]->Fill(ST,wt);

  h_TFbins_ElecLL_validation_TFbins_v2[i]->Fill(TFbins_v3,wt1);
  h_TFbins_ElecLL_validation_TFbins_v2_v1[i]->Fill(TFbins_v3,wt1);
  h_Sbins_LL_Validation_TFbins_V2[i]->Fill(searchBin,wt1);
  h_Njets_validation_TFbins_v2[i]->Fill(Njets,wt1);
  h_Nbjets_validation_TFbins_v2[i]->Fill(btags,wt1);
  h_MET_validation_TFbins_v2[i]->Fill(MET,wt1);
  h_PhotonPt_validation_TFbins_v2[i]->Fill(pho_Pt,wt1);
  h_St_validation_TFbins_v2[i]->Fill(ST,wt1);

  h_TFbins_ElecLL_validation_TFbins_v3[i]->Fill(TFbins_v4,wt2);
  h_TFbins_ElecLL_validation_TFbins_v3_v1[i]->Fill(TFbins_v4,wt2);
  h_Sbins_LL_Validation_TFbins_V3[i]->Fill(searchBin,wt2);
  h_Njets_validation_TFbins_v3[i]->Fill(Njets,wt2);
  h_Nbjets_validation_TFbins_v3[i]->Fill(btags,wt2);
  h_MET_validation_TFbins_v3[i]->Fill(MET,wt2);
  h_PhotonPt_validation_TFbins_v3[i]->Fill(pho_Pt,wt2);
  h_St_validation_TFbins_v3[i]->Fill(ST,wt2);

  
}
void AnalyzeLightBSM::FillHistogram_Kinematics_varBin(int i, int Njets, int btags, double pho_Pt, double ST, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;                                        
  h_Njets_Varbin[i]->Fill(Njets,wt);
  h_Nbjets_Varbin[i]->Fill(btags,wt);
  h_MET_Varbin[i]->Fill(MET,wt);
  h_PhotonPt_Varbin[i]->Fill(pho_Pt,wt);
  h_St_Varbin[i]->Fill(ST,wt);
  //h_HT_Varbin[i]->Fill(HT,wt);
  //cout<<h_Njets[i]->GetMean()<<endl;                                                                                                            

}
int AnalyzeLightBSM::Photons_OriginType(){
  int flag_photon = -1;
  if(!(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons) && !(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons) && hasGenPromptPhoton)
    flag_photon=0;
  else if( (*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==1)
    flag_photon=1;
  else if ((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons))
    flag_photon=2;
  else
    flag_photon = -1;
  return flag_photon;
}
// double AnalyzeLightBSM::getMVArespone(float MET, float NJets, float BTags, float HT)
// {
//   TMVA::Tools::Instance();

//   TMVA::Reader *reader1 = new TMVA::Reader("Color:Silent");
//   float met=0.0,njets=0.0,btags=0.0,ht=0.0;
//   reader1->AddVariable( "MET", &met );
//   reader1->AddVariable( "NJets", &njets );
//   reader1->AddVariable( "BTags", &btags );
//   reader1->AddVariable( "HT", &ht );
//   reader1->BookMVA( "BDT_200trees_2maxdepth method", "./TMVAClassification_BDT_200trees_2maxdepth.weights.xml" );
//   met = MET;
//   njets= NJets;
//   btags = BTags;
//   ht = HT;
//   Double_t mvaValue = reader1->EvaluateMVA( "BDT_200trees_2maxdepth method");
//   return mvaValue;
//   delete reader1;
// }

#define ANALYZETPROXYTBSM_cxx

#include "AnalyzeTProxytBSM.h"

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
//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;

//#ifdef __MAKECINT__
//#pragma link C++ class NtupleVarsTProxy+;
//#endif

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

  AnalyzeTProxytBSM ana(inputFileList, outFileName, data,sample, elec,phoID);
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  //ana.EventLoop(data,inputFileList,sample,outFileName,elec,phoID);
  Long64_t ievent=10;
  //ana.Process(ievent);
  ana.EventLoop(data,sample);
  //  ana.EventLoop(const char *,const char *);
  Tools::Instance();
  return 0;
}

//void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
void AnalyzeTProxytBSM::EventLoop(const char *data, const char *sample) {

  std::cout << "AnalyzeTProxytBSM::EventLoop() " << std::endl;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  //cout << "Analyzing dataset " << data << " " << endl;

  TString s_sample= sample;
  TString s_data=data;
    fChain->SetBranchStatus("*DeltaPhi*",0);
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  TTree* outtree = fChain->CloneTree(0);
  bool Debug=false;
  double wt=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<10;jentry++) {

    // // ==============print number of events done == == == == == == == =                                                                    
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;

    // ===============read this entry == == == == == == == == == == ==                                                                     
    //fDirector.SetReadEntry(jentry);
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;

   if(jentry<10 ) {
   std::cout<< "jentry " << jentry << " RunNum " << RunNum << std::endl;
   //std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
   // for(Long64_t ii=0; ii<GenParticles->size(); ii++){
   //   ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = GenParticles[(int)ii];
   //   std::cout <<" ii, Pt, Eta, Phi, E " << ii << GenParticles[(int)ii].Pt() << " " << GenParticles[(int)ii].Eta() << " " << GenParticles[(int)ii].Phi() << " " << GenParticles[(int)ii].E() << " pdgid, parentid, status " << GenParticles_PdgId[(int)ii] << " " << GenParticles_ParentId[(int)ii] << " " << GenParticles_Status[(int)ii]  << std::endl;
   // }
   
   std::cout << std::endl; 
   std::cout << "Electrons->size() "<< Electrons->size() << std::endl;
   for(Long64_t ii=0; ii<Electrons->size(); ii++){
     ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = Electrons[(int)ii];
     std::cout <<" ii, Pt, Eta, Phi, E " << ii << Electrons[(int)ii].Pt() 
	       << " " << Electrons[(int)ii].Eta() << " " << Electrons[(int)ii].Phi() 
	       << " " << Electrons[(int)ii].E()  
	       << " iso, mediumID " << Electrons_iso[(int)ii] << " " << Electrons_mediumID[(int)ii]
	       << std::endl;
   }

   std::cout << std::endl; 
   std::cout << "Photons->size() "<< Photons->size() << std::endl;
   for(Long64_t ii=0; ii<Photons->size(); ii++){
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mypho = Photons[(int)ii];
     std::cout <<" ii, Pt, Eta, Phi, E " << ii << Photons[(int)ii].Pt() 
	       << " " << Photons[(int)ii].Eta() << " " << Photons[(int)ii].Phi() 
	       << " " << Photons[(int)ii].E()  
	       << " mvavalueID, pfGammaIso " << Photons_mvaValuesID[(int)ii] << " " << Photons_pfGammaIso[(int)ii]
	       << std::endl;
   }

   std::cout << "================================" << std::endl;
   }
   
   wt = Weight*1000.0;
   
   h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
   myLV goodPho = getBestPhoton(0);
   if(jentry<10 ) {
     std::cout << "Good photon " << goodPho.Pt() << " " << goodPho.Eta()  << " " << goodPho.Phi() << std::endl;
   }
   // *******************  Selecting Jet objects ********************************//                                                                                       
   int minDRindx=-100,phoMatchingJetIndx=-100,hadJetID=-999,bJet1Idx=-100,nHadJets=0;
   double minDR=99999,ST=0,Ht=0;
   vector<myLV> hadJets,bjets;
   hadJets.clear();
   bjets.clear();
   TLorentzVector bestPhoton, Jets_v1;
   vector<int> jetMatchindx;
   bestPhoton.SetPtEtaPhiE(goodPho.Pt(),goodPho.Eta(),goodPho.Phi(),goodPho.E());
   bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
   if(Debug)
     cout<<"===load tree entry check2 at entry ==="<<"\t"<<jentry<<endl;
   for(int i=0;i<Jets->size();i++)
     {
       if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	 Jets_v1.SetPtEtaPhiE(Jets[i].Pt(),Jets[i].Eta(),Jets[i].Phi(),Jets[i].E());
	 double dR = bestPhoton.DeltaR(Jets_v1);//DeltaR(goodPho.Eta(),goodPho.Phi(),Jets[i].Eta(),Jets[i].Phi());//goodPho.DeltaR(Jets[i]);
	 if(dR<minDR){minDR=dR;minDRindx=i;}
       }
     }
   if(Debug)
     cout<<"===load tree entry  ==="<<"\t"<<jentry<<"\t"<<"Jets check == "<<minDR<<endl;

   for(int i=0;i<Jets->size();i++)
     { if(Debug)
	 cout<<"  = Jets.Pt()  ==  "<<Jets[i].Pt()<<"\t"<< " = Jets.Eta() == "<<Jets[i].Eta()<<endl;
       if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	 if(Debug)
	   cout<< "==== loadjets ==="<<"\t"<<i<<"\t"<<minDR<<endl;
	 if( !(minDR < 0.3 && i==minDRindx) )
	   {
	     hadJetID= (*Jets_ID)[i];
	     //if(true){                                                                                                                                                    
	     hadJets.push_back(Jets[i]);
	   }
       }
     }

   storeBTagEff(data);

   //if(hadJets.size()==0) continue;

   // ===== Seems to work =======
   // Copies entries satisfying criteria to a new tree
   // bool tighte_trgpass = true;
   // if(s_sample.Contains("data") )//&& applyTrgEff)
   //   {
      
   //     if(TriggerPass->size()!=213) continue;
   //     if((*TriggerPass)[128]==1 || (*TriggerPass)[129]==1 || (*TriggerPass)[131]==1 || (*TriggerPass)[133]==1 || (*TriggerPass)[137]==1 || (*TriggerPass)[139]==1\
   // 	  ) tighte_trgpass=true;
   //     else tighte_trgpass=false;
   //     if(tighte_trgpass==false)  continue;
   //   }
   if(goodPho.Pt()>20) h_selectBaselineYields_->Fill("Good #gamma with Pt > 20",wt);
   else continue;
   if(MET>100) h_selectBaselineYields_->Fill("MET > 100",wt);
   else continue;
   if(hadJets.size()>=2)
     h_selectBaselineYields_->Fill("Good nHadJets >= 2",wt);
   else continue;
   
   if(MET>100 && goodPho.Pt()>20 && hadJets.size()>=2){
     if(jentry < 30)
       std::cout << "Good photon " << goodPho.Pt() << " " << goodPho.Eta()  << " " << goodPho.Phi() << std::endl;
     outtree->Fill();
     if(jentry < 30)
       std::cout << "Good photon " << goodPho.Pt() << " " << goodPho.Eta()  << " " << goodPho.Phi() << std::endl;
     if(jentry<100)
       std::cout<< "Good nJets "<<hadJets.size()<<" Jets size "<<Jets->size()<<std::endl;
   }

  } // loop over entries                                                                                                                   
  // ===== this seemed to have worked =====
  // === Clones a new tree from and old tree
  //const char* sel = "(MET>200)";
  //oFile->cd();
  //outtree = fChain->CopyTree( sel );
  // =========================
}


myLV AnalyzeTProxytBSM::getBestPhoton(int pho_ID){
  //vector<TLorentzVector> goodPho;
  vector<myLV> goodPho;
  vector<int> goodPhoIndx;
  for(int iPho=0;iPho<Photons->size();iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    if(((*Photons_hasPixelSeed)[iPho]<0.001) )//&& ( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001) &&( pho_ID==0 || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42))) ) 
      {
	goodPho.push_back(Photons[iPho] );
	goodPhoIndx.push_back(iPho);
      }
  }
  
  int highPtIndx=-100;
   for(int i=0;i<goodPho.size();i++){
     if(i==0) highPtIndx=0;
     else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
   }
   
   if(highPtIndx>=0){
     bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
   }
   else bestPhotonIndxAmongPhotons = -100;
   if(highPtIndx==-100){myLV v0;return v0;}
   else return goodPho[highPtIndx];
   
}




  
Bool_t AnalyzeTProxytBSM::Process(Long64_t entry) {

  std::cout << entry << std::endl;
   fDirector.SetReadEntry(entry);
   std::cout<< "entry " << entry << " RunNum " << RunNum << std::endl;
   std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
  return 0;
}

void AnalyzeTProxytBSM::storeBTagEff(const char *data){
  // DeepCSV values:                                                                                                                                                 
  //   2016: 0.6321  2017: 0.4941  2018: 0.4184                                                                                                                \
                                                                                                                                                                     
  TString s_data = data;
  double deepCSVvalue=0;
  double lumiInfb=0, p0=0,p1=0,p2=0; 
  if(s_data.Contains("2016preVFP")){ lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV                                             
  if(s_data.Contains("2016postVFP")) { lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016                                         
  if(s_data.Contains("2017")) {lumiInfb=41.48; deepCSVvalue = 0.4506;}
  if(s_data.Contains("2018")){ lumiInfb=59.83;deepCSVvalue = 0.4168;}
  // if(s_data.Contains("signal"))lumiInfb= 137.19;

  // if(s_data.Contains("Summer16")) CSVv2WP= 0.6321;
  // if(s_data.Contains("Fall17")) CSVv2WP= 0.4941;
  // if(s_data.Contains("Autumn18")) CSVv2WP= 0.4184;
  for(unsigned ja = 0; ja < Jets->size(); ++ja){
    if(!Jets_HTMask[ja]) continue;
    int flav = abs(Jets_hadronFlavor[ja]);
    double csv = Jets_bJetTagDeepCSVBvsAll[ja]; //Jets_bJetTagDeepCSVBvsAll
    double pt = Jets[ja].Pt();
    //use abs(eta) for now                                                                                                                                           
    double eta = fabs(Jets[ja].Eta());
    if(flav==5){
      d_eff_b->Fill(pt,eta);
      if(csv > deepCSVvalue) n_eff_b->Fill(pt,eta);
    }
    else if(flav==4){
      d_eff_c->Fill(pt,eta);
      if(csv > deepCSVvalue) n_eff_c->Fill(pt,eta);
    }
    else if(flav<4 || flav==21){
      d_eff_udsg->Fill(pt,eta);
      if(csv > deepCSVvalue) n_eff_udsg->Fill(pt,eta);
    }
  }
}


//void AnalyzeTProxytBSM::Process(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
//
//
//}

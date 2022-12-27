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
using namespace TMVA;
int main(int argc, char* argv[])
{

  if (argc < 5) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" "<<"Process"<< " "<<"LostlepFlag"<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  const char *elec = argv[5];
  AnalyzeLightBSM ana(inputFileList, outFileName, data,sample, elec);
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  ana.EventLoop(data,inputFileList,sample,outFileName,elec);
  Tools::Instance();

  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec) {
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  TString Elec_flag = elec;
  bool lost_elec_flag=false;
  if(Elec_flag.Contains("Electron"))
    lost_elec_flag = true;
  else
    lost_elec_flag = false;
  
  TString s_data=data;
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
  bool higMET=false;
  if(s_data.Contains("2016")) lumiInfb=35.922;
  if(s_data.Contains("2017")) lumiInfb=41.529;
  if(s_data.Contains("2018")) lumiInfb=59.74;
  if(s_data.Contains("signal"))lumiInfb= 137.19; //since we only have 2018 dataset
  std::string s_process = sample;
  double cross_section = getCrossSection(s_process);
  cout<<cross_section<<"\t"<<"analyzed process"<<endl;
  cout<<"Event"<<"\t"<<"par-ID "<<"\t"<<"parentID"<<"\t"<<"GenMET"<<"\t"<<"MET"<<"\t"<<"pT"<<"\t"<<"Eta"<<"\t"<<"Phi"<<"\t"<<"E"<<endl;
  float met=0.0,st=0.0, njets=0, btags=0,mTPhoMET=0.0,dPhi_PhoMET=0.0,dPhi_MetJet=0.0;
  
  TMVA::Reader *reader1 = new TMVA::Reader();
  reader1->AddVariable( "MET", &met );
  reader1->AddVariable( "NhadJets", &njets );
  reader1->AddVariable( "BTags", &btags );
  reader1->AddVariable("mTPhoMET_",&mTPhoMET);
  reader1->AddVariable("dPhi_PhoMET_",&dPhi_PhoMET);
  reader1->AddVariable("dPhi_Met_Jet",&dPhi_MetJet);
  reader1->AddVariable( "ST", &st );
  reader1->BookMVA( "BDT_100trees_2maxdepth method", "/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/Susy_Analysis_2020/BDT_training/dataset_bdt_all_Equalweight_pMSSM_allvariable_NominalBasleine_100trees_2maxdepth/weights/TMVAClassification_BDT_100trees_2maxdepth.weights.xml" );
  float nsurVived=0.0;
  int searchBin=0, Tfbins=0;
  float nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0,wt_LL=0.0;

// counters for events yields after each selection/rejection
 float nEvents_Selec[100]={};
 const char* out_nEventsTags[100] ={};
 //filetag to read the TF histograms from the input root file
 const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
 //reading TF for muons
 TFile* f_muon= new TFile("TF_allin1_LLEstimation_muon_newTFbins.root");

 //reading TF for electrons
 TFile* f_elec =new TFile("TF_allin1_LLEstimation_electron_newTFbins.root");
 char* histname = new char[2000];
 TH1F* h_TF;
 if(lost_elec_flag)
   {
     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_Elec__%s_%s",sample,data);
     h_TF = (TH1F*)f_elec->Get(histname);
     cout<<"Reading TF for lost electron estimation:  "<<"\t"<<histname<<endl;
     for(int i=0; i<h_TF->GetNbinsX();i++)
       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}

   }
 else
   {
     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_FailAcep_TauHadronic_%s_%s",sample,data);
     cout<<"Reading TF for lost muon estimation:  "<<"\t"<<histname<<endl;
     h_TF = (TH1F*)f_muon->Get(histname);
     for(int i=0; i<h_TF->GetNbinsX();i++)
       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}
   }
 //  nentries=100000;
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
	cout<<"===load tree entry ==="<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // ==========================================================================================\\
      // ====================================== Calculating weight per event ======================\\
      // ==========================================================================================\\

      if(s_sample.Contains("MCMC_86_7257_1750"))
        wt = (0.002991*lumiInfb*1000)/nentries; //gluino mass = 1750 GeV                                                                                                      
      else if(s_sample.Contains("MCMC_86_7257_1400"))
        wt = (0.0284*lumiInfb*1000)/nentries; //need to check luminosity                                                                                                      
      else if(s_sample.Contains("MCMC_86_7257_1500"))
        wt = (0.0157*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1600"))
	wt = (0.00887*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1700"))
        wt = (0.00507*lumiInfb*1000)/nentries;
      else  if(s_sample.Contains("MCMC_86_7257_1800"))
        wt = (0.00293*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1900"))
        wt = (0.00171*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_2000"))
        wt = (0.00101*lumiInfb*1000)/nentries;
      else if(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))
	{
	  std::string s_process = sample;
	  double cross_section = getCrossSection(s_process);
	  //cout<<cross_section<<"\t"<<"analyzed process"<<endl;
	  wt = Weight*lumiInfb*1000.0; //(cross_section*lumiInfb*1000.0)/nentries;
	}
      //      cout<<"===load tree entry ==="<<"\t"<<jentry<<endl;

      h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
      nocut++;
      int count=0;
      int branch_size = (*GenParticles).size();
      int ele_branch=(*Electrons).size();
      int pho_branch = (*Photons).size();
      //variables to be used in getting dR
      float iGenEle_eta =99999.0, iGenEle_phi=99999.0, iRecEle_eta =99999.0, iRecEle_phi=99999.0, iRecphot_eta =99999.0, iRecphot_phi=99999.0,iGen_Wpt=99999.0,iGen_Wp=99999.0,iGen_WEta=99999.0,iGen_WPhi=99999.0;
      float dR_Ele=0, min_dR=9999, min_dR_pho=9999;
      int pdgID=0, parentId=0, status;
      int Ids[4]={11,13,15,22};
      const char* ids[4]={"e","mu","tau","pho"};
      // loop over gen particles
      int count_genElec=0,count_genMu=0,count_genTau=0, igenPho=0, count_genEl_tau=0,count_genMu_tau=0,count_haddecay=0;
      TLorentzVector genPho_W,genElec_W,genMu_W,genTau_W,genElec_Tau,genMu_Tau,genW,genNuElec_W,genNuMu_W,genNuTau_W,genElecNu_Tau,genMuNu_Tau;
      int elec_reco=0,elec_reco0_before=0,elec_reco1_before=0,muon_reco=0,elec_gen3=0,elec_gen2=0, elec_gen=0, muon_gen=0,elec_reco0=0,elec_reco1=0,evtSurvived_preselec=0,elec_reco2=0,elec2_reco=0,survived_vetohad=0,elec_reco1_CR=0,survived_elecge1=0,events_cr=0,events_sr=0,total=0,remain=0,elec_reco0_genel=0,efakepho=0,ele=0,genphomatch_after=0,genphomatch_before=0,elec_gen4=0,gentauhad2=0,lep2=0,lep=0;
      int fail_realphoton=0;
      int pass_realphoton=0;
      int fail_acceptance=0;
      int pass_acceptance=0;
      int fail_id=0;
      int pass_id=0;
      int fail_iso=0;
      int pass_iso=0;      
      int total_lost_el = 0,cr_el=0,sr_el,e_index=-1,nlep=0, NgenElec=0,leadGenPhoIdx1=0;
      double mt_ele=0,mt_pho=0,mt_ele1=0;
      TLorentzVector leadGenPho,genPho,genLep,v_lep1,v_lep2;//,v_genLep1,v_genLep2,v_genMu1,v_genMu2,v_genElec1,v_genElec2,v_genTau1,v_genTau2;
      vector<TLorentzVector> v_genPho, v_genLep,v_genElec,v_genMu,v_genTau, v_genAccepLep,v_genAccepElec,v_genAccepMu,v_genAccepTau;
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////   Making gen level collections of all particles  ////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // for (int igen=0; igen<branch_size; igen++)
      // 	{	  
      // 	  pdgID= (*GenParticles_PdgId)[igen];
      // 	  parentId=(*GenParticles_ParentId)[igen];
      // 	  // if(abs(parentId)==24)
      // 	  // cout<<"entry: "<<jentry<<" "<<pdgID<<"\t"<<parentId<<"\t"<<endl;
      // 	  if((*GenParticles)[igen].Pt()!=0){
      // 	    if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
      // 	      if(genPho.Pt() < (*GenParticles)[igen].Pt()){
      // 		leadGenPhoIdx1 = igen;
      // 		leadGenPho = ((*GenParticles)[igen]);	      
      // 	      }
      // 	      genPho = ((*GenParticles)[igen]);
      // 	      v_genPho.push_back(genPho);	      
      // 	    }	  
      // 	}
      // 	}
      
      //      h_2d_genPromptPho_VsGenPho->Fill(hasGenPromptPhoton,v_genPho.size(),wt);
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////     Identifying the decay modes for W & TTbar  //////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int igen=0; igen<branch_size; igen++)
	{
	  pdgID= (*GenParticles_PdgId)[igen];
	  parentId=(*GenParticles_ParentId)[igen];
	  if(abs(parentId)==24) // || parentId==-24)
	    {
	      if(abs(pdgID)==11) { h_GenpT[0]->Fill((*GenParticles)[igen].Pt()); h_GenEta[0]->Fill((*GenParticles)[igen].Eta()); count_genElec++;}
	      else if(abs(pdgID)==13) { h_GenpT[1]->Fill((*GenParticles)[igen].Pt()); h_GenEta[1]->Fill((*GenParticles)[igen].Eta()); count_genMu++;}
	      else if(abs(pdgID)==15){h_GenpT[2]->Fill((*GenParticles)[igen].Pt()); h_GenEta[2]->Fill((*GenParticles)[igen].Eta()); count_genTau++;}
	      else if (abs(pdgID)==1 || abs(pdgID)==2 ||abs(pdgID)==3 ||abs(pdgID)==4 || abs(pdgID)==5) count_haddecay++;
	
	      if(abs(pdgID)==11 || abs(pdgID)==13 ||abs(pdgID)==15)		
		  count++;//t<<"lepton"<<"\t"<<jentry<<endl;
	      else if (abs(pdgID)==12 ||abs(pdgID)==14||abs(pdgID)==16)
		count++;//t<<"neutrino"<<"\t"<<jentry<<endl;
	      else
		{
		  igenPho++;
		  h_counts->Fill("pho",wt);
		  //		  cout<<"pdgID"<<"\t"<<pdgID<<"\t"<<jentry<<endl;
		}	      
	      for (int i=0; i<3;i++)
		{
		  if(abs(pdgID)==Ids[i])
		    h_counts->Fill(ids[i],wt);
		}
	    }
	}
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////            checking whether tau decayed to muon/electron or hadronically   ///////////////////////////////
      //////////////////             Not fully corrected                                             //////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      int count_had_tau=0;
      if((*GenTaus).size()!=0 && GenTaus_had)
	{
      	  if((*GenElectrons).size()!=0)	    
	    {
	      genElec_Tau.SetPtEtaPhiE((*GenElectrons)[0].Pt(),(*GenElectrons)[0].Eta(),(*GenElectrons)[0].Phi(),(*GenElectrons)[0].E());
	      count_genEl_tau++;
	    }
	  else if((*GenMuons).size()!=0)
	    {
	      genMu_Tau.SetPtEtaPhiE((*GenMuons)[0].Pt(),(*GenMuons)[0].Eta(),(*GenMuons)[0].Phi(),(*GenMuons)[0].E());
	      count_genMu_tau++;
	    }
	}
      else if ((*GenTaus).size()!=0 &&!GenTaus_had)
	count_had_tau++;
      if(Debug)
      cout<<"===load tree entry check1 ==="<<"\t"<<jentry<<endl;           
      int count_genEle=0,count_recEle=0;
      // TLorentzVector genPho1,genLep1;
      // int leadGenPhoIdx=-100;
      vector<TLorentzVector> goodPho_n;
      vector<int> goodPhoIndx_n;
      for(int iPho=0;iPho<Photons->size();iPho++){
	if( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
	  goodPho_n.push_back( (*Photons)[iPho] );
	  goodPhoIndx_n.push_back(iPho);
	}
      }

      int nPhotons = goodPho_n.size();
       TLorentzVector bestPhoton=getBestPhoton();
      
      // *******************  Selecting Jet objects ********************************//
      int minDRindx=-100,phoMatchingJetIndx=-100,nHadJets=0;
      double minDR=99999,ST=0,Ht=0;
      vector<TLorentzVector> hadJets;
      vector<int> jetMatchindx;
      bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
      for(int i=0;i<Jets->size();i++)
	{
	  if( ((*Jets)[i].Pt() > 30.0) && (abs((*Jets)[i].Eta()) <= 2.4) ){
	    double dR=bestPhoton.DeltaR((*Jets)[i]);
	    if(dR<minDR){minDR=dR;minDRindx=i;}
	  }
      }
    for(int i=0;i<Jets->size();i++)
      {
	if( ((*Jets)[i].Pt() > 30.0) && (abs((*Jets)[i].Eta()) <= 2.4) ){
          if( !(minDR < 0.3 && i==minDRindx) )
	    {
	      hadJets.push_back((*Jets)[i]);
	      jetMatchindx.push_back(i);
	    }
	  
	}
      }
    if( minDR<0.3 ) {phoMatchingJetIndx=minDRindx; recoJetMatch_recoPho=true;}
    double genmindr=99999, recojetmindr=99999;
    genmindr = MinDr(bestPhoton,(*GenJets));
    recojetmindr =  minDR;// MinDr(bestPhoton,hadJets);
    if(genmindr<0.3)
      genJetMatch_recoPho = true;
    //    if(hadJets.size() == 0) continue;
    
    for(int i=0;i<hadJets.size();i++){
      if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());Ht=Ht+(hadJets[i].Pt());}
      if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
    }
    sortTLorVec(&hadJets);
    // ********* This is to account into all visible energy: Adding photon matching with Jet ******************//
    if( minDR<0.3 ) ST=ST+bestPhoton.Pt(); 
    //transverse mass of best photon+MET                                                                                  
    double mTPhoMET=sqrt(2*(bestPhoton.Pt())*MET*(1-cos(DeltaPhi(METPhi,bestPhoton.Phi()))));
    // ******* Dphi between MET and Best Photon     ******************** //                                                                 
    double dPhi_PhoMET= abs(DeltaPhi(METPhi,bestPhoton.Phi()));
    h_Njets[0]->Fill(nHadJets,wt);
    h_Nbjets[0]->Fill(BTags,wt);
    h_MET_[0]->Fill(MET,wt);
    h_PhotonPt[0]->Fill(bestPhoton.Pt(),wt);
    h_Mt_PhoMET[0]->Fill(mTPhoMET,wt);
    h_dPhi_PhoMET[0]->Fill(dPhi_PhoMET,wt);
    h_St[0]->Fill(ST,wt);
    h_HT[0]->Fill(Ht,wt);
    h_njets_vs_ST[0]->Fill(nHadJets,ST,wt);
    h_njets_vs_HT[0]->Fill(nHadJets,Ht,wt);
    h_ST_vs_ptPho[0]->Fill(ST,bestPhoton.Pt(),wt);
    if(Debug)
      cout<<"===load tree entry variable define ==="<<"\t"<<jentry<<endl;

    //    h_mvaResponse_baseline[0]->Fill(mvaValue,wt);
    //if (nHadJets<2) cout<<"wrong event"<<endl;
    nEvents_Selec[0]+=wt;
    out_nEventsTags[0] = "No cut (except skimming for backgrounds)";
    if(bestPhotonIndxAmongPhotons<0) continue;
    nEvents_Selec[1]+=wt;
    out_nEventsTags[1] = "making sure 1 #gamma";
    bool bestPhoHasPxlSeed=true;
    if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
    if( bestPhoHasPxlSeed ) continue;
    nEvents_Selec[2]+=wt;
    out_nEventsTags[2] = "#gamma with pixel seed veto";
    //    if(bestPhoton.Pt()>30) 
    if(bestPhoton.Pt()>40)                                                                                                                       
      h_selectBaselineYields_->Fill("pt>30",wt);                                                                                                
    else continue;
    out_nEventsTags[3]="#gamma pT>30";
    nEvents_Selec[3]+=wt;
    // ///////////////////  Best photon - to be prompt photon  //////////////////////////////////////
    // if(!s_data.Contains("data"))
    //   {
    // 	double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);
    // 	if(mindR_genPho>0.2)
    // 	  continue;
    // 	else
    // 	  {
    // 	    double mindrPho_genlep=getGenLep(bestPhoton);
    // 	    if(mindrPho_genlep<0.2) continue;
    // 	    h_selectBaselineYields_->Fill("prompt Photon Checks",wt);
    // 	  }	
    //   }

    h_selectBaselineYields_->Fill("photon_selec",wt);
    if (k > decade)
      cout << 10 * k << " %-3" << endl;
    if (nHadJets>=2) 
    h_selectBaselineYields_->Fill("njets>=2",wt);
    else continue;
    nEvents_Selec[4]+=wt;
    out_nEventsTags[4]="nHadJets>=2";
    if(MET>250) 
    h_selectBaselineYields_->Fill("MET>100",wt);
    else continue;
    out_nEventsTags[5]="MET>100";
    nEvents_Selec[5]+=wt;

    if(ST>300)
      h_selectBaselineYields_->Fill("ST>300",wt);
    else continue;
    out_nEventsTags[6]="ST>300";
    nEvents_Selec[6]+=wt;
    
    if(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))
      {
        if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 &&BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0) ) continue;
      }
    h_selectBaselineYields_->Fill("MetCleaning",wt);
    nEvents_Selec[8]+=wt;
    out_nEventsTags[8]="MetCleaning filters";
    if (k > decade)
      cout << 10 * k << " %-4" << endl;
    TLorentzVector Met;
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    Double_t deta_jet_pho= 0.0,deta_jet_met=0.0,deta_met_pho=0.0;
    deta_jet_pho = abs(bestPhoton.Eta()-hadJets[0].Eta());
    double mT= 0.0, dPhi_METjet1=5, dPhi_METjet2=5, dPhi_phojet1=5, dPhi_phojet2=5, dPhi_phoMET=5;
    if(hadJets.size() > 0) dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    if(hadJets.size() > 1) dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));

    dPhi_METjet1 = abs(Met.DeltaPhi(hadJets[0]));
    dPhi_METjet2 = abs(Met.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    dPhi_METjet1=3.8,dPhi_METjet2=3.8;//dphi3=3.8,dphi4=3.8,dphiG_MET=3.8;                                                                                                                                         
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    dPhi_METjet1 = abs(DeltaPhi(METPhi,(hadJets)[0].Phi()));
    dPhi_METjet2 = abs(DeltaPhi(METPhi,(hadJets)[1].Phi()));
    
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
        h_selectBaselineYields_->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    nEvents_Selec[9]+=wt;
    out_nEventsTags[9]="dPhi1 & dPhi2 >= 0.3";
    // ********************  photon jet matching : pT comparison *************************//                                                                                                                       
    if(phoMatchingJetIndx>=0 && ((*Jets)[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
    nEvents_Selec[10]+=wt;
    out_nEventsTags[10]="pT(jet)/pT #gamma>1.0";
    if(phoMatchingJetIndx<0) continue;
    nEvents_Selec[11]+=wt;
    out_nEventsTags[11]="phoMatchingJetindex";
    h_selectBaselineYields_->Fill("additional Bhumika",wt);
    if (bestPhoton.Pt()>40)
      h_selectBaselineYields_->Fill("pho-pT>30",wt);
    else continue;
    if(higMET)
      { 
	if (bestPhoton.Pt()<100)continue;
	if(MET<250) continue;                                                                                              
      }
    h_selectBaselineYields_v1->Fill("Pre-Selection",wt);
    nEvents_Selec[9]+=wt;    
    h_madHT->Fill(madHT,wt);
    h_madminPhotonDeltaR_->Fill(madMinPhotonDeltaR,wt);
    //cout<<hasGenPromptPhoton<<endl;
    h_hasGenPromptPhoton->Fill(hasGenPromptPhoton,wt);
    double mindr_Pho_genlep=getGenLep(bestPhoton);
    //h_phoPt_promptPho->Fill(
    // photon-pT distribution checks
    if(hasGenPromptPhoton)
      {
	h_phoPt_promptPho->Fill(bestPhoton.Pt(),wt);
	if(madMinPhotonDeltaR<0.5)
	  h_phoPt_promptPho_mindr_qgl0p5->Fill(bestPhoton.Pt(),wt); 
	else
	  h_phoPt_promptPho_mindr_qgh0p5->Fill(bestPhoton.Pt(),wt);
	double mindr_Pho_genlep=getGenLep(bestPhoton);
	if(mindr_Pho_genlep<0.5)
	  h_phoPt_promptPho_mindr_lepl0p5->Fill(bestPhoton.Pt(),wt);
	else
	  h_phoPt_promptPho_mindr_leph0p5->Fill(bestPhoton.Pt(),wt);
	
      }
    else
      {
	h_phoPt_NopromptPho->Fill(bestPhoton.Pt(),wt);
	
      }
    if(s_sample.Contains("TTJets-HT") && madHT<600) continue;    
    
    if(s_sample.Contains("TTJets-inc") && madHT>600) continue;
    if(!genphocheck)
      {
        genphomatch_before++;
        double mindr_Pho_genlep=getGenLep(bestPhoton);
	h_minPho_lep->Fill(mindr_Pho_genlep,wt);
        if( s_sample.Contains("TTG") )
          {
            if(!hasGenPromptPhoton)
              {		
		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); continue;}
		else
		  {
		    if(madMinPhotonDeltaR >= 0.5)
		      h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
		    if(mindr_Pho_genlep >=0.5)
		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
		  }
              }
          }//Gen prompt                                                                                                                                                                                  
        if(s_sample.Contains("WG"))
          {
            if(!hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);

                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){ h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); continue;}
		else
                  {
                    if(madMinPhotonDeltaR >= 0.5)
                      h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)
		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
                  }

              }
          }//Gen prompt    
                                                                                                                                                                               
        if(s_sample.Contains("WJets"))
          {
            if(!hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);

		if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);   continue;}
		else
                  {
                    if(madMinPhotonDeltaR >= 0.5)
                      h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)
		      h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);
                  }

              }
          }
        if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets"))
          {
            if(!hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){ h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);  continue;}
		else
                  {
                    if(madMinPhotonDeltaR >= 0.5)
                      h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)
                      h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);
                  }


              }
          }
        if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
          {
            if(!(madMinPhotonDeltaR>0.4)) continue;

          }
        if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
          {
            if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;
          }
	if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG"))))
          {
            if(!(madMinPhotonDeltaR>0.5)) continue;
          }
        if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
          {
            if(!(madMinPhotonDeltaR<=0.5)) continue;
          }
        genphomatch_after++;

      }

    h_madHT_after->Fill(madHT,wt);
    h_madminPhotonDeltaR_after->Fill(madMinPhotonDeltaR,wt);
    h_minPho_lep_after->Fill(mindr_Pho_genlep,wt);
    if(hasGenPromptPhoton)
      h_phoPt_promptPho_after->Fill(bestPhoton.Pt(),wt);
    h_selectBaselineYields_->Fill("after bkg cmp",wt);
    nEvents_Selec[7]+=wt;
    out_nEventsTags[7]="overlap bkg";
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////   Evaluating the response from BDT  /////////////////////////////////////////
    
    met = MET;
    njets = nHadJets;
    btags = BTags;
    mTPhoMET = mTPhoMET;
    dPhi_PhoMET = dPhi_phoMET;
    dPhi_MetJet =dPhi_METjet1;
    st = ST;
    Double_t mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");
    h_mvaResponse->Fill(mvaValue,wt);
    if(Debug)
      cout<<"===load tree entry before SR ==="<<"\t"<<jentry<<endl;
    
    h_Njets[1]->Fill(nHadJets,wt);
    h_Nbjets[1]->Fill(BTags,wt);
    h_MET_[1]->Fill(MET,wt);
    h_PhotonPt[1]->Fill(bestPhoton.Pt(),wt);
    h_Mt_PhoMET[1]->Fill(mTPhoMET,wt);
    h_dPhi_PhoMET[1]->Fill(dPhi_PhoMET,wt);
    h_St[1]->Fill(ST,wt);
    h_HT[1]->Fill(Ht,wt);
    h_njets_vs_ST[1]->Fill(nHadJets,ST,wt);
    h_njets_vs_HT[1]->Fill(nHadJets,Ht,wt);
    h_ST_vs_ptPho[1]->Fill(ST,bestPhoton.Pt(),wt);
    h_mvaResponse_baseline[1]->Fill(mvaValue,wt);
    h_Njets_Varbin[1]->Fill(nHadJets,wt);
    h_Nbjets_Varbin[1]->Fill(BTags,wt);
    h_MET_Varbin[1]->Fill(MET,wt);
    h_PhotonPt_Varbin[1]->Fill(bestPhoton.Pt(),wt);
    h_St_Varbin[1]->Fill(ST,wt);
    h_mindr_Pho_genlep[1]->Fill(mindr_Pho_genlep);
    if (NElectrons == 0 && NMuons == 0 ) h_selectBaselineYields_->Fill("veto electron & Muon",wt);
    {
      if(isoElectronTracks==0 && isoMuonTracks ==0 && isoPionTracks==0)
	{
	  int Sbin_prev=getBinNoV6_WithOnlyBLSelec(BTags,nHadJets); 
	  h_Sbins_LL[30]->Fill(Sbin_prev,wt);
	}
    }
    nsurVived+=wt;
    //h_selectBaselineYields_->Fill("veto minDr(#photon,lept) <0.2",wt);
    int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
    TLorentzVector genPho1,genEle1,neutr,genMu1,genTau1,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
    vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_genLep1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_recEle,v_recMu;
    v_genEle2.clear();
    v_genEle1.clear();
    v_genMu1.clear();
    v_genMu2.clear();
    v_genTau1.clear();
    v_genTau2.clear();
    
    int leadGenPhoIdx=-100, mu_index=-100;
    int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
    int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
    for(int i=0;i<GenParticles->size();i++){
      if((*GenParticles)[i].Pt()!=0){
	if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=100) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 )
	  {
	  leadGenPhoIdx = i;
	  genPho1 = ((*GenParticles)[i]);
	  v_genPho1.push_back(genPho1);
	}
	       }        
    }	  
    bool hadtau = false;       
    for(int i=0 ; i < GenElectrons->size(); i++)
      {
	if((*GenElectrons)[i].Pt()!=0)// && (*GenElectrons)[i].Pt()!=inf)
	  {	     
	    nGenEle1++;
	    genEle1 = ((*GenElectrons)[i]);
	    v_genEle2.push_back(genEle1);
	    v_genLep2.push_back(genEle1);	      
	    h_GenpT[3]->Fill(genEle1.Pt()); h_GenEta[3]->Fill(genEle1.Eta());
	  }
      }
    for(int i=0 ; i < GenMuons->size(); i++)
      {
	if((*GenMuons)[i].Pt()!=0)// && (*GenMuons)[i].Pt()!=inf)
	  {
	    nGenMu1++;
	    genMu1 = ((*GenMuons)[i]);
	    v_genMu2.push_back(genMu1);
	    v_genLep2.push_back(genMu1);
	    h_GenpT[4]->Fill(genMu1.Pt()); h_GenEta[4]->Fill(genMu1.Eta());
	  }
      }
    //cout<<"entry:"<<"\t"<<jentry<<"\t"<<nGenMu1<<"\t"<<"genMuon size"<<"\t"<<v_genMu2.size()<<endl;
    for(int i=0 ; i < GenTaus->size(); i++)
      {
	if((*GenTaus)[i].Pt()!=0) // && (*GenTaus)[i].Pt()!=inf)
	  {
	    nGenTau1++;
	    genTau1 = ((*GenTaus)[i]);
	    v_genTau2.push_back(genTau1);
	    //v_genLep2.push_back(genTau1);
	    h_GenpT[5]->Fill(genTau1.Pt()); h_GenEta[5]->Fill(genTau1.Eta());
	    if((*GenTaus_had)[i])
	      nGentau_had1++;
	    
	  }
      }

    sortTLorVec(&v_genTau2);
    sortTLorVec(&v_genMu2);
    sortTLorVec(&v_genEle2);
    
    nlep=0;
    // if(Electrons->size() != NElectrons) //==0  && (v_genEle2.size()!=0))//GenElectrons->size()!=0))
    //   cout<<"event: "<<" "<<jentry<<"\t"<<NElectrons<<"\t"<<Electrons->size()<<"\t"<<GenElectrons->size()<<"\t"<<endl;
    bool Elec_passEtacut=false,Elec_passpTcut=false,Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false;
    bool Mu_passEtacut=false,Mu_passpTcut=false, Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
    
    if(v_genEle2.size()>=1)// && Electrons->size()!=0)
      {
	// if(abs(v_genEle2[0].Eta())<2.5 || (v_genEle2.size()>1 && abs(v_genEle2[1].Eta())<2.5))
	//   {
	//     Elec_passEtacut=true;
	//     if(v_genEle2[0].Pt()>10 || (v_genEle2.size()>1 && v_genEle2[1].Pt()>10))
	//       Elec_passpTcut = true;
	//     else
	//       Elec_passpTcut = false;
	//   }
	// else
	//   Elec_passEtacut=false;
	if((v_genEle2[0].Pt()>10 && abs(v_genEle2[0].Eta())<2.5) || (v_genEle2.size()>1 && v_genEle2[1].Pt()>10 && abs(v_genEle2[1].Eta())<2.5))
	  {
	    nlep=0;
	    nelec_reco=0;
	    pass_accep_elec++;	
	    Elec_passAccep = true;
	    // if(Electrons->size()==0)
	    //   cout<<v_recEle.size()<<endl;

	    for(int i=0;i<Electrons->size();i++)
	      {
		//Elec_passAccep = false; Elec_failAccep= false; Elec_passId= false; Elec_failId= false;Elec_passIso = false;Elec_failIso = false;
		//fail_IdElec=0;	    
		if(((*Electrons)[i].Pt()>10) && abs((*Electrons)[i].Eta()) < 2.5 && nlep<1)
		  {
		    double dR = MinDr((*Electrons)[i],v_genEle2);
		    if(dR<0.2){ pass_IdElec++; Elec_passId= true;
		      if((*Electrons_passIso)[i]==1)
			{ pass_isoElec++; Elec_passIso = true;nelec_reco++; nlep++; e_index=i; recElec=(*Electrons)[i]; v_recEle.push_back((*Electrons)[i]);}
		      else //if(dR<0.2 && (*Electrons_passIso)[i]!=1)
			{
			  fail_isoEle++; Elec_passIso = false;}

		    }
		    else
		      { fail_IdElec++; Elec_passId= false;}
		    // if(dR<0.2 && (*Electrons_passIso)[i]==1)
		    //   { pass_isoElec++; Elec_passIso = true;nelec_reco++; nlep++; e_index=i; recElec=(*Electrons)[i]; v_recEle.push_back((*Electrons)[i]);}
		    // else if(dR<0.2 && (*Electrons_passIso)[i]!=1)
		    //   {
		    // 	fail_isoEle++; Elec_passIso = false;}
		    // if(nlep>1) cout<<i<<""<<e_index<<" , e pt = "<<v_lep1.Pt()<<endl;                                                   
		  }
		// else
		//   {
		// 	fail_accep_elec++;
		// 	Elec_passAccep = false;
		//   }
	      }
	    if(nelec_reco==0 && Elec_passAccep && Elec_passId && Elec_passIso)
	      cout<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<"\t"<<endl;
	    sortTLorVec(&v_recEle);
	    if(Debug && v_recEle.size()==1) // NElectrons)
	      cout<<"entry: "<<"\t"<<jentry<<" "<<nlep<<"\t"<<"reco e size"<<" "<<v_recEle.size()<<" "<<Electrons->size()<<"\t"<<NElectrons<<endl;
	  }
	else
	  {
	    // if(Electrons->size()==0)
	    //   cout<<v_recEle.size()<<endl;
	    fail_accep_elec++;
	    Elec_passAccep = false;
	    if(abs(v_genEle2[0].Eta())<2.5 || (v_genEle2.size()>1 && abs(v_genEle2[1].Eta())<2.5))
	      {
		Elec_passEtacut=true;
		if(v_genEle2[0].Pt()>10 || (v_genEle2.size()>1 && v_genEle2[1].Pt()>10))
		  Elec_passpTcut = true;
		else
		  Elec_passpTcut = false;
	      }
	    else
	      Elec_passEtacut=false;

	  }	     
      }
     
    if(v_genMu2.size()>=1)// && Muons->size()!=0)
      {
	//cout<<jentry<<"\t"<<"inside genMu loop"<<endl;
	// if(abs(v_genMu2[0].Eta())<2.5 || (v_genMu2.size()>1 && abs(v_genMu2[1].Eta())<2.5))
        //   {
        //     Mu_passEtacut=true;
        //     if(v_genMu2[0].Pt()>10 || (v_genMu2.size()>1 && v_genMu2[1].Pt()>10))
	//       Mu_passpTcut = true;
        //     else
	//       Mu_passpTcut = false;
        //   }
	// else
        //   Mu_passEtacut=false;

	if(((v_genMu2[0].Pt()>10 && abs(v_genMu2[0].Eta())<2.4) ||  (v_genMu2.size()>1 && v_genMu2[1].Pt()>10 && abs(v_genMu2[1].Eta())<2.4)))
	  {
	    nlep=0;
	    pass_accep_mu++;
	    Mu_passAccep = true;	    
	    for(int i=0; i<Muons->size() ; i++)
	      {
		if(((*Muons)[i].Pt()>10) && abs((*Muons)[i].Eta())<2.4 && nlep<1)
		  {
		    double dR = MinDr((*Muons)[i],v_genMu2);
		    if(dR<0.2){ pass_IdMu++; Mu_passId= true;
		      
			if((*Muons_passIso)[i]==1)
			  { pass_isoMu++; Mu_passIso = true;nmu_reco++;nlep++; mu_index=i; recMu=(*Muons)[i]; v_recMu.push_back((*Muons)[i]);}
		      else
			{ fail_isoMu++; Mu_passIso = false;}
		      }
		      else
			{fail_IdMu++; Mu_passId= false;} 
		  }
	      }
	    
	    sortTLorVec(&v_recMu);
	    if(Debug && v_recMu.size()==1)
	      cout<<"entry: "<<"\t"<<jentry<<" "<<"reco muon size"<<" "<<v_recMu.size()<<" "<<Muons->size()<<endl;
	  }
	else
	  {
	    fail_accep_mu++;
            Mu_passAccep = false;
	    if(abs(v_genMu2[0].Eta())<2.5 || (v_genMu2.size()>1 && abs(v_genMu2[1].Eta())<2.5))
	      {
		Mu_passEtacut=true;
		if(v_genMu2[0].Pt()>10 || (v_genMu2.size()>1 && v_genMu2[1].Pt()>10))
		  Mu_passpTcut = true;
		else
		  Mu_passpTcut = false;
	      }
	    else
	      Mu_passEtacut=false;

          }
      }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////             SR region   ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //    cout<<""<<"\t"<<jentry<<endl;
    
    // for(int i=0;i<Photons_nonPrompt->size();i++)                                                                                       
    //   {                                                                                                                                
    //     cout<<jentry<<"\t"<<"promptPhoton:   "<<(*Photons_nonPrompt).at(i)<<"\t"<<(*Photons_genMatched).at(i)<<"\t"<<(*Photons_electronFakes).at(i)<<"\t"<<hasGenPromptPhoton<<endl;                                                                                                                                            if((*Photons_genMatched).at(i)!=hasGenPromptPhoton)
    // 	cout<<jentry<<"\t"<<(*Photons_genMatched).at(i)<<"\t"<<hasGenPromptPhoton<<endl;  
    //   }    
    
   if((*Electrons).size()==0 && nelec_reco!=0) //checking if the contribution to SR with no reco e- is coming up or not?
      cout<<nelec_reco<<endl;
   float ratio =0.0, ratio1=0.0, mindr_genElecPho=-999, mindr_=-9999;

   h_2d_genPromptPho_VsGenPho->Fill(hasGenPromptPhoton,v_genPho.size(),wt);
   // int nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0;
    if(lost_elec_flag && nelec_reco == 1 && nmu_reco == 0)
      {
	h_selectBaselineYields_CR->Fill("e-CR:w/o isotrack veto",wt);
	h_selectBaselineYields_v1->Fill("e-CR:w/o isotrack veto",wt);
	nEvents_Selec[12]+=wt;
	out_nEventsTags[12]="e- CR -initial";
	if(!Elec_passAccep || !Elec_passId || !Elec_passIso) //cross check                                                       
    	  cout<<"cross check - elec SR"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;  
	if(isoMuonTracks !=0 || isoPionTracks!=0) continue;
	h_selectBaselineYields_CR->Fill("e-CR:w isotrack veto",wt);
	h_selectBaselineYields_v1->Fill("e-CR:w isotrack veto",wt);
	nEvents_Selec[13]+=wt;
	out_nEventsTags[13]="veto charged mu/pi tracks in e-CR";
	double dr2 = bestPhoton.DeltaR((*Electrons)[e_index]);
	if(dr2<=0.2) { h_selectBaselineYields_CR->Fill("e-CR: fake rates",wt); //continue;
	  h_Njets_v1[1]->Fill(nHadJets,wt);
	  h_Nbjets_v1[1]->Fill(BTags,wt);
	  h_MET_v1[1]->Fill(MET,wt);
	  h_PhotonPt_v1[1]->Fill(bestPhoton.Pt(),wt);
	  h_St_v1[1]->Fill(ST,wt);	 
	  continue;
	}
	else{ h_selectBaselineYields_CR->Fill("e-CR: real photon",wt); h_selectBaselineYields_v1->Fill("e-CR: real photon",wt);}
	out_nEventsTags[14]="dR(e,#gamma)>0.2, e-CR";
	nEvents_Selec[14]+=wt;
	ST = ST+(*Electrons)[e_index].Pt();
	st = ST;
	mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");
	double mTElecMET=sqrt(2*((*Electrons)[e_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,(*Electrons)[e_index].Phi()))));
	if(mTElecMET>100) {continue; h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);}
	else { h_selectBaselineYields_CR->Fill("e-CR:mT>100",wt);}
	
	//updated checks for prompt/non prompt checks
	if(!(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons) && !(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons) && hasGenPromptPhoton)//(*Photons_genMatched).at(bestPhotonIndxAmongPhotons))
	  {
	    h_Njets_v1[51]->Fill(nHadJets,wt);
            h_Nbjets_v1[51]->Fill(BTags,wt);
            h_MET_v1[51]->Fill(MET,wt);
            h_PhotonPt_v1[51]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[51]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[51]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[51]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[51]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	    h_hasGenPromptPhoton_v2[51]->Fill(hasGenPromptPhoton,wt);
	    h_genPromptPho[51]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	    h_mindR_genJets[51]->Fill(genmindr,wt);
	    h_mindR_recoJets[51]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
	      {
		h_Jets_chargedEmEnergyFraction_v1[51]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
		h_Jets_chargedHadronEnergyFraction_v1[51]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
		h_Jets_photonEnergyFraction_v1[51]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
		h_Jets_photonMultiplicity_v1[51]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
		h_Jets_neutralEmEnergyFraction_v1[51]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
		h_Jets_neutralHadronEnergyFraction_v1[51]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }
	    else
	      {
		h_Jets_chargedEmEnergyFraction[51]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);                                     
		h_Jets_chargedHadronEnergyFraction[51]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);                            
		h_Jets_photonEnergyFraction[51]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);                                           
		h_Jets_photonMultiplicity[51]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);                                               
		h_Jets_neutralEmEnergyFraction[51]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);                                     
		h_Jets_neutralHadronEnergyFraction[51]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);      
	      }
	    
	    //	    if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons)!=hasGenPromptPhoton)// && hasGenPromptPhoton==1)// && (*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)==1)
	    //cout<<jentry<<"\t"<<"gen Photon: "<<v_genPho1.size()<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_genMatched).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;
	  }
	else if((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons))
	  {
	    h_Njets_v1[53]->Fill(nHadJets,wt);
            h_Nbjets_v1[53]->Fill(BTags,wt);
            h_MET_v1[53]->Fill(MET,wt);
            h_PhotonPt_v1[53]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[53]->Fill(ST,wt);
            h_NonPrompt_ElectronFake[53]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[53]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[53]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
            h_hasGenPromptPhoton_v2[53]->Fill(hasGenPromptPhoton,wt);
            h_genPromptPho[53]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);

	    h_mindR_genJets[53]->Fill(genmindr,wt);
            h_mindR_recoJets[53]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction_v1[53]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction_v1[53]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction_v1[53]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity_v1[53]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction_v1[53]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction_v1[53]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
              }
            else//if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction[53]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction[53]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction[53]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity[53]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction[53]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction[53]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }

	    // if(hasGenPromptPhoton)//(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==0 && (*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)==0)
	      //cout<<jentry<<"\t"<<"FR : "<<"\t"<<Photons->size()<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_genMatched).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;
	  }
	else if((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==1)
	  {
	    h_Njets_v1[52]->Fill(nHadJets,wt);
            h_Nbjets_v1[52]->Fill(BTags,wt);
            h_MET_v1[52]->Fill(MET,wt);
            h_PhotonPt_v1[52]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[52]->Fill(ST,wt);
            h_NonPrompt_ElectronFake[52]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[52]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[52]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
            h_hasGenPromptPhoton_v2[52]->Fill(hasGenPromptPhoton,wt);
            h_genPromptPho[52]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	    
	    h_mindR_genJets[52]->Fill(genmindr,wt);
            h_mindR_recoJets[52]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction_v1[52]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction_v1[52]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction_v1[52]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity_v1[52]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction_v1[52]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction_v1[52]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
              }

            else //if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction[52]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction[52]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction[52]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity[52]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction[52]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction[52]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }

	  }
	else
	  {
	    if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons))
	      {
		h_Njets_v1[54]->Fill(nHadJets,wt);
		h_Nbjets_v1[54]->Fill(BTags,wt);
		h_MET_v1[54]->Fill(MET,wt);
		h_PhotonPt_v1[54]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[54]->Fill(ST,wt);
		h_NonPrompt_ElectronFake[54]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
		h_ElectronFake[54]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
		h_NonPrompt[54]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
		h_hasGenPromptPhoton_v2[54]->Fill(hasGenPromptPhoton,wt);
		h_genPromptPho[54]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
		h_mindR_genJets[54]->Fill(genmindr,wt);
		h_mindR_recoJets[54]->Fill(recojetmindr,wt);
		if(recojetmindr<0.3)
		  {
		    h_Jets_chargedEmEnergyFraction_v1[54]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_chargedHadronEnergyFraction_v1[54]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonEnergyFraction_v1[54]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonMultiplicity_v1[54]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
		    h_Jets_neutralEmEnergyFraction_v1[54]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_neutralHadronEnergyFraction_v1[54]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
		  }

		if(recojetmindr<0.3)
		  {
		    h_Jets_chargedEmEnergyFraction[54]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_chargedHadronEnergyFraction[54]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonEnergyFraction[54]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonMultiplicity[54]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
		    h_Jets_neutralEmEnergyFraction[54]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_neutralHadronEnergyFraction[54]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
		  }

	      }
	    else
	      {
		h_Njets_v1[55]->Fill(nHadJets,wt);
		h_Nbjets_v1[55]->Fill(BTags,wt);
		h_MET_v1[55]->Fill(MET,wt);
		h_PhotonPt_v1[55]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[55]->Fill(ST,wt);
		h_NonPrompt_ElectronFake[55]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
		h_ElectronFake[55]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
		h_NonPrompt[55]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
		h_hasGenPromptPhoton_v2[55]->Fill(hasGenPromptPhoton,wt);
		h_genPromptPho[55]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
		h_mindR_genJets[55]->Fill(genmindr,wt);
		h_mindR_recoJets[55]->Fill(recojetmindr,wt);
		if(recojetmindr<0.3)
		  {
		    h_Jets_chargedEmEnergyFraction_v1[55]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_chargedHadronEnergyFraction_v1[55]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonEnergyFraction_v1[55]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonMultiplicity_v1[55]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
		    h_Jets_neutralEmEnergyFraction_v1[55]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_neutralHadronEnergyFraction_v1[55]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
		  }		
		else //if(recojetmindr<0.3)
		  {
		    h_Jets_chargedEmEnergyFraction[55]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_chargedHadronEnergyFraction[55]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonEnergyFraction[55]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
		    h_Jets_photonMultiplicity[55]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
		    h_Jets_neutralEmEnergyFraction[55]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
		    h_Jets_neutralHadronEnergyFraction[55]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
		  }


	      }
	    //if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons))
	    //cout<<jentry<<"\t"<<"Else: "<<Photons->size()<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_genMatched).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;
	  }
	if(hasGenPromptPhoton)
	  {
	    // if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons)==0)
	    //   {
	    // for(int i=0;i<Photons_nonPrompt->size();i++)
	    //   {
	    // 	//cout<<jentry<<"\t"<<(*Photons_nonPrompt).at(i)<<"\t"<<(*Photons_genMatched).at(i)<<"\t"<<(*Photons_electronFakes).at(i)<<"\t"<<hasGenPromptPhoton<<endl;                                                                                                              		// if((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==hasGenPromptPhoton)
	    // 	//   cout<<jentry<<"\t"<<"propmt "<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;
	    //   }
	    // //if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons)==0)
	    // cout<<jentry<<"\t"<<"promptPhoton: "<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_genMatched).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<endl;
	    //   }
	    h_selectBaselineYields_CR->Fill("e-CR: gen prompt #gamma",wt);
	    h_Njets_v1[2]->Fill(nHadJets,wt);
	    h_Nbjets_v1[2]->Fill(BTags,wt);
	    h_MET_v1[2]->Fill(MET,wt);
	    h_PhotonPt_v1[2]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[2]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[2]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[2]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[2]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);

	  }
	else
	  {
	    // if((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==hasGenPromptPhoton)
	    //   cout<<jentry<<"\t"<<"\t"<<"non-prompt"<<"\t"<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<hasGenPromptPhoton<<endl;

	    h_selectBaselineYields_CR->Fill("e-CR:non-gen prompt #gamma",wt);
	    h_Njets_v1[3]->Fill(nHadJets,wt);
	    h_Nbjets_v1[3]->Fill(BTags,wt);
	    h_MET_v1[3]->Fill(MET,wt);
	    h_PhotonPt_v1[3]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[3]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[3]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_ElectronFake[3]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_NonPrompt[3]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	    // cout<<"entry: "<<"\t"<<jentry<<"\t"<<v_genPho1.size()<<"\t"<<Photons_nonPrompt->size()<<"\t"<<nPhotons<<endl;
	    // for(int i=0;i<Photons_nonPrompt->size();i++)
	    //   {
	    //	    if((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==0)
	    // if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons)==0)
            //   {
	    // 	for(int i=0;i<Photons_nonPrompt->size();i++)
	    // 	  {
	    // 	    cout<<jentry<<"\t"<<(*Photons_nonPrompt).at(i)<<"\t"<<(*Photons_genMatched).at(i)<<"\t"<<(*Photons_electronFakes).at(i)<<"\t"<<hasGenPromptPhoton<<endl;                                                                                                      		  }

	    // cout<<jentry<<"\t"<<"non-promptPhoton: "<<(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_genMatched).at(bestPhotonIndxAmongPhotons)<<"\t"<<(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons)<<endl;
	    //   }
		//}
	    // for (int igen=0; igen<branch_size; igen++)
	    //   {
	    // 	pdgID= (*GenParticles_PdgId)[igen];
	    // 	parentId=(*GenParticles_ParentId)[igen];
	    // 	cout<<igen<<"\t"<<pdgID<<"\t"<<parentId<<endl;
	    //   }
	    if(v_genPho1.size()==0)
	      {
		h_Njets_v1[47]->Fill(nHadJets,wt);
		h_Nbjets_v1[47]->Fill(BTags,wt);
		h_MET_v1[47]->Fill(MET,wt);
		h_PhotonPt_v1[47]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[47]->Fill(ST,wt);
	      }
	    else
	      {
		h_Njets_v1[48]->Fill(nHadJets,wt);
                h_Nbjets_v1[48]->Fill(BTags,wt);
                h_MET_v1[48]->Fill(MET,wt);
                h_PhotonPt_v1[48]->Fill(bestPhoton.Pt(),wt);
                h_St_v1[48]->Fill(ST,wt);

	      }

	  }
	

	// if(hasGenPromptPhoton)
	//   h_selectBaselineYields_CR->Fill("e-CR: gen prompt #gamma",wt);
	// else
	//   h_selectBaselineYields_CR->Fill("e-CR:non-gen prompt #gamma",wt);
	h_selectBaselineYields_v1->Fill("e-CR:final Event yields",wt);
	h_selectBaselineYields_CR->Fill("e-CR:final Event yields",wt);
	h_Njets_v1[4]->Fill(nHadJets,wt);
	h_Nbjets_v1[4]->Fill(BTags,wt);
	h_MET_v1[4]->Fill(MET,wt);
	h_PhotonPt_v1[4]->Fill(bestPhoton.Pt(),wt);
	h_St_v1[4]->Fill(ST,wt);

	nEvents_Selec[15]+=wt;
	out_nEventsTags[15]="mT<100 - e-CR -final";
	h_Njets[2]->Fill(nHadJets,wt);
	h_Nbjets[2]->Fill(BTags,wt);
	h_MET_[2]->Fill(MET,wt);
	h_PhotonPt[2]->Fill(bestPhoton.Pt(),wt);
	h_Mt_PhoMET[2]->Fill(mTPhoMET,wt);
	h_dPhi_PhoMET[2]->Fill(dPhi_PhoMET,wt);
	h_St[2]->Fill(ST,wt);
	h_HT[2]->Fill(Ht,wt);
	h_njets_vs_ST[2]->Fill(nHadJets,ST,wt);
	h_njets_vs_HT[2]->Fill(nHadJets,Ht,wt);
	h_ST_vs_ptPho[2]->Fill(ST,bestPhoton.Pt(),wt);
	h_mvaResponse_baseline[2]->Fill(mvaValue,wt);
	h_Njets_Varbin[2]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[2]->Fill(BTags,wt);
	h_MET_Varbin[2]->Fill(MET,wt);
	h_PhotonPt_Varbin[2]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[2]->Fill(ST,wt);
	//additional checks
	//h_mindr_Pho_RecoEle[2]->Fill(dr2,wt);	
	// mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
	// h_mindr_Pho_genElec[2]->Fill(mindr_genElecPho,wt);
	// h_mindr_Pho_genlep[2]->Fill(mindr_Pho_genlep);	
	// h_hasGenPromptPhoton_v1[2]>Fill(hasGenPromptPhoton);
	
	h_mindr_Pho_RecoEle[2]->Fill(dr2,wt);
	if(!s_data.Contains("data"))
	  {
	    if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[2]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
	    else if(v_genEle2.size()>1){
	      if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[2]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
	      else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[2]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
	    }
	  
	    ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
	    if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
	    h_ratio_bestPho_genlep_v1[2]->Fill(ratio);
	    h_ratio_bestPho_genlep_v2[2]->Fill(ratio1); 

	    mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
	    h_mindr_Pho_genElec[2]->Fill(mindr_genElecPho,wt);
	    h_mindr_Pho_genlep[2]->Fill(mindr_Pho_genlep);
	    h_hasGenPromptPhoton_v1[2]->Fill(hasGenPromptPhoton);

	  }
	int Sbin = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL[2]->Fill(Sbin, wt);	
	// Tfbins = getBinNoV1_le(BTags,nHadJets);
        // h_TFbins_ElecLL_v1[2]->Fill(Tfbins, wt);	
	wt_LL = wt*h_TF->GetBinContent(Sbin+1);
	// if(Sbin>15)
	//   cout<<Sbin<<"\t"<<h_TF->GetBinContent(Sbin+1)<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;
	if(Debug)
       	cout<<Sbin<<"\t"<<h_TF->GetBinContent(Sbin+1)<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;
	searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	//cout<<Sbin<<"\t"<<h_TF->GetBinContent(Sbin+1)<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;
	h_Sbins_LL_Validation[2]->Fill(searchBin,wt_LL);
	h_Njets[40]->Fill(nHadJets,wt_LL);
        h_Nbjets[40]->Fill(BTags,wt_LL);
        h_MET_[40]->Fill(MET,wt_LL);
	h_PhotonPt[40]->Fill(bestPhoton.Pt(),wt_LL);
        h_St[40]->Fill(ST,wt_LL);
	h_MET_Varbin[40]->Fill(MET,wt_LL);
        h_PhotonPt_Varbin[40]->Fill(bestPhoton.Pt(),wt_LL);

	Tfbins = getBinNoV1_le(BTags,nHadJets);
        h_TFbins_ElecLL_v1[2]->Fill(Tfbins,wt);
	h_TFbins_ElecLL_validation[2]->Fill(Sbin,wt_LL);	
	h_Sbins_LL[2]->Fill(searchBin,wt);
	nCR_elec+=wt; 
	h_selectBaselineYields_->Fill("elec CR",wt);	
	if(mvaValue<0.0)
	  {
	    nEvents_Selec[16]+=wt;
	    out_nEventsTags[16]="e-CR -after BDT resp";
	    h_Njets[12]->Fill(nHadJets,wt);
	    h_Nbjets[12]->Fill(BTags,wt);
	    h_MET_[12]->Fill(MET,wt);
	    h_PhotonPt[12]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[12]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[12]->Fill(dPhi_PhoMET,wt);
	    h_St[12]->Fill(ST,wt);
	    h_HT[12]->Fill(Ht,wt);
	    h_njets_vs_ST[12]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[12]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[12]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[12]->Fill(mvaValue,wt);
	    h_Njets_Varbin[12]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[12]->Fill(BTags,wt);
	    h_MET_Varbin[12]->Fill(MET,wt);
	    h_PhotonPt_Varbin[12]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[12]->Fill(ST,wt);
	    h_mindr_Pho_RecoEle[12]->Fill(dr2,wt);
	    if(!s_data.Contains("data"))
	      {
		if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[12]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
		else if(v_genEle2.size()>1){
		  if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[12]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
		  else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[12]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
		}
		ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
		if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
		h_ratio_bestPho_genlep_v1[12]->Fill(ratio);
		h_ratio_bestPho_genlep_v2[12]->Fill(ratio1);

		mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
		h_mindr_Pho_genElec[12]->Fill(mindr_genElecPho,wt);
		h_mindr_Pho_genlep[12]->Fill(mindr_Pho_genlep);
		h_hasGenPromptPhoton_v1[12]->Fill(hasGenPromptPhoton);
	      }

	    Sbin = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL[12]->Fill(Sbin, wt);
	    wt_LL = wt*h_TF->GetBinContent(Sbin+1);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[12]->Fill(searchBin,wt_LL);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[12]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[12]->Fill(Sbin,wt_LL);
	    h_Sbins_LL[12]->Fill(searchBin,wt);

	    h_Njets[41]->Fill(nHadJets,wt_LL);
	    h_Nbjets[41]->Fill(BTags,wt_LL);
	    h_MET_[41]->Fill(MET,wt_LL);
	    h_PhotonPt[41]->Fill(bestPhoton.Pt(),wt_LL);
	    h_St[41]->Fill(ST,wt_LL);



	  }
      }
    if(Debug)
      cout<<"CR elec "<<"\t"<<jentry<<endl;
    if(lost_elec_flag && nelec_reco ==0 && nmu_reco==0 )//&& (*Electrons).size()!=0)
      {
	//if(Elec_passAccep && Elec_passId && Elec_passIso)
	//cout<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
	
	///////////////////  Best photon - to be prompt photon  //////////////////////////////////////                               
	// if(!s_data.Contains("data"))
	//   {
	nEvents_Selec[17]+=wt;
	out_nEventsTags[17]="e-SR initial";
	h_selectBaselineYields_SR->Fill("e-SR - initial",wt);
	h_selectBaselineYields_v1->Fill("e-SR - initial",wt);
	if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
        nEvents_Selec[19]+=wt;
        out_nEventsTags[19]="e-SR: veto charged pion track";
	h_selectBaselineYields_SR->Fill("e-SR: iso track veto",wt);
	h_selectBaselineYields_v1->Fill("e-SR: iso track veto",wt);
        if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue;//to reject W->qq' type of event                                                                  
	h_selectBaselineYields_SR->Fill("e-SR: After rejecting W->qq events",wt);
	h_selectBaselineYields_v1->Fill("e-SR: iso track veto",wt);
        nEvents_Selec[20]+=wt;
        out_nEventsTags[20]="w-qq: e-SR";
        if(nGentau_had1>1) continue;
	h_selectBaselineYields_SR->Fill("e-SR: After rejecting tau-hadronic>1",wt);
	h_selectBaselineYields_v1->Fill("e-SR: After rejecting tau-hadronic>1",wt);
        nEvents_Selec[21]+=wt;
        out_nEventsTags[21]="More taus-had decay: e-SR";
        survived_vetohad+=wt;
        if(nGenMu1>0) continue;
	h_selectBaselineYields_SR->Fill("e-SR: GenMu1>0",wt);
	h_selectBaselineYields_v1->Fill("e-SR: GenMu1>0",wt);
        out_nEventsTags[22]="e-SR: Gen Mu>0";
        nEvents_Selec[22]+=wt;
        if(nGenEle1==0) continue;
	h_selectBaselineYields_SR->Fill("e-SR: no gen e-",wt);
	h_selectBaselineYields_v1->Fill("e-SR: no gen e-",wt);
        out_nEventsTags[23]="gen e-!=0, e-SR final";
        nEvents_Selec[23]+=wt;
        h_selectBaselineYields_->Fill("Elec SR",wt);
	//h_selectBaselineYields_v1->Fill("e-SR:final",wt);
        nSR_elec+=wt;
        survived_elecge1+=wt;
        if(v_genEle2.size()==0) {TLorentzVector v1;v_genEle2.push_back(v1);}
        sortTLorVec(&v_genEle2);
        v_genEle1=v_genEle2;
	h_selectBaselineYields_SR->Fill("before #gamma checks",wt);
	
	if(!(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons) && !(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons) && hasGenPromptPhoton)//(*Photons_genMatched).at(bestPhotonIndxAmongPhotons))          
          {
            h_Njets_v1[56]->Fill(nHadJets,wt);
            h_Nbjets_v1[56]->Fill(BTags,wt);
            h_MET_v1[56]->Fill(MET,wt);
            h_PhotonPt_v1[56]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[56]->Fill(ST,wt);
            h_NonPrompt_ElectronFake[56]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[56]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[56]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
            h_hasGenPromptPhoton_v2[56]->Fill(hasGenPromptPhoton,wt);
            h_genPromptPho[56]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	    h_mindR_genJets[56]->Fill(genmindr,wt);
            h_mindR_recoJets[56]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction_v1[56]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction_v1[56]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction_v1[56]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity_v1[56]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction_v1[56]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction_v1[56]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
              }
	    else
	      //if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction[56]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction[56]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction[56]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity[56]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction[56]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction[56]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }

             
      }
    else if((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons))
      {
	h_Njets_v1[58]->Fill(nHadJets,wt);
	h_Nbjets_v1[58]->Fill(BTags,wt);
	h_MET_v1[58]->Fill(MET,wt);
	h_PhotonPt_v1[58]->Fill(bestPhoton.Pt(),wt);
	h_St_v1[58]->Fill(ST,wt);
	h_NonPrompt_ElectronFake[58]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	h_ElectronFake[58]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	h_NonPrompt[58]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	h_hasGenPromptPhoton_v2[58]->Fill(hasGenPromptPhoton,wt);
	h_genPromptPho[58]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	h_mindR_genJets[58]->Fill(genmindr,wt);
	h_mindR_recoJets[58]->Fill(recojetmindr,wt);
	if(recojetmindr<0.3)
	  {
	    h_Jets_chargedEmEnergyFraction_v1[58]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_chargedHadronEnergyFraction_v1[58]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonEnergyFraction_v1[58]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonMultiplicity_v1[58]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
	    h_Jets_neutralEmEnergyFraction_v1[58]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_neutralHadronEnergyFraction_v1[58]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	  }
	else //if(recojetmindr<0.3)
	  {
	    h_Jets_chargedEmEnergyFraction[58]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_chargedHadronEnergyFraction[58]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonEnergyFraction[58]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonMultiplicity[58]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
	    h_Jets_neutralEmEnergyFraction[58]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_neutralHadronEnergyFraction[58]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	  }


     }

    else if((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==1)
      {
	h_Njets_v1[57]->Fill(nHadJets,wt);
	h_Nbjets_v1[57]->Fill(BTags,wt);
	h_MET_v1[57]->Fill(MET,wt);
	h_PhotonPt_v1[57]->Fill(bestPhoton.Pt(),wt);
	h_St_v1[57]->Fill(ST,wt);
	h_NonPrompt_ElectronFake[57]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	h_ElectronFake[57]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	h_NonPrompt[57]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	h_hasGenPromptPhoton_v2[57]->Fill(hasGenPromptPhoton,wt);
	h_genPromptPho[57]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	h_mindR_genJets[57]->Fill(genmindr,wt);
	h_mindR_recoJets[57]->Fill(recojetmindr,wt);
	if(recojetmindr<0.3)
	  {
	    h_Jets_chargedEmEnergyFraction_v1[57]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_chargedHadronEnergyFraction_v1[57]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonEnergyFraction_v1[57]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonMultiplicity_v1[57]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
	    h_Jets_neutralEmEnergyFraction_v1[57]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_neutralHadronEnergyFraction_v1[57]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	  }

	else //if(recojetmindr<0.3)
	  {
	    h_Jets_chargedEmEnergyFraction[57]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_chargedHadronEnergyFraction[57]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonEnergyFraction[57]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
	    h_Jets_photonMultiplicity[57]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
	    h_Jets_neutralEmEnergyFraction[57]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
	    h_Jets_neutralHadronEnergyFraction[57]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	  }


      }
    else
      {
	if((*Photons_genMatched).at(bestPhotonIndxAmongPhotons))
	  {
	    h_Njets_v1[59]->Fill(nHadJets,wt);
	    h_Nbjets_v1[59]->Fill(BTags,wt);
	    h_MET_v1[59]->Fill(MET,wt);
	    h_PhotonPt_v1[59]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[59]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[59]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_ElectronFake[59]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_NonPrompt[59]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	    h_hasGenPromptPhoton_v2[59]->Fill(hasGenPromptPhoton,wt);
	    h_genPromptPho[59]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	    h_mindR_genJets[59]->Fill(genmindr,wt);
            h_mindR_recoJets[59]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction_v1[59]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction_v1[59]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction_v1[59]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity_v1[59]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction_v1[59]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction_v1[59]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
              }
            else //if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction[59]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction[59]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction[59]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity[59]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction[59]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction[59]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }

	  }
	else
	  {
	    h_Njets_v1[60]->Fill(nHadJets,wt);
	    h_Nbjets_v1[60]->Fill(BTags,wt);
	    h_MET_v1[60]->Fill(MET,wt);
	    h_PhotonPt_v1[60]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[60]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[60]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_ElectronFake[60]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
	    h_NonPrompt[60]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);
	    h_hasGenPromptPhoton_v2[60]->Fill(hasGenPromptPhoton,wt);
	    h_genPromptPho[60]->Fill((*Photons_genMatched).at(bestPhotonIndxAmongPhotons),wt);
	    h_mindR_genJets[60]->Fill(genmindr,wt);
            h_mindR_recoJets[60]->Fill(recojetmindr,wt);
	    if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction_v1[60]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction_v1[60]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction_v1[60]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity_v1[60]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction_v1[60]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction_v1[60]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
              }
            else //if(recojetmindr<0.3)
              {
                h_Jets_chargedEmEnergyFraction[60]->Fill((*Jets_chargedEmEnergyFraction).at(minDRindx),wt);
                h_Jets_chargedHadronEnergyFraction[60]->Fill((*Jets_chargedHadronEnergyFraction).at(minDRindx),wt);
                h_Jets_photonEnergyFraction[60]->Fill((*Jets_photonEnergyFraction).at(minDRindx),wt);
                h_Jets_photonMultiplicity[60]->Fill((*Jets_photonMultiplicity).at(minDRindx),wt);
                h_Jets_neutralEmEnergyFraction[60]->Fill((*Jets_neutralEmEnergyFraction).at(minDRindx),wt);
                h_Jets_neutralHadronEnergyFraction[60]->Fill((*Jets_neutralHadronEnergyFraction).at(minDRindx),wt);
	      }


	  }
      }
	if(hasGenPromptPhoton)
	  {
	    h_Njets_v1[5]->Fill(nHadJets,wt);
	    h_Nbjets_v1[5]->Fill(BTags,wt);
	    h_MET_v1[5]->Fill(MET,wt);
	    h_PhotonPt_v1[5]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[5]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[5]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[5]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[5]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);

	    double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton); // gen information is not saved for photon
	    if(mindR_genPho<=0.2) 
	      {
		h_selectBaselineYields_SR->Fill("e-SR: dR(reco-#gamma,gen-#gamma)<0.2",wt);
		h_Njets_v1[6]->Fill(nHadJets,wt);
		h_Nbjets_v1[6]->Fill(BTags,wt);
		h_MET_v1[6]->Fill(MET,wt);
		h_PhotonPt_v1[6]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[6]->Fill(ST,wt);
	      }
	    else if(mindR_genPho>0.2)
	      {		
		if(v_genPho1.size()==0)
		  {
		    h_Njets_v1[7]->Fill(nHadJets,wt);
		    h_Nbjets_v1[7]->Fill(BTags,wt);
		    h_MET_v1[7]->Fill(MET,wt);
		    h_PhotonPt_v1[7]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[7]->Fill(ST,wt);
		  }
		else
		  {		    
                    h_Njets_v1[8]->Fill(nHadJets,wt);
                    h_Nbjets_v1[8]->Fill(BTags,wt);
                    h_MET_v1[8]->Fill(MET,wt);
                    h_PhotonPt_v1[8]->Fill(bestPhoton.Pt(),wt);
                    h_St_v1[8]->Fill(ST,wt);
		  }
		h_selectBaselineYields_SR->Fill("e-SR: !dR(reco-#gamma,gen-#gamma)<0.2",wt);
		if(v_genEle2.size()==1)
		  {
		    if(bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[0].Pt()/bestPhoton.Pt()<0.79)))
		      {
			h_Njets_v1[9]->Fill(nHadJets,wt);
			h_Nbjets_v1[9]->Fill(BTags,wt);
			h_MET_v1[9]->Fill(MET,wt);
			h_PhotonPt_v1[9]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[9]->Fill(ST,wt);
			h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR",wt);		    
			int flag=-1,pdgID=-1;
			std::vector<int> array1;
			array1 = dR_recoPho_GenParticle(bestPhoton);
			pdgID = array1[0];
			flag=array1[1];
			if(abs(pdgID)==15 && flag==1)
			  {
			    h_Njets_v1[10]->Fill(nHadJets,wt);
			    h_Nbjets_v1[10]->Fill(BTags,wt);
			    h_MET_v1[10]->Fill(MET,wt);
			    h_PhotonPt_v1[10]->Fill(bestPhoton.Pt(),wt);
			    h_St_v1[10]->Fill(ST,wt);
			    h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && dR(tau,reco #gamma)<0.2",wt);
			  }
			else if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
			  {
			    h_Njets_v1[11]->Fill(nHadJets,wt);
                            h_Nbjets_v1[11]->Fill(BTags,wt);
                            h_MET_v1[11]->Fill(MET,wt);
                            h_PhotonPt_v1[11]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[11]->Fill(ST,wt);
			  }
			else if (!flag)
			  {
			    h_Njets_v1[12]->Fill(nHadJets,wt);
                            h_Nbjets_v1[12]->Fill(BTags,wt);
                            h_MET_v1[12]->Fill(MET,wt);
                            h_PhotonPt_v1[12]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[12]->Fill(ST,wt);
			  }
			// cout<<"1e 1#gammaSRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
		      }
		    else if(bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))
		      {
			h_Njets_v1[13]->Fill(nHadJets,wt);
                        h_Nbjets_v1[13]->Fill(BTags,wt);
                        h_MET_v1[13]->Fill(MET,wt);
                        h_PhotonPt_v1[13]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[13]->Fill(ST,wt);
			h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && FR",wt);
			NphotonsVsPhotPt->Fill(nPhotons,bestPhoton.Pt(),wt);
			BestPhotPt_->Fill(goodPho_n[0].Pt(),wt);
			if(nPhotons>1)
			  BestPhotPt_1->Fill(goodPho_n[1].Pt(),wt);

			//cout<<"1e, 1gamma FRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<"\t"<<v_genEle2[0].Pt()<<endl;
			continue;
		      }
		    else
		      {
			h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && !dR(tau/q/g,reco #gamma)<0.2",wt);
			int flag=-1,pdgID=-1;
			std::vector<int> array1;
                        array1 = dR_recoPho_GenParticle(bestPhoton);
                        pdgID = array1[0];
                        flag=array1[1];
                        if(abs(pdgID)==15 && flag==1)
                          {
                            h_Njets_v1[10]->Fill(nHadJets,wt);
                            h_Nbjets_v1[10]->Fill(BTags,wt);
                            h_MET_v1[10]->Fill(MET,wt);
                            h_PhotonPt_v1[10]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[10]->Fill(ST,wt);
                            h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && dR(tau,reco #gamma)<0.2",wt);
                          }
                        else if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
                          {
                            h_Njets_v1[11]->Fill(nHadJets,wt);
                            h_Nbjets_v1[11]->Fill(BTags,wt);
                            h_MET_v1[11]->Fill(MET,wt);
                            h_PhotonPt_v1[11]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[11]->Fill(ST,wt);
                          }
                        else if (!flag)
                          {
                            h_Njets_v1[12]->Fill(nHadJets,wt);
                            h_Nbjets_v1[12]->Fill(BTags,wt);
                            h_MET_v1[12]->Fill(MET,wt);
                            h_PhotonPt_v1[12]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[12]->Fill(ST,wt);
                          }

			// h_Njets_v1[14]->Fill(nHadJets,wt);
                        // h_Nbjets_v1[14]->Fill(BTags,wt);
                        // h_MET_v1[14]->Fill(MET,wt);
                        // h_PhotonPt_v1[14]->Fill(bestPhoton.Pt(),wt);
                        // h_St_v1[14]->Fill(ST,wt);
		      }
		      
	      } // 1-electron SR
		else if(v_genEle2.size()>1)
		  {
		    if((bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[0].Pt()/bestPhoton.Pt()<0.79))) && (bestPhoton.DeltaR(v_genEle2[1])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[1].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[1].Pt()/bestPhoton.Pt()<0.79))))
		      {

			h_Njets_v1[14]->Fill(nHadJets,wt);
                        h_Nbjets_v1[14]->Fill(BTags,wt);
                        h_MET_v1[14]->Fill(MET,wt);
                        h_PhotonPt_v1[14]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[14]->Fill(ST,wt);

			h_selectBaselineYields_SR->Fill("2e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR",wt);
			//cout<<"1e 1#gammaSRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
		      }
		    else if( ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))) || ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[1].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[1].Pt()/bestPhoton.Pt()>0.79))))
		      {
			h_Njets_v1[15]->Fill(nHadJets,wt);
                        h_Nbjets_v1[15]->Fill(BTags,wt);
                        h_MET_v1[15]->Fill(MET,wt);
                        h_PhotonPt_v1[15]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[15]->Fill(ST,wt);
			h_selectBaselineYields_SR->Fill("2e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && FR",wt);
                        //cout<<"1e, 1gamma FRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<"\t"<<v_genEle2[0].Pt()<<endl;
			continue;
		      }
		    else
		      {
			h_Njets_v1[16]->Fill(nHadJets,wt);
                        h_Nbjets_v1[16]->Fill(BTags,wt);
                        h_MET_v1[16]->Fill(MET,wt);
                        h_PhotonPt_v1[16]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[16]->Fill(ST,wt);
			h_selectBaselineYields_SR->Fill("2e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && !dR(tau,reco #gamma)<0.2",wt);
		      }
		  }
	      }
      } //Checking prompt photons
	h_selectBaselineYields_SR->Fill("before !#gamma checks",wt);
	
    if(!hasGenPromptPhoton)
	  {	    
	    h_Njets_v1[17]->Fill(nHadJets,wt);
	    h_Nbjets_v1[17]->Fill(BTags,wt);
	    h_MET_v1[17]->Fill(MET,wt);
	    h_PhotonPt_v1[17]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[17]->Fill(ST,wt);
	    h_NonPrompt_ElectronFake[17]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_ElectronFake[17]->Fill((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons),wt);
            h_NonPrompt[17]->Fill((*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons),wt);

	    if(v_genPho1.size()==0)
	      {
		h_Njets_v1[49]->Fill(nHadJets,wt);
		h_Nbjets_v1[49]->Fill(BTags,wt);
		h_MET_v1[49]->Fill(MET,wt);
		h_PhotonPt_v1[49]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[49]->Fill(ST,wt);

	      }
	  
	    else
	      {
		h_Njets_v1[50]->Fill(nHadJets,wt);
		h_Nbjets_v1[50]->Fill(BTags,wt);
		h_MET_v1[50]->Fill(MET,wt);
		h_PhotonPt_v1[50]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[50]->Fill(ST,wt);
		
	      }
	    
	    if(v_genEle2.size()==1)
	      {
		if(bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[0].Pt()/bestPhoton.Pt()<0.79)))
		  { ratio=v_genEle2[0].Pt()/bestPhoton.Pt(); mindr_=bestPhoton.DeltaR(v_genEle2[0]);
		    h_selectBaselineYields_SR->Fill("e-SR:!gen-#gamma && !FR",wt);
		    h_Njets_v1[18]->Fill(nHadJets,wt);
		    h_Nbjets_v1[18]->Fill(BTags,wt);
		    h_MET_v1[18]->Fill(MET,wt);
		    h_PhotonPt_v1[18]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[18]->Fill(ST,wt);
		    int flag=-1,pdgID=-1;
		    std::vector<int> array1;
	
		    array1 = dR_recoPho_GenParticle(bestPhoton);
		    pdgID = array1[0];
		    flag=array1[1];
		    if(abs(pdgID)==15 && flag==1)
		      {
			h_Njets_v1[19]->Fill(nHadJets,wt);
			h_Nbjets_v1[19]->Fill(BTags,wt);
			h_MET_v1[19]->Fill(MET,wt);
			h_PhotonPt_v1[19]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[19]->Fill(ST,wt);
			h_selectBaselineYields_SR->Fill("e-SR:!gen-#gamma && !FR && dR(tau,reco #gamma)<0.2",wt);

		      }
		    else if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
		      {
			h_Njets_v1[20]->Fill(nHadJets,wt);
			h_Nbjets_v1[20]->Fill(BTags,wt);
			h_MET_v1[20]->Fill(MET,wt);
			h_PhotonPt_v1[20]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[20]->Fill(ST,wt);
		      }
		    else if (!flag)
		      {
			h_Njets_v1[21]->Fill(nHadJets,wt);
			h_Nbjets_v1[21]->Fill(BTags,wt);
			h_MET_v1[21]->Fill(MET,wt);
			h_PhotonPt_v1[21]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[21]->Fill(ST,wt);
		      }
		 
		      //cout<<"1e SRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
		  }
		else if((bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79)))
		  {
		    h_Njets_v1[22]->Fill(nHadJets,wt);
		    h_Nbjets_v1[22]->Fill(BTags,wt);
		    h_MET_v1[22]->Fill(MET,wt);
		    h_PhotonPt_v1[22]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[22]->Fill(ST,wt);

		    h_selectBaselineYields_SR->Fill("e-SR:!gen-#gamma && FR && !dR(tau,reco #gamma)<0.2",wt); 
		    //cout<<"1 e FR entry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
		    continue;
		  }
			
		else //if(bestPhoton.DeltaR(v_genEle2[0])<0.2 || ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))
		  {
		    int flag=-1,pdgID=-1;
		    std::vector<int> array1;
		    array1 = dR_recoPho_GenParticle(bestPhoton);
		    pdgID = array1[0];
		    flag=array1[1];
	            if(abs(pdgID)==15 && flag==1)
                      {
			h_Njets_v1[19]->Fill(nHadJets,wt);
                        h_Nbjets_v1[19]->Fill(BTags,wt);
                        h_MET_v1[19]->Fill(MET,wt);
                        h_PhotonPt_v1[19]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[19]->Fill(ST,wt);
                        h_selectBaselineYields_SR->Fill("e-SR:!gen-#gamma && !FR && dR(tau,reco #gamma)<0.2",wt);

                      }
                    else if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
                      {
                        h_Njets_v1[20]->Fill(nHadJets,wt);
                        h_Nbjets_v1[20]->Fill(BTags,wt);
                        h_MET_v1[20]->Fill(MET,wt);
                        h_PhotonPt_v1[20]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[20]->Fill(ST,wt);
                      }
                    else if (!flag)
                      {
			h_Njets_v1[21]->Fill(nHadJets,wt);
                        h_Nbjets_v1[21]->Fill(BTags,wt);
                        h_MET_v1[21]->Fill(MET,wt);
			h_PhotonPt_v1[21]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[21]->Fill(ST,wt);
                      }
		    h_selectBaselineYields_SR->Fill("e-SR:!gen-#gamma && !FR && !dR(tau,reco #gamma)<0.2",wt);
		    //cout<<"1e FRentry: !FR && !pho"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<"\t"<<v_genEle2[0].Pt()<<endl;
		    //continue;
		  }
	      }
	    else if(v_genEle2.size()>1)
	      {
		if((bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21)||  v_genEle2[0].Pt()/bestPhoton.Pt()<0.79))) && (bestPhoton.DeltaR(v_genEle2[1])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[1].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[1].Pt()/bestPhoton.Pt()<0.79))))
	      {ratio=v_genEle2[0].Pt()/bestPhoton.Pt();  mindr_=bestPhoton.DeltaR(v_genEle2[0]);h_selectBaselineYields_SR->Fill("2e-SR:!gen-#gamma && !FR",wt);
		h_Njets_v1[23]->Fill(nHadJets,wt);
		h_Nbjets_v1[23]->Fill(BTags,wt);
		h_MET_v1[23]->Fill(MET,wt);
		h_PhotonPt_v1[23]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[23]->Fill(ST,wt);

		// if(bestPhoton.DeltaR(v_genTau2[0])<0.2)
		//   h_selectBaselineYields_SR->Fill("2e-SR:!gen-#gamma && !FR && dR(tau,reco #gamma)<0.2",wt);

		//cout<<"e-2: SRentry"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
	      }
		else  if( ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))) || ((bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[1].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[1].Pt()/bestPhoton.Pt()>0.79))))
		  {
		    h_Njets_v1[24]->Fill(nHadJets,wt);
		    h_Nbjets_v1[24]->Fill(BTags,wt);
		    h_MET_v1[24]->Fill(MET,wt);
		    h_PhotonPt_v1[24]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[24]->Fill(ST,wt);

		    h_selectBaselineYields_SR->Fill("2e-SR:!gen-#gamma && FR",wt);
		    continue;
                    //cout<<"2e tau faking SRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<endl;
		  }
		else
		  {
		    //cout<<"e2: FRentry:"<<jentry<<"\t"<<hasGenPromptPhoton<<"\t"<<v_genEle2.size()<<"\t"<<nGentau_had1<<"\t"<<"\t"<<bestPhoton.Pt()<<"\t"<<bestPhoton.DeltaR(v_genEle2[0])<<"\t"<<v_genEle2[0].Pt()/bestPhoton.Pt()<<"\t"<<v_genEle2[1].Pt()/bestPhoton.Pt()<<endl;
		    h_Njets_v1[25]->Fill(nHadJets,wt);
		    h_Nbjets_v1[25]->Fill(BTags,wt);
		    h_MET_v1[25]->Fill(MET,wt);
		    h_PhotonPt_v1[25]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[25]->Fill(ST,wt);

		    h_selectBaselineYields_SR->Fill("2e-SR:!gen-#gamma && !FR && !dR(tau,reco #gamma)<0.2",wt);
		    //continue;
		  }
	      }
	  }
	nEvents_Selec[18]+=wt;
	out_nEventsTags[18]="prompt Photon Checks, e-SR";
	h_selectBaselineYields_SR->Fill("e-SR: Final yields",wt);
	h_selectBaselineYields_v1->Fill("e-SR:final",wt);
	
	if(!Elec_passAccep )
	  {
	    if(!Elec_passEtacut)
	      {
	    h_Njets[22]->Fill(nHadJets,wt);
            h_Nbjets[22]->Fill(BTags,wt);
            h_MET_[22]->Fill(MET,wt);
	    h_PhotonPt[22]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[22]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[22]->Fill(dPhi_PhoMET,wt);
	    h_St[22]->Fill(ST,wt);
            h_HT[22]->Fill(Ht,wt);

            h_njets_vs_ST[22]->Fill(nHadJets,ST,wt);
            h_njets_vs_HT[22]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[22]->Fill(ST,bestPhoton.Pt(),wt);
            h_mvaResponse_baseline[22]->Fill(mvaValue,wt);
            h_Njets_Varbin[22]->Fill(nHadJets,wt);
            h_Nbjets_Varbin[22]->Fill(BTags,wt);
            h_MET_Varbin[22]->Fill(MET,wt);
            h_PhotonPt_Varbin[22]->Fill(bestPhoton.Pt(),wt);
            h_St_Varbin[22]->Fill(ST,wt);
            int Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[22]->Fill(Sbin, wt);
            searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[22]->Fill(searchBin,wt);
            Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[22]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[22]->Fill(Sbin,wt);
            h_Sbins_LL[22]->Fill(searchBin,wt);
	    if(mvaValue>0.0)
              {
                h_Njets[26]->Fill(nHadJets,wt);
                h_Nbjets[26]->Fill(BTags,wt);
                h_MET_[26]->Fill(MET,wt);
                h_PhotonPt[26]->Fill(bestPhoton.Pt(),wt);
                h_Mt_PhoMET[26]->Fill(mTPhoMET,wt);
                h_dPhi_PhoMET[26]->Fill(dPhi_PhoMET,wt);
                h_St[26]->Fill(ST,wt);
                h_HT[26]->Fill(Ht,wt);
                h_njets_vs_ST[26]->Fill(nHadJets,ST,wt);
                h_njets_vs_HT[26]->Fill(nHadJets,Ht,wt);
                h_ST_vs_ptPho[26]->Fill(ST,bestPhoton.Pt(),wt);
                h_mvaResponse_baseline[26]->Fill(mvaValue,wt);
                h_Njets_Varbin[26]->Fill(nHadJets,wt);
                h_Nbjets_Varbin[26]->Fill(BTags,wt);
                h_MET_Varbin[26]->Fill(MET,wt);
                h_PhotonPt_Varbin[26]->Fill(bestPhoton.Pt(),wt);
                h_St_Varbin[26]->Fill(ST,wt);
                int Sbin = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL[26]->Fill(Sbin, wt);
                searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
                h_Sbins_LL_Validation[26]->Fill(searchBin,wt);
                Tfbins = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL_v1[26]->Fill(Tfbins,wt);
                h_TFbins_ElecLL_validation[26]->Fill(Sbin,wt);
                h_Sbins_LL[26]->Fill(searchBin,wt);

              }

	      }
	else if(Elec_passEtacut && !Elec_passpTcut)
	  {
	    //cout<<v_genEle2[0].Pt()<<"\t"<<abs(v_genEle2[0].Eta())<<"\t"<<v_genEle2[1].Pt()<<"\t"<<abs(v_genEle2[1].Eta())<<endl;
	    h_Njets[23]->Fill(nHadJets,wt);
            h_Nbjets[23]->Fill(BTags,wt);
            h_MET_[23]->Fill(MET,wt);
            h_PhotonPt[23]->Fill(bestPhoton.Pt(),wt);
            h_Mt_PhoMET[23]->Fill(mTPhoMET,wt);
            h_dPhi_PhoMET[23]->Fill(dPhi_PhoMET,wt);
            h_St[23]->Fill(ST,wt);
            h_HT[23]->Fill(Ht,wt);
            h_njets_vs_ST[23]->Fill(nHadJets,ST,wt);
            h_njets_vs_HT[23]->Fill(nHadJets,Ht,wt);
            h_ST_vs_ptPho[23]->Fill(ST,bestPhoton.Pt(),wt);
            h_mvaResponse_baseline[23]->Fill(mvaValue,wt);
            h_Njets_Varbin[23]->Fill(nHadJets,wt);
            h_Nbjets_Varbin[23]->Fill(BTags,wt);
            h_MET_Varbin[23]->Fill(MET,wt);
            h_PhotonPt_Varbin[23]->Fill(bestPhoton.Pt(),wt);
            h_St_Varbin[23]->Fill(ST,wt);
            int Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[23]->Fill(Sbin, wt);
            searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[23]->Fill(searchBin,wt);
            Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[23]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[23]->Fill(Sbin,wt);
            h_Sbins_LL[23]->Fill(searchBin,wt);
	    if(mvaValue>0.0)
              {
                h_Njets[27]->Fill(nHadJets,wt);
                h_Nbjets[27]->Fill(BTags,wt);
                h_MET_[27]->Fill(MET,wt);
                h_PhotonPt[27]->Fill(bestPhoton.Pt(),wt);
                h_Mt_PhoMET[27]->Fill(mTPhoMET,wt);
                h_dPhi_PhoMET[27]->Fill(dPhi_PhoMET,wt);
                h_St[27]->Fill(ST,wt);
                h_HT[27]->Fill(Ht,wt);
                h_njets_vs_ST[27]->Fill(nHadJets,ST,wt);
                h_njets_vs_HT[27]->Fill(nHadJets,Ht,wt);
                h_ST_vs_ptPho[27]->Fill(ST,bestPhoton.Pt(),wt);
                h_mvaResponse_baseline[27]->Fill(mvaValue,wt);
                h_Njets_Varbin[27]->Fill(nHadJets,wt);
                h_Nbjets_Varbin[27]->Fill(BTags,wt);
                h_MET_Varbin[27]->Fill(MET,wt);
                h_PhotonPt_Varbin[27]->Fill(bestPhoton.Pt(),wt);
                h_St_Varbin[27]->Fill(ST,wt);
                int Sbin = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL[27]->Fill(Sbin, wt);
                searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
                h_Sbins_LL_Validation[27]->Fill(searchBin,wt);
                Tfbins = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL_v1[27]->Fill(Tfbins,wt);
                h_TFbins_ElecLL_validation[27]->Fill(Sbin,wt);
                h_Sbins_LL[27]->Fill(searchBin,wt);

              }

	  }
	  }// not passign acceptance condition
	if(!Elec_passAccep)
	  {
	    // if(Debug && Elec_passId || Elec_passIso)
	    // cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
	    //	    cout<<v_genEle2[0].Pt()<<"\t"<<abs(v_genEle2[0].Eta())<<"\t"<<v_genEle2[1].Pt()<<"\t"<<abs(v_genEle2[1].Eta())<<endl;
	    nEvents_Selec[24]+=wt;
	    out_nEventsTags[24]="e-SR: fail Accep";
	    h_selectBaselineYields_->Fill("fail Accept e-SR",wt);
	    FailAccept_Elec+=wt;
	    h_Njets[4]->Fill(nHadJets,wt);
	    h_Nbjets[4]->Fill(BTags,wt);
	    h_MET_[4]->Fill(MET,wt);
	    h_PhotonPt[4]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[4]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[4]->Fill(dPhi_PhoMET,wt);
	    h_St[4]->Fill(ST,wt);
	    h_HT[4]->Fill(Ht,wt);
	    h_njets_vs_ST[4]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[4]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[4]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[4]->Fill(mvaValue,wt);
	    h_Njets_Varbin[4]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[4]->Fill(BTags,wt);
	    h_MET_Varbin[4]->Fill(MET,wt);
	    h_PhotonPt_Varbin[4]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[4]->Fill(ST,wt);
	    if(!s_data.Contains("data"))
	      {
		if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[4]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
		else if(v_genEle2.size()>1){
		  if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[4]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
		  else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[4]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
		}

		ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
		if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
		h_ratio_bestPho_genlep_v1[4]->Fill(ratio);
		h_ratio_bestPho_genlep_v2[4]->Fill(ratio1);
		mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
		h_mindr_Pho_genElec[4]->Fill(mindr_genElecPho,wt);
		h_mindr_Pho_genlep[4]->Fill(mindr_Pho_genlep);
		h_hasGenPromptPhoton_v1[4]->Fill(hasGenPromptPhoton);

	      }

	    int Sbin = getBinNoV1_le(BTags,nHadJets);	    
	    h_TFbins_ElecLL[4]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    //	    cout<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;
	    // cout<<searchBin<<"\t"<<BTags<<"\t"<<nHadJets<<"\t"<<MET<<endl;
	    h_Sbins_LL_Validation[4]->Fill(searchBin,wt);

	    Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[4]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[4]->Fill(Sbin,wt);
            h_Sbins_LL[4]->Fill(searchBin,wt);
	    if(Sbin<3)
	      {
		h_GenpT[6]->Fill(v_genEle2[0].Pt());
		h_GenEta[6]->Fill(v_genEle2[0].Eta());
		if(v_genEle2.size()>1)
		  {
		  h_GenpT[7]->Fill(v_genEle2[1].Pt());
		  h_GenEta[7]->Fill(v_genEle2[1].Eta());
		  }	      
	      }
	    else
	      {
                h_GenpT[8]->Fill(v_genEle2[0].Pt());
                h_GenEta[8]->Fill(v_genEle2[0].Eta());
                if(v_genEle2.size()>1)
                  {
		    h_GenpT[9]->Fill(v_genEle2[1].Pt());
		    h_GenEta[9]->Fill(v_genEle2[1].Eta());
                  }
              }

	    if(mvaValue>0.0)
	      {
		h_Njets[14]->Fill(nHadJets,wt);
		h_Nbjets[14]->Fill(BTags,wt);
		h_MET_[14]->Fill(MET,wt);
		h_PhotonPt[14]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[14]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[14]->Fill(dPhi_PhoMET,wt);
		h_St[14]->Fill(ST,wt);
		h_HT[14]->Fill(Ht,wt);
		h_njets_vs_ST[14]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[14]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[14]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[14]->Fill(mvaValue,wt);
		h_Njets_Varbin[14]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[14]->Fill(BTags,wt);
		h_MET_Varbin[14]->Fill(MET,wt);
		h_PhotonPt_Varbin[14]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[14]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[14]->Fill(Sbin, wt);
		nEvents_Selec[25]+=wt;
		out_nEventsTags[25]="e-SR: fail accept -BDT resp";
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[14]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[14]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[14]->Fill(Sbin,wt);
		h_Sbins_LL[14]->Fill(searchBin,wt);

	      }

	  }
	else if(Elec_passAccep && !Elec_passId) 
	  {
	    if(Debug && !Elec_passAccep|| Elec_passIso)
	      cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
	    //if((v_genEle2[0].Pt()>10 && abs(v_genEle2[0].Eta())<2.4) || (v_genEle2[1].Pt()>10 && abs(v_genEle2[1].Eta())<2.4))	      
	    //cout<<"Elec_passAccep"<<"\t"<<v_genEle2[0].Pt()<<"\t"<<abs(v_genEle2[0].Eta())<<"\t"<<v_genEle2[1].Pt()<<"\t"<<abs(v_genEle2[1].Eta())<<"\t"<<v_genEle2.size()<<endl;
	    FailId_Elec+=wt;
	    out_nEventsTags[26]="e-SR: fail ID";
	    nEvents_Selec[26]+=wt;

	    h_selectBaselineYields_->Fill("fail Id e-SR",wt);
	    h_Njets[5]->Fill(nHadJets,wt);
	    h_Nbjets[5]->Fill(BTags,wt);
	    h_MET_[5]->Fill(MET,wt);
	    h_PhotonPt[5]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[5]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[5]->Fill(dPhi_PhoMET,wt);
	    h_St[5]->Fill(ST,wt);
	    h_HT[5]->Fill(Ht,wt);
	    h_njets_vs_ST[5]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[5]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[5]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[5]->Fill(mvaValue,wt);
	    h_Njets_Varbin[5]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[5]->Fill(BTags,wt);
	    h_MET_Varbin[5]->Fill(MET,wt);
	    h_PhotonPt_Varbin[5]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[5]->Fill(ST,wt);
	    if(!s_data.Contains("data"))
              {
                if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[5]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
                else if(v_genEle2.size()>1){
                  if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[5]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
                  else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[5]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
                }

                ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
                if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
                h_ratio_bestPho_genlep_v1[5]->Fill(ratio);
                h_ratio_bestPho_genlep_v2[5]->Fill(ratio1);

                mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
                h_mindr_Pho_genElec[5]->Fill(mindr_genElecPho,wt);
                h_mindr_Pho_genlep[5]->Fill(mindr_Pho_genlep);
                h_hasGenPromptPhoton_v1[5]->Fill(hasGenPromptPhoton);

              }

	    int Sbin = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL[5]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[5]->Fill(searchBin,wt);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[5]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[5]->Fill(Sbin,wt);
            h_Sbins_LL[5]->Fill(searchBin,wt);

	    if(mvaValue>0.0)
	      {
		nEvents_Selec[27]+=wt;
		out_nEventsTags[27]="e-SR: fail ID -BDT resp";
		h_Njets[15]->Fill(nHadJets,wt);
		h_Nbjets[15]->Fill(BTags,wt);
		h_MET_[15]->Fill(MET,wt);
		h_PhotonPt[15]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[15]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[15]->Fill(dPhi_PhoMET,wt);
		h_St[15]->Fill(ST,wt);
		h_HT[15]->Fill(Ht,wt);
		h_njets_vs_ST[15]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[15]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[15]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[15]->Fill(mvaValue,wt);
		h_Njets_Varbin[15]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[15]->Fill(BTags,wt);
		h_MET_Varbin[15]->Fill(MET,wt);
		h_PhotonPt_Varbin[15]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[15]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[15]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[15]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[15]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[15]->Fill(Sbin,wt);
		h_Sbins_LL[15]->Fill(searchBin,wt);

	      }

	  }
	else if(Elec_passId && Elec_passAccep && !Elec_passIso)
	  {
	    // if(Debug && Elec_passAccep|| Elec_passId)
            //   cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
	    //	    cout<<"Elec_passAccep"<<"\t"<<v_genEle2[0].Pt()<<"\t"<<abs(v_genEle2[0].Eta())<<"\t"<<v_genEle2[1].Pt()<<"\t"<<abs(v_genEle2[1].Eta())<<"\t"<<v_genEle2.size()<<endl;
	    FailIso_Elec+=wt;
	    nEvents_Selec[28]+=wt;
	    out_nEventsTags[28]="e-SR: fail Iso";
	    h_selectBaselineYields_->Fill("pass Iso e-SR",wt);
	    h_Njets[6]->Fill(nHadJets,wt);
	    h_Nbjets[6]->Fill(BTags,wt);
	    h_MET_[6]->Fill(MET,wt);
	    h_PhotonPt[6]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[6]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[6]->Fill(dPhi_PhoMET,wt);
	    h_St[6]->Fill(ST,wt);
	    h_HT[6]->Fill(Ht,wt);
	    h_njets_vs_ST[6]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[6]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[6]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[6]->Fill(mvaValue,wt);
	    h_Njets_Varbin[6]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[6]->Fill(BTags,wt);
	    h_MET_Varbin[6]->Fill(MET,wt);
	    h_PhotonPt_Varbin[6]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[6]->Fill(ST,wt);
	    if(!s_data.Contains("data"))
              {
                if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[6]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
                else if(v_genEle2.size()>1){
                  if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[6]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
                  else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[6]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
                }

                ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
                if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
                h_ratio_bestPho_genlep_v1[6]->Fill(ratio);
                h_ratio_bestPho_genlep_v2[6]->Fill(ratio1);

                mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
                h_mindr_Pho_genElec[6]->Fill(mindr_genElecPho,wt);
                h_mindr_Pho_genlep[6]->Fill(mindr_Pho_genlep);
                h_hasGenPromptPhoton_v1[6]->Fill(hasGenPromptPhoton);

              }

	    int Sbin = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL[6]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[6]->Fill(searchBin,wt);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[6]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[6]->Fill(Sbin,wt);
            h_Sbins_LL[6]->Fill(searchBin,wt);

	    if(mvaValue>0.0)
	      {
		nEvents_Selec[29]+=wt;
		out_nEventsTags[29]="e-SR: fail ID-BDT resp";
		h_Njets[16]->Fill(nHadJets,wt);
		h_Nbjets[16]->Fill(BTags,wt);
		h_MET_[16]->Fill(MET,wt);
		h_PhotonPt[16]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[16]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[16]->Fill(dPhi_PhoMET,wt);
		h_St[16]->Fill(ST,wt);
		h_HT[16]->Fill(Ht,wt);
		h_njets_vs_ST[16]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[16]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[16]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[16]->Fill(mvaValue,wt);
		h_Njets_Varbin[16]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[16]->Fill(BTags,wt);
		h_MET_Varbin[16]->Fill(MET,wt);
		h_PhotonPt_Varbin[16]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[16]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[16]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[16]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[16]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[16]->Fill(Sbin,wt);
		h_Sbins_LL[16]->Fill(searchBin,wt);

	      }

	  }

	h_Njets[3]->Fill(nHadJets,wt);
	h_Nbjets[3]->Fill(BTags,wt);
	h_MET_[3]->Fill(MET,wt);
	h_PhotonPt[3]->Fill(bestPhoton.Pt(),wt);
	h_Mt_PhoMET[3]->Fill(mTPhoMET,wt);
	h_dPhi_PhoMET[3]->Fill(dPhi_PhoMET,wt);
	h_St[3]->Fill(ST,wt);
	h_HT[3]->Fill(Ht,wt);
	h_njets_vs_ST[3]->Fill(nHadJets,ST,wt);
	h_njets_vs_HT[3]->Fill(nHadJets,Ht,wt);
	h_ST_vs_ptPho[3]->Fill(ST,bestPhoton.Pt(),wt);
	if(!s_data.Contains("data"))
          {
            if(v_genEle2.size()==1) h2d_mindRvspT_bestPho_genlep_v1[3]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
            else if(v_genEle2.size()>1){
              if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[3]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
              else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[3]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);
            }

            ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
            if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
            h_ratio_bestPho_genlep_v1[3]->Fill(ratio);
            h_ratio_bestPho_genlep_v2[3]->Fill(ratio1);

            mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
            h_mindr_Pho_genElec[3]->Fill(mindr_genElecPho,wt);
            h_mindr_Pho_genlep[3]->Fill(mindr_Pho_genlep);
            h_hasGenPromptPhoton_v1[3]->Fill(hasGenPromptPhoton);

          }

	h_mvaResponse_baseline[3]->Fill(mvaValue,wt);
	int Sbin = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL[3]->Fill(Sbin, wt);
	//cout<<searchBin<<"\t"<<BTags<<"\t"<<nHadJets<<"\t"<<MET<<endl;

	//	cout<<"ll bins: SR"<<"\t"<<Sbin<<endl;
	// wt_LL = wt*h_TF->GetBinContent(Sbin);
	searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	//cout<<"e-sr"<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;
	h_Sbins_LL_Validation[3]->Fill(searchBin,wt);//_LL);
	Tfbins = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL_v1[3]->Fill(Tfbins,wt);
	h_TFbins_ElecLL_validation[3]->Fill(Sbin,wt);
	h_Sbins_LL[3]->Fill(searchBin,wt);

	h_Njets_Varbin[3]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[3]->Fill(BTags,wt);
	h_MET_Varbin[3]->Fill(MET,wt);
	h_PhotonPt_Varbin[3]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[3]->Fill(ST,wt);	
	if(mvaValue>0.0)
	  {
	    nEvents_Selec[30]+=wt;
	    out_nEventsTags[30]="e-SR: final with BDT";
	h_Njets[13]->Fill(nHadJets,wt);
        h_Nbjets[13]->Fill(BTags,wt);
        h_MET_[13]->Fill(MET,wt);
        h_PhotonPt[13]->Fill(bestPhoton.Pt(),wt);
        h_Mt_PhoMET[13]->Fill(mTPhoMET,wt);
        h_dPhi_PhoMET[13]->Fill(dPhi_PhoMET,wt);
        h_St[13]->Fill(ST,wt);
        h_HT[13]->Fill(Ht,wt);
        h_njets_vs_ST[13]->Fill(nHadJets,ST,wt);
        h_njets_vs_HT[13]->Fill(nHadJets,Ht,wt);
        h_ST_vs_ptPho[13]->Fill(ST,bestPhoton.Pt(),wt);
        h_mvaResponse_baseline[13]->Fill(mvaValue,wt);
	
	h_Njets_Varbin[13]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[13]->Fill(BTags,wt);
	h_MET_Varbin[13]->Fill(MET,wt);
	h_PhotonPt_Varbin[13]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[13]->Fill(ST,wt);
	Sbin = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL[13]->Fill(Sbin, wt);

        searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	h_Sbins_LL_Validation[13]->Fill(searchBin,wt);//_LL);  
	Tfbins = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL_v1[13]->Fill(Tfbins,wt);
	h_TFbins_ElecLL_validation[13]->Fill(Sbin,wt);
	h_Sbins_LL[13]->Fill(searchBin,wt);

	  }
      }
    if(Debug)
      cout<<"SR elec "<<"\t"<<jentry<<endl;
    if(!lost_elec_flag && nelec_reco ==0 && nmu_reco==1)
      {
	nEvents_Selec[31]+=wt;
	out_nEventsTags[31]="mu-CR: initial";
	if(isoElectronTracks !=0 || isoPionTracks!=0) continue;
	nEvents_Selec[32]+=wt;
	out_nEventsTags[32]="mu-CR: veto charged tracks";
	double dr2=bestPhoton.DeltaR((*Muons)[mu_index]);
        if(dr2<=0.2){
	  h_Njets_v1[26]->Fill(nHadJets,wt);
	  h_Nbjets_v1[26]->Fill(BTags,wt);
	  h_MET_v1[26]->Fill(MET,wt);
	  h_PhotonPt_v1[26]->Fill(bestPhoton.Pt(),wt);
	  h_St_v1[26]->Fill(ST,wt);
	  continue;	
	}
	nEvents_Selec[33]+=wt;
	out_nEventsTags[33]="mu-CR: dR(mu,#gamma)>0.2";
	ST = ST+(*Muons)[mu_index].Pt();
	double mTMuMET=sqrt(2*((*Muons)[mu_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,(*Muons)[mu_index].Phi()))));
	if(mTMuMET>100) continue;
	out_nEventsTags[34]="mu-CR final: mT<100";
	nEvents_Selec[34]+=wt;
	st = ST;
        mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");	
	if(hasGenPromptPhoton)
          {
            h_selectBaselineYields_CR->Fill("M-CR: gen prompt #gamma",wt);
            h_Njets_v1[27]->Fill(nHadJets,wt);
            h_Nbjets_v1[27]->Fill(BTags,wt);
            h_MET_v1[27]->Fill(MET,wt);
            h_PhotonPt_v1[27]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[27]->Fill(ST,wt);
          }
        else
          {
            h_selectBaselineYields_CR->Fill("e-CR:non-gen prompt #gamma",wt);
            h_Njets_v1[28]->Fill(nHadJets,wt);
            h_Nbjets_v1[28]->Fill(BTags,wt);
            h_MET_v1[28]->Fill(MET,wt);
            h_PhotonPt_v1[28]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[28]->Fill(ST,wt);

          }

	if(!s_data.Contains("data") && v_genEle2.size()==1)
	  {
	   h2d_mindRvspT_bestPho_genlep_v1[7]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
	  
	   ratio = v_genEle2[0].Pt()/bestPhoton.Pt();	  
	   h_ratio_bestPho_genlep_v1[7]->Fill(ratio);
	   //h_ratio_bestPho_genlep_v2[12]->Fill(ratio1);
	   mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
	   h_mindr_Pho_genElec[7]->Fill(mindr_genElecPho,wt);
	   h_mindr_Pho_genlep[7]->Fill(mindr_Pho_genlep);
	   h_hasGenPromptPhoton_v1[7]->Fill(hasGenPromptPhoton);
	  }
	if(!s_data.Contains("data"))
          {
            if(v_genEle2.size()==1){ h2d_mindRvspT_bestPho_genlep_v1[7]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
	    
	    // else if(v_genEle2.size()>1){
            //   if(bestPhoton.DeltaR(v_genEle2[0])<bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[7]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
            //   else if(bestPhoton.DeltaR(v_genEle2[0])>bestPhoton.DeltaR(v_genEle2[1])) h2d_mindRvspT_bestPho_genlep_v1[7]->Fill(bestPhoton.DeltaR(v_genEle2[1]),v_genEle2[1].Pt(),wt);       
	float ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
	//if(v_genEle2.size()>1) ratio1 = v_genEle2[1].Pt()/bestPhoton.Pt();
	h_ratio_bestPho_genlep_v1[7]->Fill(ratio);
	//h_ratio_bestPho_genlep_v2[7]->Fill(ratio1);
	float mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
	h_mindr_Pho_genElec[7]->Fill(mindr_genElecPho,wt);
	h_mindr_Pho_genlep[7]->Fill(mindr_Pho_genlep);
	h_hasGenPromptPhoton_v1[7]->Fill(hasGenPromptPhoton);
	    }
          }
      
	nCR_mu+=wt;
	h_selectBaselineYields_->Fill("Mu CR",wt); 
	h_Njets[7]->Fill(nHadJets,wt);
	h_Nbjets[7]->Fill(BTags,wt);
	h_MET_[7]->Fill(MET,wt);
	h_PhotonPt[7]->Fill(bestPhoton.Pt(),wt);
	h_Mt_PhoMET[7]->Fill(mTPhoMET,wt);
	h_dPhi_PhoMET[7]->Fill(dPhi_PhoMET,wt);
	h_St[7]->Fill(ST,wt);
	h_HT[7]->Fill(Ht,wt);
	h_njets_vs_ST[7]->Fill(nHadJets,ST,wt);
	h_njets_vs_HT[7]->Fill(nHadJets,Ht,wt);
	h_ST_vs_ptPho[7]->Fill(ST,bestPhoton.Pt(),wt);
	h_mvaResponse_baseline[7]->Fill(mvaValue,wt);
	h_Njets_Varbin[7]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[7]->Fill(BTags,wt);
	h_MET_Varbin[7]->Fill(MET,wt);
	h_PhotonPt_Varbin[7]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[7]->Fill(ST,wt);
	int Sbin = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL[7]->Fill(Sbin, wt);
	wt_LL = wt*h_TF->GetBinContent(Sbin+1);
        searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	h_Sbins_LL_Validation[7]->Fill(searchBin,wt_LL);  
	Tfbins = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL_v1[7]->Fill(Tfbins,wt);
	h_TFbins_ElecLL_validation[7]->Fill(Sbin,wt_LL);
	h_Sbins_LL[7]->Fill(searchBin,wt);
	h_Njets[42]->Fill(nHadJets,wt_LL);
        h_Nbjets[42]->Fill(BTags,wt_LL);
        h_MET_[42]->Fill(MET,wt_LL);
        h_PhotonPt[42]->Fill(bestPhoton.Pt(),wt_LL);
        h_St[42]->Fill(ST,wt_LL);
	h_Njets_v1[29]->Fill(nHadJets,wt);
	h_Nbjets_v1[29]->Fill(BTags,wt);
	h_MET_v1[29]->Fill(MET,wt);
	h_PhotonPt_v1[29]->Fill(bestPhoton.Pt(),wt);
	h_St_v1[29]->Fill(ST,wt);



	if(Debug)
	cout<<Sbin<<"\t"<<h_TF->GetBinContent(Sbin+1)<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<endl;

	if(mvaValue<0.0)
          {
	    out_nEventsTags[35]="mu-CR: BDT resp";
	    nEvents_Selec[35]+=wt;
	    h_Njets[17]->Fill(nHadJets,wt);
	    h_Nbjets[17]->Fill(BTags,wt);
	    h_MET_[17]->Fill(MET,wt);
	    h_PhotonPt[17]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[17]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[17]->Fill(dPhi_PhoMET,wt);
	    h_St[17]->Fill(ST,wt);
	    h_HT[17]->Fill(Ht,wt);
	    h_njets_vs_ST[17]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[17]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[17]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[17]->Fill(mvaValue,wt);
	    h_Njets_Varbin[17]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[17]->Fill(BTags,wt);
	    h_MET_Varbin[17]->Fill(MET,wt);
	    h_PhotonPt_Varbin[17]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[17]->Fill(ST,wt);
	    Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[17]->Fill(Sbin, wt);
	    wt_LL = wt*h_TF->GetBinContent(Sbin+1);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[17]->Fill(searchBin,wt_LL);
	    
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[17]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[17]->Fill(Sbin,wt_LL);
	    h_Sbins_LL[17]->Fill(searchBin,wt_LL);
	    h_Njets[43]->Fill(nHadJets,wt_LL);
	    h_Nbjets[43]->Fill(BTags,wt_LL);
	    h_MET_[43]->Fill(MET,wt_LL);
	    h_PhotonPt[43]->Fill(bestPhoton.Pt(),wt_LL);
	    h_St[43]->Fill(ST,wt_LL);

          }

      }
    if(Debug)
      cout<<"CR mu "<<"\t"<<jentry<<endl;

    if(!lost_elec_flag && nelec_reco ==0 && nmu_reco==0)
      {
	//nSR_mu+=wt;
	nEvents_Selec[36]+=wt;
	out_nEventsTags[36]="mu-SR: initial";
        if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
	out_nEventsTags[48]="mu-SR: veto charged tracks";
	nEvents_Selec[48]+=wt;
	// if(!s_data.Contains("data"))
        //   {
	// if(v_genPho1.size()!=0)
        //   {
        //     double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);
        //     if(mindR_genPho>=0.2)// || v_genPho1.size()==0)                                                                                  
        //       {
        //         double mindrPho_genlep=getGenLep(bestPhoton);//, true);                                                                      
        //         if(mindrPho_genlep<0.2 || (v_genMu2.size()>=1 && v_genMu2[0].Pt()/bestPhoton.Pt()>0.8 && v_genMu2[0].Pt()/bestPhoton.Pt()<1.2)||(v_genMu2.size()>=1 && v_genMu2[1].Pt()/bestPhoton.Pt()>0.8 && v_genMu2[1].Pt()/bestPhoton.Pt()<1.2)) {continue; nfakeRatePho+=wt;}
        //         h_selectBaselineYields_->Fill("prompt Photon Checks-mu",wt);
        //       }
        //   }
        // else
        //   {
        //     double mindrPho_genlep=getGenLep(bestPhoton);//, true);                                                                          
        //     if(mindrPho_genlep<0.2 || (v_genMu2.size()>=1 && v_genMu2[0].Pt()/bestPhoton.Pt()>0.8 && v_genMu2[0].Pt()/bestPhoton.Pt()<1.2) ||(v_genMu2.size()>=1 && v_genMu2[1].Pt()/bestPhoton.Pt()>0.8 && v_genMu2[1].Pt()/bestPhoton.Pt()<1.2)) continue;
	//     h_selectBaselineYields_->Fill("prompt Photon Checks-mu1",wt);

        //   }
      
            // double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);
            // if(mindR_genPho>=0.2)
            //   {
            //     double mindrPho_genlep=getGenLep(bestPhoton);//, false);
            //     if(mindrPho_genlep<0.2) {continue; nfakeRatePho+=wt;}
            //     h_selectBaselineYields_->Fill("prompt Photon Checks",wt);
            //   }
	    nEvents_Selec[37]+=wt;
	    out_nEventsTags[37]="mu-SR: prompt #gamma checks";

          // }
	  //  if(nGentau_had1>0 && nGenMu1==0 && nGenEle1==0) continue;
	if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue;//to reject W->qq' type of events
	//survived_vetohad+=wt;
	nEvents_Selec[38]+=wt;
	out_nEventsTags[38]="mu-SR: w->qq events";
	//if(nGenEle1>1) continue;
	//	if(nGenEle1==1 && nGenMu1==1) continue;
	if(nGenEle1>1) continue;
	
	nEvents_Selec[39]+=wt;
	out_nEventsTags[39]="mu-SR: gen e=0  or gen e=1 && gen mu=1";
	//l	if(nGenEle1>0) continue;
	h_selectBaselineYields_->Fill("Mu SR",wt);  
	survived_elecge1+=wt;
	if(hasGenPromptPhoton)
          {
            h_Njets_v1[30]->Fill(nHadJets,wt);
            h_Nbjets_v1[30]->Fill(BTags,wt);
            h_MET_v1[30]->Fill(MET,wt);
            h_PhotonPt_v1[30]->Fill(bestPhoton.Pt(),wt);
            h_St_v1[30]->Fill(ST,wt);

            double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);                                                                   
            if(mindR_genPho<=0.2)
              {
                h_selectBaselineYields_SR->Fill("e-SR: dR(reco-#gamma,gen-#gamma)<0.2",wt);
                h_Njets_v1[31]->Fill(nHadJets,wt);
                h_Nbjets_v1[31]->Fill(BTags,wt);
                h_MET_v1[31]->Fill(MET,wt);
                h_PhotonPt_v1[31]->Fill(bestPhoton.Pt(),wt);
                h_St_v1[31]->Fill(ST,wt);
              }
            else if(mindR_genPho>0.2)
              {
                if(v_genPho1.size()==0)
                  {
                    h_Njets_v1[32]->Fill(nHadJets,wt);
                    h_Nbjets_v1[32]->Fill(BTags,wt);
                    h_MET_v1[32]->Fill(MET,wt);
                    h_PhotonPt_v1[32]->Fill(bestPhoton.Pt(),wt);
                    h_St_v1[32]->Fill(ST,wt);
                  }
                else
                  {
                    h_Njets_v1[33]->Fill(nHadJets,wt);
                    h_Nbjets_v1[33]->Fill(BTags,wt);
                    h_MET_v1[33]->Fill(MET,wt);
                    h_PhotonPt_v1[33]->Fill(bestPhoton.Pt(),wt);
                    h_St_v1[33]->Fill(ST,wt);
                  }
                h_selectBaselineYields_SR->Fill("e-SR: !dR(reco-#gamma,gen-#gamma)<0.2",wt);
		if(v_genEle2.size()==1)
                  {
                    if(bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[0].Pt()/bestPhoton.Pt()<0.79)))
                      {
                        h_Njets_v1[34]->Fill(nHadJets,wt);
                        h_Nbjets_v1[34]->Fill(BTags,wt);
                        h_MET_v1[34]->Fill(MET,wt);
                        h_PhotonPt_v1[34]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[34]->Fill(ST,wt);
                        h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR",wt);
                        int flag=-1,pdgID=-1;
			std::vector<int> array1;
                        array1 = dR_recoPho_GenParticle(bestPhoton);
                        pdgID = array1[0];
                        flag=array1[1];
			if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
                          {
                            h_Njets_v1[35]->Fill(nHadJets,wt);
                            h_Nbjets_v1[35]->Fill(BTags,wt);
                            h_MET_v1[35]->Fill(MET,wt);
                            h_PhotonPt_v1[35]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[35]->Fill(ST,wt);
                          }
                        else if(!flag)
                          {
                            h_Njets_v1[36]->Fill(nHadJets,wt);
                            h_Nbjets_v1[36]->Fill(BTags,wt);
                            h_MET_v1[36]->Fill(MET,wt);
                            h_PhotonPt_v1[36]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[36]->Fill(ST,wt);
                          }
		      }
                    else if(bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))
                      {
                        h_Njets_v1[37]->Fill(nHadJets,wt);
                        h_Nbjets_v1[37]->Fill(BTags,wt);
                        h_MET_v1[37]->Fill(MET,wt);
                        h_PhotonPt_v1[37]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[37]->Fill(ST,wt);
                        h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && FR",wt);                            
			continue;
                      }
		    else
                      {
                        h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && !dR(tau/q/g,reco #gamma)<0.2",wt);
                        int flag=-1,pdgID=-1;
			std::vector<int> array1;
                        array1 = dR_recoPho_GenParticle(bestPhoton);
                        pdgID = array1[0];
                        flag=array1[1];
			if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
                          {
                            h_Njets_v1[35]->Fill(nHadJets,wt);
                            h_Nbjets_v1[35]->Fill(BTags,wt);
                            h_MET_v1[35]->Fill(MET,wt);
                            h_PhotonPt_v1[35]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[35]->Fill(ST,wt);
                          }
                        else if(!flag)
                          {
                            h_Njets_v1[36]->Fill(nHadJets,wt);
                            h_Nbjets_v1[36]->Fill(BTags,wt);
                            h_MET_v1[36]->Fill(MET,wt);
                            h_PhotonPt_v1[36]->Fill(bestPhoton.Pt(),wt);
                            h_St_v1[36]->Fill(ST,wt);
                          }

		      }
		  }
		if(nGentau_had1==1 && bestPhoton.DeltaR(v_genTau2[0])<0.2)
		  {
		    h_Njets_v1[38]->Fill(nHadJets,wt);
		    h_Nbjets_v1[38]->Fill(BTags,wt);
		    h_MET_v1[38]->Fill(MET,wt);
		    h_PhotonPt_v1[38]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[38]->Fill(ST,wt);
		  }
		if(nGentau_had1>1 && (bestPhoton.DeltaR(v_genTau2[0])<0.2 || bestPhoton.DeltaR(v_genTau2[1])<0.2))
		  {
		    h_Njets_v1[39]->Fill(nHadJets,wt);
                    h_Nbjets_v1[39]->Fill(BTags,wt);
                    h_MET_v1[39]->Fill(MET,wt);
                    h_PhotonPt_v1[39]->Fill(bestPhoton.Pt(),wt);
                    h_St_v1[39]->Fill(ST,wt);

		  }
	      }
	  }
	if(!hasGenPromptPhoton)
	  {
	    h_Njets_v1[40]->Fill(nHadJets,wt);
	    h_Nbjets_v1[40]->Fill(BTags,wt);
	    h_MET_v1[40]->Fill(MET,wt);
	    h_PhotonPt_v1[40]->Fill(bestPhoton.Pt(),wt);
	    h_St_v1[40]->Fill(ST,wt);


	    if(v_genEle2.size()==1)
	      {
		if(bestPhoton.DeltaR(v_genEle2[0])>0.2 || (bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()>1.21) || v_genEle2[0].Pt()/bestPhoton.Pt()<0.79)))
		  {
		    h_Njets_v1[41]->Fill(nHadJets,wt);
		    h_Nbjets_v1[41]->Fill(BTags,wt);
		    h_MET_v1[41]->Fill(MET,wt);
		    h_PhotonPt_v1[41]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[41]->Fill(ST,wt);
		    h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR",wt);
		    int flag=-1,pdgID=-1;
		    std::vector<int> array1;
		    array1 = dR_recoPho_GenParticle(bestPhoton);
		    pdgID = array1[0];
		    flag=array1[1];
		    if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
		      {
			h_Njets_v1[42]->Fill(nHadJets,wt);
			h_Nbjets_v1[42]->Fill(BTags,wt);
			h_MET_v1[42]->Fill(MET,wt);
			h_PhotonPt_v1[42]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[42]->Fill(ST,wt);
		      }
		    else if(!flag)
		      {
			h_Njets_v1[43]->Fill(nHadJets,wt);
			h_Nbjets_v1[43]->Fill(BTags,wt);
			h_MET_v1[43]->Fill(MET,wt);
			h_PhotonPt_v1[43]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[43]->Fill(ST,wt);
		      }
		  }
		else if(bestPhoton.DeltaR(v_genEle2[0])<0.2 && ((v_genEle2[0].Pt()/bestPhoton.Pt()<1.21) && v_genEle2[0].Pt()/bestPhoton.Pt()>0.79))
		  {
		    h_Njets_v1[44]->Fill(nHadJets,wt);
		    h_Nbjets_v1[44]->Fill(BTags,wt);
		    h_MET_v1[44]->Fill(MET,wt);
		    h_PhotonPt_v1[44]->Fill(bestPhoton.Pt(),wt);
		    h_St_v1[44]->Fill(ST,wt);
		    h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && FR",wt);
		    continue;
		  }
		else
		  {
		    h_selectBaselineYields_SR->Fill("e-SR:!dR(reco-#gamma,gen-#gamma)<0.2 && !FR && !dR(tau/q/g,reco #gamma)<0.2",wt);
		    int flag=-1,pdgID=-1;
		    std::vector<int> array1;
		    array1 = dR_recoPho_GenParticle(bestPhoton);
		    pdgID = array1[0];
		    flag=array1[1];
		    if((abs(pdgID)<6 || abs(pdgID)==21) && flag==1)
                      {
                        h_Njets_v1[42]->Fill(nHadJets,wt);
                        h_Nbjets_v1[42]->Fill(BTags,wt);
                        h_MET_v1[42]->Fill(MET,wt);
                        h_PhotonPt_v1[42]->Fill(bestPhoton.Pt(),wt);
                        h_St_v1[42]->Fill(ST,wt);
                      }
                    else if(!flag)
                      {
                        h_Njets_v1[43]->Fill(nHadJets,wt);
			h_Nbjets_v1[43]->Fill(BTags,wt);
			h_MET_v1[43]->Fill(MET,wt);
                        h_PhotonPt_v1[43]->Fill(bestPhoton.Pt(),wt);
			h_St_v1[43]->Fill(ST,wt);
                      }
		  }
	      }
	    if(nGentau_had1==1 && bestPhoton.DeltaR(v_genTau2[0])<0.2)
	      {
		h_Njets_v1[45]->Fill(nHadJets,wt);
		h_Nbjets_v1[45]->Fill(BTags,wt);
		h_MET_v1[45]->Fill(MET,wt);
		h_PhotonPt_v1[45]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[45]->Fill(ST,wt);
	      }
	    if(nGentau_had1>1 && (bestPhoton.DeltaR(v_genTau2[0])<0.2 || bestPhoton.DeltaR(v_genTau2[1])<0.2))
	      {
		h_Njets_v1[46]->Fill(nHadJets,wt);
		h_Nbjets_v1[46]->Fill(BTags,wt);
		h_MET_v1[46]->Fill(MET,wt);
		h_PhotonPt_v1[46]->Fill(bestPhoton.Pt(),wt);
		h_St_v1[46]->Fill(ST,wt);

	      }
	  }
      
	//if(nGentau_had1>1) continue;
	
	//if(v_genMu2.size()==0) {TLorentzVector v1;v_genMu2.push_back(v1);}
	//if(nGenEle1>0) continue;
	//if(v_genMu2.size()==0)//nGentau_had1>1 || (v_genMu2.size()==0 && nGentau_had1>=1))                                                 	
	//if(v_genMu2.size()!=0 && nGentau_had1>=1) cout<<"entry1:"<<"\t"<<jentry<<"\t"<<v_genMu2.size()<<"\t"<<nGentau_had1<<endl;
	nEvents_Selec[40]+=wt;
	out_nEventsTags[40]="mu-SR: gen e=0";
        //if(nGenMu1==0) continue;
	nSR_mu+=wt;
	sortTLorVec(&v_genMu2);

	if(!s_data.Contains("data"))
	  {
            if(v_genEle2.size()==1){ h2d_mindRvspT_bestPho_genlep_v1[8]->Fill(bestPhoton.DeltaR(v_genEle2[0]),v_genEle2[0].Pt()/bestPhoton.Pt(),wt);
    float ratio = v_genEle2[0].Pt()/bestPhoton.Pt();
    h_ratio_bestPho_genlep_v1[8]->Fill(ratio);
    float mindr_genElecPho = MinDr(bestPhoton,v_genEle2);
    h_mindr_Pho_genElec[8]->Fill(mindr_genElecPho,wt);
    h_mindr_Pho_genlep[8]->Fill(mindr_Pho_genlep);
    h_hasGenPromptPhoton_v1[8]->Fill(hasGenPromptPhoton);
	    }
	  }
//v_genEle1=v_genMu2;
	h_Njets[8]->Fill(nHadJets,wt);
	h_Nbjets[8]->Fill(BTags,wt);
	h_MET_[8]->Fill(MET,wt);
	h_PhotonPt[8]->Fill(bestPhoton.Pt(),wt);
	h_Mt_PhoMET[8]->Fill(mTPhoMET,wt);
	h_dPhi_PhoMET[8]->Fill(dPhi_PhoMET,wt);
	h_St[8]->Fill(ST,wt);
	h_HT[8]->Fill(Ht,wt);
	h_njets_vs_ST[8]->Fill(nHadJets,ST,wt);
	h_njets_vs_HT[8]->Fill(nHadJets,Ht,wt);
	h_ST_vs_ptPho[8]->Fill(ST,bestPhoton.Pt(),wt);
	h_mvaResponse_baseline[8]->Fill(mvaValue,wt);
	h_Njets_Varbin[8]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[8]->Fill(BTags,wt);
	h_MET_Varbin[8]->Fill(MET,wt);
	h_PhotonPt_Varbin[8]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[8]->Fill(ST,wt);
	int Sbin = getBinNoV1_le(BTags,nHadJets);
	h_TFbins_ElecLL[8]->Fill(Sbin, wt);
	searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	h_Sbins_LL_Validation[8]->Fill(searchBin,wt);
	Tfbins = getBinNoV1_le(BTags,nHadJets);
        h_TFbins_ElecLL_v1[8]->Fill(Tfbins,wt);
        h_TFbins_ElecLL_validation[8]->Fill(Sbin,wt);
        h_Sbins_LL[8]->Fill(searchBin,wt);
	
	if(mvaValue>0.0)
          {
	    nEvents_Selec[41]+=wt;
	    out_nEventsTags[41]="mu-SR: BDt resp";

            h_Njets[18]->Fill(nHadJets,wt);
            h_Nbjets[18]->Fill(BTags,wt);
            h_MET_[18]->Fill(MET,wt);
            h_PhotonPt[18]->Fill(bestPhoton.Pt(),wt);
            h_Mt_PhoMET[18]->Fill(mTPhoMET,wt);
            h_dPhi_PhoMET[18]->Fill(dPhi_PhoMET,wt);
            h_St[18]->Fill(ST,wt);
            h_HT[18]->Fill(Ht,wt);
            h_njets_vs_ST[18]->Fill(nHadJets,ST,wt);
            h_njets_vs_HT[18]->Fill(nHadJets,Ht,wt);
            h_ST_vs_ptPho[18]->Fill(ST,bestPhoton.Pt(),wt);
            h_mvaResponse_baseline[18]->Fill(mvaValue,wt);
	    h_Njets_Varbin[18]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[18]->Fill(BTags,wt);
	    h_MET_Varbin[18]->Fill(MET,wt);
	    h_PhotonPt_Varbin[18]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[18]->Fill(ST,wt);
	    Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[18]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[18]->Fill(searchBin,wt);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[18]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[18]->Fill(Sbin,wt);
	    h_Sbins_LL[18]->Fill(searchBin,wt);


          }
	if(!Mu_passAccep){
	  if(!Mu_passEtacut) // && v_genMu2.size()>0 && nGentau_had1<=1)
          {
	    if(!(nGentau_had1>1 || (v_genMu2.size()==0 && nGentau_had1>=1)))
	      {
            h_Njets[24]->Fill(nHadJets,wt);
            h_Nbjets[24]->Fill(BTags,wt);
            h_MET_[24]->Fill(MET,wt);
            h_PhotonPt[24]->Fill(bestPhoton.Pt(),wt);
            h_Mt_PhoMET[24]->Fill(mTPhoMET,wt);
            h_dPhi_PhoMET[24]->Fill(dPhi_PhoMET,wt);
            h_St[24]->Fill(ST,wt);
            h_HT[24]->Fill(Ht,wt);
            h_njets_vs_ST[24]->Fill(nHadJets,ST,wt);
            h_njets_vs_HT[24]->Fill(nHadJets,Ht,wt);
            h_ST_vs_ptPho[24]->Fill(ST,bestPhoton.Pt(),wt);
            h_mvaResponse_baseline[24]->Fill(mvaValue,wt);
            h_Njets_Varbin[24]->Fill(nHadJets,wt);
            h_Nbjets_Varbin[24]->Fill(BTags,wt);
            h_MET_Varbin[24]->Fill(MET,wt);
            h_PhotonPt_Varbin[24]->Fill(bestPhoton.Pt(),wt);
            h_St_Varbin[24]->Fill(ST,wt);
            int Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[24]->Fill(Sbin, wt);
            searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[24]->Fill(searchBin,wt);
            Tfbins = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL_v1[24]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[24]->Fill(Sbin,wt);
            h_Sbins_LL[24]->Fill(searchBin,wt);
	    if(mvaValue>0.0)
	      {
		h_Njets[28]->Fill(nHadJets,wt);
		h_Nbjets[28]->Fill(BTags,wt);
		h_MET_[28]->Fill(MET,wt);
		h_PhotonPt[28]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[28]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[28]->Fill(dPhi_PhoMET,wt);
		h_St[28]->Fill(ST,wt);
		h_HT[28]->Fill(Ht,wt);
		h_njets_vs_ST[28]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[28]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[28]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[28]->Fill(mvaValue,wt);
		h_Njets_Varbin[28]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[28]->Fill(BTags,wt);
		h_MET_Varbin[28]->Fill(MET,wt);
		h_PhotonPt_Varbin[28]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[28]->Fill(ST,wt);
		int Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[28]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[28]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[28]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[28]->Fill(Sbin,wt);
		h_Sbins_LL[28]->Fill(searchBin,wt);

	      }
	      }
	    if(nGentau_had1>1 || (v_genMu2.size()==0 && nGentau_had1>=1))
	      {
		//cout<<"entry:"<<"\t"<<jentry<<"\t"<<v_genMu2.size()<<"\t"<<nGentau_had1<<endl;
		h_Njets[48]->Fill(nHadJets,wt);
		h_Nbjets[48]->Fill(BTags,wt);
		h_MET_[48]->Fill(MET,wt);
		h_PhotonPt[48]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[48]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[48]->Fill(dPhi_PhoMET,wt);
		h_St[48]->Fill(ST,wt);
		h_HT[48]->Fill(Ht,wt);
		h_njets_vs_ST[48]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[48]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[48]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[48]->Fill(mvaValue,wt);
		h_Njets_Varbin[48]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[48]->Fill(BTags,wt);
		h_MET_Varbin[48]->Fill(MET,wt);
		h_PhotonPt_Varbin[48]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[48]->Fill(ST,wt);
		int Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[48]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[48]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[48]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[48]->Fill(Sbin,wt);
		h_Sbins_LL[48]->Fill(searchBin,wt);
		if(mvaValue>0.0)
		  {
		    h_Njets[49]->Fill(nHadJets,wt);
		    h_Nbjets[49]->Fill(BTags,wt);
		    h_MET_[49]->Fill(MET,wt);
		    h_PhotonPt[49]->Fill(bestPhoton.Pt(),wt);
		    h_Mt_PhoMET[49]->Fill(mTPhoMET,wt);
		    h_dPhi_PhoMET[49]->Fill(dPhi_PhoMET,wt);
		    h_St[49]->Fill(ST,wt);
		    h_HT[49]->Fill(Ht,wt);
		    h_njets_vs_ST[49]->Fill(nHadJets,ST,wt);
		    h_njets_vs_HT[49]->Fill(nHadJets,Ht,wt);
		    h_ST_vs_ptPho[49]->Fill(ST,bestPhoton.Pt(),wt);
		    h_mvaResponse_baseline[49]->Fill(mvaValue,wt);
		    h_Njets_Varbin[49]->Fill(nHadJets,wt);
		    h_Nbjets_Varbin[49]->Fill(BTags,wt);
		    h_MET_Varbin[49]->Fill(MET,wt);
		    h_PhotonPt_Varbin[49]->Fill(bestPhoton.Pt(),wt);
		    h_St_Varbin[49]->Fill(ST,wt);
		    int Sbin = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL[49]->Fill(Sbin, wt);
		    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		    h_Sbins_LL_Validation[49]->Fill(searchBin,wt);
		    Tfbins = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL_v1[49]->Fill(Tfbins,wt);
		    h_TFbins_ElecLL_validation[49]->Fill(Sbin,wt);
		    h_Sbins_LL[49]->Fill(searchBin,wt);

		  }

	      }
          }
	  else if(Mu_passEtacut && !Mu_passpTcut)// && v_genMu2.size()>0 && nGentau_had1<=1) // && !nGentau_had1)
          {
            h_Njets[25]->Fill(nHadJets,wt);
            h_Nbjets[25]->Fill(BTags,wt);
            h_MET_[25]->Fill(MET,wt);
            h_PhotonPt[25]->Fill(bestPhoton.Pt(),wt);
            h_Mt_PhoMET[25]->Fill(mTPhoMET,wt);
            h_dPhi_PhoMET[25]->Fill(dPhi_PhoMET,wt);
            h_St[25]->Fill(ST,wt);
            h_HT[25]->Fill(Ht,wt);
            h_njets_vs_ST[25]->Fill(nHadJets,ST,wt);
            h_njets_vs_HT[25]->Fill(nHadJets,Ht,wt);
            h_ST_vs_ptPho[25]->Fill(ST,bestPhoton.Pt(),wt);
            h_mvaResponse_baseline[25]->Fill(mvaValue,wt);
            h_Njets_Varbin[25]->Fill(nHadJets,wt);
            h_Nbjets_Varbin[25]->Fill(BTags,wt);
            h_MET_Varbin[25]->Fill(MET,wt);
            h_PhotonPt_Varbin[25]->Fill(bestPhoton.Pt(),wt);
            h_St_Varbin[25]->Fill(ST,wt);
            int Sbin = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL[25]->Fill(Sbin, wt);
            searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[25]->Fill(searchBin,wt);
            Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[25]->Fill(Tfbins,wt);
            h_TFbins_ElecLL_validation[25]->Fill(Sbin,wt);
            h_Sbins_LL[25]->Fill(searchBin,wt);
	    
	    if(mvaValue>0.0)
              {
                h_Njets[29]->Fill(nHadJets,wt);
                h_Nbjets[29]->Fill(BTags,wt);
                h_MET_[29]->Fill(MET,wt);
                h_PhotonPt[29]->Fill(bestPhoton.Pt(),wt);
                h_Mt_PhoMET[29]->Fill(mTPhoMET,wt);
                h_dPhi_PhoMET[29]->Fill(dPhi_PhoMET,wt);
                h_St[29]->Fill(ST,wt);
                h_HT[29]->Fill(Ht,wt);
                h_njets_vs_ST[29]->Fill(nHadJets,ST,wt);
                h_njets_vs_HT[29]->Fill(nHadJets,Ht,wt);
                h_ST_vs_ptPho[29]->Fill(ST,bestPhoton.Pt(),wt);
                h_mvaResponse_baseline[29]->Fill(mvaValue,wt);
                h_Njets_Varbin[29]->Fill(nHadJets,wt);
                h_Nbjets_Varbin[29]->Fill(BTags,wt);
                h_MET_Varbin[29]->Fill(MET,wt);
                h_PhotonPt_Varbin[29]->Fill(bestPhoton.Pt(),wt);
                h_St_Varbin[29]->Fill(ST,wt);
                int Sbin = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL[29]->Fill(Sbin, wt);
                searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
                h_Sbins_LL_Validation[29]->Fill(searchBin,wt);
                Tfbins = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL_v1[29]->Fill(Tfbins,wt);
                h_TFbins_ElecLL_validation[29]->Fill(Sbin,wt);
                h_Sbins_LL[29]->Fill(searchBin,wt);

              }

	  }
	}
	if(!Mu_passAccep) //&& !Mu_passId && ! Mu_passIso)
	  {
	    nEvents_Selec[42]+=wt;
	    out_nEventsTags[42]="mu-SR: failed accept";//
	    //if(((v_genMu2[0].Pt()>10 && abs(v_genMu2[0].Eta())<2.4) || (v_genMu2.size()==2 && v_genMu2[1].Pt()>10 && abs(v_genMu2[1].Eta())<2.4)))
	    // cout<<"entry:"<<"\t"<<jentry<<"\t"<<"Mu_passAccep"<<"\t"<<v_genMu2[0].Pt()<<"\t"<<abs(v_genMu2[0].Eta())<<"\t"<<v_genMu2.size()<<"\t"<<v_genTau2.size()<<"\t"<<nGentau_had1<<endl;//"\t"<<v_genMu2[1].Pt()<<"\t"<<abs(v_genMu[1].Eta())<<endl;
	     // if(v_genMu2.size()>1)
	     //   cout<<"size ==2"<<"\t"<<v_genMu2[1].Pt()<<"\t"<<abs(v_genMu2[1].Eta())<<"\t"<<v_genMu2.size()<<"\t"<<v_genTau2.size()<<endl;
	    h_selectBaselineYields_->Fill("fail Accept mu-SR",wt);	    
	    FailAccept_Mu+=wt;
	    h_Njets[9]->Fill(nHadJets,wt);
	    h_Nbjets[9]->Fill(BTags,wt);
	    h_MET_[9]->Fill(MET,wt);
	    h_PhotonPt[9]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[9]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[9]->Fill(dPhi_PhoMET,wt);
	    h_St[9]->Fill(ST,wt);
	    h_HT[9]->Fill(Ht,wt);
	    h_njets_vs_ST[9]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[9]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[9]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[9]->Fill(mvaValue,wt);
	    h_Njets_Varbin[9]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[9]->Fill(BTags,wt);
	    h_MET_Varbin[9]->Fill(MET,wt);
	    h_PhotonPt_Varbin[9]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[9]->Fill(ST,wt);
	    Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[9]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
            h_Sbins_LL_Validation[9]->Fill(searchBin,wt);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[9]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[9]->Fill(Sbin,wt);
	    h_Sbins_LL[9]->Fill(searchBin,wt);

	    if(nGentau_had1>=1 && v_genMu2.size()==0)
	      {
		h_Njets[44]->Fill(nHadJets,wt);
		h_Nbjets[44]->Fill(BTags,wt);
		h_MET_[44]->Fill(MET,wt);
		h_PhotonPt[44]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[44]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[44]->Fill(dPhi_PhoMET,wt);
		h_St[44]->Fill(ST,wt);
		h_HT[44]->Fill(Ht,wt);
		h_njets_vs_ST[44]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[44]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[44]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[44]->Fill(mvaValue,wt);
		h_Njets_Varbin[44]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[44]->Fill(BTags,wt);
		h_MET_Varbin[44]->Fill(MET,wt);
		h_PhotonPt_Varbin[44]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[44]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[44]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[44]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[44]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[44]->Fill(Sbin,wt);
		h_Sbins_LL[44]->Fill(searchBin,wt);
	      }
	    else //if (!nGentau_had1)
              {
                h_Njets[46]->Fill(nHadJets,wt);
                h_Nbjets[46]->Fill(BTags,wt);
                h_MET_[46]->Fill(MET,wt);
                h_PhotonPt[46]->Fill(bestPhoton.Pt(),wt);
                h_Mt_PhoMET[46]->Fill(mTPhoMET,wt);
                h_dPhi_PhoMET[46]->Fill(dPhi_PhoMET,wt);
                h_St[46]->Fill(ST,wt);
                h_HT[46]->Fill(Ht,wt);
                h_njets_vs_ST[46]->Fill(nHadJets,ST,wt);
                h_njets_vs_HT[46]->Fill(nHadJets,Ht,wt);
                h_ST_vs_ptPho[46]->Fill(ST,bestPhoton.Pt(),wt);
                h_mvaResponse_baseline[46]->Fill(mvaValue,wt);
                h_Njets_Varbin[46]->Fill(nHadJets,wt);
                h_Nbjets_Varbin[46]->Fill(BTags,wt);
                h_MET_Varbin[46]->Fill(MET,wt);
                h_PhotonPt_Varbin[46]->Fill(bestPhoton.Pt(),wt);
                h_St_Varbin[46]->Fill(ST,wt);
                Sbin = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL[46]->Fill(Sbin, wt);
                searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
                h_Sbins_LL_Validation[46]->Fill(searchBin,wt);
                Tfbins = getBinNoV1_le(BTags,nHadJets);
                h_TFbins_ElecLL_v1[46]->Fill(Tfbins,wt);
                h_TFbins_ElecLL_validation[46]->Fill(Sbin,wt);
                h_Sbins_LL[46]->Fill(searchBin,wt);
	      }
		if(Sbin<3)
		  {
		    if(v_genMu2.size()==1)
		      {
			h_GenpT[10]->Fill(v_genMu2[0].Pt());
			h_GenEta[10]->Fill(v_genMu2[0].Eta());
		      }
		    else if(v_genMu2.size()>1)
		      {
			h_GenpT[11]->Fill(v_genMu2[1].Pt());
			h_GenEta[11]->Fill(v_genMu2[1].Eta());
		      }
		  }
		else
		  {
		    if(v_genMu2.size()==1)
		      {
		    h_GenpT[12]->Fill(v_genMu2[0].Pt());
		    h_GenEta[12]->Fill(v_genMu2[0].Eta());
		      }
		    if(v_genMu2.size()>1)
		      {
			h_GenpT[13]->Fill(v_genMu2[1].Pt());
			h_GenEta[13]->Fill(v_genMu2[1].Eta());
		      }
		  }


              
	      


	    if(mvaValue>0.0)
	      {
		nEvents_Selec[43]+=wt;
		out_nEventsTags[43]="mu-SR: faild accept - BDT resp";

		h_Njets[19]->Fill(nHadJets,wt);
		h_Nbjets[19]->Fill(BTags,wt);
		h_MET_[19]->Fill(MET,wt);
		h_PhotonPt[19]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[19]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[19]->Fill(dPhi_PhoMET,wt);
		h_St[19]->Fill(ST,wt);
		h_HT[19]->Fill(Ht,wt);
		h_njets_vs_ST[19]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[19]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[19]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[19]->Fill(mvaValue,wt);
		h_Njets_Varbin[19]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[19]->Fill(BTags,wt);
		h_MET_Varbin[19]->Fill(MET,wt);
		h_PhotonPt_Varbin[19]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[19]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[19]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[19]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[19]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[19]->Fill(Sbin,wt);
		h_Sbins_LL[19]->Fill(searchBin,wt);
		
		if(nGentau_had1>1)
		  {
		    h_Njets[45]->Fill(nHadJets,wt);
		    h_Nbjets[45]->Fill(BTags,wt);
		    h_MET_[45]->Fill(MET,wt);
		    h_PhotonPt[45]->Fill(bestPhoton.Pt(),wt);
		    h_Mt_PhoMET[45]->Fill(mTPhoMET,wt);
		    h_dPhi_PhoMET[45]->Fill(dPhi_PhoMET,wt);
		    h_St[45]->Fill(ST,wt);
		    h_HT[45]->Fill(Ht,wt);
		    h_njets_vs_ST[45]->Fill(nHadJets,ST,wt);
		    h_njets_vs_HT[45]->Fill(nHadJets,Ht,wt);
		    h_ST_vs_ptPho[45]->Fill(ST,bestPhoton.Pt(),wt);
		    h_mvaResponse_baseline[45]->Fill(mvaValue,wt);
		    h_Njets_Varbin[45]->Fill(nHadJets,wt);
		    h_Nbjets_Varbin[45]->Fill(BTags,wt);
		    h_MET_Varbin[45]->Fill(MET,wt);
		    h_PhotonPt_Varbin[45]->Fill(bestPhoton.Pt(),wt);
		    h_St_Varbin[45]->Fill(ST,wt);
		    Sbin = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL[45]->Fill(Sbin, wt);
		    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		    h_Sbins_LL_Validation[45]->Fill(searchBin,wt);
		    Tfbins = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL_v1[45]->Fill(Tfbins,wt);
		    h_TFbins_ElecLL_validation[45]->Fill(Sbin,wt);
		    h_Sbins_LL[45]->Fill(searchBin,wt);

		  }
		else //if (!nGentau_had1)
		  {
		    h_Njets[47]->Fill(nHadJets,wt);
		    h_Nbjets[47]->Fill(BTags,wt);
		    h_MET_[47]->Fill(MET,wt);
		    h_PhotonPt[47]->Fill(bestPhoton.Pt(),wt);
		    h_Mt_PhoMET[47]->Fill(mTPhoMET,wt);
		    h_dPhi_PhoMET[47]->Fill(dPhi_PhoMET,wt);
		    h_St[47]->Fill(ST,wt);
		    h_HT[47]->Fill(Ht,wt);
		    h_njets_vs_ST[47]->Fill(nHadJets,ST,wt);
		    h_njets_vs_HT[47]->Fill(nHadJets,Ht,wt);
		    h_ST_vs_ptPho[47]->Fill(ST,bestPhoton.Pt(),wt);
		    h_mvaResponse_baseline[47]->Fill(mvaValue,wt);
		    h_Njets_Varbin[47]->Fill(nHadJets,wt);
		    h_Nbjets_Varbin[47]->Fill(BTags,wt);
		    h_MET_Varbin[47]->Fill(MET,wt);
		    h_PhotonPt_Varbin[47]->Fill(bestPhoton.Pt(),wt);
		    h_St_Varbin[47]->Fill(ST,wt);
		    Sbin = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL[47]->Fill(Sbin, wt);
		    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		    h_Sbins_LL_Validation[47]->Fill(searchBin,wt);
		    Tfbins = getBinNoV1_le(BTags,nHadJets);
		    h_TFbins_ElecLL_v1[47]->Fill(Tfbins,wt);
		    h_TFbins_ElecLL_validation[47]->Fill(Sbin,wt);
		    h_Sbins_LL[47]->Fill(searchBin,wt);

		  }

	      }

	    if(Debug)
	      {  cout<<"SR mu accept "<<"\t"<<jentry<<endl;}
	  }
	else if(Mu_passAccep && !Mu_passId)//!pass_IdMu && )
	  {
	    // cout<<"Mu_passAccep"<<"\t"<<v_genMu2[0].Pt()<<"\t"<<abs(v_genMu2[0].Eta())<<"\t"<<pass_accep_mu<<endl;
            // if(v_genMu2.size()>1)
            //   cout<<"second"<<"\t"<<v_genMu2[1].Pt()<<"\t"<<abs(v_genMu[1].Eta())<<"\t"<<v_genMu2.size()<<endl;
	    out_nEventsTags[44]="mu-SR: fail ID";

	    nEvents_Selec[44]+=wt;
	    FailId_Mu+=wt; h_selectBaselineYields_->Fill("fail Id mu-SR",wt);
	    h_Njets[10]->Fill(nHadJets,wt);
	    h_Nbjets[10]->Fill(BTags,wt);
	    h_MET_[10]->Fill(MET,wt);
	    h_PhotonPt[10]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[10]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[10]->Fill(dPhi_PhoMET,wt);
	    h_St[10]->Fill(ST,wt);
	    h_HT[10]->Fill(Ht,wt);
	    h_njets_vs_ST[10]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[10]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[10]->Fill(ST,bestPhoton.Pt(),wt);	    
	    h_Njets_Varbin[10]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[10]->Fill(BTags,wt);
	    h_MET_Varbin[10]->Fill(MET,wt);
	    h_PhotonPt_Varbin[10]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[10]->Fill(ST,wt);
	    Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[10]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[10]->Fill(searchBin,wt);
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[10]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[10]->Fill(Sbin,wt);
	    h_Sbins_LL[10]->Fill(searchBin,wt);


	    if(mvaValue>0.0)
	      {
		out_nEventsTags[45]="mu-SR: fail ID- BDT resp";
		nEvents_Selec[45]+=wt;
		h_Njets[20]->Fill(nHadJets,wt);
		h_Nbjets[20]->Fill(BTags,wt);
		h_MET_[20]->Fill(MET,wt);
		h_PhotonPt[20]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[20]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[20]->Fill(dPhi_PhoMET,wt);
		h_St[20]->Fill(ST,wt);
		h_HT[20]->Fill(Ht,wt);
		h_njets_vs_ST[20]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[20]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[20]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[20]->Fill(mvaValue,wt);
		h_Njets_Varbin[20]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[20]->Fill(BTags,wt);
		h_MET_Varbin[20]->Fill(MET,wt);
		h_PhotonPt_Varbin[20]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[20]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[20]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[20]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[20]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[20]->Fill(Sbin,wt);
		h_Sbins_LL[20]->Fill(searchBin,wt);



	      }

	    h_mvaResponse_baseline[10]->Fill(mvaValue,wt);
	  }
	else if(Mu_passId && Mu_passAccep && !Mu_passIso) 
	  {
	    out_nEventsTags[46]="mu-SR: fail Iso";
	    nEvents_Selec[46]+=wt;
	    FailIso_Mu+=wt;
	    h_selectBaselineYields_->Fill("fail Iso mu-SR",wt); 
	    h_Njets[11]->Fill(nHadJets,wt);
	    h_Nbjets[11]->Fill(BTags,wt);
	    h_MET_[11]->Fill(MET,wt);
	    h_PhotonPt[11]->Fill(bestPhoton.Pt(),wt);
	    h_Mt_PhoMET[11]->Fill(mTPhoMET,wt);
	    h_dPhi_PhoMET[11]->Fill(dPhi_PhoMET,wt);
	    h_St[11]->Fill(ST,wt);
	    h_HT[11]->Fill(Ht,wt);
	    h_njets_vs_ST[11]->Fill(nHadJets,ST,wt);
	    h_njets_vs_HT[11]->Fill(nHadJets,Ht,wt);
	    h_ST_vs_ptPho[11]->Fill(ST,bestPhoton.Pt(),wt);
	    h_mvaResponse_baseline[11]->Fill(mvaValue,wt);
	    h_Njets_Varbin[11]->Fill(nHadJets,wt);
	    h_Nbjets_Varbin[11]->Fill(BTags,wt);
	    h_MET_Varbin[11]->Fill(MET,wt);
	    h_PhotonPt_Varbin[11]->Fill(bestPhoton.Pt(),wt);
	    h_St_Varbin[11]->Fill(ST,wt);
	    Sbin = getBinNoV1_le(BTags,nHadJets);
            h_TFbins_ElecLL[11]->Fill(Sbin, wt);
	    searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
	    h_Sbins_LL_Validation[11]->Fill(searchBin,wt);
	    
	    Tfbins = getBinNoV1_le(BTags,nHadJets);
	    h_TFbins_ElecLL_v1[11]->Fill(Tfbins,wt);
	    h_TFbins_ElecLL_validation[11]->Fill(Sbin,wt);
	    h_Sbins_LL[11]->Fill(searchBin,wt);


	    if(mvaValue>0.0)
	      {
		nEvents_Selec[47]+=wt;
		out_nEventsTags[47]="mu-SR: fail Iso - BDT resp";
		h_Njets[21]->Fill(nHadJets,wt);
		h_Nbjets[21]->Fill(BTags,wt);
		h_MET_[21]->Fill(MET,wt);
		h_PhotonPt[21]->Fill(bestPhoton.Pt(),wt);
		h_Mt_PhoMET[21]->Fill(mTPhoMET,wt);
		h_dPhi_PhoMET[21]->Fill(dPhi_PhoMET,wt);
		h_St[21]->Fill(ST,wt);
		h_HT[21]->Fill(Ht,wt);
		h_njets_vs_ST[21]->Fill(nHadJets,ST,wt);
		h_njets_vs_HT[21]->Fill(nHadJets,Ht,wt);
		h_ST_vs_ptPho[21]->Fill(ST,bestPhoton.Pt(),wt);
		h_mvaResponse_baseline[21]->Fill(mvaValue,wt);
		h_Njets_Varbin[21]->Fill(nHadJets,wt);
		h_Nbjets_Varbin[21]->Fill(BTags,wt);
		// h_MET_Varbin[21]->Fill(MET,wt);
		// h_PhotonPt_Varbin[21]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[21]->Fill(ST,wt);
		Sbin = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL[21]->Fill(Sbin, wt);
		searchBin = getBinNoV6_WithOnlyBLSelec(BTags,nHadJets);
		h_Sbins_LL_Validation[21]->Fill(searchBin,wt);
		Tfbins = getBinNoV1_le(BTags,nHadJets);
		h_TFbins_ElecLL_v1[21]->Fill(Tfbins,wt);
		h_TFbins_ElecLL_validation[21]->Fill(Sbin,wt);
		h_Sbins_LL[21]->Fill(searchBin,wt);
		h_MET_Varbin[21]->Fill(MET,wt_LL);
                h_PhotonPt_Varbin[21]->Fill(bestPhoton.Pt(),wt_LL);



	      }

	  }
    if(Debug)
      cout<<"SR mu "<<"\t"<<jentry<<endl;
      }

    // if (NElectrons == 0 && NMuons == 0 )
    //   {
    // 	if(isoElectronTracks==0 && isoMuonTracks ==0 && isoPionTracks==0)
    // 	  {
    // 	    h_selectBaselineYields_->Fill("Iso track",wt);
    // 	    h_mvaResponse_baseline[1]->Fill(mvaValue,wt);
    // 	    if(Debug) 
    // 	      cout<<"===load tree entry SR region ==="<<"\t"<<jentry<<"\t"<<count_genElec<<"t"<<count_genMu<<"\t"<<count_genTau<<endl;
    //     h_Njets[1]->Fill(nHadJets,wt);
    // 	h_Nbjets[1]->Fill(BTags,wt);
    // 	h_MET_[1]->Fill(MET,wt);
    // 	h_PhotonPt[1]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[1]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[1]->Fill(dPhi_PhoMET,wt);
    // 	h_St[1]->Fill(ST,wt);
    // 	h_HT[1]->Fill(Ht,wt);
    // 	h_njets_vs_ST[1]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[1]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[1]->Fill(ST,bestPhoton.Pt(),wt);
    // 	if(mvaValue>0.0)
    // 	  {
    // 	    h_Njets[18]->Fill(nHadJets,wt);
    // 	    h_Nbjets[18]->Fill(BTags,wt);
    // 	    h_MET_[18]->Fill(MET,wt);
    // 	    h_PhotonPt[18]->Fill(bestPhoton.Pt(),wt);
    // 	    h_Mt_PhoMET[18]->Fill(mTPhoMET,wt);
    // 	    h_dPhi_PhoMET[18]->Fill(dPhi_PhoMET,wt);
    // 	    h_St[18]->Fill(ST,wt);
    // 	    h_HT[18]->Fill(Ht,wt);
    // 	    h_njets_vs_ST[18]->Fill(nHadJets,ST,wt);
    // 	    h_njets_vs_HT[18]->Fill(nHadJets,Ht,wt);
    // 	    h_ST_vs_ptPho[18]->Fill(ST,bestPhoton.Pt(),wt);

    // 	  }
	  
    // 	if(mvaValue>0.3)
    //       {
    //         h_Njets[33]->Fill(nHadJets,wt);
    //         h_Nbjets[33]->Fill(BTags,wt);
    //         h_MET_[33]->Fill(MET,wt);
    //         h_PhotonPt[33]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[33]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[33]->Fill(dPhi_PhoMET,wt);
    //         h_St[33]->Fill(ST,wt);
    //         h_HT[33]->Fill(Ht,wt);
    //         h_njets_vs_ST[33]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[33]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[33]->Fill(ST,bestPhoton.Pt(),wt);

    //       }

    // 	if(count_genElec>0)
    // 	  {
    // 	    if(Debug)
    //           cout<<"===load tree entry SR region ==="<<"\t"<<jentry<<endl;
    // 	    if(mvaValue>0.0)
    // 	      {
    // 		h_Njets[20]->Fill(nHadJets,wt);
    // 		h_Nbjets[20]->Fill(BTags,wt);
    // 		h_MET_[20]->Fill(MET,wt);
    // 		h_PhotonPt[20]->Fill(bestPhoton.Pt(),wt);
    // 		h_Mt_PhoMET[20]->Fill(mTPhoMET,wt);
    // 		h_dPhi_PhoMET[20]->Fill(dPhi_PhoMET,wt);
    // 		h_St[20]->Fill(ST,wt);
    // 		h_HT[20]->Fill(Ht,wt);
    // 		h_njets_vs_ST[20]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[20]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[20]->Fill(ST,bestPhoton.Pt(),wt);

    // 	      }

    // 	    if(mvaValue>0.3)
    // 	      {
    // 		h_Njets[35]->Fill(nHadJets,wt);
    // 		h_Nbjets[35]->Fill(BTags,wt);
    // 		h_MET_[35]->Fill(MET,wt);
    // 		h_PhotonPt[35]->Fill(bestPhoton.Pt(),wt);
    // 		h_Mt_PhoMET[35]->Fill(mTPhoMET,wt);
    // 		h_dPhi_PhoMET[35]->Fill(dPhi_PhoMET,wt);
    // 		h_St[35]->Fill(ST,wt);
    // 		h_HT[35]->Fill(Ht,wt);
    // 		h_njets_vs_ST[35]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[35]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[35]->Fill(ST,bestPhoton.Pt(),wt);

    // 	      }

    // 	    h_Njets[3]->Fill(nHadJets,wt);
    // 	    h_Nbjets[3]->Fill(BTags,wt);
    // 	    h_MET_[3]->Fill(MET,wt);
    // 	    h_PhotonPt[3]->Fill(bestPhoton.Pt(),wt);
    // 	    h_Mt_PhoMET[3]->Fill(mTPhoMET,wt);
    // 	    h_dPhi_PhoMET[3]->Fill(dPhi_PhoMET,wt);
    // 	    h_St[3]->Fill(ST,wt);
    // 	    h_HT[3]->Fill(Ht,wt);
    // 	    h_njets_vs_ST[3]->Fill(nHadJets,ST,wt);
    // 	    h_njets_vs_HT[3]->Fill(nHadJets,Ht,wt);
    // 	    h_ST_vs_ptPho[3]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[3]->Fill(mvaValue,wt);

    // 	  }
    // 	else if(count_genMu>0)
    //       {
    // 	    if(Debug)
    //           cout<<"===load tree entry SR muon region ==="<<"\t"<<jentry<<"\t"<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;

    // 	    if(mvaValue>0.0)
    //           {
    // 		if(Debug)
    // 		  cout<<"===load tree entry SR muon 0region ==="<<"\t"<<jentry<<"\t"<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;

    //             h_Njets[21]->Fill(nHadJets,wt);
    //             h_Nbjets[21]->Fill(BTags,wt);
    //             h_MET_[21]->Fill(MET,wt);
    //             h_PhotonPt[21]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[21]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[21]->Fill(dPhi_PhoMET,wt);
    //             h_St[21]->Fill(ST,wt);
    //             h_HT[21]->Fill(Ht,wt);
    //             h_njets_vs_ST[21]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[21]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[21]->Fill(ST,bestPhoton.Pt(),wt);

    //           }

    //         if(mvaValue>0.3)
    //           {
    //             h_Njets[36]->Fill(nHadJets,wt);
    //             h_Nbjets[36]->Fill(BTags,wt);
    //             h_MET_[36]->Fill(MET,wt);
    //             h_PhotonPt[36]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[36 ]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[36]->Fill(dPhi_PhoMET,wt);
    //             h_St[36]->Fill(ST,wt);
    //             h_HT[36]->Fill(Ht,wt);
    //             h_njets_vs_ST[36]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[36]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[36]->Fill(ST,bestPhoton.Pt(),wt);
    // 		if(Debug)
    //               cout<<"===load tree entry SR muon  0.3region ==="<<"\t"<<jentry<<"\t"<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;
    //           }

    //         h_Njets[4]->Fill(nHadJets,wt);
    //         h_Nbjets[4]->Fill(BTags,wt);
    //         h_MET_[4]->Fill(MET,wt);
    //         h_PhotonPt[4]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[4]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[4]->Fill(dPhi_PhoMET,wt);
    //         h_St[4]->Fill(ST,wt);
    //         h_HT[4]->Fill(Ht,wt);
    //         h_njets_vs_ST[4]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[4]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[4]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[4]->Fill(mvaValue,wt);
    // 	    if(Debug)
    //           cout<<"===load tree entry SR muon 1region ==="<<"\t"<<jentry<<"\t"<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;

    //       }
    // 	else if(count_genTau>0)
    //       {

    // 	    if(mvaValue>0.0)
    //           {
    //             h_Njets[22]->Fill(nHadJets,wt);
    //             h_Nbjets[22]->Fill(BTags,wt);
    // 		h_MET_[22]->Fill(MET,wt);
    //             h_PhotonPt[22]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[22]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[22]->Fill(dPhi_PhoMET,wt);
    //             h_St[22]->Fill(ST,wt);
    //             h_HT[22]->Fill(Ht,wt);
    // 		h_njets_vs_ST[22]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[22]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[22]->Fill(ST,bestPhoton.Pt(),wt);

    //           }

    //         if(mvaValue>0.3)
    //           {
    //             h_Njets[37]->Fill(nHadJets,wt);
    // 		h_Nbjets[37]->Fill(BTags,wt);
    // 		h_MET_[37]->Fill(MET,wt);
    //             h_PhotonPt[37]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[37]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[37]->Fill(dPhi_PhoMET,wt);
    //             h_St[37]->Fill(ST,wt);
    //             h_HT[37]->Fill(Ht,wt);
    //             h_njets_vs_ST[37]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[37]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[37]->Fill(ST,bestPhoton.Pt(),wt);

    //           }
    //         h_Njets[5]->Fill(nHadJets,wt);
    //         h_Nbjets[5]->Fill(BTags,wt);
    //         h_MET_[5]->Fill(MET,wt);
    //         h_PhotonPt[5]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[5]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[5]->Fill(dPhi_PhoMET,wt);
    //         h_St[5]->Fill(ST,wt);
    //         h_HT[5]->Fill(Ht,wt);
    //         h_njets_vs_ST[5]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[5]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[5]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[5]->Fill(mvaValue,wt);

    // 	    if(count_genEl_tau>0)
    // 	      {
    // 		h_Njets[13]->Fill(nHadJets,wt);
    // 		h_Nbjets[13]->Fill(BTags,wt);
    // 		h_MET_[13]->Fill(MET,wt);
    // 		h_PhotonPt[13]->Fill(bestPhoton.Pt(),wt);
    // 		h_Mt_PhoMET[13]->Fill(mTPhoMET,wt);
    // 		h_dPhi_PhoMET[13]->Fill(dPhi_PhoMET,wt);
    // 		h_St[13]->Fill(ST,wt);
    // 		h_HT[13]->Fill(Ht,wt);
    // 		h_njets_vs_ST[13]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[13]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[13]->Fill(ST,bestPhoton.Pt(),wt);
    // 		h_mvaResponse_baseline[13]->Fill(mvaValue,wt);
		
    // 		if(mvaValue>0.0)
    // 		  {
    // 		    h_Njets[25]->Fill(nHadJets,wt);
    // 		    h_Nbjets[25]->Fill(BTags,wt);
    // 		    h_MET_[25]->Fill(MET,wt);
    // 		    h_PhotonPt[25]->Fill(bestPhoton.Pt(),wt);
    // 		    h_Mt_PhoMET[25]->Fill(mTPhoMET,wt);
    // 		    h_dPhi_PhoMET[25]->Fill(dPhi_PhoMET,wt);
    // 		    h_St[25]->Fill(ST,wt);
    // 		    h_HT[25]->Fill(Ht,wt);
    // 		    h_njets_vs_ST[25]->Fill(nHadJets,ST,wt);
    // 		    h_njets_vs_HT[25]->Fill(nHadJets,Ht,wt);
    // 		    h_ST_vs_ptPho[25]->Fill(ST,bestPhoton.Pt(),wt);

    // 		  }

    // 		if(mvaValue>0.3)
    // 		  {
    // 		    h_Njets[40]->Fill(nHadJets,wt);
    // 		    h_Nbjets[40]->Fill(BTags,wt);
    // 		    h_MET_[40]->Fill(MET,wt);
    // 		    h_PhotonPt[40]->Fill(bestPhoton.Pt(),wt);
    // 		    h_Mt_PhoMET[40]->Fill(mTPhoMET,wt);
    // 		    h_dPhi_PhoMET[40]->Fill(dPhi_PhoMET,wt);
    // 		    h_St[40]->Fill(ST,wt);
    // 		    h_HT[40]->Fill(Ht,wt);
    // 		    h_njets_vs_ST[40]->Fill(nHadJets,ST,wt);
    // 		    h_njets_vs_HT[40]->Fill(nHadJets,Ht,wt);
    // 		    h_ST_vs_ptPho[40]->Fill(ST,bestPhoton.Pt(),wt);

    // 		  }


    // 	      }
    // 	    else if(count_genMu_tau>0)
    //           {
    //             h_Njets[14]->Fill(nHadJets,wt);
    //             h_Nbjets[14]->Fill(BTags,wt);
    //             h_MET_[14]->Fill(MET,wt);
    //             h_PhotonPt[14]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[14]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[14]->Fill(dPhi_PhoMET,wt);
    //             h_St[14]->Fill(ST,wt);
    //             h_HT[14]->Fill(Ht,wt);
    //             h_njets_vs_ST[14]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[14]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[14]->Fill(ST,bestPhoton.Pt(),wt);
    // 		h_mvaResponse_baseline[14]->Fill(mvaValue,wt);

    // 		if(mvaValue>0.0)
    //               {
    //                 h_Njets[26]->Fill(nHadJets,wt);
    //                 h_Nbjets[26]->Fill(BTags,wt);
    //                 h_MET_[26]->Fill(MET,wt);
    //                 h_PhotonPt[26]->Fill(bestPhoton.Pt(),wt);
    //                 h_Mt_PhoMET[26]->Fill(mTPhoMET,wt);
    //                 h_dPhi_PhoMET[26]->Fill(dPhi_PhoMET,wt);
    //                 h_St[26]->Fill(ST,wt);
    //                 h_HT[26]->Fill(Ht,wt);
    //                 h_njets_vs_ST[26]->Fill(nHadJets,ST,wt);
    //                 h_njets_vs_HT[26]->Fill(nHadJets,Ht,wt);
    //                 h_ST_vs_ptPho[26]->Fill(ST,bestPhoton.Pt(),wt);

    //               }

    //             if(mvaValue>0.3)
    //               {
    //                 h_Njets[41]->Fill(nHadJets,wt);
    //                 h_Nbjets[41]->Fill(BTags,wt);
    //                 h_MET_[41]->Fill(MET,wt);
    //                 h_PhotonPt[41]->Fill(bestPhoton.Pt(),wt);
    //                 h_Mt_PhoMET[41]->Fill(mTPhoMET,wt);
    //                 h_dPhi_PhoMET[41]->Fill(dPhi_PhoMET,wt);
    //                 h_St[41]->Fill(ST,wt);
    //                 h_HT[41]->Fill(Ht,wt);
    //                 h_njets_vs_ST[41]->Fill(nHadJets,ST,wt);
    //                 h_njets_vs_HT[41]->Fill(nHadJets,Ht,wt);
    //                 h_ST_vs_ptPho[41]->Fill(ST,bestPhoton.Pt(),wt);

    //               }

    //           }

    // 	    else
    // 	      {
    // 		h_Njets[15]->Fill(nHadJets,wt);
    //             h_Nbjets[15]->Fill(BTags,wt);
    //             h_MET_[15]->Fill(MET,wt);
    //             h_PhotonPt[15]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[15]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[15]->Fill(dPhi_PhoMET,wt);
    //             h_St[15]->Fill(ST,wt);
    //             h_HT[15]->Fill(Ht,wt);
    //             h_njets_vs_ST[15]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[15]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[15]->Fill(ST,bestPhoton.Pt(),wt);
    //             h_mvaResponse_baseline[15]->Fill(mvaValue,wt);


    // 		if(mvaValue>0.0)
    //               {
    //                 h_Njets[27]->Fill(nHadJets,wt);
    //                 h_Nbjets[27]->Fill(BTags,wt);
    //                 h_MET_[27]->Fill(MET,wt);
    //                 h_PhotonPt[27]->Fill(bestPhoton.Pt(),wt);
    //                 h_Mt_PhoMET[27]->Fill(mTPhoMET,wt);
    //                 h_dPhi_PhoMET[27]->Fill(dPhi_PhoMET,wt);
    //                 h_St[27]->Fill(ST,wt);
    //                 h_HT[27]->Fill(Ht,wt);
    //                 h_njets_vs_ST[27]->Fill(nHadJets,ST,wt);
    //                 h_njets_vs_HT[27]->Fill(nHadJets,Ht,wt);
    //                 h_ST_vs_ptPho[27]->Fill(ST,bestPhoton.Pt(),wt);

    //               }

    //             if(mvaValue>0.3)
    //               {
    //                 h_Njets[42]->Fill(nHadJets,wt);
    //                 h_Nbjets[42]->Fill(BTags,wt);
    //                 h_MET_[42]->Fill(MET,wt);
    //                 h_PhotonPt[42]->Fill(bestPhoton.Pt(),wt);
    //                 h_Mt_PhoMET[42]->Fill(mTPhoMET,wt);
    //                 h_dPhi_PhoMET[42]->Fill(dPhi_PhoMET,wt);
    //                 h_St[42]->Fill(ST,wt);
    //                 h_HT[42]->Fill(Ht,wt);
    //                 h_njets_vs_ST[42]->Fill(nHadJets,ST,wt);
    //                 h_njets_vs_HT[42]->Fill(nHadJets,Ht,wt);
    //                 h_ST_vs_ptPho[42]->Fill(ST,bestPhoton.Pt(),wt);

    //               }

    // 	      }

    //       }

    // 	else if(count_haddecay>0 && (count_genTau==0 || count_genElec==0 || count_genMu==0) && s_sample.Contains("TTG") )
    // 	  {
    // 	    h_Njets[6]->Fill(nHadJets,wt);
    //         h_Nbjets[6]->Fill(BTags,wt);
    //         h_MET_[6]->Fill(MET,wt);
    //         h_PhotonPt[6]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[6]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[6]->Fill(dPhi_PhoMET,wt);
    //         h_St[6]->Fill(ST,wt);
    //         h_HT[6]->Fill(Ht,wt);
    //         h_njets_vs_ST[6]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[6]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[6]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[6]->Fill(mvaValue,wt);
    // 	    if(mvaValue>0.0)
    //           {
    //             h_Njets[23]->Fill(nHadJets,wt);
    //             h_Nbjets[23]->Fill(BTags,wt);
    //             h_MET_[23]->Fill(MET,wt);
    //             h_PhotonPt[23]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[23]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[23]->Fill(dPhi_PhoMET,wt);
    //             h_St[23]->Fill(ST,wt);
    //             h_HT[23]->Fill(Ht,wt);
    //             h_njets_vs_ST[23]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[23]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[23]->Fill(ST,bestPhoton.Pt(),wt);

    //           }

    //         if(mvaValue>0.3)
    //           {
    //             h_Njets[38]->Fill(nHadJets,wt);
    //             h_Nbjets[38]->Fill(BTags,wt);
    //             h_MET_[38]->Fill(MET,wt);
    //             h_PhotonPt[38]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[38]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[38]->Fill(dPhi_PhoMET,wt);
    //             h_St[38]->Fill(ST,wt);
    //             h_HT[38]->Fill(Ht,wt);
    //             h_njets_vs_ST[38]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[38]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[38]->Fill(ST,bestPhoton.Pt(),wt);

    // 	      }


    // 	  }
    // 	else if(igenPho>0 && s_sample.Contains("WG"))
    // 	       {
    // 		 h_Njets[6]->Fill(nHadJets,wt);
    // 		 h_Nbjets[6]->Fill(BTags,wt);
    // 		 h_MET_[6]->Fill(MET,wt);
    // 		 h_PhotonPt[6]->Fill(bestPhoton.Pt(),wt);
    // 		 h_Mt_PhoMET[6]->Fill(mTPhoMET,wt);
    // 		 h_dPhi_PhoMET[6]->Fill(dPhi_PhoMET,wt);
    // 		 h_St[6]->Fill(ST,wt);
    // 		 h_HT[6]->Fill(Ht,wt);
    // 		 h_njets_vs_ST[6]->Fill(nHadJets,ST,wt);
    // 		 h_njets_vs_HT[6]->Fill(nHadJets,Ht,wt);
    // 		 h_ST_vs_ptPho[6]->Fill(ST,bestPhoton.Pt(),wt);
    // 		 if(mvaValue>0.0)
    // 		   {
    // 		     h_Njets[23]->Fill(nHadJets,wt);
    // 		     h_Nbjets[23]->Fill(BTags,wt);
    // 		     h_MET_[23]->Fill(MET,wt);
    // 		     h_PhotonPt[23]->Fill(bestPhoton.Pt(),wt);
    // 		     h_Mt_PhoMET[23]->Fill(mTPhoMET,wt);
    // 		     h_dPhi_PhoMET[23]->Fill(dPhi_PhoMET,wt);
    // 		     h_St[23]->Fill(ST,wt);
    // 		     h_HT[23]->Fill(Ht,wt);
    // 		     h_njets_vs_ST[23]->Fill(nHadJets,ST,wt);
    // 		     h_njets_vs_HT[23]->Fill(nHadJets,Ht,wt);
    // 		     h_ST_vs_ptPho[23]->Fill(ST,bestPhoton.Pt(),wt);

    // 		   }

    // 		 if(mvaValue>0.3)
    // 		   {
    // 		     h_Njets[38]->Fill(nHadJets,wt);
    // 		     h_Nbjets[38]->Fill(BTags,wt);
    // 		     h_MET_[38]->Fill(MET,wt);
    // 		     h_PhotonPt[38]->Fill(bestPhoton.Pt(),wt);
    // 		     h_Mt_PhoMET[38]->Fill(mTPhoMET,wt);
    // 		     h_dPhi_PhoMET[38]->Fill(dPhi_PhoMET,wt);
    // 		     h_St[38]->Fill(ST,wt);
    // 		     h_HT[38]->Fill(Ht,wt);
    // 		     h_njets_vs_ST[38]->Fill(nHadJets,ST,wt);
    // 		     h_njets_vs_HT[38]->Fill(nHadJets,Ht,wt);
    // 		     h_ST_vs_ptPho[38]->Fill(ST,bestPhoton.Pt(),wt);
    // 		   }
    // 	       }
    // 	else
    //       {


    // 	    h_Njets[11]->Fill(nHadJets,wt);
    //         h_Nbjets[11]->Fill(BTags,wt);
    //         h_MET_[11]->Fill(MET,wt);
    //         h_PhotonPt[11]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[11]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[11]->Fill(dPhi_PhoMET,wt);
    //         h_St[11]->Fill(ST,wt);
    //         h_HT[11]->Fill(Ht,wt);
    //         h_njets_vs_ST[11]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[11]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[11]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[11]->Fill(mvaValue,wt);
	    
    // 	    if(mvaValue>0.0)
    //           {
    //             h_Njets[24]->Fill(nHadJets,wt);
    //             h_Nbjets[24]->Fill(BTags,wt);
    //             h_MET_[24]->Fill(MET,wt);
    //             h_PhotonPt[24]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[24]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[24]->Fill(dPhi_PhoMET,wt);
    //             h_St[24]->Fill(ST,wt);
    //             h_HT[24]->Fill(Ht,wt);
    //             h_njets_vs_ST[24]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[24]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[24]->Fill(ST,bestPhoton.Pt(),wt);

    //           }

    //         if(mvaValue>0.3)
    //           {
    //             h_Njets[39]->Fill(nHadJets,wt);
    //             h_Nbjets[39]->Fill(BTags,wt);
    //             h_MET_[39]->Fill(MET,wt);
    //             h_PhotonPt[39]->Fill(bestPhoton.Pt(),wt);
    //             h_Mt_PhoMET[39]->Fill(mTPhoMET,wt);
    //             h_dPhi_PhoMET[39]->Fill(dPhi_PhoMET,wt);
    //             h_St[39]->Fill(ST,wt);
    //             h_HT[39]->Fill(Ht,wt);
    //             h_njets_vs_ST[39]->Fill(nHadJets,ST,wt);
    //             h_njets_vs_HT[39]->Fill(nHadJets,Ht,wt);
    //             h_ST_vs_ptPho[39]->Fill(ST,bestPhoton.Pt(),wt);
    // 	      }
    //       }

    // 	  }
	
    //   }
    // else if (NElectrons == 1 || NMuons == 1 )
    //   {
    // 	if(Debug)
    // 	  cout<<"===load tree entry variable CR ==="<<"\t"<<jentry<<endl;

    // 	double dr2=bestPhoton.DeltaR((*Electrons)[e_index]);
    // 	if(dr2<=0.2) continue;

    //     h_Njets[2]->Fill(nHadJets,wt);
    //     h_Nbjets[2]->Fill(BTags,wt);
    //     h_MET_[2]->Fill(MET,wt);
    //     h_PhotonPt[2]->Fill(bestPhoton.Pt(),wt);
    //     h_Mt_PhoMET[2]->Fill(mTPhoMET,wt);
    //     h_dPhi_PhoMET[2]->Fill(dPhi_PhoMET,wt);
    //     //h_St[2]->Fill(ST,wt);
    //     h_HT[2]->Fill(Ht,wt);
    //     h_njets_vs_ST[2]->Fill(nHadJets,ST,wt);
    //     h_njets_vs_HT[2]->Fill(nHadJets,Ht,wt);
    //     h_ST_vs_ptPho[2]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[2]->Fill(mvaValue,wt);	
    // 	if(mvaValue>0.0)
    //       {
    //         h_Njets[19]->Fill(nHadJets,wt);
    //         h_Nbjets[19]->Fill(BTags,wt);
    //         h_MET_[19]->Fill(MET,wt);
    //         h_PhotonPt[19]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[19]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[19]->Fill(dPhi_PhoMET,wt);
    //         //h_St[19]->Fill(ST,wt);
    //         h_HT[19]->Fill(Ht,wt);
    //         h_njets_vs_ST[19]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[19]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[19]->Fill(ST,bestPhoton.Pt(),wt);

    //       }

    //     if(mvaValue>0.3)
    //       {
    //         h_Njets[33]->Fill(nHadJets,wt);
    //         h_Nbjets[33]->Fill(BTags,wt);
    //         h_MET_[33]->Fill(MET,wt);
    //         h_PhotonPt[33]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[33]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[33]->Fill(dPhi_PhoMET,wt);
    //         //h_St[33]->Fill(ST,wt);
    //         h_HT[33]->Fill(Ht,wt);
    //         h_njets_vs_ST[33]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[33]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[33]->Fill(ST,bestPhoton.Pt(),wt);

    //       }

    // 	if(count_genElec>0)
    //       {
	    
    // 	    //	    ST = ST+(*v_recEle)[0].Pt();
    // 	    h_Njets[7]->Fill(nHadJets,wt);
    // 	h_Nbjets[7]->Fill(BTags,wt);
    // 	h_MET_[7]->Fill(MET,wt);
    // 	h_PhotonPt[7]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[7]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[7]->Fill(dPhi_PhoMET,wt);
    // 	h_St[7]->Fill(ST,wt);
    // 	h_HT[7]->Fill(Ht,wt);
    // 	h_njets_vs_ST[7]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[7]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[7]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[7]->Fill(mvaValue,wt);
    // 	if(mvaValue>0.0)
    // 	  {
    // 	    h_Njets[28]->Fill(nHadJets,wt);
    // 	    h_Nbjets[28]->Fill(BTags,wt);
    // 	    h_MET_[28]->Fill(MET,wt);
    // 	    h_PhotonPt[28]->Fill(bestPhoton.Pt(),wt);
    // 	    h_Mt_PhoMET[28]->Fill(mTPhoMET,wt);
    // 	    h_dPhi_PhoMET[28]->Fill(dPhi_PhoMET,wt);
    // 	    h_St[28]->Fill(ST,wt);
    // 	    h_HT[28]->Fill(Ht,wt);
    // 	    h_njets_vs_ST[28]->Fill(nHadJets,ST,wt);
    // 	    h_njets_vs_HT[28]->Fill(nHadJets,Ht,wt);
    // 	    h_ST_vs_ptPho[28]->Fill(ST,bestPhoton.Pt(),wt);
    // 	  }
    // 	if(mvaValue>0.3)
    // 	  {
    // 	    h_Njets[43]->Fill(nHadJets,wt);
    // 	    h_Nbjets[43]->Fill(BTags,wt);
    // 	    h_MET_[43]->Fill(MET,wt);
    // 	    h_PhotonPt[43]->Fill(bestPhoton.Pt(),wt);
    // 	    h_Mt_PhoMET[43]->Fill(mTPhoMET,wt);
    // 	    h_dPhi_PhoMET[43]->Fill(dPhi_PhoMET,wt);
    // 	    h_St[43]->Fill(ST,wt);
    // 	    h_HT[43]->Fill(Ht,wt);
    // 	    h_njets_vs_ST[43]->Fill(nHadJets,ST,wt);
    // 	    h_njets_vs_HT[43]->Fill(nHadJets,Ht,wt);
    // 	    h_ST_vs_ptPho[43]->Fill(ST,bestPhoton.Pt(),wt);
    // 	  }
    //   }
    // else if(count_genMu>0)
    //   {
    // 	//ST = ST+(*v_recMu)[0].Pt(); 
    // 	h_Njets[8]->Fill(nHadJets,wt);
    // 	h_Nbjets[8]->Fill(BTags,wt);
    // 	h_MET_[8]->Fill(MET,wt);
    // 	h_PhotonPt[8]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[8]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[8]->Fill(dPhi_PhoMET,wt);
    // 	h_St[8]->Fill(ST,wt);
    // 	h_HT[8]->Fill(Ht,wt);
    // 	h_njets_vs_ST[8]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[8]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[8]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[8]->Fill(mvaValue,wt);
    // 	if(mvaValue>0.0)
    //       {
    //         h_Njets[29]->Fill(nHadJets,wt);
    //         h_Nbjets[29]->Fill(BTags,wt);
    //         h_MET_[29]->Fill(MET,wt);
    //         h_PhotonPt[29]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[29]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[29]->Fill(dPhi_PhoMET,wt);
    //         h_St[29]->Fill(ST,wt);
    //         h_HT[29]->Fill(Ht,wt);
    //         h_njets_vs_ST[29]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[29]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[29]->Fill(ST,bestPhoton.Pt(),wt);
    //       }
    //     if(mvaValue>0.3)
    //       {
    //         h_Njets[44]->Fill(nHadJets,wt);
    //         h_Nbjets[44]->Fill(BTags,wt);
    //         h_MET_[44]->Fill(MET,wt);
    //         h_PhotonPt[44]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[44]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[44]->Fill(dPhi_PhoMET,wt);
    //         h_St[44]->Fill(ST,wt);
    //         h_HT[44]->Fill(Ht,wt);
    //         h_njets_vs_ST[44]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[44]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[44]->Fill(ST,bestPhoton.Pt(),wt);
    //       }

    //   }
    // else if(count_genTau>0)
    //   {
    // 	//ST = ST+(*Taus)[0].Pt();
	
    // 	h_Njets[9]->Fill(nHadJets,wt);
    // 	h_Nbjets[9]->Fill(BTags,wt);
    // 	h_MET_[9]->Fill(MET,wt);
    // 	h_PhotonPt[9]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[9]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[9]->Fill(dPhi_PhoMET,wt);
    // 	h_St[9]->Fill(ST,wt);
    // 	h_HT[9]->Fill(Ht,wt);
    // 	h_njets_vs_ST[9]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[9]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[9]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[9]->Fill(mvaValue,wt);

    // 	if(mvaValue>0.0)
    //       {
    //         h_Njets[30]->Fill(nHadJets,wt);
    //         h_Nbjets[30]->Fill(BTags,wt);
    //         h_MET_[30]->Fill(MET,wt);
    //         h_PhotonPt[30]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[30]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[30]->Fill(dPhi_PhoMET,wt);
    //         h_St[30]->Fill(ST,wt);
    //         h_HT[30]->Fill(Ht,wt);
    //         h_njets_vs_ST[30]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[30]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[30]->Fill(ST,bestPhoton.Pt(),wt);
    //       }
    //     if(mvaValue>0.3)
    //       {
    //         h_Njets[45]->Fill(nHadJets,wt);
    //         h_Nbjets[45]->Fill(BTags,wt);
    //         h_MET_[45]->Fill(MET,wt);
    //         h_PhotonPt[45]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[45]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[45]->Fill(dPhi_PhoMET,wt);
    //         h_St[45]->Fill(ST,wt);
    //         h_HT[45]->Fill(Ht,wt);
    //         h_njets_vs_ST[45]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[45]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[45]->Fill(ST,bestPhoton.Pt(),wt);
    //       }

    //   }
    // else if(count_haddecay>0 && (count_genTau==0 || count_genElec==0 || count_genMu==0) && s_sample.Contains("TTG") )
    //   {
    // 	ST = ST+genTau_W.Pt();
    // 	h_Njets[10]->Fill(nHadJets,wt);
    // 	h_Nbjets[10]->Fill(BTags,wt);
    // 	h_MET_[10]->Fill(MET,wt);
    // 	h_PhotonPt[10]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[10]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[10]->Fill(dPhi_PhoMET,wt);
    // 	h_St[10]->Fill(ST,wt);
    // 	h_HT[10]->Fill(Ht,wt);
    // 	h_njets_vs_ST[10]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[10]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[10]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[10]->Fill(mvaValue,wt);
    // 	if(mvaValue>0.0)
    //       {
    //         h_Njets[31]->Fill(nHadJets,wt);
    //         h_Nbjets[31]->Fill(BTags,wt);
    //         h_MET_[31]->Fill(MET,wt);
    //         h_PhotonPt[31]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[31]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[31]->Fill(dPhi_PhoMET,wt);
    //         h_St[31]->Fill(ST,wt);
    //         h_HT[31]->Fill(Ht,wt);
    //         h_njets_vs_ST[31]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[31]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[31]->Fill(ST,bestPhoton.Pt(),wt);
    //       }
    //     if(mvaValue>0.3)
    //       {
    //         h_Njets[46]->Fill(nHadJets,wt);
    //         h_Nbjets[46]->Fill(BTags,wt);
    //         h_MET_[46]->Fill(MET,wt);
    //         h_PhotonPt[46]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[46]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[46]->Fill(dPhi_PhoMET,wt);
    //         h_St[46]->Fill(ST,wt);
    //         h_HT[46]->Fill(Ht,wt);
    //         h_njets_vs_ST[46]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[46]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[46]->Fill(ST,bestPhoton.Pt(),wt);
    //       }


    //   }

    // else if(igenPho>0 && s_sample.Contains("WG"))
    //   {	    h_Njets[10]->Fill(nHadJets,wt);
    // 	    h_Nbjets[10]->Fill(BTags,wt);
    // 	    h_MET_[10]->Fill(MET,wt);
    // 	    h_PhotonPt[10]->Fill(bestPhoton.Pt(),wt);
    // 	    h_Mt_PhoMET[10]->Fill(mTPhoMET,wt);
    // 	    h_dPhi_PhoMET[10]->Fill(dPhi_PhoMET,wt);
    // 	    h_St[10]->Fill(ST,wt);
    // 	    h_HT[10]->Fill(Ht,wt);
    // 	    h_njets_vs_ST[10]->Fill(nHadJets,ST,wt);
    // 	    h_njets_vs_HT[10]->Fill(nHadJets,Ht,wt);
    // 	    h_ST_vs_ptPho[10]->Fill(ST,bestPhoton.Pt(),wt);
    // 	    h_mvaResponse_baseline[10]->Fill(mvaValue,wt);
    // 	    if(mvaValue >0.0)
    // 	      {
    // 		h_Njets[31]->Fill(nHadJets,wt);
    // 		h_Nbjets[31]->Fill(BTags,wt);
    // 		h_MET_[31]->Fill(MET,wt);
    // 		h_PhotonPt[31]->Fill(bestPhoton.Pt(),wt);
    // 		h_Mt_PhoMET[31]->Fill(mTPhoMET,wt);
    // 		h_dPhi_PhoMET[31]->Fill(dPhi_PhoMET,wt);
    // 		h_St[31]->Fill(ST,wt);
    // 		h_HT[31]->Fill(Ht,wt);
    // 		h_njets_vs_ST[31]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[31]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[31]->Fill(ST,bestPhoton.Pt(),wt);
    // 	      }
    // 	    if(mvaValue>0.3)
    // 	      {
    // 		h_Njets[46]->Fill(nHadJets,wt);
    // 		h_Nbjets[46]->Fill(BTags,wt);
    // 		h_MET_[46]->Fill(MET,wt);
    // 		h_PhotonPt[46]->Fill(bestPhoton.Pt(),wt);
    // 		h_Mt_PhoMET[46]->Fill(mTPhoMET,wt);
    // 		h_dPhi_PhoMET[46]->Fill(dPhi_PhoMET,wt);
    // 		h_St[46]->Fill(ST,wt);
    // 		h_HT[46]->Fill(Ht,wt);
    // 		h_njets_vs_ST[46]->Fill(nHadJets,ST,wt);
    // 		h_njets_vs_HT[46]->Fill(nHadJets,Ht,wt);
    // 		h_ST_vs_ptPho[46]->Fill(ST,bestPhoton.Pt(),wt);
    // 	      }


    //   }

    // else
    //   {
    // 	//	cout<<count_genTau<<"\t"<<count_genElec<<"\t"<<count_genMu<<endl;
    // 	h_Njets[12]->Fill(nHadJets,wt);
    // 	h_Nbjets[12]->Fill(BTags,wt);
    // 	h_MET_[12]->Fill(MET,wt);
    // 	h_PhotonPt[12]->Fill(bestPhoton.Pt(),wt);
    // 	h_Mt_PhoMET[12]->Fill(mTPhoMET,wt);
    // 	h_dPhi_PhoMET[12]->Fill(dPhi_PhoMET,wt);
    // 	h_St[12]->Fill(ST,wt);
    // 	h_HT[12]->Fill(Ht,wt);
    // 	h_njets_vs_ST[12]->Fill(nHadJets,ST,wt);
    // 	h_njets_vs_HT[12]->Fill(nHadJets,Ht,wt);
    // 	h_ST_vs_ptPho[12]->Fill(ST,bestPhoton.Pt(),wt);
    // 	h_mvaResponse_baseline[12]->Fill(mvaValue,wt);
    // 	if(mvaValue>0.0)
    // 	  {
    //         h_Njets[32]->Fill(nHadJets,wt);
    //         h_Nbjets[32]->Fill(BTags,wt);
    //         h_MET_[32]->Fill(MET,wt);
    //         h_PhotonPt[32]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[32]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[32]->Fill(dPhi_PhoMET,wt);
    //         h_St[32]->Fill(ST,wt);
    //         h_HT[32]->Fill(Ht,wt);
    //         h_njets_vs_ST[32]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[32]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[32]->Fill(ST,bestPhoton.Pt(),wt);
    //       }
    //     if(mvaValue>0.3)
    //       {
    //         h_Njets[47]->Fill(nHadJets,wt);
    //         h_Nbjets[47]->Fill(BTags,wt);
    //         h_MET_[47]->Fill(MET,wt);
    //         h_PhotonPt[47]->Fill(bestPhoton.Pt(),wt);
    //         h_Mt_PhoMET[47]->Fill(mTPhoMET,wt);
    //         h_dPhi_PhoMET[47]->Fill(dPhi_PhoMET,wt);
    //         h_St[47]->Fill(ST,wt);
    //         h_HT[47]->Fill(Ht,wt);
    //         h_njets_vs_ST[47]->Fill(nHadJets,ST,wt);
    //         h_njets_vs_HT[47]->Fill(nHadJets,Ht,wt);
    //         h_ST_vs_ptPho[47]->Fill(ST,bestPhoton.Pt(),wt);
    //       }
    // 	// h_St[2]->Fill(ST,wt);
    // 	// if(mvaValue>0.0) h_St[19]->Fill(ST,wt);
    // 	// if(mvaValue>0.3) h_St[33]->Fill(ST,wt); 
    //   }
    // 	h_St[2]->Fill(ST,wt);
    //     if(mvaValue>0.0) h_St[19]->Fill(ST,wt);
    //     if(mvaValue>0.3) h_St[33]->Fill(ST,wt);

    //   }
    // 
    if(Debug)
      cout<<"filling the branches in tree"<<endl;
    //    if(Debug)
    //      cout<<"===load tree entry check2 ==="<<"\t"<<jentry<<endl;

 
    if (k > decade)
      cout<<"endl"<<endl;
    }//loop over entries
   //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
   cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<endl;
  cout << "============   Electron     ===============" <<endl;

  cout<<"applied weights "<<" "<<wt<<endl;
  cout<<"CR electron :==  "<<"\t"<<nCR_elec<<"\t"<<endl; 
  cout<<"SR electron :==  "<<"\t"<<nSR_elec<<"\t"<<endl;
  cout<<"SR e- : Fail acceptance "<<"\t"<<FailAccept_Elec<<"\t"<<endl;
  cout<<"SR e- : Fail Id " <<"\t"<<FailId_Elec<<"\t"<<endl;
  cout<< "SR e- : Fail Iso" <<"\t"<<FailIso_Elec<<"\t"<<endl;

  cout << "============   Muon     ===============" <<endl;

  cout<<"applied weights "<<" "<<wt<<endl;
  cout<<"CR muon :==  "<<"\t"<<nCR_mu<<"\t"<<endl;
  cout<<"SR muon :==  "<<"\t"<<nSR_mu<<"\t"<<endl;
  cout<<"SR mu : Fail acceptance "<<"\t"<<FailAccept_Mu<<"\t"<<endl;
  cout<<"SR mu : Fail Id " <<"\t"<<FailId_Mu<<"\t"<<endl;
  cout<< "SR mu : Fail Iso" <<"\t"<<FailIso_Mu<<"\t"<<wt*FailIso_Mu<<endl;

  for (int i =0;i<49;i++)
    {
      cout<<out_nEventsTags[i]<<" :"<<"\t"<<"=====: "<<"\t"<<nEvents_Selec[i]<<endl;
      
      cout<<" "<<"\t"<<endl;
      
    }

  for (int i =0;i<49;i++)
    {
      cout<<nEvents_Selec[i]<<endl;
    }  
  // outfile->cd();
  // skim_tree->Write();
  // outfile->Close();
  // cout<<"outFile: "<<outfileName<<" written!!"<<endl;

 
}

int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nbjets, int nHadJets)
{
  
  int sBin=-100,m_i=0;
  if(nbjets==0 ){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=7;}
    else if(nHadJets>=7)               { sBin=13;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=19;}
    else if(nHadJets==5 || nHadJets==6){ sBin=25;}
    else if(nHadJets>=7)               { sBin=31;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<99.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
        break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
        break; }
    }
  }
  else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
      if(METLowEdge_v3_1[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
        break;}
      else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
        break; }
    }
  }

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
int AnalyzeLightBSM::getBinNoV7_le(int nbjets, int nHadJets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(MET<150)sBin=1; else if (MET>150) sBin=2;}
    else if(nHadJets==3)     { if(MET<150)sBin=3; else if (MET>150) sBin=4;}
    else if(nHadJets==4)     { if(MET<150)sBin=5; else if (MET>150) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=7; else if (MET>150) sBin=8;}
    else if(nHadJets>=7)   { if(MET<150)sBin=9; else if (MET>150) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(MET<150)sBin=11; else if (MET>150) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=13; else if (MET>150) sBin=14;}
    else if(nHadJets>=7)   { if(MET<150)sBin=15; else if (MET>150) sBin=16;}
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV1_le(int nbjets, int nHadJets){
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
int AnalyzeLightBSM::getBinNoV7(int nHadJets){
  int sBin=-100,m_i=0;
  if(BTags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=6;}
    else if(nHadJets>=7)               { sBin=11;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=16;}
    else if(nHadJets==5 || nHadJets==6){ sBin=21;}
    else if(nHadJets>=7)               { sBin=26;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge1.size()-1;i++){
      if(METLowEdge1[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge1[i] && MET < METLowEdge1[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge1[METLowEdge1.size()-1])  { sBin = 6         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge2.size()-1;i++){
      if(METLowEdge2[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge2[i] && MET < METLowEdge2[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge2[METLowEdge2.size()-1])  { sBin = sBin+5   ;break; }
    }
  }
  return sBin;
}
TLorentzVector AnalyzeLightBSM::getBestPhoton(){
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
  for(int iPho=0;iPho<Photons->size();iPho++){
    if( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
      goodPho.push_back( (*Photons)[iPho] );
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



double AnalyzeLightBSM::getGendRLepPho(){//MC only
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  // vector<TLorentzVector> goodPho;
  // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
  //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back((*Photons)[iPhoton]);
  // }
  // if(goodPho.size()!=0) 
  genPho1 =getBestPhoton();
  
  for(int i=0;i<GenParticles->size();i++){
     if((*GenParticles)[i].Pt()!=0){
       if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	 if(genLep1.Pt() < ((*GenParticles)[i]).Pt()) genLep1 = ((*GenParticles)[i]);
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
   for (int igen=0; igen<(*GenParticles).size(); igen++)
     {
       if((*GenParticles)[igen].Pt()!=0){
	 if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
	   genPho = ((*GenParticles)[igen]);
	   v_genPho.push_back(genPho);
	 }
       }
     }
   return MinDr(bestPhoton,v_genPho);   
 }
double AnalyzeLightBSM::getGendRElecPho(){//MC only                                                                                             
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  genPho1 =getBestPhoton();
  for(int i=0;i<GenParticles->size();i++){
    if((*GenParticles)[i].Pt()!=0){
      if( abs((*GenParticles_PdgId)[i])==11 && (abs((*GenParticles_ParentId)[i])<=25 )&& (abs((*GenParticles_ParentId)[i])!=15)){
	if(genLep1.Pt() < ((*GenParticles)[i]).Pt()) genLep1 = ((*GenParticles)[i]);
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
      for(int i=0 ; i < GenElectrons->size(); i++)
	{
	  if((*GenElectrons)[i].Pt()!=0)
	    {
	      genEle1 = ((*GenElectrons)[i]);
	      v_genLep2.push_back(genEle1);
	    }

	}
    }
  else
    {
      for(int i=0 ; i < GenMuons->size(); i++)
	{
	  if((*GenMuons)[i].Pt()!=0)
	    {
	      genMu1 = ((*GenMuons)[i]);
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
      for(int i=0 ; i < GenElectrons->size(); i++)
        {
          if((*GenElectrons)[i].Pt()!=0)
            {
              genEle1 = ((*GenElectrons)[i]);
              v_genLep2.push_back(genEle1);
            }

        }
  //   }
  // else
  //   {
      for(int i=0 ; i < GenMuons->size(); i++)
        {
          if((*GenMuons)[i].Pt()!=0)
            {
              genMu1 = ((*GenMuons)[i]);
              v_genLep2.push_back(genMu1);
            }
        }
      //  }
  return MinDr(bestPhoton,v_genLep2);
}

double AnalyzeLightBSM::getGenRecodRLep(TLorentzVector recoLep1){//MC only                                                                                                
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  for(int i=0;i<GenParticles->size();i++){
    if((*GenParticles)[i].Pt()!=0){
      if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	if(genLep1.Pt() < ((*GenParticles)[i]).Pt()) genLep1 = ((*GenParticles)[i]);
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
  for(int i=0;i<GenParticles->size();i++){
    if((*GenParticles)[i].Pt()!=0){
      //if(abs((*GenParticles_PdgId)[i])==22 || abs((*GenParticles_ParentId)[i])==22) continue;
      if(recoPho.DeltaR((*GenParticles)[i])<0.2 )
	{
	  //cout<<"You got to be kidding me"<<endl;	  
	  h_parID_matchRecoPho->Fill((*GenParticles_PdgId)[i],wt);
	  h_parentID_vsPartId_matchPho->Fill((*GenParticles_PdgId)[i],(*GenParticles_ParentId)[i],wt);
	  if(abs((*GenParticles_PdgId)[i])==22)
	    {
	      h_phopT_BeamRemenant->Fill(recoPho.Pt(),wt);
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

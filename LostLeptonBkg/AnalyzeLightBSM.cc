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
  if(s_data.Contains("2016")) lumiInfb=35.922;
  if(s_data.Contains("2017")) lumiInfb=41.529;
  if(s_data.Contains("2018")) lumiInfb=59.74;
  if(s_data.Contains("signal"))lumiInfb= 137.19; //since we only have 2018 dataset
  std::string s_process = sample;
  double cross_section = getCrossSection(s_process);
  cout<<cross_section<<"\t"<<"analyzed process"<<endl;
  //



nentries = 1000;
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
float nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0;
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
      //      TLorentzVector 
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
      for (int igen=0; igen<branch_size; igen++)
      	{	  
      	  pdgID= (*GenParticles_PdgId)[igen];
      	  parentId=(*GenParticles_ParentId)[igen];
	  if((*GenParticles)[igen].Pt()!=0){
	    if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
	      if(genPho.Pt() < (*GenParticles)[igen].Pt()){
		leadGenPhoIdx1 = igen;
		leadGenPho = ((*GenParticles)[igen]);	      
	      }
	      genPho = ((*GenParticles)[igen]);
	      v_genPho.push_back(genPho);	      
	    }	  
	}
	}	  	  	 
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////     Identifying the decay modes for W & TTbar  //////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int igen=0; igen<branch_size; igen++)
	{
	  pdgID= (*GenParticles_PdgId)[igen];
	  parentId=(*GenParticles_ParentId)[igen];
	  if(abs(parentId)==24) // || parentId==-24)
	    {
	      if(abs(pdgID)==11) count_genElec++;
	      else if(abs(pdgID)==13) count_genMu++;
	      else if(abs(pdgID)==15) count_genTau++;
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
       TLorentzVector bestPhoton=getBestPhoton();
      
      // *******************  Selecting Jet objects ********************************//
      int minDRindx=-100,phoMatchingJetIndx=-100,nHadJets=0;
      double minDR=99999,ST=0,Ht=0;
      vector<TLorentzVector> hadJets;
                                                            
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
            hadJets.push_back((*Jets)[i]);}
      }
    if( minDR<0.3 ) phoMatchingJetIndx=minDRindx;
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
    if(bestPhotonIndxAmongPhotons<0) continue;
    bool bestPhoHasPxlSeed=true;
    if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
    if( bestPhoHasPxlSeed ) continue;
    //    if(bestPhoton.Pt()>30) 
    if(bestPhoton.Pt()>30)                                                                                                                       
      h_selectBaselineYields_->Fill("pt>30",wt);                                                                                                
    else continue;

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
    if(MET>100) 
    h_selectBaselineYields_->Fill("MET>100",wt);
    else continue;
    if(ST>300)
      h_selectBaselineYields_->Fill("ST>300",wt);
    else continue;
    if( s_sample.Contains("TTJets-HT") && madHT<600) continue;
    if( s_sample.Contains("TTJets-inc") && madHT>600) continue;
    if(!genphocheck)
      {
        genphomatch_before++;
        double mindr_Pho_genlep=getGenLep(bestPhoton);
        if( s_sample.Contains("TTG") )
          {
            if(!hasGenPromptPhoton)
              {
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
                if(!(madMinPhotonDeltaR >= 0.3 && mindr_Pho_genlep >=0.5 ))continue;
              }
          }//Gen prompt                                                                                                                                                                                  
        if(s_sample.Contains("WG"))
          {
            if(!hasGenPromptPhoton)
              {
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
                if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 ))continue;
              }
          }//Gen prompt    
                                                                                                                                                                               
        if(s_sample.Contains("WJets"))
          {
            if(!hasGenPromptPhoton)
              {
		if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)) continue;
              }
          }
        if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17"))
          {
            if(!hasGenPromptPhoton)
              {
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
              }
            else if(hasGenPromptPhoton)
              {
                if(!(madMinPhotonDeltaR < 0.3 || mindr_Pho_genlep < 0.5)) continue;
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

    h_selectBaselineYields_->Fill("after bkg cmp",wt);
    //MEt cleaning filters
    if(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))
      {
	if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0) ) continue;
	//	if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 &&	     BadPFMuonFilter && NVtx > 0) ) continue;

        //if(!(PrimaryVertexFilter ==1 && globalSuperTightHalo2016Filter == 1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter && ecalBadCalibReducedExtraFilter && NVtx>0 && eeBadScFilter)) continue;
	// if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadPFMuonFilter && NVtx > 0) ) continue;

      }
    h_selectBaselineYields_->Fill("MetCleaning",wt);
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

    //  ***************** dPhi between MEt & two elading jets : mainly affecting the faek MET from QCD ************************ //
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
	h_selectBaselineYields_->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    
    // ********************  photon jet matching : pT comparison *************************//
    if(phoMatchingJetIndx>=0 && ((*Jets)[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
    if(phoMatchingJetIndx<0) continue;
    h_selectBaselineYields_->Fill("additional Bhumika",wt);
    if (bestPhoton.Pt()>=30)
      h_selectBaselineYields_->Fill("pho-pT>30",wt);
    else continue;
    //cout<<"===load tree entry ==="<<"\t"<<jentry<<endl;
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

    nsurVived+=wt;
    //h_selectBaselineYields_->Fill("veto minDr(#photon,lept) <0.2",wt);
    int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
    TLorentzVector genPho1,genEle1,neutr,genMu1,genTau1,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
    vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_genLep1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_recEle,v_recMu;
    int leadGenPhoIdx=-100, mu_index=-100;
    int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
    int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
    for(int i=0;i<GenParticles->size();i++){
      if((*GenParticles)[i].Pt()!=0){
	if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=25) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 ){
	  leadGenPhoIdx = i;
	  genPho1 = ((*GenParticles)[i]);
	  v_genPho1.push_back(genPho1);
	}
      }        
    }	  
    bool hadtau = false;       
    for(int i=0 ; i < GenElectrons->size(); i++)
      {
	if((*GenElectrons)[i].Pt()!=0)
	  {	     
	    nGenEle1++;
	    genEle1 = ((*GenElectrons)[i]);
	    v_genEle2.push_back(genEle1);
	    v_genLep2.push_back(genEle1);	      
	  }
      }
    for(int i=0 ; i < GenMuons->size(); i++)
      {
	if((*GenMuons)[i].Pt()!=0)
	  {
	    nGenMu1++;
	    genMu1 = ((*GenMuons)[i]);
	    v_genMu2.push_back(genMu1);
	    v_genLep2.push_back(genMu1);
	  }
      }
       
    for(int i=0 ; i < GenTaus->size(); i++)
      {
	if((*GenTaus)[i].Pt()!=0)
	  {
	    nGenTau1++;
	    genTau1 = ((*GenTaus)[i]);
	    v_genTau2.push_back(genTau1);
	    v_genLep2.push_back(genTau1);
	    if((*GenTaus_had)[i])
	      nGentau_had1++;
	  }
      }

    
    nlep=0;
    // if(Electrons->size() != NElectrons) //==0  && (v_genEle2.size()!=0))//GenElectrons->size()!=0))
    //   cout<<"event: "<<" "<<jentry<<"\t"<<NElectrons<<"\t"<<Electrons->size()<<"\t"<<GenElectrons->size()<<"\t"<<endl;
    bool Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false;
    bool Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
    
    if(v_genEle2.size()>=1 && Electrons->size()!=0)
      {
	if((v_genEle2[0].Pt()>10 && abs(v_genEle2[0].Eta())<2.5) || (v_genEle2[1].Pt()>10 && abs(v_genEle2[1].Eta())<2.5))
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
		    if(dR<0.2){ pass_IdElec++; Elec_passId= true;}
		    else
		      { fail_IdElec++; Elec_passId= false;}
		    if(dR<0.2 && (*Electrons_passIso)[i]==1)
		      { pass_isoElec++; Elec_passIso = true;nelec_reco++; nlep++; e_index=i; recElec=(*Electrons)[i]; v_recEle.push_back((*Electrons)[i]);}
		    else if(dR<0.2 && (*Electrons_passIso)[i]!=1)
		      {
			fail_isoEle++; Elec_passIso = false;}
		    if(nlep>1) cout<<i<<""<<e_index<<" , e pt = "<<v_lep1.Pt()<<endl;                                                   
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
	  }	     
      }
    if(v_genMu2.size()>=1 && Muons->size()!=0)
      {
	if(((v_genMu2[0].Pt()>10 && abs(v_genMu2[0].Eta())<2.5) || (v_genMu2[1].Pt()>10 && abs(v_genMu2[1].Eta())<2.5)))
	  {
	    nlep=0;
	    pass_accep_mu++;
	    Mu_passAccep = true;	    
	    for(int i=0; i<Muons->size() ; i++)
	      {
		if(((*Muons)[i].Pt()>10) && abs((*Muons)[i].Eta())<2.5 && nlep<1)
		  {
		    double dR = MinDr((*Muons)[i],v_genMu2);
		    if(dR<0.2){ pass_IdMu++; Mu_passId= true;}
		    else
		      {fail_IdMu++; Mu_passId= false;} 
		    if(dR<0.2 && (*Muons_passIso)[i]==1) 
		      { pass_isoMu++; Mu_passIso = true;nmu_reco++;nlep++; mu_index=i; recMu=(*Muons)[i]; v_recMu.push_back((*Muons)[i]);}
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
          }
      }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////             SR region   ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //    cout<<""<<"\t"<<jentry<<endl;

   if((*Electrons).size()==0 && nelec_reco!=0) //checking if the contribution to SR with no reco e- is coming up or not?
      cout<<nelec_reco<<endl;
    


   // int nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0;
    if(lost_elec_flag && nelec_reco ==1 && nmu_reco==0)
      {
	if(!Elec_passAccep || !Elec_passId || !Elec_passIso) //cross check                                                       
    	  cout<<"cross check - elec SR"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;  
	if(isoMuonTracks !=0 || isoPionTracks!=0) continue;
	double dr2=bestPhoton.DeltaR((*Electrons)[e_index]);
	if(dr2<=0.2) continue;
	ST = ST+(*Electrons)[e_index].Pt();
	st = ST;
	mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");
	double mTElecMET=sqrt(2*((*Electrons)[e_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,(*Electrons)[e_index].Phi()))));
	if(mTElecMET>100) continue;
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
	int Sbin = getBinNoV7_le(BTags,nHadJets);
	//cout<<"ll bins"<<"\t"<<Sbin<<endl;
	
	h_TFbins_ElecLL[2]->Fill(Sbin, wt);

	nCR_elec+=wt; 
	h_selectBaselineYields_->Fill("elec CR",wt);	
	if(mvaValue<0.0)
	  {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
	    h_TFbins_ElecLL[12]->Fill(Sbin, wt);



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
	    double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);
	    if(mindR_genPho>=0.2)
	      {
		double mindrPho_genlep=getGenLep(bestPhoton);//, true);
		if(mindrPho_genlep<0.2) {continue; nfakeRatePho+=wt;}
		h_selectBaselineYields_->Fill("prompt Photon Checks",wt);
	      }	    
	    //	  }

	if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
	if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue;//to reject W->qq' type of events
	if(nGentau_had1>1) continue;                                                                                                   
	survived_vetohad+=wt;
	if(nGenMu1>0) continue;
	if(nGenEle1==0) continue;
	h_selectBaselineYields_->Fill("Elec SR",wt);
	nSR_elec+=wt;
	survived_elecge1+=wt;
	if(v_genEle2.size()==0) {TLorentzVector v1;v_genEle2.push_back(v1);}
	sortTLorVec(&v_genEle2);
	v_genEle1=v_genEle2;
	if(!Elec_passAccep)
	  {
	    // if(Debug && Elec_passId || Elec_passIso)
	    // cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
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
	    int Sbin = getBinNoV7_le(BTags,nHadJets);
	    h_TFbins_ElecLL[4]->Fill(Sbin, wt);


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
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[14]->Fill(Sbin, wt);


	      }

	  }
	else if(!Elec_passId)
	  {
	    if(Debug && !Elec_passAccep|| Elec_passIso)
	      cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;
	    FailId_Elec+=wt;
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
	    int Sbin = getBinNoV7_le(BTags,nHadJets);
	    h_TFbins_ElecLL[5]->Fill(Sbin, wt);


	    if(mvaValue>0.0)
	      {
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
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[15]->Fill(Sbin, wt);


	      }

	  }
	else if(!Elec_passIso)
	  {
	    // if(Debug && Elec_passAccep|| Elec_passId)
            //   cout<<jentry<<"\t"<<Elec_passAccep<<"\t"<<Elec_passId<<"\t"<<Elec_passIso<<endl;

	    FailIso_Elec+=wt;
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
	    int Sbin = getBinNoV7_le(BTags,nHadJets);
	    h_TFbins_ElecLL[6]->Fill(Sbin, wt);


	    if(mvaValue>0.0)
	      {
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
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[16]->Fill(Sbin, wt);


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
	h_mvaResponse_baseline[3]->Fill(mvaValue,wt);
	int Sbin = getBinNoV7_le(BTags,nHadJets);
	h_TFbins_ElecLL[3]->Fill(Sbin, wt);
	//	cout<<"ll bins: SR"<<"\t"<<Sbin<<endl;

	h_Njets_Varbin[3]->Fill(nHadJets,wt);
	h_Nbjets_Varbin[3]->Fill(BTags,wt);
	h_MET_Varbin[3]->Fill(MET,wt);
	h_PhotonPt_Varbin[3]->Fill(bestPhoton.Pt(),wt);
	h_St_Varbin[3]->Fill(ST,wt);


	if(mvaValue>0.0)
	  {
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
	Sbin = getBinNoV7_le(BTags,nHadJets);
	h_TFbins_ElecLL[13]->Fill(Sbin, wt);


	  }
      }
      
    if(Debug)
      cout<<"SR elec "<<"\t"<<jentry<<endl;
    if(!lost_elec_flag && nelec_reco ==0 && nmu_reco==1)
      {
	if(isoElectronTracks !=0 || isoPionTracks!=0) continue;
	double dr2=bestPhoton.DeltaR((*Muons)[mu_index]);
        if(dr2<=0.2) continue;
	ST = ST+(*Muons)[mu_index].Pt();
	double mTMuMET=sqrt(2*((*Muons)[mu_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,(*Muons)[mu_index].Phi()))));
	if(mTMuMET>100) continue;
	st = ST;
        mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");

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
	int Sbin = getBinNoV7_le(BTags,nHadJets);
	h_TFbins_ElecLL[7]->Fill(Sbin, wt);


	if(mvaValue<0.0)
          {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
            h_TFbins_ElecLL[17]->Fill(Sbin, wt);

          }

      }
    if(Debug)
      cout<<"CR mu "<<"\t"<<jentry<<endl;

    if(!lost_elec_flag && nelec_reco ==0 && nmu_reco==0)
      {
	//nSR_mu+=wt;
        if(isoElectronTracks!=0 || isoMuonTracks !=0 || isoPionTracks!=0) continue;
	// if(!s_data.Contains("data"))
        //   {
            double mindR_genPho = getdR_GenPho_RecoPho(bestPhoton);
            if(mindR_genPho>=0.2)
              {
                double mindrPho_genlep=getGenLep(bestPhoton);//, false);
                if(mindrPho_genlep<0.2) {continue; nfakeRatePho+=wt;}
                h_selectBaselineYields_->Fill("prompt Photon Checks",wt);
              }
          // }
	  //  if(nGentau_had1>0 && nGenMu1==0 && nGenEle1==0) continue;
	if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue;//to reject W->qq' type of events
	//survived_vetohad+=wt;
	if(!(nGenEle1==0 || (nGenEle1==1 && nGenMu1==1))) continue;
	//	if(nGenEle1>0) continue;
	h_selectBaselineYields_->Fill("Mu SR",wt);  
	survived_elecge1+=wt;
	//if(nGentau_had1>1) continue;

	if(v_genMu2.size()==0) {TLorentzVector v1;v_genMu2.push_back(v1);}
	if(nGenEle1>0) continue;
        //if(nGenMu1==0) continue;

	nSR_mu+=wt;
	sortTLorVec(&v_genMu2);
	v_genEle1=v_genMu2;
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
	int Sbin = getBinNoV7_le(BTags,nHadJets);
	h_TFbins_ElecLL[8]->Fill(Sbin, wt);


	if(mvaValue>0.0)
          {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
            h_TFbins_ElecLL[18]->Fill(Sbin, wt);


          }

	if(!Mu_passAccep)
	  {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
            h_TFbins_ElecLL[9]->Fill(Sbin, wt);


	    if(mvaValue>0.0)
	      {
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
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[19]->Fill(Sbin, wt);


	      }

	    if(Debug)
	      {  cout<<"SR mu accept "<<"\t"<<jentry<<endl;}
	  }
	else if(!pass_IdMu)
	  {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
            h_TFbins_ElecLL[10]->Fill(Sbin, wt);


	    if(mvaValue>0.0)
	      {
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
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[20]->Fill(Sbin, wt);


	      }

	    h_mvaResponse_baseline[10]->Fill(mvaValue,wt);
	  }
	else if(!Mu_passIso) 
	  {
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
	    Sbin = getBinNoV7_le(BTags,nHadJets);
            h_TFbins_ElecLL[11]->Fill(Sbin, wt);
	    if(mvaValue>0.0)
	      {
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
		h_MET_Varbin[21]->Fill(MET,wt);
		h_PhotonPt_Varbin[21]->Fill(bestPhoton.Pt(),wt);
		h_St_Varbin[21]->Fill(ST,wt);
		Sbin = getBinNoV7_le(BTags,nHadJets);
		h_TFbins_ElecLL[21]->Fill(Sbin, wt);


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
      if(METLowEdge_v3[i]<199.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
	break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7         ;
	break; }
    }
  }
  else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
      if(METLowEdge_v3_1[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
	break;}
      else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
	break; }
    }
  }

  else if(sBin==37){
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
      // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }

    }
  }

  else if(sBin==44){
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
      // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }

    }
  }
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

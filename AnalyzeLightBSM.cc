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

  if (argc < 4) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" "<<"N2-Mass"<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  AnalyzeLightBSM ana(inputFileList, outFileName, data,sample);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList,sample,outFileName);
  Tools::Instance();

  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

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
  if(s_data.Contains("signal"))lumiInfb= 59.74; //since we only have 2018 dataset
  std::string s_process = sample;
  double cross_section = getCrossSection(s_process);
  cout<<cross_section<<"\t"<<"analyzed process"<<endl;
  //    nentries = 100;
  cout<<"Event"<<"\t"<<"par-ID "<<"\t"<<"parentID"<<"\t"<<"GenMET"<<"\t"<<"MET"<<"\t"<<"pT"<<"\t"<<"Eta"<<"\t"<<"Phi"<<"\t"<<"E"<<endl;
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
      int count_genElec=0,count_genMu=0,count_genTau=0, igenPho=0, count_genEl_tau=0,count_genMu_tau=0;
      //      TLorentzVector 
      TLorentzVector genPho_W,genElec_W,genMu_W,genTau_W,genElec_Tau,genMu_Tau,genW,genNuElec_W,genNuMu_W,genNuTau_W,genElecNu_Tau,genMuNu_Tau;
      // for (int igen=0; igen<branch_size; igen++)
      // 	{
      // 	  pdgID= (*GenParticles_PdgId)[igen];
      // 	  parentId=(*GenParticles_ParentId)[igen];

      // 	  cout<<jentry<<"\t"<<pdgID<<"\t"<<parentId<<"\t"<<GenMET<<"\t"<<MET<<"\t"<<(*GenParticles)[igen].Pt()<<"\t"<<(*GenParticles)[igen].Eta\
      // 	    ()<<"\t"<<(*GenParticles)[igen].Phi()<<"\t"<<(*GenParticles)[igen].E()<<endl;
      // 	}

      for (int igen=0; igen<branch_size; igen++)
	{
	  pdgID= (*GenParticles_PdgId)[igen];
	  parentId=(*GenParticles_ParentId)[igen];
	  if(abs(pdgID)==24)
	    genW.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==11 && abs(parentId)==24)
	    genElec_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==13 && abs(parentId)==24)
            genMu_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==15 && abs(parentId)==24)
            genTau_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());	 
	  else if(abs(pdgID)==12 && abs(parentId)==24)
            genNuElec_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==14 && abs(parentId)==24)
            genNuMu_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==16 && abs(parentId)==24)
            genNuTau_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  else if(abs(pdgID)==22 && abs(parentId)==24)
            genPho_W.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  // else
	  //   {
	  //     cout<<pdgID<<"\t"<<parentId<<endl;
	  //   }
	  if(abs(parentId)==15)
	    {
	  //     if(abs(pdgID)==11)
	  // 	{
	  // 	  genElec_Tau.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  // 	  count_genEl_tau++;
	  // 	}
	  //     else if(abs(pdgID)==13)
	  // 	{
	  // 	genMu_Tau.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	  // 	count_genMu_tau++;
	  // 	}
	  //     else
		if(abs(pdgID)==12)
                genElecNu_Tau.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());
	      else if(abs(pdgID)==14)
                genMuNu_Tau.SetPtEtaPhiE((*GenParticles)[igen].Pt(),(*GenParticles)[igen].Eta(),(*GenParticles)[igen].Phi(),(*GenParticles)[igen].E());

	    }

	  if(abs(parentId)==24) // || parentId==-24)
	    {
	      if(abs(pdgID)==11) count_genElec++;
	      else if(abs(pdgID)==13) count_genMu++;
	      else if(abs(pdgID)==15) count_genTau++;
	      if(abs(pdgID)==11 || abs(pdgID)==13 ||abs(pdgID)==15)		
		{count++;//t<<"lepton"<<"\t"<<jentry<<endl;
		}
	      else if (abs(pdgID)==12 ||abs(pdgID)==14||abs(pdgID)==16)
		count++;//t<<"neutrino"<<"\t"<<jentry<<endl;
	      else
		{
		  igenPho++;
		  h_counts->Fill("pho",wt);
		  //		  cout<<"pdgID"<<"\t"<<pdgID<<"\t"<<jentry<<endl;
		}
	      
	      // 	dR_Ele++;
	      // else
	      // 	cout<<pdgID<<endl;
	      // if(abs(pdgID)==11 || abs(pdgID)==13 ||abs(pdgID)==15) //||(pdgID==-11 || pdgID==-13||pdgID==-15))
	      // 	//{
	      // 	  h_Lep_pT->Fill((*GenParticles)[igen].Pt(),wt);
		  
	      // 	  //double mindR_lep_reco = getGenRecodRLep(recoLep1);
	      // 	  //}
	      // if(abs(pdgID)==11||abs(pdgID)==12) //||(pdgID==-11||pdgID==12))
	      // 	//{
	      // 	  h_2d_elec_nu_pT->Fill((*GenParticles)[igen].Pt(),wt);
		  
	      // //}
	      // if(abs(pdgID)==12|| abs(pdgID)==14 || abs(pdgID)==16) //||pdgID==16)
	      // 	h_2d_nu_MET_pT->Fill((*GenParticles)[igen].Pt(),MET,wt);
	      for (int i=0; i<3;i++)
		{
		  if(abs(pdgID)==Ids[i])
		    h_counts->Fill(ids[i],wt);
		}
	    }
	}
      //      cout<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;
      //checking the reoncstruction efficiency
      // for (int ireco=0;ireco<(*Electrons).size();ireco++)
      // 	{
	  
      // 	}
      if((*GenTaus).size()!=0 && GenTaus_had)
	{
      	  if((*GenElectrons).size()!=0)
	    
	    {
	      //cout<<"(*GenTaus).size()"<<"\t"<<(*GenTaus).size()<<endl;
	      genElec_Tau.SetPtEtaPhiE((*GenElectrons)[0].Pt(),(*GenElectrons)[0].Eta(),(*GenElectrons)[0].Phi(),(*GenElectrons)[0].E());
	      count_genEl_tau++;
	    }
	  else if((*GenMuons).size()!=0)
	    {
	      genMu_Tau.SetPtEtaPhiE((*GenMuons)[0].Pt(),(*GenMuons)[0].Eta(),(*GenMuons)[0].Phi(),(*GenMuons)[0].E());
	      count_genMu_tau++;
	    }
	}

      h_GenpTvsEta[0]->Fill(genW.Pt(),genW.Eta(),wt);
      h_GenEtavsPhi[0]->Fill(genW.Eta(),genW.Phi(),wt);
      h_GenEnvsEta[0]->Fill(genW.E(),genW.Eta(),wt);
      h_GenPtvsPhi[0]->Fill(genW.Pt(),genW.Phi(),wt);
      h_GenMETvsGenpT[0]->Fill(GenMET,genW.Pt(),wt);
      
      h_2d_pT[0]->Fill(genW.Pt(),genElec_W.Pt(),wt);
      h_2d_Eta[0]->Fill(genW.Eta(),genElec_W.Eta(),wt);
      h_2d_Phi[0]->Fill(genW.Phi(),genElec_W.Phi(),wt);
      h_2d_Et[0]->Fill(genW.Et(),genElec_W.Et(),wt);
      
      h_2d_pT[1]->Fill(genW.Pt(),genMu_W.Pt(),wt);
      h_2d_Eta[1]->Fill(genW.Eta(),genMu_W.Eta(),wt);
      h_2d_Phi[1]->Fill(genW.Phi(),genMu_W.Phi(),wt);
      h_2d_Et[1]->Fill(genW.Et(),genMu_W.Et(),wt);

      h_2d_pT[2]->Fill(genW.Pt(),genTau_W.Pt(),wt);
      h_2d_Eta[2]->Fill(genW.Eta(),genTau_W.Eta(),wt);
      h_2d_Phi[2]->Fill(genW.Phi(),genTau_W.Phi(),wt);
      h_2d_Et[2]->Fill(genW.Et(),genTau_W.Et(),wt);

      h_2d_pT[3]->Fill(genTau_W.Pt(),genElec_Tau.Pt(),wt);
      h_2d_Eta[3]->Fill(genTau_W.Eta(),genElec_Tau.Eta(),wt);
      h_2d_Phi[3]->Fill(genTau_W.Phi(),genElec_Tau.Phi(),wt);
      h_2d_Et[3]->Fill(genTau_W.Et(),genElec_Tau.Et(),wt);

      h_2d_pT[4]->Fill(genTau_W.Pt(),genMu_Tau.Pt(),wt);
      h_2d_Eta[4]->Fill(genTau_W.Eta(),genMu_Tau.Eta(),wt);
      h_2d_Phi[4]->Fill(genTau_W.Phi(),genMu_Tau.Phi(),wt);
      h_2d_Et[4]->Fill(genTau_W.Et(),genMu_Tau.Et(),wt);


      h_GenpT[0]->Fill(genW.Pt(),wt);
      h_GenEta[0]->Fill(genW.Eta(),wt);
      h_GenPhi[0]->Fill(genW.Phi(),wt);
      h_GenEn[0]->Fill(genW.E(),wt);
      h_GenMET[0]->Fill(GenMET,wt);
      h_RecoMET[0]->Fill(MET,wt);

      h_GenpT[8]->Fill(genElec_W.Pt(),wt);
      h_GenEta[8]->Fill(genElec_W.Eta(),wt);
      h_GenPhi[8]->Fill(genElec_W.Phi(),wt);
      h_GenEn[8]->Fill(genElec_W.E(),wt);

      h_GenpT[11]->Fill(genMu_W.Pt(),wt);
      h_GenEta[11]->Fill(genMu_W.Eta(),wt);
      h_GenPhi[11]->Fill(genMu_W.Phi(),wt);
      h_GenEn[11]->Fill(genMu_W.E(),wt);
      h_GenpT[14]->Fill(genTau_W.Pt(),wt);
      h_GenEta[14]->Fill(genTau_W.Eta(),wt);
      h_GenPhi[14]->Fill(genTau_W.Phi(),wt);
      h_GenEn[14]->Fill(genTau_W.E(),wt);
      h_GenpT[19]->Fill(genElec_Tau.Pt(),wt);
      h_GenEta[19]->Fill(genElec_Tau.Eta(),wt);
      h_GenPhi[19]->Fill(genElec_Tau.Phi(),wt);
      h_GenEn[19]->Fill(genElec_Tau.E(),wt);

      
      h_GenpTvsEta[8]->Fill(genElec_W.Pt(),genElec_W.Eta(),wt);
      h_GenEtavsPhi[8]->Fill(genElec_W.Eta(),genElec_W.Phi(),wt);
      h_GenEnvsEta[8]->Fill(genElec_W.E(),genElec_W.Eta(),wt);
      h_GenPtvsPhi[8]->Fill(genElec_W.Pt(),genElec_W.Phi(),wt);
      h_GenMETvsGenpT[8]->Fill(GenMET,genElec_W.Pt(),wt);

      h_GenpTvsEta[11]->Fill(genMu_W.Pt(),genMu_W.Eta(),wt);
      h_GenEtavsPhi[11]->Fill(genMu_W.Eta(),genMu_W.Phi(),wt);
      h_GenEnvsEta[11]->Fill(genMu_W.E(),genMu_W.Eta(),wt);
      h_GenPtvsPhi[11]->Fill(genMu_W.Pt(),genMu_W.Phi(),wt);
      h_GenMETvsGenpT[11]->Fill(GenMET,genMu_W.Pt(),wt);

      h_GenpTvsEta[14]->Fill(genTau_W.Pt(),genTau_W.Eta(),wt);
      h_GenEtavsPhi[14]->Fill(genTau_W.Eta(),genTau_W.Phi(),wt);
      h_GenEnvsEta[14]->Fill(genTau_W.E(),genTau_W.Eta(),wt);
      h_GenPtvsPhi[14]->Fill(genTau_W.Pt(),genTau_W.Phi(),wt);
      h_GenMETvsGenpT[14]->Fill(GenMET,genTau_W.Pt(),wt);


      
      
      int count_genEle=0,count_recEle=0;
      TLorentzVector genPho1,genLep1;
      int leadGenPhoIdx=-100;
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

    //if (nHadJets<2) cout<<"wrong event"<<endl;
    if(bestPhotonIndxAmongPhotons<0) continue;
    bool bestPhoHasPxlSeed=true;
    if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
    if( bestPhoHasPxlSeed ) continue;
    h_selectBaselineYields_->Fill("photon_selec",wt);
    if (k > decade)
      cout << 10 * k << " %-3" << endl;
    if(bestPhoton.Pt()>30) 
    h_selectBaselineYields_->Fill("pt>30",wt);
    else continue;

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

    // if(isoElectronTracks==0 && isoMuonTracks ==0 && isoPionTracks==0)
    //   {
    //     h_selectBaselineYields_->Fill("Iso track",wt);
    //   }
    // else continue;

    if (bestPhoton.Pt()>=30)
      h_selectBaselineYields_->Fill("pho-pT>30",wt);
    else continue;

    // h_Njets[7]->Fill(nHadJets,wt);
    // h_Nbjets[7]->Fill(BTags,wt);
    // h_MET_[7]->Fill(MET,wt);
    // h_PhotonPt[7]->Fill(bestPhoton.Pt(),wt);
    // h_Mt_PhoMET[7]->Fill(mTPhoMET,wt);
    // h_dPhi_PhoMET[7]->Fill(dPhi_PhoMET,wt);
    // h_St[7]->Fill(ST,wt);
    // h_HT[7]->Fill(Ht,wt);
    // h_njets_vs_ST[7]->Fill(nHadJets,ST,wt);
    // h_njets_vs_HT[7]->Fill(nHadJets,Ht,wt);
    // h_ST_vs_ptPho[7]->Fill(ST,bestPhoton.Pt(),wt);


    if (NElectrons == 0 && NMuons == 0 )
      {
	if(isoElectronTracks==0 && isoMuonTracks ==0 && isoPionTracks==0)
	  {
	    h_selectBaselineYields_->Fill("Iso track",wt);
	  
	//	else continue;
	//cout<<count_genElec<<"\t"<<count_genMu<<"\t"<<count_genTau<<endl;
	    h_GenpT[1]->Fill(genW.Pt(),wt);
	    h_GenEta[1]->Fill(genW.Eta(),wt);
	    h_GenPhi[1]->Fill(genW.Phi(),wt);
	    h_GenEn[1]->Fill(genW.E(),wt);
	    h_GenMET[1]->Fill(GenMET,wt);
	    h_RecoMET[1]->Fill(MET,wt);
	    
	    h_GenpTvsEta[1]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[1]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[1]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[1]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[1]->Fill(GenMET,genW.Pt(),wt);

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
	if(count_genElec>0)
	  {
	    
	    h_GenpTvsEta[9]->Fill(genElec_W.Pt(),genElec_W.Eta(),wt);
	    h_GenEtavsPhi[9]->Fill(genElec_W.Eta(),genElec_W.Phi(),wt);
	    h_GenEnvsEta[9]->Fill(genElec_W.E(),genElec_W.Eta(),wt);
	    h_GenPtvsPhi[9]->Fill(genElec_W.Pt(),genElec_W.Phi(),wt);
	    h_GenMETvsGenpT[9]->Fill(GenMET,genElec_W.Pt(),wt);

	    h_GenpTvsEta[3]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[3]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[3]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[3]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[3]->Fill(GenMET,genW.Pt(),wt);

	    h_GenpT[9]->Fill(genElec_W.Pt(),wt);
	    h_GenEta[9]->Fill(genElec_W.Eta(),wt);
	    h_GenPhi[9]->Fill(genElec_W.Phi(),wt);
	    h_GenEn[9]->Fill(genElec_W.E(),wt);

	    h_GenpT[3]->Fill(genW.Pt(),wt);
            h_GenEta[3]->Fill(genW.Eta(),wt);
            h_GenPhi[3]->Fill(genW.Phi(),wt);
            h_GenEn[3]->Fill(genW.E(),wt);
            h_GenMET[3]->Fill(GenMET,wt);
	    h_RecoMET[3]->Fill(MET,wt);


	    // for (int igen=0; igen<branch_size; igen++)
            //   {
            //     pdgID= (*GenParticles_PdgId)[igen];
            //     parentId=(*GenParticles_ParentId)[igen];

	    // cout<<jentry<<"\t"<<pdgID<<"\t"<<parentId<<"\t"<<GenMET<<"\t"<<MET<<"\t"<<(*GenParticles)[igen].Pt()<<"\t"<<(*GenParticles)[igen].Eta()<<"\t"<<(*GenParticles)[igen].Phi()<<"\t"<<(*GenParticles)[igen].E()<<endl;
	    //   }

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

	  }
	else if(count_genMu>0)
          {
	    h_GenpTvsEta[12]->Fill(genMu_W.Pt(),genMu_W.Eta(),wt);
            h_GenEtavsPhi[12]->Fill(genMu_W.Eta(),genMu_W.Phi(),wt);
            h_GenEnvsEta[12]->Fill(genMu_W.E(),genMu_W.Eta(),wt);
            h_GenPtvsPhi[12]->Fill(genMu_W.Pt(),genMu_W.Phi(),wt);
            h_GenMETvsGenpT[12]->Fill(GenMET,genMu_W.Pt(),wt);

	    h_GenpTvsEta[4]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[4]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[4]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[4]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[4]->Fill(GenMET,genW.Pt(),wt);

	    h_GenpT[12]->Fill(genMu_W.Pt(),wt);
	    h_GenEta[12]->Fill(genMu_W.Eta(),wt);
	    h_GenPhi[12]->Fill(genMu_W.Phi(),wt);
	    h_GenEn[12]->Fill(genMu_W.E(),wt);

	    h_GenpT[4]->Fill(genW.Pt(),wt);
            h_GenEta[4]->Fill(genW.Eta(),wt);
            h_GenPhi[4]->Fill(genW.Phi(),wt);
            h_GenEn[4]->Fill(genW.E(),wt);
            h_GenMET[4]->Fill(GenMET,wt);
	    h_RecoMET[4]->Fill(MET,wt);


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

          }
	else if(count_genTau>0)
          {

            h_GenpTvsEta[15]->Fill(genTau_W.Pt(),genTau_W.Eta(),wt);
            h_GenEtavsPhi[15]->Fill(genTau_W.Eta(),genTau_W.Phi(),wt);
            h_GenEnvsEta[15]->Fill(genTau_W.E(),genTau_W.Eta(),wt);
            h_GenPtvsPhi[15]->Fill(genTau_W.Pt(),genTau_W.Phi(),wt);
            h_GenMETvsGenpT[15]->Fill(GenMET,genTau_W.Pt(),wt);

	    h_GenpTvsEta[5]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[5]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[5]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[5]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[5]->Fill(GenMET,genW.Pt(),wt);

	    h_GenpT[15]->Fill(genTau_W.Pt(),wt);
	    h_GenEta[15]->Fill(genTau_W.Eta(),wt);
	    h_GenPhi[15]->Fill(genTau_W.Phi(),wt);
	    h_GenEn[15]->Fill(genTau_W.E(),wt);
	    h_GenpT[18]->Fill(genElec_Tau.Pt(),wt);
	    h_GenEta[18]->Fill(genElec_Tau.Eta(),wt);
	    h_GenPhi[18]->Fill(genElec_Tau.Phi(),wt);
	    h_GenEn[18]->Fill(genElec_Tau.E(),wt);
	    

	    h_GenpT[5]->Fill(genW.Pt(),wt);
            h_GenEta[5]->Fill(genW.Eta(),wt);
            h_GenPhi[5]->Fill(genW.Phi(),wt);
            h_GenEn[5]->Fill(genW.E(),wt);
            h_GenMET[5]->Fill(GenMET,wt);
	    h_RecoMET[5]->Fill(MET,wt);


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
	    
	    if(count_genEl_tau>0)
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

	      }
	    if(count_genMu_tau>0)
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

              }


          }

	else if(igenPho>0)
	  {
	    h_GenpTvsEta[6]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[6]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[6]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[6]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[6]->Fill(GenMET,genW.Pt(),wt);

	    //cout<<count_genTau<<"\t"<<count_genElec<<"\t"<<count_genMu<<endl;
	    h_GenpT[6]->Fill(genW.Pt(),wt);
            h_GenEta[6]->Fill(genW.Eta(),wt);
            h_GenPhi[6]->Fill(genW.Phi(),wt);
            h_GenEn[6]->Fill(genW.E(),wt);
            h_GenMET[6]->Fill(GenMET,wt);
	    h_RecoMET[6]->Fill(MET,wt);


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
	    // for (int igen=0; igen<branch_size; igen++)
            //   {
            //     pdgID= (*GenParticles_PdgId)[igen];
            //     parentId=(*GenParticles_ParentId)[igen];
	    // 	cout<<jentry<<"\t"<<pdgID<<"\t"<<parentId<<"\t"<<GenMET<<"\t"<<MET<<"\t"<<(*GenParticles)[igen].Pt()<<"\t"<<(*GenParticles)[igen].Eta()<<"\t"<<(*GenParticles)[igen].Phi()<<"\t"<<(*GenParticles)[igen].E()<<endl;
	    // 	//                cout<<jentry<<"\t"<<pdgID<<"\t"<<parentId<<"\t"<<GenMET<<"\t"<<MET<<endl;
            //   }

	  }
	else
          {
	    h_GenpTvsEta[7]->Fill(genW.Pt(),genW.Eta(),wt);
	    h_GenEtavsPhi[7]->Fill(genW.Eta(),genW.Phi(),wt);
	    h_GenEnvsEta[7]->Fill(genW.E(),genW.Eta(),wt);
	    h_GenPtvsPhi[7]->Fill(genW.Pt(),genW.Phi(),wt);
	    h_GenMETvsGenpT[7]->Fill(GenMET,genW.Pt(),wt);

	    h_GenpT[7]->Fill(genW.Pt(),wt);
            h_GenEta[7]->Fill(genW.Eta(),wt);
            h_GenPhi[7]->Fill(genW.Phi(),wt);
            h_GenEn[7]->Fill(genW.E(),wt);
            h_GenMET[7]->Fill(GenMET,wt);
	    h_RecoMET[7]->Fill(MET,wt);



	    // for (int igen=0; igen<branch_size; igen++)
	    //   {
	    // 	pdgID= (*GenParticles_PdgId)[igen];
	    // 	parentId=(*GenParticles_ParentId)[igen];
	    // 	cout<<jentry<<"\t"<<pdgID<<"\t"<<parentId<<"\t"<<GenMET<<"\t"<<MET<<"\t"<<(*GenParticles)[igen].Pt()<<"\t"<<(*GenParticles)[igen].Eta()<<"\t"<<(*GenParticles)[igen].Phi()<<"\t"<<(*GenParticles)[igen].E()<<endl;
	    //   }
	    //  cout<<count_genTau<<"\t"<<count_genElec<<"\t"<<count_genMu<<endl;
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

          }

	  }
	
      }
    else // if(NElectrons == 1 || NMuons == 1 )
      {

	h_GenpT[2]->Fill(genW.Pt(),wt);
	h_GenEta[2]->Fill(genW.Eta(),wt);
	h_GenPhi[2]->Fill(genW.Phi(),wt);
	h_GenEn[2]->Fill(genW.E(),wt);
	h_GenMET[2]->Fill(GenMET,wt);
	h_RecoMET[2]->Fill(MET,wt);

	h_GenpTvsEta[2]->Fill(genW.Pt(),genW.Eta(),wt);
	h_GenEtavsPhi[2]->Fill(genW.Eta(),genW.Phi(),wt);
	h_GenEnvsEta[2]->Fill(genW.E(),genW.Eta(),wt);
	h_GenPtvsPhi[2]->Fill(genW.Pt(),genW.Phi(),wt);
	h_GenMETvsGenpT[2]->Fill(GenMET,genW.Pt(),wt);

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
	if(count_genElec>0)
          {
	    h_GenpTvsEta[10]->Fill(genElec_W.Pt(),genElec_W.Eta(),wt);
            h_GenEtavsPhi[10]->Fill(genElec_W.Eta(),genElec_W.Phi(),wt);
            h_GenEnvsEta[10]->Fill(genElec_W.E(),genElec_W.Eta(),wt);
            h_GenPtvsPhi[10]->Fill(genElec_W.Pt(),genElec_W.Phi(),wt);
            h_GenMETvsGenpT[10]->Fill(GenMET,genElec_W.Pt(),wt);

	    h_GenpT[10]->Fill(genElec_W.Pt(),wt);
            h_GenEta[10]->Fill(genElec_W.Eta(),wt);
            h_GenPhi[10]->Fill(genElec_W.Phi(),wt);
            h_GenEn[10]->Fill(genElec_W.E(),wt);

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

      }
    else if(count_genMu>0)
      {

	h_GenpTvsEta[13]->Fill(genMu_W.Pt(),genMu_W.Eta(),wt);
	h_GenEtavsPhi[13]->Fill(genMu_W.Eta(),genMu_W.Phi(),wt);
	h_GenEnvsEta[13]->Fill(genMu_W.E(),genMu_W.Eta(),wt);
	h_GenPtvsPhi[13]->Fill(genMu_W.Pt(),genMu_W.Phi(),wt);
	h_GenMETvsGenpT[13]->Fill(GenMET,genMu_W.Pt(),wt);
	h_GenpT[13]->Fill(genMu_W.Pt(),wt);
	h_GenEta[13]->Fill(genMu_W.Eta(),wt);
	h_GenPhi[13]->Fill(genMu_W.Phi(),wt);
	h_GenEn[13]->Fill(genMu_W.E(),wt);

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

      }
    else if(count_genTau>0)
      {

	h_GenpTvsEta[16]->Fill(genTau_W.Pt(),genTau_W.Eta(),wt);
	h_GenEtavsPhi[16]->Fill(genTau_W.Eta(),genTau_W.Phi(),wt);
	h_GenEnvsEta[16]->Fill(genTau_W.E(),genTau_W.Eta(),wt);
	h_GenPtvsPhi[16]->Fill(genTau_W.Pt(),genTau_W.Phi(),wt);
	h_GenMETvsGenpT[16]->Fill(GenMET,genTau_W.Pt(),wt);

	h_GenpT[16]->Fill(genTau_W.Pt(),wt);
	h_GenEta[16]->Fill(genTau_W.Eta(),wt);
	h_GenPhi[16]->Fill(genTau_W.Phi(),wt);
	h_GenEn[16]->Fill(genTau_W.E(),wt);
	
	h_GenpT[17]->Fill(genElec_Tau.Pt(),wt);
        h_GenEta[17]->Fill(genElec_Tau.Eta(),wt);
        h_GenPhi[17]->Fill(genElec_Tau.Phi(),wt);
        h_GenEn[17]->Fill(genElec_Tau.E(),wt);
	
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

      }
    else if(igenPho>0)
      {	    h_Njets[10]->Fill(nHadJets,wt);
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

      }

    else
      {
	//	cout<<count_genTau<<"\t"<<count_genElec<<"\t"<<count_genMu<<endl;
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

      }

      }

    if(Debug)
      cout<<"filling the branches in tree"<<endl;

 
    if (k > decade)
      cout<<"endl"<<endl;
    }//loop over entries
  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
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
double AnalyzeLightBSM::getGenLep(TLorentzVector bestPhoton){                                     
  vector<TLorentzVector> v_genLep2;
  TLorentzVector genMu1, genEle1;
  for(int i=0 ; i < GenElectrons->size(); i++)
    {
      if((*GenElectrons)[i].Pt()!=0)
        {
          genEle1 = ((*GenElectrons)[i]);
          v_genLep2.push_back(genEle1);
        }

    }
  for(int i=0 ; i < GenMuons->size(); i++)
    {
      if((*GenMuons)[i].Pt()!=0)
        {
          genMu1 = ((*GenMuons)[i]);
          v_genLep2.push_back(genMu1);
        }
    }
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

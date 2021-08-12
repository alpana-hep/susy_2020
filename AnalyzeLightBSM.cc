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

using namespace std;

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

  ana.EventLoop(data,inputFileList,sample);

  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  //Weight 
  //  NumEvents = nentries;
  //cout<<Weight<<endl;
  //Weight = CrossSection/NumEvents;    
  //  wt = Weight*1000.0*35.9;
  TString s_data=data;
  TString s_sample= sample;
  Long64_t nbytes = 0, nb = 0; 
  int decade = 0;
  int count=0;
  double lumiInfb=137;
  if(s_data.Contains("2016")) lumiInfb=35.9;
  if(s_data.Contains("2017")) lumiInfb=41.59;
  if(s_data.Contains("2018")) lumiInfb=59.74;
  if(s_data.Contains("signal"))lumiInfb= 137;
  int count_QCD=0;
  cout<<s_sample<<nentries<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      // if(s_data.Contains("2016")) lumiInfb=35.9;
      // if(s_data.Contains("2017")) lumiInfb=41.59;
      // if(s_data.Contains("2018")) lumiInfb=59.74;
      if(s_data.Contains("signal"))//lumiInfb= 1000*137;
	wt = (0.165*137)/nentries;
      else
	wt = Weight*lumiInfb*1000.0;
      // if(s_sample.Contains("QCD_Jets"))
      // 	{
      // if(madMinPhotonDeltaR<0.4)
      // 	{cout<<madMinPhotonDeltaR<<endl; count_QCD++;}// continue;}
      // 	}
      //checking for G+jets, QCD+jets sample
      // if(s_sample.Contains("GJets_DR") && madMinPhotonDeltaR<0.4) continue;//{ //cout<<"exiting the event"<<endl; continue;}
      // else if(s_sample.Contains("QCD_Jets")&& madMinPhotonDeltaR>0.4) { cout<<madMinPhotonDeltaR<<endl; }
      // else
      // 	{
	  // cout<<"enteringg the event"<<endl;
	  // cout<<Weight<<endl;
      // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
    

    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //my code starts from here .... 19Aug2020
    int branch_size = (*GenParticles).size();
    int ele_branch=(*Electrons).size();
    int pho_branch = (*Photons).size();
    float iGenEle_eta =99999.0, iGenEle_phi=99999.0, iRecEle_eta =99999.0, iRecEle_phi=99999.0, iRecphot_eta =99999.0, iRecphot_phi=99999.0,iGen_Wpt=99999.0,iGen_Wp=99999.0,iGen_WEta=99999.0,iGen_WPhi=99999.0;
    float dR_Ele=0, min_dR=9999, min_dR_pho=9999;
    int count_genEle=0,count_recEle=0;
    
    //cout<<branch_size<<"\t"<<"ele_size"<<"\t"<<ele_branch<<endl;
    if(NMuons ==0 && NElectrons==0)
      {
    	if(isoElectronTracks==0 && isoMuonTracks ==0)// continue;// && isoPionTracks==0)//veto isolated leptons                                                                                      
	  {
	    
	    vector<TLorentzVector> goodPho;
	    for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
	      if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
	    }
	    
	    if(goodPho.size()==0) {count++;continue;}
	    
	    TLorentzVector bestPhoton=goodPho[0];
	    if(bestPhoton.Pt()<20) continue;// cout<< "notPt cut"<<endl;
	    int minDRindx=-100,phoMatchingJetIndx=-100,nHadJets=0;
	    double minDR=99999,ST=0;
	    vector<TLorentzVector> hadJets;

	    for(int i=0;i<Jets->size();i++)
	      {
		if( ((*Jets)[i].Pt() > 30.0) && (abs((*Jets)[i].Eta()) <= 2.4) ){
		  double dR=bestPhoton.DeltaR((*Jets)[i]);
		  if(dR<minDR){minDR=dR;minDRindx=i;}
		  if( !(minDR < 0.3 && i==minDRindx) )
		    hadJets.push_back((*Jets)[i]);}
	      }
	    if( minDR<0.3 ) phoMatchingJetIndx=minDRindx;
	    for(int i=0;i<hadJets.size();i++){
	      if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());}
	      if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
	    }
	    if( minDR<0.3 ) ST=ST+bestPhoton.Pt();
	    if (nHadJets>=2) 
	      {
		if(MET>100) 
		  {
		    if(s_sample.Contains("GJets_DR"))
		      {
			if(madMinPhotonDeltaR<0.4) continue;
			else
			  {
			    h_NJets->Fill(nHadJets,wt);
			    h_NbJets->Fill(BTags,wt);
			    if(MET<1170)
			      h_MeT->Fill(MET,wt);
			    else
			      h_MeT->SetBinContent(h_MeT->GetNbinsX(),wt);
			    if(bestPhoton.Pt()<1000)
			      h_PhoPt->Fill(bestPhoton.Pt(),wt);
			    else
			      h_PhoPt->SetBinContent(h_PhoPt->GetNbinsX(),wt);
			    h_check_PhoPt->Fill(bestPhoton.Pt(),wt);
			    h_HT->Fill(ST,wt);

			  }
		      }
		    else if(s_sample.Contains("QCD_Jets"))
		      {
			if (madMinPhotonDeltaR>0.4) continue;
			else
			  {
			    //  cout<<"success"<<endl;
			    h_NJets->Fill(nHadJets,wt);
                            h_NbJets->Fill(BTags,wt);
                            if(MET<1170)
                              h_MeT->Fill(MET,wt);
                            else
                              h_MeT->SetBinContent(h_MeT->GetNbinsX(),wt);
                            if(bestPhoton.Pt()<1000)
                              h_PhoPt->Fill(bestPhoton.Pt(),wt);
                            else
                              h_PhoPt->SetBinContent(h_PhoPt->GetNbinsX(),wt);
                            h_check_PhoPt->Fill(bestPhoton.Pt(),wt);
                            h_HT->Fill(ST,wt);
			  }
		      }
		    else if(s_sample.Contains("temp"))
		      {
			
			h_NJets->Fill(nHadJets,wt);
			h_NbJets->Fill(BTags,wt);
			if(MET<1170)
			  h_MeT->Fill(MET,wt);
			else
			  h_MeT->SetBinContent(h_MeT->GetNbinsX(),wt);
			if(bestPhoton.Pt()<1000)
			  h_PhoPt->Fill(bestPhoton.Pt(),wt);
			else
			  h_PhoPt->SetBinContent(h_PhoPt->GetNbinsX(),wt);
			h_check_PhoPt->Fill(bestPhoton.Pt(),wt);
			h_HT->Fill(ST,wt);
		      }
		    // if(MET>250)
		    //   {
			
		    //   }
		    // int sBin7 = getBinNoV7(nHadJets);
		    // h_SBins_v7_CD->Fill(sBin7,wt);
		  }
	      }
	  }
      }
    }
  //    }
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
  TLorentzVector v1;
  vector<TLorentzVector> goodPho;
  for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
    if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
  }

  if(goodPho.size()==0) return v1;
  sortTLorVec(&goodPho);
  return goodPho[0];
  
}




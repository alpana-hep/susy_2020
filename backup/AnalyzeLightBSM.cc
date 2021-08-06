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

  if (argc < 3) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" "<<"N2-Mass"<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  AnalyzeLightBSM ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList) {
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
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  int count=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      wt = (0.165*137)/nentries;
      //wt = Weight*1000.0;//*35.9;
      //           cout<<Weight<<endl;
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
    h_NJets_inclusive->Fill(Jets->size(),wt);
    h_NbJets_inclusive->Fill(BTags,1000*Weight);
    //veto events with muon and lepton 
    if(NMuons ==0 && NElectrons==0)
      {
    	if(isoElectronTracks==0 && isoMuonTracks ==0)// && isoPionTracks==0)//veto isolated leptons
    	  {
	    h_NJets_exclusive->Fill(Jets->size(),wt);
	    h_NbJets_exclusive->Fill(BTags,wt);
	    h_Met_exclusive->Fill(MET,wt);
	   
	    //cout<<BTags<<endl;
	    // if(MET>100)
	    //   {
	    vector<TLorentzVector> goodPho;
	    for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
	      if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
	    }

	    if(goodPho.size()==0) {count++;continue;}
	    //else cout<<goodPho.size()<<endl;
	    TLorentzVector bestPhoton=goodPho[0];//getBestPhoton();
	    //if (bestPhoton.Pt()==0) {count++;continue;}
	      
	    h_PhoPt_exclusive->Fill(bestPhoton.Pt(),wt);
	    h_PhoPhi_exclusive->Fill(bestPhoton.Phi(),wt);
	    
		//cout<<bestPhoton.Pt()<<endl;
		// if(bestPhoton.Pt()>100)
		//   {
		//selection//  for best jets
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

	    h_NJets_LepIsoVeto->Fill(nHadJets,wt);
            h_NbJets_LepIsoVeto->Fill(BTags,wt);
            h_Met_LepIsoVeto->Fill(MET,wt);
	    h_PhoPt_LepIsoVeto->Fill(bestPhoton.Pt(),wt);
	    h_PhoPhi_LepIsoVeto->Fill(bestPhoton.Phi(),wt);
	    h_HT_LepIsoVeto->Fill(ST,wt);
	    
	    //applying jet cut
	    if(nHadJets>=2)
	      {
		h_NJets_aJetcut->Fill(nHadJets,wt);
		h_NbJets_aJetcut->Fill(BTags,wt);
		h_Met_aJetcut->Fill(MET,wt);
		h_PhoPt_aJetcut->Fill(bestPhoton.Pt(),wt);
		h_PhoPhi_aJetcut->Fill(bestPhoton.Phi(),wt);
		h_HT_aJetcut->Fill(ST,wt);
		if(MET>100)
		  {
		    h_NJets_aMetcut->Fill(nHadJets,wt);
		    h_NbJets_aMetcut->Fill(BTags,wt);
		    h_Met_aMetcut->Fill(MET,wt);
		    h_PhoPt_aMetcut->Fill(bestPhoton.Pt(),wt);
		    h_PhoPhi_aMetcut->Fill(bestPhoton.Phi(),wt);
		    h_HT_aMetcut->Fill(ST,wt);
		    if(bestPhoton.Pt()>10)
		      {
			h_NJets_PhoPt->Fill(nHadJets,wt);
			h_NbJets_PhoPt->Fill(BTags,wt);
			h_Met_PhoPt->Fill(MET,wt);
			h_PhoPt_PhoPt->Fill(bestPhoton.Pt(),wt);
			h_PhoPhi_PhoPt->Fill(bestPhoton.Phi(),wt);
			h_HT_PhoPt->Fill(ST,wt);
		      }
		  }
		
	      }
	    if(MET>100)
	      {
		h_NJets_ExMET->Fill(nHadJets,wt);
                h_NbJets_ExMET->Fill(BTags,wt);
                h_Met_ExMET->Fill(MET,wt);
                h_PhoPt_ExMET->Fill(bestPhoton.Pt(),wt);
                h_PhoPhi_ExMET->Fill(bestPhoton.Phi(),wt);
		h_HT_ExMET->Fill(ST,wt);
	      }
	    //   {
	    // 	h_NJets_aJetcut->Fill(nHadJets);
            //     h_NbJets_aJetcut->Fill(BTags);
            //     h_Met_aJetcut->Fill(MET);
            //     h_PhoPt_aJetcut->Fill(bestPhoton.Pt());
            //     h_PhoPhi_aJetcut->Fill(bestPhoton.Phi());

	    

		    /////cout<<nHadJets<<endl;
		    // if( (nHadJets >= 2) )
		    //   {
		    // 	//St selection
		    // 	if(bestPhoton.Pt()>190 && ST >500)
		    // 	  {
	    // h_photonPt->Fill(bestPhoton.Pt());//,wt);
	    // h_njets->Fill(NJets);//,wt);
	    // h_MET->Fill(MET);//,wt);
	    // h_bjets->Fill(BTags);//,wt);
		    //}
		 	// else if(bestPhoton.Pt()>100 && ST >800)
			//   {
			   //  h_photonPt->Fill(bestPhoton.Pt(),wt);
                        //     h_njets->Fill(NJets,wt);
                        //     h_MET->Fill(MET,wt);
                        //     h_bjets->Fill(BTags,wt);
                        // // 
	  }
			// else
			//   continue;
	//}
	  // 	  }
	  //     }
	  // }
      }
    }   // loop over entries
  //h_photonPt->Draw();
  //  h->Write();
  cout<<count<<endl;
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




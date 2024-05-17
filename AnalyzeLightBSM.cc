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

   if (argc <5 ) {
     cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"Which pho_ID"<<endl;
     return -1;
   }
   const char *inputFileList = argv[1];
   const char *outFileName   = argv[2];
   const char *data          = argv[3];
   const char *sample=argv[4];
   const char *phoID = argv[5];
   //TString pho_ID = phoID;





   AnalyzeLightBSM ana(inputFileList, outFileName, data,sample, phoID);
   cout << "dataset " << data << " " << endl;
   cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
   ana.EventLoop(data,inputFileList,sample,outFileName,phoID);
   Tools::Instance();
   return 0;
 }

 void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName,  const char* phoID) {
   cout<<"inside event loop"<<endl;
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nentries " << nentries << endl;
   cout << "Analyzing dataset " << data << " " << endl;
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
   bool FR_elec_flag=false;
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
   bool higMET=false, highdphi=true;
   double deepCSVvalue=0,p0=0,p1=0,p2=0;
   bool applyTrgEff=true;
   bool applyHEMveto=true;
   bool applyL1TrigFire_cal=true;
   bool applyL1TrigFire_prob=true;
   bool applyPUwt=true;
   bool apply_pixelveto=false;
   bool check_flag_wBL =false;
   if(s_sample.Contains("UL")){
   if(s_data.Contains("2016preVFP")){ lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV
   if(s_data.Contains("2016postVFP")) { lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016

   if(s_data.Contains("2017")) {lumiInfb=41.48; deepCSVvalue = 0.4506;}
   if(s_data.Contains("2018")){ lumiInfb=59.83;deepCSVvalue = 0.4168;}
   if(s_data.Contains("signal"))lumiInfb= 137.19;

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
   }
   if(s_sample.Contains("data") && s_sample.Contains("UL"))
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
   int coutt=0;
   double  count_nelec=0.0;
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
       if(Debug)         cout<<"===load tree entry ==="<<jentry<<endl;
       wt =1.0;
       h_selectBaselineYields_v1->Fill("no cuts",wt);
       if(s_Process.Contains("2018.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2018.WGJets_MonoPhoton_PtG-130UL")|| s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-40to130UL") ||s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-130UL") || s_Process.Contains("2017.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2017.WGJets_MonoPhoton_PtG-130UL") ||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-130UL")||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-40to130UL"))//|| s_Process.Contains("2018.TTGJets_incUL")) //(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018") || s_sample.Contains("WGJets"))
       	{	  
 	  if(jentry==0)
       	    cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_process<<endl;
       	  wt = cross_section*lumiInfb*1000.0;//)/nentries; // Weight*lumiInfb*1000.0; //(cross_section*lumiInfb*1000.0)/nentries;
 	  if(Debug &&jentry<10) cout<<"Alpana urgent check "<<wt<<"\t"<<cross_section<<"\t"<<lumiInfb<<endl;
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
	}
      bool tighte_trgpass=true;
      bool tightmu_trgpass=true;
      bool elec=false, muon=false, elecmu=false;
      if(NElectrons==1 && NMuons==1){ elec=false; muon=false;elecmu=true;}
      else if(NElectrons==2 && NMuons==0){ elec=true; muon=false; elecmu=false;}
      else if(NElectrons==0 && NMuons==2){ elec=false; muon=true;elecmu=false;}
      else { elec=false; muon=false;elecmu=false;}
      if(s_sample.Contains("data") && applyTrgEff)
        {
          wt=1;

          if(TriggerPass->size()!=213) continue;
	  if(elec && !muon && !elecmu && NElectrons==2){

	  if((*TriggerPass)[30]==1 || (*TriggerPass)[31]==1 || (*TriggerPass)[32]==1 || (*TriggerPass)[33]==1 || (*TriggerPass)[45]==1 || (*TriggerPass)[47]==1 || (*TriggerPass)[42]==1 || (*TriggerPass)[46]==1 || (*TriggerPass)[51]==1 || (*TriggerPass)[52]==1 || (*TriggerPass)[41]==1 || (*TriggerPass)[22]==1 || (*TriggerPass)[24]==1)
	    tighte_trgpass=true;
	  else tighte_trgpass=false;
	  }
	  
	  else	  if(NMuons==2 && !elec && muon && !elecmu){
	    if((*TriggerPass)[65]==1 || (*TriggerPass)[66]==1 || (*TriggerPass)[67]==1 || (*TriggerPass)[68]==1 || (*TriggerPass)[69]==1 || (*TriggerPass)[70]==1 || (*TriggerPass)[71]==1 || (*TriggerPass)[72]==1 || (*TriggerPass)[73]==1 || (*TriggerPass)[98]==1 || (*TriggerPass)[101]==1 || (*TriggerPass)[102]==1 || (*TriggerPass)[28]==1 || (*TriggerPass)[29]==1 || (*TriggerPass)[81]==1 || (*TriggerPass)[82]==1 || (*TriggerPass)[85]==1 || (*TriggerPass)[86]==1 )
	      tightmu_trgpass=true;
	    else tightmu_trgpass=false;
	    
	  }

	  else if (NMuons==1 && NElectrons ==1 &&  !elec && !muon && elecmu){

	    if((*TriggerPass)[65]==1 || (*TriggerPass)[66]==1 || (*TriggerPass)[67]==1 || (*TriggerPass)[68]==1 || (*TriggerPass)[69]==1 || (*TriggerPass)[70]==1  || (*TriggerPass)[71]==1 || (*TriggerPass)[72]==1 || (*TriggerPass)[73]==1 || (*TriggerPass)[98]==1 || (*TriggerPass)[101]==1 || (*TriggerPass)[102]==1  || (*TriggerPass)[81]==1 || (*TriggerPass)[82]==1 || (*TriggerPass)[85]==1 || (*TriggerPass)[86]==1 ) // muon triggers
              tighte_trgpass=true;
	    else    if((*TriggerPass)[30]==1 || (*TriggerPass)[31]==1 || (*TriggerPass)[32]==1 || (*TriggerPass)[33]==1 || (*TriggerPass)[45]==1 || (*TriggerPass)[47]==1 || (*TriggerPass)[42]==1 || (*TriggerPass)[46]==1 || (*TriggerPass)[51]==1 || (*TriggerPass)[52]==1 || (*TriggerPass)[41]==1 || (*TriggerPass)[22]==1 || (*TriggerPass)[24]==1) // electron triggers
	      tightmu_trgpass=true;
	    else tightmu_trgpass=false;

	  }
          // if((*TriggerPass)[46]==1 && s_data.Contains("2016"))
          //   tighte_trgpass=true;
          // else if((*TriggerPass)[51]==1 && !s_data.Contains("2016"))
          //   {
          //     if((*TriggerPrescales)[51]!=1)     cout<<jentry<<" :  TriggerPrescales = "<<(*TriggerPrescales)[51]<<endl;
          //     tighte_trgpass=true;
          //   }
	  else { tightmu_trgpass=false;	  tighte_trgpass=false;}
	}      
      if((s_data.Contains("2016") || s_data.Contains("2017")) && applyL1TrigFire_prob){ //L1 prefire issue
	wt= wt*NonPrefiringProb;
      }
     
      if(!s_sample.Contains("data") && applyPUwt)
	wt = wt*puWeight;

      
      // if(!s_sample.Contains("data") && !s_sample.Contains("signal") && applyTrgEff )
      // 	{
      // 	  wt = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
      // 	}
      if(Debug)
	cout<<"alpa"<<"\t"<<jentry<<"\t"<<endl;
      if(jentry<10)
	{    
	  cout<<"Event "<<jentry<<" event weight"<<"\t"<<wt<<"\t"<<cross_section<<"\t"<<Weight<<"lumiInfb  "<<lumiInfb<<endl;
	  cout<<"Event "<<jentry<<"\t"<<wt<<"\t"<<"MET  "<<MET<<"\t"<<(((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2)<<endl;
	}
           //      h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
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
      	HLTElectronObjects_v1.clear();
      	TAPElectronTracks_v1.clear();
      	if(Debug)
      	  cout<<jentry<<"\t"<<"After setting pointer  lorentz vectors"<<endl;
      	Electrons_v1 = getLorentzVector(Electrons_,Electrons_fCoordinates_fPt,Electrons_fCoordinates_fEta,Electrons_fCoordinates_fPhi,Electrons_fCoordinates_fE);
      	Muons_v1 = getLorentzVector(Muons_,Muons_fCoordinates_fPt,Muons_fCoordinates_fEta,Muons_fCoordinates_fPhi,Muons_fCoordinates_fE);
      	Photons_v1 = getLorentzVector(Photons_,Photons_fCoordinates_fPt,Photons_fCoordinates_fEta,Photons_fCoordinates_fPhi,Photons_fCoordinates_fE);
      	Jets_v1 = getLorentzVector(Jets_,Jets_fCoordinates_fPt,Jets_fCoordinates_fEta,Jets_fCoordinates_fPhi,Jets_fCoordinates_fE);
      	// if(s_sample.Contains("data"))
      	//    HLTElectronObjects_v1 = getLorentzVector(HLTElectronObjects_,HLTElectronObjects_fCoordinates_fPt,HLTElectronObjects_fCoordinates_fEta,HLTElectronObjects_fCoordinates_fPhi,HLTElectronObjects_fCoordinates_fE);
      	TAPElectronTracks_v1 = getLorentzVector(TAPElectronTracks_,TAPElectronTracks_fCoordinates_fPt,TAPElectronTracks_fCoordinates_fEta,TAPElectronTracks_fCoordinates_fPhi,TAPElectronTracks_fCoordinates_fE);
	
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
      int ele_branch=Electrons_v1.size();
      if(Debug)
        cout<<"Electrons "<<"\t"<<Electrons_<<"\t"<<Muons_<<"\t"<<Jets_<<endl;
      //      int total_lost_el = 0,cr_el=0,sr_el;
      int pho_branch = Photons_v1.size();
      if(Debug)
    	        cout<<"Photons "<<"\t"<<Photons_<<endl;
      //      h_selectBaselineYields_v1->Fill("no cuts",wt);
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
     ///////////////////////////////////   Making gen level collections of all leptons  ////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(Debug)
      cout<<"===load tree entry check1 - after filling gen level at entry ==="<<"\t"<<jentry<<endl;           
           
      int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
      TLorentzVector genPho1,genEle1,genMu1,genTau1,recElec, recMu;//,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
      vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_genLep1,v_recEle,v_recMu;
      v_genEle2.clear();
      v_genEle1.clear();
      v_genMu1.clear();
      v_genMu2.clear();
      v_genTau1.clear();
      v_genTau2.clear();
      bool HEMaffected=false, extracondition=false;
      if(s_data.Contains("2018") && applyHEMveto){
        for(int i=0; i<Electrons_v1.size();i++){
	  if((*Electrons_passIso)[i]==1) {
          if(Electrons_v1[i].Pt() >30 && Electrons_v1[i].Eta() > -3.0 && Electrons_v1[i].Eta() < -1.4 && Electrons_v1[i].Phi() > -1.57 && Electrons_v1[i].Phi() < -0.87) {HEMaffected = true; break;}
	  }
        }
        for(int i=0; i<Jets_v1.size();i++){
          if(Jets_v1[i].Pt() > 30 && Jets_v1[i].Eta() > -3.2 && Jets_v1[i].Eta() < -1.2 && Jets_v1[i].Phi() > -1.77 && Jets_v1[i].Phi() < -0.67 && DeltaPhi(Jets_v1[i].Pt(),METPhi)<0.5) {HEMaffected = true; break;}
        }
	for(int i=0; i<Muons_v1.size();i++)
	  {
	    if((*Muons_passIso)[i]==1)
	      {
		if(Muons_v1[i].Pt() >30 && Muons_v1[i].Eta() > -3.0 && Muons_v1[i].Eta() < -1.4 && Muons_v1[i].Phi() > -1.57 && Muons_v1[i].Phi() < -0.87) {HEMaffected = true; break;}
	      }
	  }
        if(HEMaffected == true) continue;
      }

      // if(s_data.Contains("2018") && applyHEMveto){
      // 	for(int i=0; i<Electrons->size();i++)
      // 	  if((*Electrons)[i].Pt() >30 && (*Electrons)[i].Eta() > -3.0 && (*Electrons)[i].Eta() < -1.4 && (*Electrons)[i].Phi() > -1.57 && (*Electrons)[i].Phi() < -0.87)   {HEMaffected = true; break;}
      // 	for(int i=0; i<Jets->size();i++)
      // 	  if((*Jets)[i].Pt() > 30 && (*Jets)[i].Eta() > -3.2 && (*Jets)[i].Eta() < -1.2 && (*Jets)[i].Phi() > -1.77 && (*Jets)[i].Phi() < -0.67 && DeltaPhi((*Jets)[i].Pt(),METPhi)<0.5) {HEMaffected = true; break;}
      // 	if(HEMaffected == true) continue;
      // }
      h_selectBaselineYields_v1->Fill("HEM veto ",wt);
      //recommended MET filtser for UL - taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_should
      //check nvtx filter - less efficiency for 2016 - https://indico.cern.ch/event/1057110/#27-met-filters-performance-stu
      if(PFCaloMETRatio >=  5) continue;
      if (s_data.Contains("2017") || s_data.Contains("2018"))
	if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 && NVtx>0)) continue;
      if (s_data.Contains("2016")){
	if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 &&EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 )) continue;
      }
      
      if(MET/CaloMET > 2.0) continue;
      
      h_selectBaselineYields_v1->Fill("MET filters ",wt);
      int hasEle=0,hasEle2=0, hasPho=0, hasPho_px=0,hastagEle=0,hasprobeEle=0;
      TLorentzVector bestPhoton=getBestPhoton(pho_ID);
      TLorentzVector bestPhoton_woFullID = getPhoton_withoutfullID();
      bool bestPhoHasPxlSeed=true, noPhoHasPxlSeed=true;
      if(bestPhotonIndxAmongPhotons>=0){
        if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
	if( bestPhoHasPxlSeed ) continue;
      }
      //   if(!bestPhoHasPxlSeed && bestPhoton.Pt()>40)
      //     {
      //       hasPho=1;hasPho_px=0;
      //     }

      //   else if(bestPhoHasPxlSeed && bestPhoton.Pt()>40) {hasPho_px=true;hasPho=0;}// npho_px++;}                                                                   
      //   else
      //     {
      //       hasPho_px=0;
      //       hasPho=0;
      //     }
      // }

      double mindr_Pho_genlep=getGenLep(bestPhoton);
      
      if(!s_sample.Contains("data")){
	if(Debug)
	  cout<<"===loading tree entry === "<<jentry<<"\t"<<"  ===after preselection  ==="<<endl;
	h_madminPhotonDeltaR_beforeStitching->Fill(madMinPhotonDeltaR,wt);
	h_phoPt_promptPho->Fill(mindr_Pho_genlep,wt);
	//remove the overlapped events from different MC samples as discussed in these slides - https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf                                                                                                                                 
	if((s_sample.Contains("TTJets_HT")||s_sample.Contains("TTJets-HT")) && madHT<600) continue;
	if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept")) && madHT>600) continue;
	if(!genphocheck)      {        genphomatch_before++;
	  //double mindr_Pho_genlep=getGenLep(bestPhoton);                                                                                                                  
	  if( s_sample.Contains("TTG") )
	    {  if(!hasGenPromptPhoton)              {
                h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;          }
	      else if(hasGenPromptPhoton)              {
		h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
		if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                            
		  if(madMinPhotonDeltaR<0.5)                  h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
                  if(mindr_Pho_genlep<0.5)                  h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);
                  continue;}
                else
                  {                 if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);}
	      }
	    }
	  if(s_sample.Contains("WG"))
	    {            if(!hasGenPromptPhoton){                                 h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;              }
	      else if(hasGenPromptPhoton)              {        h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
		if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                             
                  if(madMinPhotonDeltaR<0.5)                    h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
                  if(mindr_Pho_genlep<0.5)                    h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);
                  continue;}
                else
                  {           if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);
                  }}}
	  if(s_sample.Contains("WJets"))
	    { if(!hasGenPromptPhoton){h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
		if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}
	      else if(hasGenPromptPhoton){  h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                             
                  continue;}
                else{
		  if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
		  if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);
                }}}

	  if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets") || s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept"))
	    {            if(!hasGenPromptPhoton){
                h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}
	      else if(hasGenPromptPhoton){
                h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){// h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                            
                  continue;}
                else{
		  if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
		  if(mindr_Pho_genlep >=0.5)                      h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);}                          }}
	  if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
	    {            if(!(madMinPhotonDeltaR>0.4)) continue;}

	  if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
	    {            if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;}
	  if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG")) || s_sample.Contains("ZNuNuGJets")))          {
	    if(!(madMinPhotonDeltaR>0.5)) continue;}

	  if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
	    {            if(!(madMinPhotonDeltaR<=0.5)) continue;}
	  if(hasGenPromptPhoton && (s_sample.Contains("DYJets")))
	    {
	      if(!(madMinPhotonDeltaR<=0.5 || mindr_Pho_genlep<=0.5) ) continue;

	    }
	  if(hasGenPromptPhoton && (s_sample.Contains("ZLLG")))
            {
	      if(!(madMinPhotonDeltaR>0.5 && mindr_Pho_genlep>0.5)) continue;//if(!(madMinPhotonDeltaR<=0.5 || mindr_Pho_genlep<=0.5) ) continue;

            }

	  genphomatch_after++;
	}
	h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt);
	h_phoPt_promptPho_after->Fill(mindr_Pho_genlep,wt);
      }
      h_selectBaselineYields_v1->Fill("overlap removal",wt);

      if(s_sample.Contains("data") && !s_sample.Contains("UL") && applyTrgEff){
	if(tighte_trgpass==false || tightmu_trgpass==false)  continue;
      }
      h_selectBaselineYields_v1->Fill("tight trigger",wt);

      //veto muons
      // if(NMuons==0) h_selectBaselineYields_v1->Fill("0 #mu",wt);
      // else continue;
      // //apply iso track veto                                                                                                                                          
      if(isoPionTracks==0)
	{
	  h_selectBaselineYields_v1->Fill("Iso  pion track",wt);
	}
      else continue;

      int leadGenPhoIdx=-100;//, mu_index=-100;
      int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
      int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
      if(Debug)
    	cout<<"=== load entry "<<jentry<<"\t"<<"GenTaus_  "<<GenTaus_<<"\t"<<"had tau decay"<<"\t"<<nGentau_had1<<endl;

      bool hasEMobj_pho=false, hasEMobj_pho_px=false,hasEMobj_elec=false; //flags - true if a good photon or good electron is present in the event 
      TLorentzVector emobj_pho, emobj_pho_px,emobj_elec;
      // int  hasPho=0;
      // int hasPho_px=0;
      // int hasEle=0;
      double mt_ele=0.0;
      int total_lost_el = 0,cr_el=0,sr_el,e_index=-1,nlep=0, NgenElec=0,leadGenPhoIdx1=0,lep1=0,lep3=0,zwin=0;
      bool Elec_passEtacut=false,Elec_passpTcut=false,Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false, genElec_passEtacut=false,genElec_passpTcut=false,genElec_passAccep = false, genElec_passId= false, genElec_passIso = false,genMu_passEtacut=false,genMu_passpTcut=false,genMu_passAccep = false, genMu_passId= false, genMu_passIso = false;
      bool Mu_passEtacut=false,Mu_passpTcut=false, Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
      int fail_isoElec=0;
      nlep=0; nelec_reco=0;
      sortTLorVec(&Electrons_v1);
      sortTLorVec(&Muons_v1);
      int e_index2=-1, two_e=0, one_e=0, one_e2=0,two_e2=0,two_e3=0, mu_index=-1,mu_index2=-1;
      //    TLorentzVector //lep2;
      //      TLorentzVector tmp,tmp2;
      TLorentzVector genzvec,zvec,muon1,muon2,elec1,elec2,vlep1,vlep2,vgenl1,vgenl2, v_leadLept;
      double zmass= -10,nmuons=0,nelectrons=0,dR1=999,dR2=999;
      if(bestPhoton.Pt()>40)
	h_selectBaselineYields_v1->Fill("bestPhoton.Pt()>40",wt);
      else continue;
      if(Debug)
	cout<<"Entrye: "<<jentry<<" bestPhotonpT >40 "<<endl;
      if(s_sample.Contains("ZLLGJets") || s_sample.Contains("DYJetsToLL") || s_sample.Contains("data")  || s_sample.Contains("WGJets") || s_sample.Contains("TTGJets") || s_sample.Contains("TTJets") ){
	if(NElectrons==2 || NElectrons==1){
	  for(int i=0 ; i<Electrons_v1.size() ; i++)
	    {
	      if((*Electrons_passIso)[i]==1 )
		{
		  nlep++;
		  if(nlep==1) e_index=i;
		  if(nlep==2) e_index2=i;
		}
	    }
	}
	if(NMuons==2 || NMuons==1){
	  for(int i=0 ; i<Muons_v1.size() ; i++)
 	    {
	      if((*Muons_passIso)[i]==1 )
		{
		  nlep++;
		  if(nlep==1)  e_index=i;
		  if(nlep==2)  e_index2=i;
		}
	    }
	}

      }
      
      if(Debug)
	cout<<"Entrye: "<<jentry<<" 1 or 2 e/mu "<<endl;

      if(s_sample.Contains("ZLLGJets") || s_sample.Contains("DYJetsToLL") || s_sample.Contains("data")  || s_sample.Contains("WGJets") || s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")){
      
	if(jentry<3) cout<<"Processing as Z(G)->LL(G) Decay"<<endl;
	if(Debug)
	  cout<<"Entrye: "<<jentry<<" inside CR loop "<<NElectrons<<"\t"<<NMuons<<"\t"<<(*Muons_passIso)[0]<<"\t"<<nlep<<"\t"<<mu_index<<endl;

	if(NElectrons==2 && (s_sample.Contains("Electron") || s_sample.Contains("EGamma"))){
	  //  if(!(elec && !muon) )continue;
	  //h_selectbaselineyields_v1
	  if(Debug)
	    cout<<"Entrye: "<<jentry<<" inside CR loop "<<endl;
	  h_selectBaselineYields_v1->Fill("neve with 2e",wt);
	   if(! tighte_trgpass) continue;
	   h_selectBaselineYields_v1->Fill("neve with 2e - not passing trigger",wt);
	   if(Debug)
	     cout<<"Entrye: "<<jentry<<" tighepass "<<endl;
	  if( (NMuons!=0) ) continue;//exactly 0 RECO muons and 2 RECO e      
	  h_selectBaselineYields_v1->Fill("neve with 2e - #mu veto",wt);
	  if(isoMuonTracks!=0) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2e - iso #mu tracks veto",wt);
	  if((*Electrons_passIso)[e_index]!=1 || (*Electrons_passIso)[e_index2]!=1) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2e - After iso e",wt);
	  if( (*Electrons_charge)[e_index] == (*Electrons_charge)[e_index2] ) continue;//opp sign same flav(OSSF)                 
	  h_selectBaselineYields_v1->Fill("neve with 2e - opp sign flavour",wt);
	  if(e_index!=1 &&  e_index!=0)cout<<jentry<<" : "<<e_index<<" : "<<e_index2<<endl;               
	  elec1=Electrons_v1[e_index];
	  elec2=Electrons_v1[e_index2];
	  if(elec1.Pt()<20 || elec2.Pt()<20) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2e - both e pT>=20",wt);
	  dR1=bestPhoton.DeltaR(Electrons_v1[e_index]);
	  dR2=bestPhoton.DeltaR(Electrons_v1[e_index2]);
	  if(!(dR1>0.2 && dR2>0.2)) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2e - dR(e,#gamma)>0.2",wt);
	  //	  bool eMatchedG=check_eMatchedtoGamma(bestPhoton);                                                                                                                   
	  // if(eMatchedG) continue;                                                                                                                                        
	  nelectrons=2;
	  //zvec =  ((*Electrons)[0] + (*Electrons)[1]);
	  zvec =  (elec1 + elec2);
	  zmass = zvec.M();
	  h_InvariantMass_noselec->Fill(zmass,wt);
	  nmuons=0;
	  vlep1=elec1;
	  vlep2=elec2;
      	  lep1++;
	  if(Debug)
	    cout<<"Entrye: "<<jentry<<" just last line in  CR loop "<<endl;

	}
	else if(NMuons==2 && s_sample.Contains("Muon")) {
	  //if(!(!elec && muon) )continue;       
	  //       if(!(!elec && muon && tightmu_trgpass) )continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu",wt);
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" 2 muons in  CR loop "<<endl;


	  if(!tightmu_trgpass )continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu - not passing trigger",wt);
	  //       if(isoMuonTracks!=0) continue;
	  if(NElectrons!=0 ) continue;//exactly 2 RECO muons and 0 RECO e                                                           
	  h_selectBaselineYields_v1->Fill("neve with 2mu - e veto",wt);
	  if(isoElectronTracks!=0) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu - iso e track veto",wt);
	  if((*Muons_passIso)[e_index]!=1 || (*Muons_passIso)[e_index2]!=1) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu - After iso e",wt);
	  if( (*Muons_charge)[e_index] == (*Muons_charge)[e_index2] ) continue;//opp sign same flav(OSSF)                                        
	  h_selectBaselineYields_v1->Fill("neve with 2mu - opp sign flavour",wt);
	  muon1=Muons_v1[e_index];
	  muon2=Muons_v1[e_index2];
	  if(muon1.Pt()<20 || muon2.Pt()<20) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu - both e pT>20",wt);
	  dR1=bestPhoton.DeltaR(Muons_v1[e_index]);
	  dR2=bestPhoton.DeltaR(Muons_v1[e_index2]);
	  //	  bool muMatchedG=check_muMatchedtoGamma(bestPhoton);                                                                                               
	  if(!(dR1>0.2 && dR2>0.2)) continue;
	  h_selectBaselineYields_v1->Fill("neve with 2mu - dR(e, #gamm)>0.2",wt);
	  nmuons=2;
	  //zvec =  ((*Muons)[0] + (*Muons)[1]);
	  zvec =  muon1 + muon2;
	  zmass = zvec.M();
	  h_InvariantMass_noselec->Fill(zmass,wt);
	  nelectrons=0;
	  vlep1=muon1;
	  vlep2=muon2;
	  lep2++;
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" just last line in muons  CR loop "<<endl;

	}

	else if(NMuons==1 && NElectrons==1 && s_sample.Contains("Muon")){
	  //if(!(tighte_trgpass || tightmu_trgpass)) continue;                                                    
	  if(! tightmu_trgpass) continue;
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" in the ne==1 and nmu==1"<<"\t"<<e_index<<"\t"<<mu_index<<endl;
	  //exactly 1 RECO muons and 1 RECO e                                                           
	  h_selectBaselineYields_v1->Fill("neve with 1e, 1#mu",wt);
	  //	  if(mu_index==-1 or e_index==-1) continue;
	  if((*Electrons_passIso)[e_index]!=1 || (*Muons_passIso)[e_index]!=1) continue;
	  h_selectBaselineYields_v1->Fill("neve with 1e, 1#mu, not the first index",wt);
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" in the ne==1 and nmu==1 v2 "<<endl;

	  if( (*Electrons_charge)[e_index] == (*Muons_charge)[e_index] ) continue;//opp sign same flav(OSSF)                                        
	  h_selectBaselineYields_v1->Fill("neve with 1e, 1#mu, opp sign flavour",wt);
	  elec1=Electrons_v1[e_index];
	  muon1=Muons_v1[e_index];
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" in the ne==1 and nmu==1 v1"<<endl;

	  if(elec1.Pt()<20 || muon1.Pt()<20) continue;
	  h_selectBaselineYields_v1->Fill("neve with 1e, 1#mu, both lep pT>=20",wt);
	  dR1=bestPhoton.DeltaR(Electrons_v1[e_index]);
	  dR2=bestPhoton.DeltaR(Muons_v1[e_index]);
	  //	  bool muMatchedG=check_muMatchedtoGamma(bestPhoton);                                                                                                                
	  //       if(muMatchedG) continue;                                                                                                                                  
	  if(!(dR1>0.2 && dR2>0.2)) continue;
	  h_selectBaselineYields_v1->Fill("neve with 1e, 1#mu, dR(e/#mu, #gamma)>0.2",wt);
	  //zvec =  ((*Muons)[e_index] + (*Muons)[1]);
	  zvec =  muon1 + elec1;
	  zmass = zvec.M();
	  //       cout<<jentry<<" : "<<zmass<<endl;
	  nelectrons=1;
	  nmuons=1;
	  vlep1=elec1;
	  vlep2=muon1;
	  lep3++;
	  if(Debug)
            cout<<"Entrye: "<<jentry<<" just last line in  1 e and 1 mu CR loop "<<endl;


	}
	else continue;
	h_selectBaselineYields_v1->Fill("neve with no lep",wt);
      
      //if(s_sample.Contains("elemu")){
	//   if(!(nelectrons==1 && nmuons==1)) continue;
	// }
	// else{
	if(!(nelectrons==2 || nmuons==2 || (nelectrons==1 && nmuons==1))) continue;
	  //}
	h_selectBaselineYields_v1->Fill("neve with no lep1",wt);
	lep++;
    
	// if(nelectrons==2)  {  h2_ElectronEta_Phi->Fill(elec1.Eta(),elec1.Phi(),wt); h2_leadElectronEta_Phi->Fill(elec2.Eta(),elec2.Phi(),wt);}
	// if(nmuons==2){    h2_ElectronEta_Phi->Fill(muon1.Eta(),muon1.Phi(),wt); h2_leadElectronEta_Phi->Fill(muon2.Eta(),muon2.Phi(),wt);}
	if(zvec.M() < 80 || zvec.M() > 100) continue;
	zwin++;
	h_selectBaselineYields_v1->Fill("80<= Minv <=100",wt);
	//	h_zpt->Fill(zvec.Pt(),wt);
	//	h_invariantMass[2]->Fill(invariantmass,wt); 
	if(MET>200) continue;
	h_selectBaselineYields_v1->Fill("MET<200",wt);
	if(Debug)
	  cout<<"Entrye: "<<jentry<<" just last line in  final CR "<<endl;
      }

      TLorentzVector metvec, metstar;
      metstar = metvec;
      metvec.SetPtEtaPhiE(MET,0.0,METPhi,0.0);
      metstar = metvec + zvec;
      int genzvec_id=-1;
      TLorentzVector tmp;
      if(s_sample.Contains("ZNuNuGJets") || s_sample.Contains("ZJetsToNuNu_HT") || s_sample.Contains("ZGTo2NuG")){
	if(jentry<3) cout<<"Processing as Z(G)->NuNu(G) Decay"<<endl;
	h_selectBaselineYields_v1->Fill("Events in SR",wt);
	if(isoElectronTracks!=0 || isoMuonTracks!=0 || isoPionTracks!=0) continue;
	h_selectBaselineYields_v1->Fill("SR: veto iso e/#mu/#pi tracks",wt);
	if(NElectrons!=0 || NMuons!=0) continue;
	h_selectBaselineYields_v1->Fill("SR: veto e/#mu",wt);
	bool genpho=true;
	//  else cout<<NElectrons<<endl;                                                                                                                                 
	for(int i=0;i<GenParticles_v1.size();i++){
	  if( abs((*GenParticles_PdgId)[i])==23)
	    {
	      genpho=false;
	      genzvec_id=i;
	      tmp = (GenParticles_v1[i]);
	      //cout<<"Parent id : "<<abs((*GenParticles_ParentId)[i])<<" status : "<<(*GenParticles_Status)[i]<<"\t"<<tmp.M()<<endl;                                                  
	      if(tmp.M()>80 && tmp.M()<100)// && tmp.Pt()>200)                                                                                                           
		{
		  genzvec=tmp;
		  //		  h_ZpT[3]->Fill(genzvec.Pt(), wt);
		  // h_zpt->Fill(genzvec.Pt(),wt);
		}
	    }
	}
	//if(!genpho) continue; // check it later what it is doing 
      }
      
      if( metstar.Pt() > 200 )    
	{    h_selectBaselineYields_v1->Fill("MET* > 200",wt);
	}
      else  continue;

      if(Debug)
	cout<<"Entrye: "<<jentry<<" after met star cut "<<endl;


      MET = metstar.Pt();
      /////======Pixel Power supply issue ===============//////////////                                                                                                    // if(apply_pixelveto){
      // 	if(s_data.Contains("2017") && ((bestEMObj.Eta() > 0.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 2.7) || (tagEMObj.Eta() > 0.0 && tagEMObj.Eta() < 1.5 && tagEMObj.Phi() > 2.7))) continue;
      // 	if(s_data.Contains("2018") && ((bestEMObj.Eta() > 0.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 0.4 && bestEMObj.Phi() < 0.8) || (tagEMObj.Eta() > 0.0 && tagEMObj.Eta() < 1.5 && tagEMObj.Phi() > 0.4 && tagEMObj.Phi() < 0.8))) continue;
      // }

      // h_selectBaselineYields_v1->Fill("after pixel veto first",wt);

           
      // *******************  Selecting Jet objects ********************************//
      int minDRindx=-100,minDR2indx=-100, minDRindxl1 =-100, minDRindxl2=-100,phoMatchingJetIndx=-100,hadJetID=-999,bJet1Idx=-100,nHadJets=0, nbjets=0,photonMatchingJetIndx=-100;
      double qmulti=-1,leadjet_qmulti=-1,leadjet_Pt=-1,leadbjet_tag=-1 ,leadbjet_tag1=-1,leadbjet_tag2=-1,leadbjet_tag3=-1,qmulti1=-1,qmulti2=-1,qmulti3=-1;
      //double qmulti=0,leadjet_qmulti=-1,leadjet_Pt=-1,leadbjet_tag=-1;
      double minDR=99999,minDR2=9999,ST=0,Ht=0,remJetPt=0,minDRl1=0,minDRl2=0;
       vector<TLorentzVector> hadJets,bjets,remJets;
      hadJets.clear();
      bjets.clear();
      remJets.clear();
      vector<int> hadJets_hadronFlavor;
      vector<bool> hadJets_HTMask;
      vector<double> hadJets_bJetTagDeepCSVBvsAll;
      vector<TLorentzVector> nonbjets;
      vector<int> jetMatchindx;
      TLorentzVector Jet_matched;
      TLorentzVector leadJet;
      bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
      if(Debug)
        cout<<"===load tree entry check2 at entry ==="<<"\t"<<jentry<<endl;
      for(int i=0;i<Jets_;i++){
	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) )
	  {
            
	      double dR=bestPhoton.DeltaR(Jets_v1[i]);
	      //double dR2=tagEMObj.DeltaR(Jets_v1[i]);
	      //if(dR2<minDR2){minDR2=dR2;minDR2indx=i;}
	      if(dR<minDR){minDR=dR;minDRindx=i;}
	      if(nmuons==2){
		if( muon1.DeltaR(Jets_v1[i]) < minDRl1 ){ minDRl1 = muon1.DeltaR(Jets_v1[i]); minDRindxl1=i;}
		if( muon2.DeltaR(Jets_v1[i]) < minDRl2 ){ minDRl2 = muon2.DeltaR(Jets_v1[i]); minDRindxl2=i;}        
	      }
	      else if(nelectrons==2){ 
		if( elec1.DeltaR(Jets_v1[i]) < minDRl1 ){ minDRl1 = elec1.DeltaR(Jets_v1[i]); minDRindxl1=i;}
		if( elec2.DeltaR(Jets_v1[i]) < minDRl2 ){ minDRl2 = elec2.DeltaR(Jets_v1[i]); minDRindxl2=i;}
	      }
	      else if(nelectrons==1 && nmuons==1){
		if( elec1.DeltaR(Jets_v1[i]) < minDRl1 ){ minDRl1 = elec1.DeltaR(Jets_v1[i]); minDRindxl1=i;}
		if( muon1.DeltaR(Jets_v1[i]) < minDRl2 ){ minDRl2 = muon1.DeltaR(Jets_v1[i]); minDRindxl2=i;}
	      }
       
	    
	  }
      }
      if(Debug)
	cout<<"Entrye: "<<jentry<<"minDR condition jets selection "<<endl;


      // if( minDR > 0.3 )  minDRindx   = -100;
      // if( minDRl1 > 0.3) minDRindxl1 = -100;
      // if( minDRl2 > 0.3) minDRindxl2 = -100;

      int q=0;
      for(int i=0;i<Jets_;i++){
	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
	  //q++;
	  //if((*Jets_ID)[i]) h2_mindr_jetmatchphoratio->Fill(minDR,Jets_v1[minDRindx].Pt()/bestEMObj.Pt(),wt);
	  if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
	    //	    	  if( !(minDR < 0.3 && i==minDRindx) ){
	    //	    if((i!=minDRindx) && ){//&& (i!=minDR2indx) ){
	    if( !(minDR < 0.3 && i==minDRindx) && !(minDRl1 < 0.3 && i==minDRindxl1) && !(minDRl2 < 0.3 && i==minDRindxl2) ){
	      hadJetID=(*Jets_ID)[i];
	      if(hadJetID)
		{
		  q++;
		  hadJets.push_back(Jets_v1[i]);
		  hadJets_hadronFlavor.push_back((*Jets_hadronFlavor)[i]);
		  hadJets_HTMask.push_back((*Jets_HTMask)[i]);
		  hadJets_bJetTagDeepCSVBvsAll.push_back((*Jets_bJetTagDeepCSVBvsAll)[i]);
		  if(q==1) leadjet_qmulti=(*Jets_chargedMultiplicity)[q];
		  if(q==1) leadjet_Pt=Jets_v1[q].Pt();
		  if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
		    bjets.push_back(Jets_v1[i]); bJet1Idx = i;}
		}
	    }
	  else if(minDR<0.3 && i==minDRindx)
	    Jet_matched = Jets_v1[minDRindx];
	  }  
	}
      }
      //    cout<<jentry<<" ====== break ========"<<endl;
      if( minDR<0.3 )
	{
	  phoMatchingJetIndx=minDRindx;
	  leadbjet_tag1=(*Jets_bJetTagDeepCSVBvsAll)[phoMatchingJetIndx];
	  qmulti1=(*Jets_chargedMultiplicity)[phoMatchingJetIndx];
	}
      if( minDRl1<0.3 )
	{
	  leadbjet_tag2=(*Jets_bJetTagDeepCSVBvsAll)[minDRindxl1];
	  qmulti2=(*Jets_chargedMultiplicity)[minDRindxl1];
	}
      if( minDRl2<0.3 )
	{
	  leadbjet_tag3=(*Jets_bJetTagDeepCSVBvsAll)[minDRindxl2];
	  qmulti3=(*Jets_chargedMultiplicity)[minDRindxl2];
	}

      if(Debug)
	cout<<"Entrye: "<<jentry<<" Jets remvoing matching with others "<<endl;


      // if( photonMatchingJetIndx>=0 ){
      // 	if( (Jets_v1[photonMatchingJetIndx].Pt()) > 1.1*(bestEMObj.Pt()) ){
      // 	  if( ((Jets_v1[photonMatchingJetIndx] - bestEMObj).Pt())>30){
      // 	    hadJets.push_back( Jets_v1[photonMatchingJetIndx] - bestEMObj );
      // 	    remJetPt=(Jets_v1[photonMatchingJetIndx] - bestEMObj).Pt();
      // 	    remJets.push_back( Jets_v1[photonMatchingJetIndx] - bestEMObj );
      // 	  }
      // 	}
      // }


      for(int i=0;i<hadJets.size();i++)
	{
	  if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());}
	  if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
	}
    
      for(int i=0;i<bjets.size();i++){
	if( (abs(bjets[i].Eta()) < 2.4) ){nbjets++;}
      }
      if( minDR<0.3 ){
	ST=ST+bestPhoton.Pt();
      }

      // if( minDRl1<0.3 ){
      //   ST=ST+vlep1.Pt();
      // }
      // if( minDRl2<0.3 ){
      //   ST=ST+vlep2.Pt();
      // }


      if(phoMatchingJetIndx>=0 && (Jets_v1[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
      if(phoMatchingJetIndx<0) continue;     
      h_selectBaselineYields_v1->Fill("Jet matching photon pT cut",wt);
      // if( minDR2<0.3 ){
      //   ST=ST+tagEMObj.Pt();
      // }

     
      sortTLorVec(&hadJets);
      //      if(minDRindx<0) continue;
    if(hadJets.size()==0) continue;
    if(Debug)
      cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;
    BTags = nbjets;
     double genmindr=99999, recojetmindr=99999;
    if(Debug)
      cout<<"===load tree entry === "<<"\t"<<jentry<<" hadJets.size()  "<< hadJets.size()<<endl;    
    if(hadJets.size()<BTags)
      cout<<"Entry: "<<jentry<<"\t"<<"njets "<<nHadJets<<"\t"<<" hadJets.size() "<<hadJets.size()<<endl;
    // // ********* This is to account into all visible energy: Adding photon matching with Jet ******************//

    TLorentzVector Met; 
    Met.SetPtEtaPhiE(metstar.Pt(),0,metstar.Phi(),0);
    double mTPhoMET = sqrt(2*(bestPhoton.Pt())*metstar.Pt()*(1-cos(DeltaPhi(metstar.Phi(),bestPhoton.Phi()))));
    double dPhi_PhoMET= abs(metstar.DeltaPhi(bestPhoton));//DeltaPhi(METPhi,bestPhoton.Phi()));
     
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" mTPhoMET "<<mTPhoMET<<" METPhi  "<<bestPhoton.Eta()<<" Photons_ "<<Photons_<<" hadJets.size() " <<"\t"<<hadJets.size()<<" NElectrons "<<NElectrons<< "  Jets  "<<Jets_<<"\t"<< "Jets size  "<<Jets_v1.size()<<endl;
    if(jentry<10)
      cout<<"madMinPhotonDeltaR  "<<madMinPhotonDeltaR<<endl;
    //    h_madminPhotonDeltaR_noSelection->Fill(madMinPhotonDeltaR,wt);
			    		    
    // TLorentzVector Met;
    // Met.SetPtEtaPhiE(MET,0,METPhi,0);
    Double_t deta_jet_pho= 0.0,deta_jet_met=0.0,deta_met_pho=0.0;
    if(Debug)
      cout<<"METPhi "<<METPhi<<endl;
    double mT= 0.0, dPhi_METjet1=5, dPhi_METjet2=5, dPhi_phojet1=5, dPhi_phojet2=5, dPhi_phoMET=5;
    dPhi_METjet1 = abs(metstar.DeltaPhi(hadJets[0]));
    dPhi_METjet2 = abs(metstar.DeltaPhi(hadJets[1]));
     
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" hadJets.size() "<<hadJets.size()<<endl;
    // h_selectBaselineYields_->Fill("all",wt);
    
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" bestPhoHasPxlSeed "<<bestPhoHasPxlSeed<<endl;
    
    //cout<<"Event : " <<jentry <<"\t"<<(*Photons_mvaValuesID)[bestPhotonIndxAmongPhotons] << " (*Photons_mvaValuesID)[iPho]"<< "\t"<<bestPhotonIndxAmongPhotons<< " Photons _ "<<Photons_ <<"\t Photons_mvaValuesID->size() "<<Photons_mvaValuesID->size()<<"\t"<<Photons_v1.size()<<"\t"<<Photons_fullID->size()<<"\t"<<Photons_cutBasedID->size()<<endl;
    if (nHadJets>=2)h_selectBaselineYields_v1->Fill("njets>=2",wt);
    else continue;
    // if(MET<200) 
    //   h_selectBaselineYields_v1->Fill("MET<200",wt);
    // else continue;
    if(ST>300)
      h_selectBaselineYields_v1->Fill("ST>300",wt);
    else continue;
    // //recommended MET filtser for UL - taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_should
    // //check nvtx filter - less efficiency for 2016 - https://indico.cern.ch/event/1057110/#27-met-filters-performance-stu
    if(highdphi){
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
        h_selectBaselineYields_v1->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    }
    bool process=true, cr=false, sr=false; 
    if(s_sample.Contains("ZNuNuGJets") || s_sample.Contains("ZJetsToNuNu") || s_sample.Contains("ZGTo2NuG")) 
      {    
	if(highdphi)  process = process &&  ST>300 && metstar.Pt()>200 && nHadJets >=2 && dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 && bestPhoton.Pt() > 40 ;
	if(!highdphi)  process = process &&  ST>300 && metstar.Pt()>200 && nHadJets >=2 && (!(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3)) && bestPhoton.Pt() > 40 ;
	// if(!highdphi)  process = process &&  ST>300 && metstar.Pt()>200 && nHadJets >=2 &&  bestPhoton.Pt() > 100 ;

	sr=true;
      } 
    if(s_sample.Contains("ZLLGJets") || s_sample.Contains("DYJetsToLL") || s_sample.Contains("data")   || s_sample.Contains("WGJets") || s_sample.Contains("TTGJets") || s_sample.Contains("TTJets"))
      {
	if(highdphi) process = process &&  ST>300 && metstar.Pt()>200 && nHadJets >=2 && dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 && bestPhoton.Pt() > 40; 
	if(!highdphi)  process = process &&  ST>300 && metstar.Pt()>200 && nHadJets >=2 && (!(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3)) && bestPhoton.Pt() > 40 ;
	cr=true;
      }
    if(!(process)) continue;

    
    FillHistogram_Kinematics(1,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
    //    MET = metstar.Pt();

    // if(print_)
    //   {
    // 	cout<<"======================== hadJet collections before cleaning ================="<<endl;
    // 	for(int i=0;i<Jets->size();i++){
    // 	  if( ((*Jets)[i].Pt() > 30) && (abs((*Jets)[i].Eta()) <= 2.4) ){
    // 	    cout<<jentry<<" : "<<i<<" : JetID = "<<(*Jets_ID)[i]<<" Jet Pt = "<<(*Jets)[i].Pt()<<" Jet Eta = "<<(*Jets)[i].Eta()<<" Jet Phi = "<<(*Jets)[i].Phi()<<endl;
    // 	    cout<<jentry<<" : "<<i<<" , Jet pt = "<<(*Jets)[i].Pt()<<" dR(Jet,pho) = "<<bestPhoton.DeltaR((*Jets)[i])<<" dR(Jet,l1) = "<<vlep1.DeltaR((*Jets)[i])<<" dR(Jet,l2) = "<<vlep2.DeltaR((*Jets)[i])<<endl;       
    // 	  }}
    // 	cout<<"======================== had Ak4Jet collections after cleaning================="<<endl;
    // 	for(int i=0;i<hadJets.size();i++)
    // 	  cout<<jentry<<" : "<<i<<" , Jet pt = "<<hadJets[i].Pt()<<" dR(Jet,pho) = "<<bestPhoton.DeltaR(hadJets[i])<<" dR(Jet,l1) = "<<vlep1.DeltaR(hadJets[i])<<" dR(Jet,l2) = "<<vlep2.DeltaR(hadJets[i])<<endl;
	
    // 	cout<<"======================== lepton collections ================="<<endl;
    // 	cout<<jentry<<" , lepton 1 pt = "<<vlep1.Pt()<<" , lepton 2 pt = "<<vlep2.Pt()<<" , Z Pt = "<<zvec.Pt()<<endl;
    // 	cout<<jentry<<" , lepton 1 Eta = "<<vlep1.Eta()<<" , lepton 2 Eta = "<<vlep2.Eta()<<" , Z Eta = "<<zvec.Eta()<<endl;
    // 	cout<<jentry<<" , lepton 1 Phi = "<<vlep1.Phi()<<" , lepton 2 Phi = "<<vlep2.Phi()<<" , Z Phi = "<<zvec.Phi()<<endl;
    // 	cout<<jentry<<" dr b/w photon and l1 = "<<dR1<<" , dr b/w photon and l2 = "<<dR2<<" , dr b/w l1 and l2 = "<<vlep1.DeltaR(vlep2)<<endl;

    // 	cout<<"======================== Triggers collections ================="<<endl;
    // 	if(s_data.Contains("data"))
    // 	  for(int i=0; i<TriggerPass->size(); i++)
    // 	    if((*TriggerPass)[i]==1) cout<<jentry<<" : Trg id : "<<i<<endl;//" , trigName : "<<trigNames1[i]<<endl; \                

    // 	cout<<"======================== MET collections ================="<<endl;
    // 	cout<<jentry<<" , MET = "<<MET<<" , METPhi = "<<METPhi<<" vec(zvec + MET).Pt = "<<metstar.Pt()<<" , dPhi(MET,Pho) = "<<DeltaPhi(metstar.Phi(),bestPhoton.Phi())<<" , dPhi(MET,Jet1) = "<<DeltaPhi(metstar.Phi(),(hadJets)[0].Phi())<<" , dPhi(MET,Jet2) ="<<DeltaPhi(metstar.Phi(),(hadJets)[1].Phi())<<endl;

    // 	cout<<"======================== Photon collections ================="<<endl;
    //     cout<<jentry<<" , Photon Pt = "<<bestPhoton.Pt()<<" , Photon Phi = "<<bestPhoton.Phi()<<" , Photon Eta = "<<bestPhoton.Eta()<<endl;


    // 	cout<<"======================== had Ak8Jet collections ================="<<endl;
    // 	for(int i=0;i<JetsAK8->size();i++)
    // 	  if( ((*JetsAK8)[i].Pt() > 200.0) && (abs((*JetsAK8)[i].Eta()) <= 2.4) ){
    // 	    cout<<jentry<<" : "<<i<<" , Ak8JetID = "<<(*JetsAK8_ID)[i]<<" , Ak8Jet Pt = "<<(*JetsAK8)[i].Pt()<<" , Ak8Jet Eta = "<<(*JetsAK8)[i].Eta()<<" , Ak8Jet Phi = "<<(*JetsAK8)[i].Phi()<<endl;
    // 	    cout<<jentry<<" : "<<i<<" ,Ak8Jet pt = "<<(*JetsAK8)[i].Pt()<<" dR(ak8Jet,pho) = "<<bestPhoton.DeltaR((*JetsAK8)[i])<<" dR(ak8Jet,l1) = "<<vlep1.DeltaR((*JetsAK8)[i])<<" dR(Ak8Jet,l2) = "<<vlep2.DeltaR((*JetsAK8)[i])<<endl;   
    // 	    cout<<"nhadJets = "<<nHadJets<<" & Ak8 Jet Mass = "<<(*JetsAK8_softDropMass)[i]<<endl;
    // 	  }
    //   }
     	 //    h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt);
    if(Debug)
      cout<<" just before photon identification - prompt/non-prompt ========   ===="<<jentry<<endl;

    ////////////////////////////////////////////////////////////////////////////////////
    //  ====  ========== ======== 1 Electron CR + 0 photon  =====================  ///
    ///////////////////////////////////////////////////////////////////////////////////
    
    bool elec_CR= false; 
    bool pho_SR= false; 
    elec_CR = cr;
    pho_SR= sr;
    if(cr){//bestEMobj && bestEMObjIsEle && !hasPho ){
      // if(isoMuonTracks !=0 || isoPionTracks!=0) continue; // veto muon/pions from 1 electron CR
      // if(NMuons!=0 || NElectrons==0) continue;     
      // if(Debug)
      // 	cout<<"in Electron CR ====== "<<jentry<<endl;
      //      mt_ele=sqrt(2*bestEMObj.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestEMObj.Phi()))));
      //double mTElecMET=sqrt(2*(Electrons_v1[e_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,Electrons_v1[e_index].Phi()))));
      //if(mt_ele>100) { continue;}//h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);continue;} // remove signal contamination 
      elec_CR = true;      
      int searchBin = getBinNoV6_WithOnlyBLSelec(nHadJets,BTags,metstar.Pt());
      if(nHadJets>=7 && Debug)
	cout<<"Entry: "<<jentry<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<bestPhoton.Pt()<<"\t"<<wt<<endl;
      //if (bestEMObj.Pt()<=40) continue;
      h_selectBaselineYields_v1->Fill("Elec CR ",wt);
      //cout<<" Debugging CR "<<"\t"<<hasPho<<"\t"<<Electrons_v1[e_index].Pt()<<"\t"<<bestEMObj.Pt()<<"\t"<<MET<<"\t"<<nHadJets<<"\t"<<ST<<"\t"<<bestPhoton.Pt()<<"\t"<<hasPho_px<<"\t"<<isoPionTracks<<"\t"<<isoMuonTracks<<"\t"<<NMuons<<"\t"<<NElectrons<<"\t"<<isoElectronTracks<<"\t"<<(*Electrons_passIso)[e_index]<<"\t"<<mt_ele<<"\t"<<mTElecMET<<endl;
      //double wt_LL = wt*h_TF->GetBinContent(TFbins+1);
      //FillTFBins_Valid(6, nHadJets, BTags,wt_LL, wt_LL1, wt_LL2, bestPhoton.Pt(),ST);//
      // h_TFbins_ElecLL_validation[6]->Fill(TFbins+1,wt_LL);
      // h_TFbins_ElecLL_validation_v1[6]->Fill(TFbins+1,wt_LL);
      // h_Sbins_LL_Validation[6]->Fill(searchBin,wt_LL);
      //      cout<<"TFbins " <<TFbins<<"\t"<<h_TF->GetBinContent(TFbins+1)<<"\t"<<wt_LL<<endl;
      if(elec_CR && !elecmu) {    
	//	h_invariantMass_noCut[2]->Fill(zmassinvariantmass,wt);
		// if(!(invariantmass>=80 && invariantmass<=100) )continue; 
	// h_selectBaselineYields_v1->Fill("Elec CR - inv mass cut ",wt);

	FillHistogram_Kinematics(2,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	FillHistogram_Kinematics_varBin(2,metstar.Pt(),nHadJets, BTags, bestPhoton.Pt(),ST,qmulti,wt);
	h_invariantMass[2]->Fill(zmass,wt);
	h_ZpT[2]->Fill(zvec.Pt(), wt);
      }
      if(elec_CR){
	if(!elec && muon && !elecmu){
	  h_muon1_pT->Fill(v_lep1.Pt(),wt);
	  h_muon1_Eta->Fill(v_lep1.Eta(),wt);
	  h_muon1_Phi->Fill(v_lep1.Phi(),wt);
	  h_muon1_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
	  h_muon1_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
	  h_muon1_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);

	  h_muon2_pT->Fill(v_lep1.Pt(),wt);
          h_muon2_Eta->Fill(v_lep1.Eta(),wt);
          h_muon2_Phi->Fill(v_lep1.Phi(),wt);
          h_muon2_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
          h_muon2_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
          h_muon2_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);
	  FillHistogram_Kinematics(4,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	  if(BTags==0)
	    FR_nbtagBins[4]->Fill(1,wt);
	  else if(BTags>=1)
	    FR_nbtagBins[4]->Fill(2,wt);

	}
	if(elec && !muon && !elecmu){
	  h_electron1_pT->Fill(v_lep1.Pt(),wt);
          h_electron1_Eta->Fill(v_lep1.Eta(),wt);
          h_electron1_Phi->Fill(v_lep1.Phi(),wt);
          h_electron1_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
          h_electron1_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
          h_electron1_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);

          h_electron2_pT->Fill(v_lep1.Pt(),wt);
          h_electron2_Eta->Fill(v_lep1.Eta(),wt);
          h_electron2_Phi->Fill(v_lep1.Phi(),wt);
          h_electron2_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
          h_electron2_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
          h_electron2_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);
	  FillHistogram_Kinematics(5,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti,leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);

	  if(BTags==0)
	    FR_nbtagBins[5]->Fill(1,wt);
	  else if(BTags>=1)
	    FR_nbtagBins[5]->Fill(2,wt);

	}

	if(!elec && !muon && elecmu){
	  h_emu_electron1_pT->Fill(v_lep1.Pt(),wt);
          h_emu_electron1_Eta->Fill(v_lep1.Eta(),wt);
          h_emu_electron1_Phi->Fill(v_lep1.Phi(),wt);
          h_emu_electron1_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
          h_emu_electron1_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
          h_emu_electron1_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);

          h_emu_muon1_pT->Fill(v_lep1.Pt(),wt);
          h_emu_muon1_Eta->Fill(v_lep1.Eta(),wt);
          h_emu_muon1_Phi->Fill(v_lep1.Phi(),wt);
          h_emu_muon1_EtaVsPhi->Fill(v_lep1.Eta(), v_lep1.Phi(),wt);
          h_emu_muon1_PtVsPhi->Fill(v_lep1.Pt(), v_lep1.Phi(),wt);
          h_emu_muon1_PtVsEta->Fill(v_lep1.Pt(), v_lep1.Eta(),wt);
	  FillHistogram_Kinematics(6,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti,leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);

	  if(BTags==0)
	    FR_nbtagBins[6]->Fill(1,wt);
	  else if(BTags>=1)
	    FR_nbtagBins[6]->Fill(2,wt);

	}
	// h_tagEle_pT_Elec_CR->Fill(tagEMObj.Pt(),wt);
	// h_tagEle_Eta_Elec_CR->Fill(tagEMObj.Eta(),wt);
	// h_tagEle_Phi_Elec_CR->Fill(tagEMObj.Phi(),wt);
	// h_trackEle_pT_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),wt);
	// h_trackEle_Eta_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Eta(),wt);
	// h_trackEle_Phi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Phi(),wt);
	// h_tagEle_EtaVsPhi_Elec_CR->Fill(tagEMObj.Eta(),tagEMObj.Phi(),wt);
	// h_tagEle_PtVsPhi_Elec_CR->Fill(tagEMObj.Pt(),tagEMObj.Phi(),wt);
	// h_tagEle_PtVsPhi_Elec_CR->Fill(tagEMObj.Pt(),tagEMObj.Eta(),wt);
	
	// h_tagEle_EtaVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Eta(),TAPElectronTracks_v1[track_idx].Phi(),wt);
	// h_tagEle_PtVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Phi(),wt);
	// h_tagEle_PtVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Eta(),wt);
	//h_ZpT[2]->Fill(zvec.Pt(),wt);
	// h2_minDr_tagEle_pt[2]->Fill();
	// h_minDr_bestphovsTagEle[2]->Fill();
	if(BTags==0)
	  FR_nbtagBins[2]->Fill(1,wt);
	else if(BTags>=1)
	  FR_nbtagBins[2]->Fill(2,wt);
	// h_invariantMass[2]->Fill(invariantmass,wt);
      
      }
    }


         ////////////////////////////////////////////////////////////////////////////////////
    //  ====  ========== ======== 1 photon SR =====================  ///
    /////////////////////////////////////////////////////////////////////////////////// 
    bool fakeElectron=false, proc=false, proc1=false;
    if(sr){//hasPho && bestEMobj && !bestEMObjIsEle ){ //&& !s_sample.Contains("data")){
      // if(isoMuonTracks !=0 || isoPionTracks!=0 ) continue; // veto muon/pions/electrons 
      //if(NMuons!=0 ) continue;
      // h_selectBaselineYields_v1->Fill("Pho SR ",wt);
      // if(!s_sample.Contains("data")){
      // 	if(GenMuons_v1.size()==0 && GenElectrons_v1.size()==0 && GenTaus_v1.size()==0) continue; //reject w->qq events or all hadronic decay of TTbar                                                    
      // 	if(GenMuons_v1.size()>0) continue;
      // 	if (GenElectrons_v1.size()==0 ) continue;
	
      // //  // if(nGenEle1!=2) continue;
      
      // // 	//	if(((s_sample.Contains("WGJets") || s_sample.Contains("WJets"))&& nGentau_had1>0) ||( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")) && nGentau_had1>1)) continue;
      // 	if(v_genz.size()==0) continue;      
      // 	proc = v_genz.size()==1;
      // 	proc1 = v_genz.size()>1;
      // 	fakeElectron =(proc && (bestPhoton.DeltaR(v_genz[0])<0.2 && bestPhoton.Pt()/v_genz[0].Pt()<1.2 && bestPhoton.Pt()/v_genz[0].Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genz[0])<0.2 && bestPhoton.Pt()/v_genz[0].Pt()<1.2 && bestPhoton.Pt()/v_genz[0].Pt()>0.8)|| (v_genz.size()>1 && (bestPhoton.DeltaR(v_genz[1])<0.2 && bestPhoton.Pt()/v_genz[1].Pt()<1.2 && bestPhoton.Pt()/v_genz[1].Pt()>0.8))));
      // 	fakeElectron= true;
      // 	if(fakeElectron) pho_SR= true;
      // 	else continue;
      // 	h_selectBaselineYields_v1->Fill("Pho SR fake photon condition ",wt);

      // }
      pho_SR= true;
      //if(bestEMObj.Pt()<=40) continue;
      //cout<<" Debugging SR "<<"\t"<<hasPho<<"\t"<<"\t"<<bestEMObj.Pt()<<"\t"<<MET<<"\t"<<nHadJets<<"\t"<<ST<<"\t"<<bestPhoton.Pt()<<"\t"<<hasPho_px<<"\t"<<isoPionTracks<<"\t"<<isoMuonTracks<<"\t"<<NMuons<<"\t"<<NElectrons<<"\t"<<isoElectronTracks<<"\t"<<mt_ele<<"\t"<<endl;

      //if(ntracks>1) continue; //check this condition
      h_selectBaselineYields_v1->Fill("Pho SR nTrack conditions ",wt);

      if(pho_SR) {
	//h_invariantMass_noCut[3]->Fill(invariantmass,wt);
        // if(!(invariantmass>=80 && invariantmass<=100) )continue;

	// h_selectBaselineYields_v1->Fill("Pho SR inv mass cut",wt);

	FillHistogram_Kinematics(3,metstar,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestPhoton.Eta(),bestPhoton.Phi(),bestPhoton.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	FillHistogram_Kinematics_varBin(3,metstar.Pt(),nHadJets, BTags, bestPhoton.Pt(),ST,qmulti,wt);
	h_ZpT[3]->Fill(genzvec.Pt(), wt);
	// h_tagEle_pT_Pho_SR->Fill(tagEMObj.Pt(),wt);
        // h_tagEle_Eta_Pho_SR->Fill(tagEMObj.Eta(),wt);
        // h_tagEle_Phi_Pho_SR->Fill(tagEMObj.Phi(),wt);
        // h_trackEle_pT_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),wt);
	// h_trackEle_Eta_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Eta(),wt);
        // h_trackEle_Phi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Phi(),wt);
        // h_tagEle_EtaVsPhi_Pho_SR->Fill(tagEMObj.Eta(),tagEMObj.Phi(),wt);
	// h_tagEle_PtVsPhi_Pho_SR->Fill(tagEMObj.Pt(),tagEMObj.Phi(),wt);
        // h_tagEle_PtVsPhi_Pho_SR->Fill(tagEMObj.Pt(),tagEMObj.Eta(),wt);

        // h_tagEle_EtaVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Eta(),TAPElectronTracks_v1[track_idx].Phi(),wt);
        // h_tagEle_PtVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Phi(),wt);
        // h_tagEle_PtVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Eta(),wt);
	//h_ZpT[3]->Fill(zvec.Pt(),wt);
        if(BTags==0)
          FR_nbtagBins[3]->Fill(1,wt);
        else if(BTags>=1)
          FR_nbtagBins[3]->Fill(2,wt);
        //h_invariantMass[3]->Fill(invariantmass,wt);


      }
    }
			    
     
     
    nsurVived+=wt;
     }			     //loop over entries

     
  if(Debug)
    cout<<"filling the branches in tree"<<endl;
   //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
  cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;
  cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<"\t"<< count_nelec<<endl;
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
int AnalyzeLightBSM::getBinNo_v0FR(double pho_pt,double qmulti, double minDRindx){
  int sBin=0,m_i=1,sBin1=0,n_i=0; 
  for(int i=0;i<BestPhotonPtBinLowEdge.size()-1;i++){
    if(BestPhotonPtBinLowEdge[i]<99.99) continue; 
    //if(qmulti>=10) continue;
    if(i!=0)
      m_i++;  
    if(pho_pt >= BestPhotonPtBinLowEdge[i] && pho_pt < BestPhotonPtBinLowEdge[i+1])
      {
	sBin = sBin+((m_i-1)*4);

	break;

      }
    else if(pho_pt >= BestPhotonPtBinLowEdge[BestPhotonPtBinLowEdge.size()-1])
      {
	sBin = 40;
	break;
      }
  }

  if(sBin%4==0)// || sBin==4 || sBin==8 || sBin==12 || sBin==16 || sBin==20 )
    {
      //cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<qmulti<<" , sBin = "<<sBin1<<endl;

      for(int i=0;i<QMultLowedge.size()-1;i++){
	n_i++;
	if(qmulti>=QMultLowedge[i] && qmulti<QMultLowedge[i+1]) {sBin1=sBin+n_i; break;}
	else if(qmulti>=QMultLowedge[QMultLowedge.size()-1]){sBin1=sBin1+(QMultLowedge.size()-1); break;}
	//else sBin1=-999;
      }
    }
  
  //   if((qmulti>=7 && qmulti<10) && pho_pt>=100 && pho_pt <120)
  // if(pho_pt>=200 && pho_pt <230)
  //   cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<qmulti<<" , sBin = "<<sBin1<<endl;
  return sBin1;  
}
int AnalyzeLightBSM::getBinNo_v1FR(double pho_pt, int nHadJets){
  int sBin=0,m_i=1,sBin1=0,n_i=0; 
  for(int i=0;i<BestPhotonPtBinLowEdge.size()-1;i++){
    if(BestPhotonPtBinLowEdge[i]<99.99) continue; 
    if (i!=0)
      m_i++;  
    if(pho_pt >= BestPhotonPtBinLowEdge[i] && pho_pt < BestPhotonPtBinLowEdge[i+1])
      {
	sBin = sBin+((m_i-1)*3);

	break;

      }
    else if(pho_pt >= BestPhotonPtBinLowEdge[BestPhotonPtBinLowEdge.size()-1])
      {
	sBin = 30;
	break;
      }
  }

  if(sBin%3==0)// || sBin==4 || sBin==8 || sBin==12 || sBin==16 || sBin==20 )
    {
      for(int i=0;i<nJetsLowedge.size()-1;i++){
	n_i++;
	if(nHadJets>=nJetsLowedge[i] && nHadJets<nJetsLowedge[i+1]) {sBin1=sBin+n_i; break;}
	else if(nHadJets>=nJetsLowedge[nJetsLowedge.size()-1]){sBin1=sBin+(nJetsLowedge.size()-1); break;}
	//else sBin1=-999;
      }
    }
  
  //   if((nHadJets>=7 && nHadJets<10) && pho_pt>=100 && pho_pt <120)
  // if(pho_pt>=200 && pho_pt <230)
  //   cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<nHadJets<<" , sBin = "<<sBin1<<endl;
  return sBin1;  
}

int AnalyzeLightBSM::getBinNo_v2FR(double qmulti, int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets>=2 && nHadJets<=4)      { sBin=1;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=5;}
    else if(nHadJets>=7)   { sBin=9;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      { sBin=13;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=17;}
    else if(nHadJets>=7)   { sBin=21;}
  }
  
  for(int i=0;i<QMultLowedge.size()-1;i++){
    m_i++;
    if(qmulti>=QMultLowedge[i] && qmulti<QMultLowedge[i+1]) {sBin=sBin+m_i; break;}
    else if(qmulti>=QMultLowedge[QMultLowedge.size()-1]){sBin=sBin+(QMultLowedge.size()-1); break;}
    //else sBin1=-999;                                                                                                                                            
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nHadJets,int nbjets,double metstar)
{

  MET = metstar;
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
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;
  //  cout<<"\t"<<(*Photons_mvaValuesID).size()<<endl;
  for(int iPho=0;iPho<Photons_;iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    //    cout<<(*Photons_hasPixelSeed)[iPho] <<"\t (*Photons_hasPixelSeed)[iPho]  "<<" (*Photons_fullID)[iPho]"<<"\t"<<(*Photons_fullID)[iPho]<<"\t"<<Photons_v1[iPho].Eta()<<endl;
    if((*Photons_fullID)[iPho] && (*Photons_hasPixelSeed)[iPho]<0.001)//abs(Photons_v1[iPho].Eta())<2.4 && ((*Photons_fullID)[iPho] && pho_ID==0 && (*Photons_hasPixelSeed)[iPho]<0.001))// && ( ((*Photons_fullID)[iPho] && pho_ID==0) || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42)))
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

void AnalyzeLightBSM::FillHistogram_Kinematics(int i, TLorentzVector metstar, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST,float Eta , float Phi, double E, float METPhi,double qmulti,double leadJets_qmulti, double leadJet_Pt, double leadbjet_tag, double mindR, TLorentzVector Jet_matched, vector<TLorentzVector>hadJets, TLorentzVector leadJet,int nvrtx,double mindr_Pho_genElec, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;
  MET = metstar.Pt();
  h_Njets[i]->Fill(Njets,wt);
  h_Nbjets[i]->Fill(btags,wt);
  h_MET_[i]->Fill(MET,wt);
  h_PhotonPt[i]->Fill(pho_Pt,wt);
  h_Mt_PhoMET[i]->Fill(mt_phoMET,wt);
  h_dPhi_PhoMET[i]->Fill(dPhi,wt);
  h_St[i]->Fill(ST,wt);
  h_HT[i]->Fill(HT,wt);
  h_Photon_Eta[i]->Fill(Eta,wt);
  h_Photon_Phi[i]->Fill(Phi,wt);
  
  h_Photon_E[i]->Fill(E,wt);
  h_MET_Phi[i]->Fill(METPhi,wt);
  h_qmulti_1[i]->Fill(qmulti,wt);
  h_qmultiVsEmobjPT[i]->Fill(pho_Pt,qmulti,wt);
  h_qmultiVsnJets[i]->Fill(Njets,qmulti,wt);
  h_ST_vs_EMObjPt[i]->Fill( pho_Pt,ST,wt);
  h_Emobj_PtvsEta[i]->Fill(pho_Pt,Eta,wt);
  h_Emobj_PtvsPhi[i]->Fill(pho_Pt,Phi,wt);
  h_Emobj_EtavsPhi[i]->Fill(Eta,Phi,wt);
  //  cout<<"just after eta vs phi  "<<endl;
  h_nvrtx[i]->Fill(nvrtx,wt);
  h_minDR_Jets_EMObject[i]->Fill(mindR,wt);
  h_Phi_matchedJet[i]->Fill(Jet_matched.Phi(),wt);
  h_Pt_matchedJet[i]->Fill(Jet_matched.Pt(),wt);
  h_Eta_matchedJet[i]->Fill(Jet_matched.Eta(),wt);
  h_EtavsPhi_matchedJet[i]->Fill(Jet_matched.Eta(),Jet_matched.Phi(),wt);
  h_PtvsEta_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Eta(),wt);
  h_PtvsPhi_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Phi(),wt);
  h_minDR_Jets_vs_Em_Pt[i]->Fill(Jet_matched.Pt(),mindR,wt);
  h_btaggervalue_vs_qmulti[i]->Fill(qmulti,leadbjet_tag,wt);
  h_btaggervalue_vs_minDR_Jets[i]->Fill(mindR,leadbjet_tag,wt);
  h_minDR_Jets_vsqMulti[i]->Fill(qmulti,mindR,wt);
  h_HT5HT[i]->Fill(HT5/HT,wt);
  h_Emobje_pt_vs_Jet_Pt[i]->Fill(pho_Pt,Jet_matched.Pt(),wt);
  for(int j=0; j<hadJets.size();j++){
    if (j>=4) continue;
    h_Phi_leadJet[j][i]->Fill(hadJets[j].Phi(),wt);
    h_Pt_leadJet[j][i]->Fill(hadJets[j].Pt(),wt);
    h_Eta_leadJet[j][i]->Fill(hadJets[j].Eta(),wt);
    h_EtavsPhi_leadJet[j][i]->Fill(hadJets[j].Eta(),hadJets[j].Phi(),wt);
    h_PtvsPhi_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Phi(),wt);
    h_PtvsEta_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Eta(),wt);
    // TLorentzVector Met;
    // Met.SetPtEtaPhiE(MET,0,METPhi,0);
    double dPhi_METjet = abs(metstar.DeltaPhi(hadJets[j])); 
    h_HT5HT_vsdPhi_METJet[j][i]->Fill(dPhi_METjet,HT5/HT,wt);
    h_dPhi_METJet[j][i]->Fill(dPhi_METjet,wt);
    h_dPhi_METJet_vsMET[j][i]->Fill(MET,dPhi_METjet,wt);
  }

    
  h_nBjets_vs_qmulti[i]->Fill(btags,qmulti,wt);
  h_qmultiVs_MET[i]->Fill(MET,qmulti,wt);
  h_qmultiVs_ST[i]->Fill(ST,qmulti,wt);
  h_leadJets_qmulti[i]->Fill(leadJets_qmulti,wt);
  h_leadJet_Pt[i]->Fill(leadJet_Pt,wt);
  h_leadbjet_tag[i]->Fill(leadbjet_tag,wt);
  h_leadbjet_tag_vs_leadQmulti[i]->Fill(leadbjet_tag,leadJets_qmulti,wt);
  h_leadbjet_tag_vsQmulti[i]->Fill(leadbjet_tag,qmulti,wt);
  h_leadjet_ptvsqmulti[i]->Fill(leadJet_Pt,qmulti,wt);

  int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags,MET);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  int TFbins_v5 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
  int TFbins_v6 = getBinNo_v1FR(pho_Pt,Njets);
  int TFbins_v7 = getBinNo_v2FR(qmulti,Njets,btags);

  h_TFbins_LL_v4[i]->Fill(TFbins_v4,wt);
  h_TFbins_LL_v2[i]->Fill(TFbins_v2,wt);
  
  h_Sbins_LL[i]->Fill(searchBin,wt);
  h_TFbins_LL_v3[i]->Fill(TFbins_v3,wt);
  h_TFbins_LL_v5[i]->Fill(TFbins_v5,wt);
  h_TFbins_LL_v6[i]->Fill(TFbins_v6,wt);
  h_TFbins_LL_v7[i]->Fill(TFbins_v7,wt);
  

  
  //cout<<"Alps "<<i<<"\t"<<"TFbins_v1  "<<TFbins_v1<<" TFbins_v2 "<< TFbins_v2<<"\t MET "<<MET<<"\t NJets  "<<Njets<<"  btags "<<btags<<"\t"<<h_TFbins_LL_v1[i]->GetBinContent(TFbins_v1)<<"\t"<<wt<<"\t"<<searchBin<<endl;
  //cout<<h_Njets[i]->GetMean()<<endl;

}
void AnalyzeLightBSM::FillTFBins_Valid(int i, TLorentzVector metstar, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST,float Eta , float Phi, double E, float METPhi,double qmulti,double leadJets_qmulti, double leadJet_Pt, double leadbjet_tag, double mindR, TLorentzVector Jet_matched, vector<TLorentzVector>hadJets, TLorentzVector leadJet,int nvrtx,double mindr_Pho_genElec, double wt,double wt1, double wt2){
  //cout<<"enetring the validation loopjust after EM phi  "<<endl;
  MET = metstar.Pt();
  h_Mt_PhoMET_validation[i]->Fill(mt_phoMET,wt);
  h_dPhi_PhoMET_validation[i]->Fill(dPhi,wt);
  h_St_validation[i]->Fill(ST,wt);
  h_Photon_Eta_validation[i]->Fill(Eta,wt);
  h_Photon_Phi_validation[i]->Fill(Phi,wt);
  //  cout<<"just after EM phi  "<<endl;

  h_Photon_E_validation[i]->Fill(E,wt);
  h_MET_Phi_validation[i]->Fill(METPhi,wt);
  h_qmulti_1_validation[i]->Fill(qmulti,wt);
  h_nvrtx_validation[i]->Fill(nvrtx,wt);
  //cout<<"just after n vertex "<<endl;

  h_minDR_Jets_EMObject_validation[i]->Fill(mindR,wt);
  h_Phi_matchedJet_validation[i]->Fill(Jet_matched.Phi(),wt);
  h_Pt_matchedJet_validation[i]->Fill(Jet_matched.Pt(),wt);
  h_Eta_matchedJet_validation[i]->Fill(Jet_matched.Eta(),wt);
  //  cout<<"just after matched jet  "<<endl;

  h_HT5HT_validation[i]->Fill(HT5/HT,wt);
  for(int j=0; j<hadJets.size();j++){
    if (j>=4) continue;
    h_Phi_leadJet_validation[j][i]->Fill(hadJets[j].Phi(),wt);
    h_Pt_leadJet_validation[j][i]->Fill(hadJets[j].Pt(),wt);
    h_Eta_leadJet_validation[j][i]->Fill(hadJets[j].Eta(),wt);
    // TLorentzVector Met;
    // Met.SetPtEtaPhiE(MET,0,METPhi,0);
    double dPhi_METjet = abs(metstar.DeltaPhi(hadJets[j]));
    h_dPhi_METJet_validation[j][i]->Fill(dPhi_METjet,wt);
  }
  //cout<<"just after jet loop  "<<endl;

  h_leadJets_qmulti_validation[i]->Fill(leadJets_qmulti,wt);
  h_leadJet_Pt_validation[i]->Fill(leadJet_Pt,wt);
  h_leadbjet_tag_validation[i]->Fill(leadbjet_tag,wt);
  //cout<<"just after lead jet  "<<endl;

  //old ones
  int TFbins_v2 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
  //int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags,MET);
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
void AnalyzeLightBSM::FillHistogram_Kinematics_varBin(int i, double metstar,int Njets, int btags, double pho_Pt, double ST,double qmulti, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;                                        
  MET = metstar;
  h_Njets_Varbin[i]->Fill(Njets,wt);
  h_Nbjets_Varbin[i]->Fill(btags,wt);
  h_MET_Varbin[i]->Fill(MET,wt);
  h_PhotonPt_Varbin[i]->Fill(pho_Pt,wt);
  h_St_Varbin[i]->Fill(ST,wt);
  h_qmulti_Varbin[i]->Fill(qmulti,wt);
    h_qmultiVs_phopT_Varbin[i]->Fill(pho_Pt,qmulti,wt);
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

TLorentzVector AnalyzeLightBSM::getPhoton_withoutfullID(){
  vector<TLorentzVector> nogoodPho;
  vector<int> nogoodPhoIndx;
  for(int iPho=0;iPho<Photons_;iPho++){
    if( !(*Photons_fullID)[iPho] ) 
      {
	if(abs(Photons_v1[iPho].Eta())<2.4){
	nogoodPho.push_back( Photons_v1[iPho] );
	nogoodPhoIndx.push_back(iPho);
	}
      }
  }

  int highPtIndx=-100;
  for(int i=0;i<nogoodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (nogoodPho[highPtIndx].Pt()) < (nogoodPho[i].Pt()) ){highPtIndx=i;}
  }
  int eIndxAmongPhotons=-100;
  if(highPtIndx>=0){
    eIndxAmongPhotons = nogoodPhoIndx[highPtIndx];
  }
  else eIndxAmongPhotons = -100;

  if(highPtIndx==-100){TLorentzVector v0;return v0;}
  else return nogoodPho[highPtIndx];
}


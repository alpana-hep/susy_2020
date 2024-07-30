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
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
using namespace TMVA;
int main(int argc, char* argv[])
{

  if (argc < 4) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<<  " "<<"Which pho_ID"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  //  const char *elec = argv[5];
  //  const char *phoID = argv[5];
  //TString pho_ID = phoID;

  AnalyzeLightBSM ana(inputFileList, outFileName, data,sample);
  cout << "dataset " << data << " " << endl;
  // cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  //  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  cout<<inputFileList <<"\t"<<outFileName<<"\t"<<data<<"\t"<<sample<<endl;
  ana.EventLoop(data,inputFileList,sample,outFileName);
  Tools::Instance();
  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName) {
  cout<<"inside event loop"<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
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
  bool higMET=false, highdphi=false;
  double deepCSVvalue=0,p0=0,p1=0,p2=0;
  bool applyTrgEff=false;
  if(s_data.Contains("2016preVFP")){ lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV
  if(s_data.Contains("2016postVFP")) { lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016

  if(s_data.Contains("2017")) {lumiInfb=41.48; deepCSVvalue = 0.4506;}
  if(s_data.Contains("2018")){ lumiInfb=59.83;deepCSVvalue = 0.4168;}
  if(s_data.Contains("signal"))lumiInfb= 137.19; 
  if(s_sample.Contains("signal") && s_data.Contains("Run2")) {deepCSVvalue = 0.4168; lumiInfb= 137.19;} //taking 2018 deepCSV medium working point for signals because it has more efficiency as compared to other eras.
  
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
	if(Debug)
	  cout<<"alpa"<<"\t"<<jentry<<"\t"<<endl;
	// do something with tree branches here

      }
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
  return sBin;
}

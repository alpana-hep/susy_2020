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

  // TFile *outfile;                                            
  // outfile = TFile::Open(outFileName, "RECREATE");  
  // if (fChain == 0) return;
  // Long64_t nentries = fChain->GetEntriesFast();
  // fChain->SetBranchStatus("*",1);
  // TTree* outtree1 = fChain->CopyTree("MET>100" );
  // unsigned int output_entries1 = outtree1->GetEntries();
  // cout << "Output tree has entries: " << output_entries1 << endl
  //      << "Reduction factor of: " << double(nentries)/double(output_entries1) << endl;

  
 //  //  Long64_t nentries = fChain->GetEntriesFast();
 //  cout<<"original entries  "<<nentries<<endl;
 //  cout << "nentries " << nentries << endl;
 //  cout << "Analyzing dataset " << data << " " << endl;
 //  TString Elec_flag = elec;
 //  TString pho_ID_str = phoID;
  
 //  int pho_ID=-1;
 //  if(pho_ID_str.Contains("loose"))
 //    pho_ID=0;
 //  else if (pho_ID_str.Contains("medium"))
 //    pho_ID=1;
 //  else if (pho_ID_str.Contains("tight"))    
 //    pho_ID=2;
 //  else if (pho_ID_str.Contains("mva_wp90"))
 //    pho_ID=3;

 //  else if (pho_ID_str.Contains("mva_wp80"))
 //    pho_ID=4;
 //  cout<<"which photon ID "<<"\t"<<"case"<<"\t"<<pho_ID<<endl;
 //  // TString pho_D = phoID;
 //  bool lost_elec_flag=false;
 //  if(Elec_flag.Contains("Electron"))
 //    lost_elec_flag = true;
 //  else
 //    lost_elec_flag = false;
  
 //  float counter=0.0;
 //  TString s_data=data;
 // char* s_cross = new char[100];
 //  sprintf(s_cross,"%s.%s",data,sample);
 //  TString s_sample= sample;
  

 //  cout<<"Calculating the weight for process "<<s_cross<<endl;
 //  std::string s_process = s_cross;
 //  TString s_Process = s_process;
 //  double cross_section = getCrossSection(s_process);
 //  cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_cross<<endl;

  //Skimming

  const char* sel = "(MET>100)";

  cout << "Skimming with selection : "<<sel<<endl;

  //--------------------------------------------------                                                                                                                                                             
  // input and output file                                                                                                                                                                                         
  //--------------------------------------------------                                                                                                                                                             

  //const char* infilename  = Form("WGJets_MonoPhoton_PtG-40to130_RA2AnalysisTree.root");                                                                                                                          
  ifstream infile(inputFileList, ifstream::in);
  std::string infilename;
  const char* outfilename = Form(outFileName);//"output1.root");
  
  TFile *outfile;
  bool fileopen = false;
  //--------------------------------------------------                                                                                                                                                             
  // cout stuff                                                                                                                                                                                                    
  //--------------------------------------------------                                                                                                                                                             
  TChain *chain = new TChain("TreeMaker2/PreSelection");
  TTree*outtree;
TTree* FinalTree;
 int loop_count=0;
  while(1) {
    infile >> infilename;
    if(!infile.good()) break;
    cout << "Reading in : " << infilename << endl;
    cout << "Writing to : " << outfilename << endl;
    cout << "Selection : " << sel << endl;

  //--------------------------------------------------                                                                                                                                                           
  // read input file, write to output files                                                                                                                                                                      
  //--------------------------------------------------                                                                                                                                                           

  //long long max_tree_size = 5000000000LL; // 5GB                                                                                                                                                               
  long long max_tree_size = 3000000000LL; // 3.0GB                                                                                                                                                               
  TTree::SetMaxTreeSize(max_tree_size);
  loop_count++;
  chain->Add(infilename.c_str());

  unsigned int input_entries = chain->GetEntries();
  cout << "Input tree has entries: " << input_entries << endl;

  //-------------------                                                                                                                                                                                          
  // skim                                                                                                                                                                                                        
  //-------------------                                                                                                                                                                                         
  TFile *outfile1;
  if(!fileopen) {
    outfile = TFile::Open(outfilename, "RECREATE");
  }
  else {
    outfile1 = TFile::Open(outfilename, "UPDATE");
    assert( outfile1 != 0 );
  }
  //outfile->cd();                                                                                                                                                                                               
  outtree = chain->CopyTree( sel );
  unsigned int output_entries = outtree->GetEntries();
  cout << "File - "<<loop_count<< "\t"<<"Output tree has entries: " << output_entries << endl
       << "Reduction factor of: " << double(input_entries)/double(output_entries) << endl;
  if(!fileopen) {
    outtree->Write();
    outfile->Close();
    fileopen = true;
  }
  else {
  if(loop_count==10){
    outtree->Write();
    outfile->Close();
    }
  }
  //FinalTree = outtree->Clone(0);
  //cout<< "FInal tree entries  "<<FinalTree->GetEntries()<<endl;

  }
  //unsigned int  output_entries = outtree->GetEntries();
  
  //cout << "Output tree has entries: " << output_entries << endl;



  //  newtree->CopyEntries(oldtree);
  //  TTree* outtree = chain->CopyTree( sel );
  // TTree* originalTree;
  // originalTree->CopyEntries(fChain);//->CloneTree();
  //originalTree->CopyEntries(fChain->GetTree());
  //cout<<"originalTree "<< originalTree->GetEntries()<<endl;
  //  TTree* selectedTree = (TTree*)originalTree->CopyTree("MET>100");
  // TTree* originalTree = (TTree*)fChain->GetTree();
  // TTree* selectedTree = originalTree->CopyTree("MET>100");

  //selectedTree->Write();
  //  if(Debug)
  //   cout<<"filling the branches in tree"<<endl;
  //  //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
  // cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;

  //  cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<endl;
  // cout << "============   Electron     ===============" <<endl;

  // cout<<"applied weights "<<" "<<wt<<endl;
  // cout<<"CR electron :==  "<<"\t"<<nCR_elec<<"\t"<<endl; 
  // cout<<"SR electron :==  "<<"\t"<<nSR_elec<<"\t"<<endl;
  // cout<<"SR e- : Fail acceptance "<<"\t"<<FailAccept_Elec<<"\t"<<endl;
  // cout<<"SR e- : Fail Id " <<"\t"<<FailId_Elec<<"\t"<<endl;
  // cout<< "SR e- : Fail Iso" <<"\t"<<FailIso_Elec<<"\t"<<endl;

  // cout << "============   Muon     ===============" <<endl;

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
				// sspq


// int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nbjets, int nHadJets)
// {
  
//   int sBin=-100,m_i=0;
//   if(nbjets==0 ){
//     if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
//     else if(nHadJets==5 || nHadJets==6){ sBin=7;}
//     else if(nHadJets>=7)               { sBin=13;}
//   }
//   else{
//     if(nHadJets>=2 && nHadJets<=4)     { sBin=19;}
//     else if(nHadJets==5 || nHadJets==6){ sBin=25;}
//     else if(nHadJets>=7)               { sBin=31;}
//   }
//   if(sBin==0){
//     for(int i=0;i<METLowEdge_v3.size()-1;i++){
//       if(METLowEdge_v3[i]<99.99) continue;
//       int sBin1=sBin;
//       m_i++;
//       if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
//         break; }
//       else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
//         break; }
//     }
//   }
//   else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
//     int sBin1=sBin;
//     for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
//       if(METLowEdge_v3_1[i]<99.99) continue;
//       m_i++;
//       if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
//         break;}
//       else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
//         break; }
//     }
//   }

//   // else if(sBin==37){
//   //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
//   //     if(METLowEdge_v3[i]<99.99) continue;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
//   //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }                                                                                                                    

//   //   }
//   // }

//   // else if(sBin==44){
//   //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
//   //     if(METLowEdge_v3[i]<99.99) continue;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
//   //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }                                                                                                                    

// //     }
// //   }
// // -
//   // int sBin=-100,m_i=0;
//   // if(nbjets==0 ){
//   //   if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
//   //   else if(nHadJets==5 || nHadJets==6){ sBin=8;}
//   //   else if(nHadJets>=7)               { sBin=15;}
//   // }
//   // else{
//   //   if(nHadJets>=2 && nHadJets<=4)     { sBin=22;}
//   //   else if(nHadJets==5 || nHadJets==6){ sBin=29;}
//   //   else if(nHadJets>=7)               { sBin=36;}
//   // }
//   // if(sBin==0){
//   //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
//   //     if(METLowEdge_v3[i]<99.99) continue;
//   //     int sBin1=sBin;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
//   // 	break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
//   // 	break; }
//   //   }
//   // }
//   // else if(sBin==8 || sBin==15 || sBin==22 || sBin==29 || sBin==36){
//   //   int sBin1=sBin;
//   //   for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
//   //     if(METLowEdge_v3_1[i]<99.99) continue;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
//   // 	break;}
//   //     else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
//   // 	break; }
//   //   }
//   // }

//   // else if(sBin==37){
//   //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
//   //     if(METLowEdge_v3[i]<99.99) continue;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }

//   // }
//   // }

//   // else if(sBin==44){
//   //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
//   //     if(METLowEdge_v3[i]<99.99) continue;
//   //     m_i++;
//   //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
//   //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
//   //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }

//   //   }
//   //}
//   return sBin;
// }
// int AnalyzeLightBSM::getBinNoV7_le(int nbjets, int nHadJets){
//   int sBin=-100,m_i=0;
//   if(nbjets==0){
//     if(nHadJets==2) { if(MET<150)sBin=1; else if (MET>150) sBin=2;}
//     else if(nHadJets==3)     { if(MET<150)sBin=3; else if (MET>150) sBin=4;}
//     else if(nHadJets==4)     { if(MET<150)sBin=5; else if (MET>150) sBin=6;}
//     else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=7; else if (MET>150) sBin=8;}
//     else if(nHadJets>=7)   { if(MET<150)sBin=9; else if (MET>150) sBin=10;}
//   }
//   else{
//     if(nHadJets>=2 && nHadJets<=4)      {if(MET<150)sBin=11; else if (MET>150) sBin=12;}
//     else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=13; else if (MET>150) sBin=14;}
//     else if(nHadJets>=7)   { if(MET<150)sBin=15; else if (MET>150) sBin=16;}
//   }
//   return sBin;
// }

// int AnalyzeLightBSM::getBinNoV1_le(int nbjets, int nHadJets){
//   int sBin=-100,m_i=0;
//   if(nbjets==0){
//     if(nHadJets==2)     { sBin=1;}
//     else if(nHadJets==3)     { sBin=2;}
//     else if(nHadJets==4)     { sBin=3;}
//     else if((nHadJets==5 || nHadJets==6)){ sBin=4;}
//     else if(nHadJets>=7)   { sBin=5;}
//   }
//   else{
//     if(nHadJets>=2 && nHadJets<=4)      { sBin=6;}
//     else if((nHadJets==5 || nHadJets==6)){ sBin=7;}
//     else if(nHadJets>=7)   { sBin=8;}
//   }
//   return sBin;
// }
// int AnalyzeLightBSM::getBinNoV7(int nHadJets, int btags){
//   int sBin=-100,m_i=0;
//   if(btags==0){
//     if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
//     else if(nHadJets==5 || nHadJets==6){ sBin=6;}
//     else if(nHadJets>=7)               { sBin=11;}
//   }
//   else{
//     if(nHadJets>=2 && nHadJets<=4)     { sBin=16;}
//     else if(nHadJets==5 || nHadJets==6){ sBin=21;}
//     else if(nHadJets>=7)               { sBin=26;}
//   }
//   if(sBin==0){
//     for(int i=0;i<METLowEdge1.size()-1;i++){
//       if(METLowEdge1[i]<99.99) continue;
//       m_i++;
//       if(MET >= METLowEdge1[i] && MET < METLowEdge1[i+1]){ sBin = sBin+m_i;break; }
//       else if(MET >= METLowEdge1[METLowEdge1.size()-1])  { sBin = 6         ;break; }
//     }
//   }
//   else{
//     for(int i=0;i<METLowEdge2.size()-1;i++){
//       if(METLowEdge2[i]<99.99) continue;
//       m_i++;
//       if(MET >= METLowEdge2[i] && MET < METLowEdge2[i+1]){ sBin = sBin+m_i;break; }
//       else if(MET >= METLowEdge2[METLowEdge2.size()-1])  { sBin = sBin+5   ;break; }
//     }
//   }
//   return sBin;
// }
// TLorentzVector AnalyzeLightBSM::getBestPhoton(int pho_ID){
//   // TLorentzVector v1;
//   // vector<TLorentzVector> goodPho;
//   // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
//   //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
//   // }

//   // if(goodPho.size()==0) return v1;
//   // sortTLorVec(&goodPho);
//   // return goodPho[0];
//   vector<TLorentzVector> goodPho;
//   vector<int> goodPhoIndx;
//   for(int iPho=0;iPho<Photons_;iPho++){
//     //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
//     if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001) &&( pho_ID==0 || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42))) ) 
//       {
//       goodPho.push_back(Photons_v1[iPho] );
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

// // TLorentzVector AnalyzeLightBSM::getBestPhoton_tightID(){
// //   vector<TLorentzVector> goodPho;
// //   vector<int> goodPhoIndx;
// //   for(int iPho=0;iPho<Photons->size();iPho++){
// //     if((*Photons_hadTowOverEM)[iPho]<0.02148 && (*Photons_sigmaIetaIeta)[iPho]<0.00996 &&  (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
// //       goodPho.push_back( (*Photons)[iPho] );
// //       goodPhoIndx.push_back(iPho);
// //     }
// //   }

// //   int highPtIndx=-100;
// //   for(int i=0;i<goodPho.size();i++){
// //     if(i==0) highPtIndx=0;
// //     else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
// //   }

// //   if(highPtIndx>=0){
// //     bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
// //   }
// //   else bestPhotonIndxAmongPhotons = -100;
// //   if(highPtIndx==-100){TLorentzVector v0;return v0;}
// //   else return goodPho[highPtIndx];
// // }



// double AnalyzeLightBSM::getGendRLepPho(int pho_ID){//MC only
//   TLorentzVector genPho1,genLep1;
//   int leadGenPhoIdx=-100;
//   // vector<TLorentzVector> goodPho;
//   // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
//   //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back((*Photons)[iPhoton]);
//   // }
//   // if(goodPho.size()!=0) 
//   genPho1 =getBestPhoton(pho_ID);
  
//   for(int i=0;i<GenParticles_;i++){
//      if(GenParticles_v1[i].Pt()!=0){
//        if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
// 	 if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
//        }
//      }
//   }//for
//   if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
//   else return 1000.0;
// }

//  double AnalyzeLightBSM::getdR_GenPho_RecoPho(TLorentzVector bestPhoton)
//  {
   
//    TLorentzVector genPho;
//    int leadGenPhoIdx=-100;
//    int minDR = 9999;
//    vector<TLorentzVector> v_genPho;
//    for (int igen=0; igen<GenParticles_; igen++)
//      {
//        if(GenParticles_v1[igen].Pt()!=0){
// 	 if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
// 	   genPho = (GenParticles_v1[igen]);
// 	   v_genPho.push_back(genPho);
// 	 }
//        }
//      }
//    return MinDr(bestPhoton,v_genPho);   
//  }
// double AnalyzeLightBSM::getGendRElecPho(int pho_ID){//MC only                                                                                             
//   TLorentzVector genPho1,genLep1;
//   int leadGenPhoIdx=-100;
//   genPho1 =getBestPhoton(pho_ID);
//   for(int i=0;i<GenParticles_;i++){
//     if(GenParticles_v1[i].Pt()!=0){
//       if( abs((*GenParticles_PdgId)[i])==11 && (abs((*GenParticles_ParentId)[i])<=25 )&& (abs((*GenParticles_ParentId)[i])!=15)){
// 	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
//       }
//     }
//   }//for                                                                                                                                       
//   if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
//   else return 1000.0;

// }
// double AnalyzeLightBSM::getGenLep1(TLorentzVector bestPhoton, bool flag){                                     
//   vector<TLorentzVector> v_genLep2;
//   TLorentzVector genMu1, genEle1;
//   if(flag)
//     {
//       for(int i=0 ; i < GenElectrons_; i++)
// 	{
// 	  if(GenElectrons_v1[i].Pt()!=0)
// 	    {
// 	      genEle1 = (GenElectrons_v1[i]);
// 	      v_genLep2.push_back(genEle1);
// 	    }

// 	}
//     }
//   else
//     {
//       for(int i=0 ; i < GenMuons_; i++)
// 	{
// 	  if(GenMuons_v1[i].Pt()!=0)
// 	    {
// 	      genMu1 = (GenMuons_v1[i]);
// 	      v_genLep2.push_back(genMu1);
// 	    }
// 	}
//     }
//   return MinDr(bestPhoton,v_genLep2);
// }


// double AnalyzeLightBSM::getGenLep(TLorentzVector bestPhoton){
//   vector<TLorentzVector> v_genLep2;
//   TLorentzVector genMu1, genEle1;
//   // if(flag)
//   //   {
//       for(int i=0 ; i < GenElectrons_; i++)
//         {
//           if(GenElectrons_v1[i].Pt()!=0)
//             {
//               genEle1 = (GenElectrons_v1[i]);
//               v_genLep2.push_back(genEle1);
//             }

//         }
//   //   }
//   // else
//   //   {
//       for(int i=0 ; i < GenMuons_; i++)
//         {
//           if(GenMuons_v1[i].Pt()!=0)
//             {
//               genMu1 = (GenMuons_v1[i]);
//               v_genLep2.push_back(genMu1);
//             }
//         }
//       //  }
//   return MinDr(bestPhoton,v_genLep2);
// }

// double AnalyzeLightBSM::getGenRecodRLep(TLorentzVector recoLep1){//MC only                                                                                                
//   TLorentzVector genPho1,genLep1;
//   int leadGenPhoIdx=-100;
//   for(int i=0;i<GenParticles_;i++){
//     if(GenParticles_v1[i].Pt()!=0){
//       if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
// 	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
//       }
//     }
//   }//for                                                                                                                                          
//   if(recoLep1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(recoLep1);
//   else return 1000.0;
// }

// std::vector<int> AnalyzeLightBSM::dR_recoPho_GenParticle(TLorentzVector recoPho){//MC only                                                                     
//   TLorentzVector genPho1,genLep1;
//   int leadGenPhoIdx=-100, pdgID =-99999, parentID=-99999, count=0;
//   int flag = 0;
//   std::vector<int> array;//={};
//   for(int i=0;i<GenParticles_;i++){
//     if(GenParticles_v1[i].Pt()!=0){
//       //if(abs((*GenParticles_PdgId)[i])==22 || abs((*GenParticles_ParentId)[i])==22) continue;
//       if(recoPho.DeltaR(GenParticles_v1[i])<0.2 )
// 	{
// 	  //cout<<"You got to be kidding me"<<endl;	  
// 	  //h_parID_matchRecoPho->Fill((*GenParticles_PdgId)[i],wt);
// 	  //h_parentID_vsPartId_matchPho->Fill((*GenParticles_PdgId)[i],(*GenParticles_ParentId)[i],wt);
// 	  if(abs((*GenParticles_PdgId)[i])==22)
// 	    {
// 	      //h_phopT_BeamRemenant->Fill(recoPho.Pt(),wt);
// 	      continue;
// 	    }
// 	  if(abs((*GenParticles_ParentId)[i])==22) continue;	  
// 	  pdgID =(*GenParticles_PdgId)[i];
// 	  count++;
// 	  flag = 1;	  
// 	  //cout<<"You got to be kidding me"<<"\t"<<flag<<"\t"<<pdgID<<"\t"<<(*GenParticles_ParentId)[i]<<"\t"<<(*GenParticles_Status)[i]<<endl;
// 	}
//     }
//   }
//   array.push_back(pdgID);
//   array.push_back(flag);
//   //cout<<"after kidding me"<<"\t"<<array[0]<<"\t"<<array[1]<<endl;
//   return array;//pdgID, flag;
 
// }
// // double AnalyzeLightBSM::getMVArespone(float MET, float NJets, float BTags, float HT)
// // {
// //   TMVA::Tools::Instance();

// //   TMVA::Reader *reader1 = new TMVA::Reader("Color:Silent");
// //   float met=0.0,njets=0.0,btags=0.0,ht=0.0;
// //   reader1->AddVariable( "MET", &met );
// //   reader1->AddVariable( "NJets", &njets );
// //   reader1->AddVariable( "BTags", &btags );
// //   reader1->AddVariable( "HT", &ht );
// //   reader1->BookMVA( "BDT_200trees_2maxdepth method", "./TMVAClassification_BDT_200trees_2maxdepth.weights.xml" );
// //   met = MET;
// //   njets= NJets;
// //   btags = BTags;
// //   ht = HT;
// //   Double_t mvaValue = reader1->EvaluateMVA( "BDT_200trees_2maxdepth method");
// //   return mvaValue;
// //   delete reader1;
// // }

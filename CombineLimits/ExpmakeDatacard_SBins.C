#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"

using namespace std;
string getfname(const char *fname1){string fname=fname1;fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();return fname;}
void setLastBinAsOverFlow(TH1D*);

void makeDatacard_SBins(double mGl,double mNLSP,TString sigFile,string histname1,string histname, string modelname){//[1000],char histname[1000]){
  char name[1000];
  int jmax=6;//no. of backgrounds
  int nSig=1;//1 signal
  int nFiles=nSig+jmax;
  TFile *f[nFiles];
  char* fnames = new char[1000];
  sprintf(fnames,"FullRun2_WGJets_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[1] = new TFile(fnames);
  sprintf(fnames,"FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[2] = new TFile(fnames);
  sprintf(fnames,"FullRun2_TTJets_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[3] = new TFile(fnames);
  sprintf(fnames,"FullRun2_TTGJets_inc_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[4] = new TFile(fnames);
  sprintf(fnames,"FullRun2_ZNuNu_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[5] = new TFile(fnames);
  sprintf(fnames,"FullRun2_GJets_QCD_PhoIdloose_phopt40_BL_BDTwith13variables_%s.root",modelname.c_str());
  f[6] = new TFile(fnames);



  f[0] = new TFile(sigFile);
  // f[1] = new TFile("FullRun2_WGJets_PhoIdloose_phopt40.root");
  // f[2] = new TFile("FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40.root");
  // f[3] = new TFile("FullRun2_TTJets_PhoIdloose_phopt40.root");
  // f[4] = new TFile("FullRun2_TTGJets_inc_PhoIdloose_phopt40.root");
  // f[5] = new TFile("FullRun2_ZNuNu_PhoIdloose_phopt40.root");
  //  f[6] = new TFile("FullRun2_GJets_QCD_PhoIdloose_phopt40.root");

  //  char histname[100]="MET_R1";
  //  char histname[100]="METvBin_EW";
  TH1D *hist1[nFiles];
  //  double min_binLowedge=99.9999;
  std::ofstream outf;
  vector<int> observation;
  for(int i=0;i<nFiles;i++){
    hist1[i] = (TH1D*)f[i]->FindObjectAny(histname.c_str()); //h_Sbins_LL_newSbins_v7_Pred_SR
    //setLastBinAsOverFlow(hist1[i]);
  }

  TH1D *hist=(TH1D*)f[0]->FindObjectAny(histname.c_str());
  int imax=hist->GetNbinsX();
  for(int i=1;i<=imax;i++){observation.push_back(1);}
  for(int i=1;i<=imax;i++){
    cout<<"MEt bin --> "<<i<<" bin content --> "<<hist->GetBinContent(i)<<endl;
  }
  // cout<<"-----"<<getfname(f[j]->GetName())<<"-----"<<endl;}
  cout<<"file name : "<<getfname(f[0]->GetName())<<endl;
  for(int i=1;i<=imax;i++){
    bool flag_write=true;
    int counter=0;
    for(int j=0;j<jmax+nSig;j++){
      if(hist1[j]->GetBinContent(i) <= 0) {   cout<<"file name : "<<getfname(f[j]->GetName())<<"\t"<<i<<endl;
	counter++;
      }      
    }
    //    if(i==13) cout<<"counter "<<counter<<"\t"<<flag_write<<endl;
    if (counter ==jmax+nSig) { flag_write=false; cout<<"total N files has 0 bin content "<<i<<"\t"<<counter<<endl;}
    // string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    if(flag_write){
    string name2="dataCards/"+getfname(f[0]->GetName())+"_"+histname1+"/"+getfname(f[0]->GetName())+"_"+histname+"_bin"+to_string(i)+".txt";
    //cout<<name2<<endl;
    sprintf(name,"%s",name2.c_str());
    //    cout<<name<<endl;
    outf.open(name,ios::out);
    outf<<"# - - - - - - - - - - - - - - - - - - -"<<endl<< 
      "# Datacard for mGl= "<<mGl<<" mNLSP= "<<mNLSP<<endl<<
      "# - - - - - - - - - - - - - - - - - - - "<<endl<<
      "imax 1 number of channels"<<endl<<
      "jmax *  number of backgrounds('*' = automatic)"<<endl<<
      "kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl<<
      "------------"<<endl<<
      "bin "<<histname<<i<<endl<<
      "observation "<<observation[i-1]<<endl<<
      "------------"<<endl<<
      "bin ";
    for(int j=0;j<jmax+nSig;j++){outf<<histname<<i<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<getfname(f[j]->GetName())<<" ";}
    outf<<endl<<
      "process ";
    for(int j=0;j<jmax+nSig;j++){outf<<j<<" ";}
    outf<<endl<<
      "rate ";
    for(int j=0;j<jmax+nSig;j++){
      if(hist1[j]->GetBinContent(i) >= 0) outf<<hist1[j]->GetBinContent(i)<<" ";
      else outf<<"0 ";
    }
    outf<<endl<<"------------"<<endl;

    for(int j1=0;j1<jmax+nSig;j1++){
      outf<<getfname(f[j1]->GetName())<<"bin"<<i<<" lnN ";
      for(int j2=0;j2<jmax+nSig;j2++){
	if(j1==j2){
	  //  outf<<"1.20 ";
	  if((hist1[j2]->GetBinContent(i))>0.00000001) outf<<1.2 + (hist1[j2]->GetBinError(i))/(hist1[j2]->GetBinContent(i))<<" ";
	  else outf<<"1.20 ";
	}
	else{
	  outf<<"- ";
	}
      }
      outf<<endl;
    }
  
    outf.close();
    //delete hist2;
    }
  }
}


void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);

}

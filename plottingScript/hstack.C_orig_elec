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
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"


const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[5] = {"Failed Id","Failed Iso","Failed acceptance","(0e,0#mu)"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[5]={"Lost Electron","Lost Mu","Taus","Unidentified-Photon","Unidentified"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
int line_color[9] = {9,7,45,kGray+2,kMagenta,kRed,kBlue+2,kMagenta,kCyan}; //{kBlack, kBlue, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
vector<int> col={9,7,45,kViolet,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
void setLastBinAsOverFlow(TH1F*);
std::map<int,int> rebin_;
std::map<int,float> x_max;
std::map<int,float> x_min;
void rebin_init_layer_pi(){
  rebin_.insert(pair<int, int>(20,2));
  rebin_.insert(pair<int, int>(50,2));
  rebin_.insert(pair<int, int>(100,2));
  rebin_.insert(pair<int, int>(200,4));
    
  x_max.insert(pair<int, int>(20,2500.0));
  x_max.insert(pair<int, int>(50,6000.0));
  x_max.insert(pair<int, int>(100,12000.0));
  x_max.insert(pair<int, int>(200,22000.0));
  
  x_min.insert(pair<int, int>(20,0.0));
  x_min.insert(pair<int, int>(50,0.0));
  x_min.insert(pair<int, int>(100,0.0));
  x_min.insert(pair<int, int>(200,0.0));
  x_min.insert(pair<int, int>(300,0.0));
  
}
//void decorate(THStack*,int);
void decorate(TH1F*,int );
// void decorate(THStack *hs,int i){
//   hs->GetXaxis()->SetLabelSize(.05);
//   hs->GetYaxis()->SetLabelSize(.05);
//   hs->GetXaxis()->SetTitleSize(0.05);
//   hs->GetYaxis()->SetTitleSize(0.05);
  
//   gStyle->SetOptStat(0);
// }
void decorate(TH1F* hist,int i){
  hist->SetLineColor(col[i]);
  //  if(i<nBG) {
    hist->SetFillColor(col[i]);
    // hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    //}
}
void setLastBinAsOverFlow(TH1F* h_hist){
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

void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", float energy=-1, char const *leg_head="",
		     bool flag=true, bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", int xrange=-1, char const *title1=""){  
  

  TCanvas *c = new TCanvas(tag_name, tag_name, 800, 700);
  c->SetRightMargin(0.009);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.025);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.06);
   c->SetBottomMargin(0.12);

   if(flag)
     TString legend_text[5]={"Lost Electron","Lost Mu","Taus","Unidentified-Photon","Unidentified"};

    gStyle->SetOptStat(1111111);
  //   gStyle->SetOptStat(0);
  //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  gStyle->SetOptStat(0);
  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","");
  cout<<xrange<<endl;
  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.54,0.63,0.9,0.9);  
  legend->SetTextSize(0.04);
  char* lhead = new char[100];
  sprintf(lhead,"#bf{%s}",title1);
  legend->SetHeader(lhead);

  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.025;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;
  float integral = hist.at(0)->Integral()+hist.at(1)->Integral()+hist.at(2)->Integral();
  //  cout<<"SR muontron : " <<"\t"<<hist.at(3)->Integral()<<"\t"<<integral<<endl;
  legend->AddEntry(hist.at(0),"Lost Electron","");
  //    leg_entry[(int)hist.size()]->SetTextColor(hist.at(i)->GetLineColor());

  for(int i = 0; i < (int)hist.size(); i++) {
    //    normalize = true;
    //    if(i>=3) continue;
    if(normalize) {
      	hist.at(i)->Scale(1.0/integral);
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
    
    hist.at(i)->SetLineWidth(line_width[i]);
    
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    
    hist.at(i)->GetXaxis()->SetTitleSize(0.04);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.030);
    // hist.at(i)->GetXaxis()->SetRangeUser(0,xrange);
    hist.at(i)->GetYaxis()->SetLabelSize(0.030);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.45);
    if(DoRebin) {
     hist.at(i)->Rebin(2);
      //hist.at(i)->Rebin(1);
    }
    hist.at(i)->GetXaxis()->SetRangeUser(0,xrange);

    /* hist.at(i)->GetXaxis()->SetRangeUser(x_min[energy],x_max[energy]); */
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //    setLastBinAsOverFlow(hist.at(i));
    
    

}
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-3;
  //  if(ymax<=10) ymax=10;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.1,ymax);
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.01,ymax);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
	// hs_var->Add(hist.at(i));
	// hs_var->SetMinimum(0.01);
	//,ymax*10.0);
	hist.at(i)->GetXaxis()->SetTitleOffset(1.0);
	hist.at(i)->GetXaxis()->SetTitle(title); hist.at(i)->GetYaxis()->SetTitle("Events");//hs_var->SetTitle(0);
	setLastBinAsOverFlow(hist.at(i));

	// hs_var->GetYaxis()->SetTitleOffset(.90);


      }
    if(i<3)
      hs_var->Add(hist.at(i));
    hs_var->SetMinimum(0.0001);
    hs_var->SetMaximum(ymax*100);
    //    hs_var->GetXaxis()->SetTitle(title);
    // 	if(!i) hist.at(i)->Draw("hist");
    // else   hist.at(i)->Draw("hist sames");
    
  }
  

  
hs_var->Draw("BAR HIST");
 hs_var->Draw("HIST");

 // hist.at(3)->Scale(1.0/hist.at(3)->Integral());
 // hist.at(3)->Draw(" sames");
 
 hs_var->GetXaxis()->SetRangeUser(0,xrange);
 hs_var->GetXaxis()->SetTitle(title);
hs_var->GetXaxis()->SetTitleOffset(1.0);
  gPad->Modified(); gPad->Update();
  hs_var->GetXaxis()->SetTitle(title); hs_var->GetYaxis()->SetTitle("Events");
  hs_var->SetTitle(0);
   hs_var->GetYaxis()->SetTitleOffset(1.2);
hs_var->GetXaxis()->SetTitleSize(00.05);
    hs_var->GetXaxis()->SetLabelSize(0.04);
    hs_var->GetYaxis()->SetLabelSize(0.04);
    hs_var->GetYaxis()->SetTitleSize(00.05);
    hs_var->GetYaxis()->SetTitleOffset(1.2);


    //hs_var->Draw("HIST");

    legend->Draw();



  if(log_flag) {
      gPad->SetLogy();
    }
  if(logx)
    gPad->SetLogx();

  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);
  textOnTop->DrawLatexNDC(0.12,0.96,"CMS #it{#bf{Preliminary}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.7,0.96,en_lat);

 

  gPad->Modified();
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"./Results/hstack_elec/%s_stacked.png",tag_name);
    c->SaveAs(canvas_name);
    
    sprintf(canvas_name,"./Results/hstack_elec/%s_stacked.pdf",tag_name);//./plots/%s.pdf",tag_name);
    c->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=10,nBG=6;                                                                                                                                                              
TFile *f[nfiles];


void hstack()
{
  char* hname = new char[200];
  char* hist_name  = new char[200];
  char* hist_name1 = new char[200];
  char* hist_name2 = new char[200];
  char* hist_name3 = new char[200];
  char* hist_name4 = new char[200];
  char* hist_name5 = new char[200];
  char* hist_name6 = new char[200];

  char* full_path = new char[2000];
  char* full_path1 = new char[2000];
  char* full_path2 = new char[2000];
  char* path2 = new char[2000];
  char* title= new char[2000];
  //string filetag;//=new char[20000];                                                                                                                                                                   
  char* full_path3 = new char[2000];
  char* full_path4 = new char[2000];
  char* full_path5 = new char[2000];
  char* full_path6 = new char[2000];
  char* full_path7 = new char[2000];
  char* full_path8 = new char[2000];
  char* full_path9 = new char[2000];
  char* full_path10 = new char[2000];
  char* full_path11= new char[2000];
  char *leg_head = new char[200];
  sprintf(path2,"./Results");
  // f[0] = new TFile("Out_Autumn18_WGJets_LL_estimation_v18_Jul22.root");//Autumn18_lowPhotpT_WGJets_MonoPhoton_v18.root");
  // f[1] = new TFile("Out_Fall17_WGJets_LL_estimation_v18_Jul22.root");
  // f[2] = new TFile("Out_Summer16v3_WGJets_LL_estimation_v18_Jul22.root");
 f[0] = new TFile("Out_Autumn18_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[1] = new TFile("Out_Fall17_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[2] = new TFile("Out_Summer16v3_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[3]= new TFile("Out_FullRun2_TTGJets_LL_estimation_v18_Electron_Aug26.root");

  f[4] = new TFile("Out_Autumn18_WGJets_LL_estimation_v18_Electron_Aug26.root");
  f[5] = new TFile("Out_Fall17_WGJets_LL_estimation_v18_Electron_Aug26.root");
  f[6] = new TFile("Out_Summer16v3_WGJets_LL_estimation_v18_Electron_Aug26.root");

  f[7]= new TFile("Out_FullRun2_WGJets_LL_estimation_v18_Electron_Aug26.root");

  //Out_Autumn18_TTGJets_LL_estimation_v18_wBDTcut_STcorr.root
  THStack *hs_var=new THStack("var_Stack","");

  //  const char* baseline[11] = {"nocut","bkg_comp","veto-lep","veto-iso","Phot-pt_20","NhadJetsCut","dPhi_Met","Met_100","Met_250","st_300_Met250","pt_st_Met_250"};
  //  const char* baseline[9]= {"Nocut","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250","nocut"};
  //const char* baseline1[11]= {"Nocut","nocut_sam","nocut","Met_100","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250"};
  //const char *baseline1[25]={"Nocut","photon_smuon","Phot_pT_20","nHadJets_2","MET_100","ST_300","bkg_comp","Met_cleaning","lept_veto","veto_chargedTracks","dPhi_MET","jet_pT_Pho_pT","MET_250","pho_pt_100","Final","Pho_pT_30","HT_1TeV_Met250","HT_1TeV_Met250_pt_100","HT_15TeV_Met100","HT_15TeV_Met250","HT_15TeV_Met250_pt_100","HT_175TeV_Met100","nocut_sam","basic_sam"};
  // const char *baseline1[3]={"Nocut","SignalRegion","ControlRegion"};//
  vector<string> baseline = {"Nocut", "PreSmuontion","Electron_CR","Electron_SR","FailAcep_ElectronSR","FailId_ElectronSR","FailIso_ElectronSR","Mu_CR","Mu_SR","FailAcep_MuSR","FailId_MuSR","FailIso_MuSR","Electron_CR_BDTcut1","Electron_SR_BDTcut1","FailAcep_ElectronSR_BDTcut1","FailId_ElectronSR_BDTcut1","FailIso_ElectronSR_BDTcut1","Mu_CR_BDTcut1","Mu_SR_BDTcut1","FailAcep_MuSR_BDTcut1","FailId_MuSR_BDTcut1","FailIso_MuSR_BDTcut1","Electron_CR_BDTcut2","Electron_SR_BDTcut2","FailAcep_ElectronSR_BDTcut2","FailId_ElectronSR_BDTcut2","FailIso_ElectronSR_BDTcut2","Mu_CR_BDTcut2","Mu_SR_BDTcut2","FailAcep_MuSR_BDTcut2","FailId_MuSR_BDTcut2","FailIso_MuSR_BDTcut2"};

     vector<string> baseline1 = {"FailId_ElecSR","FailIso_ElecSR","FailAcep_ElecSR"};
  //   vector<string> baseline1 ={"FailId_ElecSR_BDTcut1","FailIso_ElecSR_BDTcut1","FailAcep_ElecSR_BDTcut1"};//,"FailId_ElecSR_BDTcut1","FailIso_ElecSR_BDTcut1"};//,"Elec_SR"};
  
  //   vector<string> baseline1 = {"FailId_MuSR_BDTcut1","FailIso_MuSR_BDTcut1","FailAcep_MuSR_BDTcut1"};
  //      vector<string> baseline1 ={"FailId_MuSR","FailIso_MuSR","FailAcep_MuSR"};
   // const char *baseline[48]={"Nocut","SignalRegion","ControlRegion","lostElectron_SR","lostMu_SR","lostTau_SR","else_pho_SR","lostElectron_CR","lostMu_CR","lostTau_CR","else_pho_CR","else_SR","else_CR","lostElectronTau_SR","lostMuTau_SR","Tau_hadronic_SR","lostElectron_SR_Accept","lostElectron_SR_ident","SignalRegion_BDTcut1","ControlRegion_BDTcut1","lostElectron_SR_BDTcut1","lostMu_SR_BDTcut1","lostTau_SR_BDTcut1","else_pho_SR_BDTcut1","else_SR_BDTcut1","lostElectronTau_SR_BDTcut1","lostMuTau_SR_BDTcut1","Tau_hadronic_SR_BDTcut1","lostElectron_CR_BDTcut1","lostMu_CR_BDTcut1","lostTau_CR_BDTcut1","else_pho_CR_BDTcut1","else_CR_BDTcut1","SignalRegion_BDTcut2","ControlRegion_BDTcut2","lostElectron_SR_BDTcut2","lostMu_SR_BDTcut2","lostTau_SR_BDTcut2","else_pho_SR_BDTcut2","else_SR_BDTcut2","lostElectronTau_SR_BDTcut2","lostMuTau_SR_BDTcut2","Tau_hadronic_SR_BDTcut2","lostElectron_CR_BDTcut2","lostMu_CR_BDTcut2","lostTau_CR_BDTcut2","else_pho_CR_BDTcut2","else_CR_BDTcut2"};
   //   const char *baseline1[11]={"Nocut","lostElectron_SR","lostMu_SR","lostTau_SR","else_pho_SR","else_SR","SignalRegion","lostElectron_SR_iso","lostElectron_SR_Accept","lostElectron_SR_ident"}; 
   //  const char *baseline1[11]={"Nocut","lostElectron_SR_BDTcut2","lostMu_SR_BDTcut2","lostTau_SR_BDTcut2","else_pho_SR_BDTcut2","else_SR_BDTcut2","SignalRegion_BDTcut2","lostElectron_SR_iso_BDTcut2","lostElectron_SR_Accept_BDTcut2","lostElectron_SR_ident"};
  //   const char *baseline1[11]={"Nocut","lostElectron_SR_BDTcut1","lostMu_SR_BDTcut1","lostTau_SR_BDTcut1","else_pho_SR_BDTcut1","else_SR_BDTcut1","SignalRegion_BDTcut1","lostElectron_SR_iso_BDTcut1","lostElectron_SR_Accept_BDTcut1","lostElectron_SR_ident"};
  const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
  float energyy[8]={59.74,41.529,35.922,137.19,59.74,41.529,35.922,137.19};
  bool flag=false;
   /* vector<TH1F*> hist_list_Njets; */
   /*    vector<TH1F*> hist_list_Bjets; */
   /*    vector<TH1F*> hist_list_MET; */
   /*    vector<TH1F*> hist_list_PhoPt; */
   /*    //vector<TH1F*> hist_list_Mt;                                                                                                                                    */
   /*    vector<TH1F*> hist_list_ST; */
   /*    vector<TH1F*> hist_list_HT; */

  for(int i_file=0; i_file<8;i_file++)
    {
      vector<TH1F*> hist_list_Njets;
      vector<TH1F*> hist_list_Bjets;
      vector<TH1F*> hist_list_MET;
      vector<TH1F*> hist_list_PhoPt;
      //vector<TH1F*> hist_list_Mt;
      vector<TH1F*> hist_list_ST;
      vector<TH1F*> hist_list_HT;
      THStack *hs_Njets=new THStack("NJets_Stack","");
      int nsel=0.0;
      // if(i_file<4)
      // 	nsel=5;
      // else
      // 	nsel=6;
      for(int i_cut=0; i_cut<baseline1.size();i_cut++)
	{
	  //if(i_cut>0 && i_cut!=6)continue;
	  // if(i_file<6 )
	  //{
	      sprintf(hist_name,"h_NhadJets_%s",baseline1[i_cut].c_str());
	      sprintf(hist_name1,"h_NBJets_%s",baseline1[i_cut].c_str());
	      sprintf(hist_name2,"h_MET_%s",baseline1[i_cut].c_str());
	      sprintf(hist_name3,"h_PhoPt_%s",baseline1[i_cut].c_str());
	      //sprintf(hist_name4,"h_Mt_phoMET_%s",baseline1[i_cut].c_str());
	      sprintf(hist_name5,"h_HT_%s",baseline1[i_cut].c_str());
	      sprintf(hist_name6,"h_St_%s",baseline1[i_cut].c_str());
	      //}
	  // else
	  //   {
	  //     sprintf(hist_name,"h_NhadJets_%s",baseline[i_cut]);
	  //     sprintf(hist_name1,"h_NBJets_%s",baseline[i_cut]);
	  //     sprintf(hist_name2,"h_MET_%s",baseline[i_cut]);
	  //     sprintf(hist_name3,"h_PhoPt_%s",baseline[i_cut]);
	  //     sprintf(hist_name4,"h_Mt_phoMET_%s",baseline[i_cut]);
	  //     sprintf(hist_name5,"h_dPhi_phoMet_%s",baseline[i_cut]);
	  //     sprintf(hist_name6,"h_St_%s",baseline[i_cut]);
	  //   }
	  TH1F* hist_Njets_temp = (TH1F*)f[i_file]->Get(hist_name);
	  TH1F* hist_Bjets_temp = (TH1F*)f[i_file]->Get(hist_name1);
	  TH1F* hist_MET_temp = (TH1F*)f[i_file]->Get(hist_name2);
	  TH1F* hist_phopt_temp = (TH1F*)f[i_file]->Get(hist_name3);
	  //	  TH1F* hist_Mt_temp = (TH1F*)f[i_file]->Get(hist_name4);
	  TH1F* hist_Ht_temp = (TH1F*)f[i_file]->Get(hist_name5);
	  TH1F* hist_ST_temp = (TH1F*)f[i_file]->Get(hist_name6);
	  
	  cout<<hist_Njets_temp->Integral()<<"\t"<<endl;
	  hs_Njets->Add(hist_Njets_temp);
//  decorate(hs_var,i);
 	  hist_MET_temp->Rebin(2);
	  //	  hist_phopt_temp->Rebin(2);
	  hist_Ht_temp->Rebin(2);
	  hist_ST_temp->Rebin(2);
	  hist_Njets_temp->GetXaxis()->SetTitle("N_{jets}");
	  //sprintf(title,"%s",f[i_file]);
	  //hist_Njets_temp->SetTitle("");
	  //	   hist_Bjets_temp->
	  hist_Bjets_temp->GetXaxis()->SetTitle("N_{B-jets}");
	  hist_MET_temp->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
	  hist_phopt_temp->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
	  //	  hist_Mt_temp->GetXaxis()->SetTitle("M_{T}^{#gamma & MET}");
	  hist_Ht_temp->GetXaxis()->SetTitle("HT");
	  //	  hist_Ht_temp->GetXaxis()->SetRangeUser(0,3000);
	  hist_ST_temp->GetXaxis()->SetTitle("EnSum for jets");
	  hist_ST_temp->GetXaxis()->SetRangeUser(0,3000);
	  // if(i_file==0)
	  //   {
	  //     hist_Njets_temp->GetXaxis()->SetRangeUser(0,10);
	  //     hist_Bjets_temp->GetXaxis()->SetRangeUser(0,10);
	  //   }
	  // else if(i_file==1)
          //   {
              hist_Njets_temp->GetXaxis()->SetRangeUser(0,16);
              hist_Bjets_temp->GetXaxis()->SetRangeUser(0,10);
          //   }
	  // else
	  //   {
	  //     hist_Njets_temp->GetXaxis()->SetRangeUser(0,10);
          //     hist_Bjets_temp->GetXaxis()->SetRangeUser(0,6);

	  //
	     
	  hist_MET_temp->GetXaxis()->SetRangeUser(0,1500);
	  hist_phopt_temp->GetXaxis()->SetRangeUser(0,1000);
	  decorate(hist_Njets_temp,i_cut);
          decorate(hist_Bjets_temp,i_cut);
          decorate(hist_MET_temp,i_cut);
          decorate(hist_phopt_temp,i_cut);
          decorate(hist_Ht_temp,i_cut);
          decorate(hist_ST_temp,i_cut);
          hist_Njets_temp->GetXaxis()->SetRangeUser(0,12);
           hist_Bjets_temp->GetXaxis()->SetRangeUser(0,8);
              hist_Ht_temp->GetXaxis()->SetTitle("HT");
          hist_Ht_temp->GetXaxis()->SetRangeUser(0,3000);
          hist_ST_temp->GetXaxis()->SetTitle("EnSum for jets");
          hist_ST_temp->GetXaxis()->SetRangeUser(0,3000);
          hist_MET_temp->GetXaxis()->SetRangeUser(0,1500);
          hist_phopt_temp->GetXaxis()->SetRangeUser(0,1000);

	  hist_list_Njets.push_back(hist_Njets_temp);
	  hist_list_Bjets.push_back(hist_Bjets_temp);
	  hist_list_MET.push_back(hist_MET_temp);
	  hist_list_PhoPt.push_back(hist_phopt_temp);
	  hist_list_HT.push_back(hist_Ht_temp);
	  hist_list_ST.push_back(hist_ST_temp);
	  // decorate(hist_Njets_temp,i_cut);
	  // decorate(hist_Bjets_temp,i_cut);
	  // decorate(hist_MET_temp,i_cut);
	  // decorate(hist_phopt_temp,i_cut);
	  // decorate(hist_Ht_temp,i_cut);
	  // decorate(hist_ST_temp,i_cut);
	  // hist_Njets_temp->GetXaxis()->SetRangeUser(0,12);
	  //  hist_Bjets_temp->GetXaxis()->SetRangeUser(0,8);
	  //     hist_Ht_temp->GetXaxis()->SetTitle("HT");
          // hist_Ht_temp->GetXaxis()->SetRangeUser(0,3000);
          // hist_ST_temp->GetXaxis()->SetTitle("EnSum for jets");
          // hist_ST_temp->GetXaxis()->SetRangeUser(0,3000);
	  // hist_MET_temp->GetXaxis()->SetRangeUser(0,2000);
          // hist_phopt_temp->GetXaxis()->SetRangeUser(0,1500);


	  //hist_list_dPhi.push_back(hist_dPhi_temp);
	}
      // hs_var->SetMinimum(yMin);
      // hs_var->SetMaximum(yMax);
      //      decorate(hs_var,i_file);
// decorate(hs_Njets,i_file);
///hs_Njets->Draw("BAR HIST");
      // }
      //path to save the png file
      float energy=energyy[i_file];
      //path to save the file
      /* sprintf(full_path6,"Overlay_dPhi_photon_MET_%s",filetag[i_file]); */
      /* generate_1Dplot(hist_list_dPhi,full_path6,energy,leg_head,false,true,false,true); */
      // if(i_file>4)
      // 	flag = true;
      sprintf(full_path,"NJets_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_Njets,full_path,energy,leg_head,flag,false,true,false,true,"N_{Jets}",16,filetag[i_file]);
      sprintf(full_path1,"NBJets_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_Bjets,full_path1,energy,leg_head,flag,false,true,false,true,"N_{B-Jets}",8,filetag[i_file]);
      sprintf(full_path2,"MET_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_MET,full_path2,energy,leg_head,flag,false,true,true,true,"p_{T}^{Miss} [GeV]",1500,filetag[i_file]);
      sprintf(full_path3,"PhoPt_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_PhoPt,full_path3,energy,leg_head,flag,false,true,true,true,"p_{T}^{#gamma}",1500,filetag[i_file]);
      sprintf(full_path4,"Ht_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_HT,full_path4,energy,leg_head,flag,false,true,true,true,"HT",3000,filetag[i_file]);
      sprintf(full_path5,"ST_%s_ElectronSR_CR",filetag[i_file]);
      generate_1Dplot(hist_list_ST,full_path5,energy,leg_head,flag,false,true,true,true,"ST",3000,filetag[i_file]);
	

      }
      

}






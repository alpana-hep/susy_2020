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

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=13,nBG=6;    //Specify no. of files
TFile *f[nfiles];
bool savePlots=1;
bool isPaper=0;
// bool v_17=false, v_12=true, v_12_vinay= false;
bool v_17=true, v_12=false, v_12_vinay= false;
// bool v_17=false, v_12=false, v_12_vinay=true;

//int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kMagenta+2,kBlue,kCyan,kRed,kBlue+2,kMagenta};  //Specify Colors b's
////int col[11]={kTeal+9,kGreen,kYellow,kOrange,kPink+1,kPink-2,kBlue,kCyan,kRed,kBlue+2,kMagenta};  //Specify Colors b's
//int col[11]={kPink-2,kTeal+9,kGreen,kYellow,kOrange,kBlue,kCyan,kRed,kBlue+2,kMagenta,kPink+1};  //Specify Colors b's
vector<int> col={kBlack, kBlue, kRed, kGreen+2, kMagenta, kBlue+2, kOrange,kGreen+3,kTeal+9 };//{kPink+1,kTeal+9,kYellow,kGray,kOrange,kCyan,kBlue,kRed,kBlue+2,kMagenta,kCyan};  //Specify Colors b's

TString legend_text[11] = {"No cuts","1-photon","pt>30","lep-veto","isotrk-veto","Njets>=2","MET>100","ST>300","falana","ST>300","Pho-pt>100"};

TCanvas *c_cA=new TCanvas("kinVar","plot of a kin var",1500,900);

string getfname(const char *fname1){string fname=fname1;fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();fname.pop_back();return fname;}
void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
// TLegend *legend1=new TLegend(0.4, 0.75,  0.87, 0.88);
// TLegend *legend2=new TLegend(0.38, 0.6,  0.85, 0.72);

//TLegend *legend1=new TLegend(0.5660881,0.5702076,0.8491322,0.6996337);
TLegend *legend1=new TLegend(0.5553672,0.7299145,0.85247,0.8706337);
//TLegend *legend2=new TLegend(0.2803738,0.7350427,0.8190921,0.8669109);
TLegend *legend2=new TLegend(0.15,0.7350427,0.55190921,0.8669109);

void setLastBinAsOverFlow(TH1D*);
TH1D* setMyRange(TH1D*,double,double);

void plotKinStack_1(){
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
 f[6] = new TFile("out_pMSSM_MCMC_86_7257_Mg1400.root");
  f[7] = new TFile("out_pMSSM_MCMC_86_7257_Mg1500.root");
  f[8] = new TFile("out_pMSSM_MCMC_86_7257_Mg1600.root");
  f[9] = new TFile("out_pMSSM_MCMC_86_7257_Mg1700.root");
  f[10] = new TFile("out_pMSSM_MCMC_86_7257_Mg1800.root");
  f[11] = new TFile("out_pMSSM_MCMC_86_7257_Mg1900.root");
  f[12] = new TFile("out_pMSSM_MCMC_86_7257_Mg2000.root");
  //const char* baseline1[9]= {"Nocut","nphotons_1","Phot_pT_30","lep_veto","iso_trk","nHadJets_2","MET_100","ST_300","dPhi_Met"};
   const char* baseline1[9]= {"Nocut","nphotons_1","Phot_pT_30","lep_veto","iso_trk","nHadJets_2","MET_100","ST_300","Met_250"};
 const char* xlabel_[5]={"Sum of P_{T}^{Jets} & P_{T}^{#gamma} (GeV)","p_{T}^{miss} (GeV)","p_{T}^{#gamma} (GeV)","N_{ jets}","N_{ b-jets}"};
 int rebin[5]={4,4,1,1,1};
 // int yMin_[5]={0.001,
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  double yMin=0.001,yMax = 10000;
  double xMin=0.0,xMax = 1000;
   f[0] = new TFile("Out_TTGJets_v18.root");
   f[1] = new TFile("Out_TTjets_v18.root");
   //f[1]= new TFile("Out_TTJets.root");
   f[2]= new TFile("Out_ZJets_Gamma_v18.root");
   f[5] = new TFile("Out_GJets_DR_combined_v18.root");
   f[3] = new TFile("Out_WGJets_v18.root");
   f[4] = new TFile("Out_WJets_v18.root");

  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","MET Stacked");
  //TH1D *h_R;
  TH1D *h_MET_R[nfiles];
  for(int i=0;i<nfiles;i++){
    sprintf(name,"hist_file%i",i);
    h_MET_R[i]=new TH1D(name,name,21,0.5,21.5);
  }
  vector<double> Bcnt;
  double intLumi=0.0;
  TLatex tl1;
 for(int i_file=6; i_file<13;i_file++)
    {
      vector<TH1F*> hist_list_Njets;
      vector<TH1F*> hist_list_Bjets;
      vector<TH1F*> hist_list_MET;
      vector<TH1F*> hist_list_PhoPt;
  vector<TH1F*> hist_list_ST;
      vector<TH1F*> hist_list_HT;
      for(int i_cut=0; i_cut<8;i_cut++)
        {
          sprintf(hist_name,"h_NhadJets_%s",baseline1[i_cut]);
	  sprintf(hist_name1,"h_NBJets_%s",baseline1[i_cut]);
	  sprintf(hist_name2,"h_MET_%s",baseline1[i_cut]);
	  sprintf(hist_name3,"h_PhoPt_%s",baseline1[i_cut]);
	  sprintf(hist_name5,"h_HT_%s",baseline1[i_cut]);
	  sprintf(hist_name6,"h_St_%s",baseline1[i_cut]);
	  TH1F* hist_Njets_temp = (TH1F*)f[i_file]->Get(hist_name);
          TH1F* hist_Bjets_temp = (TH1F*)f[i_file]->Get(hist_name1);
          TH1F* hist_MET_temp = (TH1F*)f[i_file]->Get(hist_name2);
	  TH1F* hist_phopt_temp = (TH1F*)f[i_file]->Get(hist_name3);
	  TH1F* hist_Ht_temp = (TH1F*)f[i_file]->Get(hist_name5);
          TH1F* hist_ST_temp = (TH1F*)f[i_file]->Get(hist_name6);
	  hist_list_Njets.push_back(hist_Njets_temp);
          hist_list_Bjets.push_back(hist_Bjets_temp);
          hist_list_MET.push_back(hist_MET_temp);
          hist_list_PhoPt.push_back(hist_phopt_temp);
          hist_list_HT.push_back(hist_Ht_temp);
          hist_list_ST.push_back(hist_ST_temp);

	}
    }
 
    h_MET->Rebin(rebin);
    //    h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    //    h_MET->SetMinimum(100);
    //    decorate(h_MET,i,f[i]->GetName());
    setLastBinAsOverFlow(h_MET);

    h_MET = setMyRange(h_MET,xMin,xMax);
    //    h_MET->GetXaxis()->SetRangeUser(xMin,xMax);
    
    if(i<=(nBG-1))  hs_var->Add(h_MET);

    if(i==nBG-1) {
      c_cA->cd();
      hs_var->Draw("BAR HIST");
      hs_var->Draw("HIST");
      hs_var->SetMinimum(yMin);
      hs_var->SetMaximum(yMax);
      //      decorate(hs_var,i,f[i]->GetName());
      //setLastBinAsOverFlow(h_MET);
      if(xMin > -10000 && xMax < 10000) hs_var->GetXaxis()->SetRangeUser(xMin-0.1,0.1+xMax);
    }
    if(i>=nBG){ 
      c_cA->cd(); 
      h_MET->SetMarkerStyle(20);
      h_MET->SetMarkerColor(col[i]);
      h_MET->SetLineColor(col[i]);
      h_MET->SetLineWidth(3);
      if(i>=10)  h_MET->SetLineStyle(2);
      h_MET->Draw("hist same");
      //      h_MET->GetYaxis()->SetRangeUser(0.5,20000);
      //      h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    }
    drawlegend(h_MET,i,f[i]->GetName());
    if(i==nfiles-1){ 
      hs_var->GetXaxis()->SetTitleOffset(1.0);
      hs_var->GetXaxis()->SetTitle(xLabel); hs_var->GetYaxis()->SetTitle("Events");hs_var->SetTitle(0);
      hs_var->GetYaxis()->SetTitleOffset(.90);
      TString temp=h_MET->GetName(),temp2;
      if(temp.Contains("nHadJets") || temp.Contains("nBTags")){
	gPad->SetTickx(0);
	hs_var->GetXaxis()->SetLabelSize(0.08);
	for(int i=1;i<=h_MET->GetNbinsX();i++){
	  temp2 = to_string(i-1);
	  if(i%2==0 && temp.Contains("nHadJets")) continue;
	    hs_var->GetXaxis()->SetBinLabel(i,temp2);
	}
	//      cout<<hist->GetName()<<endl;
      }
    }
    
  }

  legend1->SetFillStyle(0); legend2->SetFillStyle(0);  
  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  legend2->SetBorderSize(0);
  legend2->SetMargin(0.12);
  c_cA->cd(); gPad->SetLogy();legend1->Draw();
  c_cA->cd(); gPad->SetLogy();legend2->Draw();
  //  gPad->RedrawAxis();
  //  hs_var->GetXaxis()->SetTitle(xLabel);
 
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  if(isPaper) textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Supplementary}}");
  else textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.7,0.91,name2);
  TLatex Tl;
  Tl.SetTextSize(0.04);
  // if(isPaper) Tl.DrawLatexNDC(0.48,0.91,"#bf{arXiv:xxxx.xxxxx}");

  if(varName == "mindPhi1dPhi2"){
    TLine *line1=new TLine( 0.3,0.11,  0.3,yMax);
    line1->Draw();
    line1->SetLineStyle(2);
    TArrow *arrow1 = new TArrow(0.3,100,1.2,100,0.01,"|>");
    arrow1->Draw();
    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatex(0.33,140,"#bf{Signal Region}");
    //  Tl.SetTextSize(0.04);
    //  Tl.DrawLatexNDC(0.48,0.91,"#bf{arXiv:xxxx.xxxxx}");
  }
  // //  c_cB->SaveAs("searchBins.png");
  //------------------------------------
  c_cA->cd(); c_cA->SetGridx(0); c_cA->SetGridy(0);
  if(varName=="h_Sbins_v6_withOnlyBL_Selec_Met100"){
    // TLine *line1V7=new TLine( 7.0,0.1,  7.0,100000);
    // TLine *line2V7=new TLine(13.0,0.1, 13.0,100000);
    // TLine *line3V7=new TLine(19.0,0.1, 19.0,100000);
    // TLine *line4V7=new TLine(25.0,0.1, 25.0,100000);
    // TLine *line5V7=new TLine(31.0,0.1, 31.0,100000);
    TLine *line1V7=new TLine( 7.0,0.1,  7.0,100000);
    TLine *line2V7=new TLine(13.0,0.1, 13.0,100000);
    TLine *line3V7=new TLine(19.0,0.1, 19.0,100000);
    TLine *line4V7=new TLine(25.0,0.1, 25.0,100000);
    TLine *line5V7=new TLine(31.0,0.1, 31.0,100000);
    //    TLine *line6V7=new TLine(31.5,0.1, 31.5,10000);
    

    c_cA->cd(); c_cA->SetGridx(0); c_cA->SetGridy(0);
    line1V7->Draw();      line2V7->Draw();  line3V7->Draw();
    line4V7->Draw();      line5V7->Draw(); //line6V7->Draw();

    TArrow *arrow1 = new TArrow( 0.0,100000, 7.0,100000,0.01,"<|>");
    TArrow *arrow2 = new TArrow( 7.0,100000,13.0,100000,0.01,"<|>");
    TArrow *arrow3 = new TArrow(13.0,100000,19.0,100000,0.01,"<|>");
    TArrow *arrow4 = new TArrow(19.0,100000, 25.0,100000,0.01,"<|>");
    TArrow *arrow5 = new TArrow(25.0,100000, 31.0,100000,0.01,"<|>");
    TArrow *arrow6 = new TArrow(31.0,100000, 38.0,100000,0.01,"<|>");

    arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
    arrow4->Draw(); arrow5->Draw(); arrow6->Draw();

    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatex(3.5,10000,"N^{ 0}_{ 2-4}");
    Tl.DrawLatex(10.5,10000,"N^{ 0}_{ 5-6}");
    Tl.DrawLatex(15.5,10000,"N^{ 0}_{ #geq7}");
    Tl.DrawLatex(21.5,10000,"N^{ #geq1}_{ 2-4}");
    Tl.DrawLatex(26.5,10000,"N^{ #geq1}_{ 5-6}");
    Tl.DrawLatex(33.5,10000,"N^{ #geq1}_{ #geq7}");
  }
  //------------------------------------
  if(savePlots){
    TString saveName = "Results/supp_Sim_"+varName;
    TString modelName = f[6]->GetName();
    if(modelName.Contains("T5bbbb")) modelName = "T5bbbbZG";
    else if(modelName.Contains("T5qqqq")) modelName = "T5qqqqHG";
    else if(modelName.Contains("T5tttt")) modelName = "T5ttttZG";
    else if(modelName.Contains("T6tt")) modelName = "T6ttZG";
    saveName = saveName+"_MET100GeV.png";
    //saveName = saveName+"_"+modelName+".png";
    c_cA->SaveAs(saveName);
  }
}

void decorate(THStack *hs,int i,const char* fname){
  //  hs->SetMinimum(0.5);
  //hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.05);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  //  gPad->Update();
  gStyle->SetOptStat(0);
}
void decorate(TH1D* hist,int i,const char* fname){
  hist->SetLineColor(col[i]);
  if(i<nBG) {
    hist->SetFillColor(col[i]);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
  }
  else hist->SetLineWidth(2);
  hist->SetTitle(0);
  hist->GetXaxis()->SetLabelSize(.06);
  hist->GetYaxis()->SetLabelSize(.06);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.06);
  // drawlegend(hist,i,fname);
  //  gPad->Update();
  setLastBinAsOverFlow(hist);
  gStyle->SetOptStat(0);
  
  //Hlist.Add(hist);
}

void drawlegend(TH1D *hist,int i,const char* fname){
  gStyle->SetLegendBorderSize(0);
 
  TString lName=fname;
  
  if(lName.Contains("ZGZJ")){lName="Z(#nu#bar{#nu}) + #gamma";}
  //  else if(lName.Contains("ZJets")){lName="Z(#nu#bar{#nu}) + jets";}
  else if(lName.Contains("DYJetsToLL")){lName="DY(l^{+}l^{-})";}
  else if(lName.Contains("WJets")){lName="W(l#nu) + jets";}
  else if(lName.Contains("RareProcess")){}
  else if(lName.Contains("TTjets")){lName="t #bar{t}";}
  else if(lName.Contains("WGJets")){lName="W(l#nu) + #gamma";}
  else if(lName.Contains("ZJets")){lName="Z(#nu#bar{#nu}) + #gamma";}
  else if(lName.Contains("TTGJets")){lName="t #bar{t} + #gamma";}
  //  else if(lName.Contains("QCD")){lName="QCD";}
  else if(lName.Contains("GJets")){lName="#gamma + jets";}
  else if(lName.Contains("Run2016")){lName="Data";}
  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="T5bbbbZg 1600, 150";}
  else if(lName.Contains("T5bbbbZg_1600_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5bbbbZG_1800_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZG_1800_1750")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}

  else if(lName.Contains("T5bbbbZg_1800_150")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5bbbbZg_1800_1750")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}
  else if(lName.Contains("T5bbbbZg_1800_1000")){lName = "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)";}
  else if(lName.Contains("T5qqqqHg_1800_150")){lName = "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5qqqqHg_1800_1750")){lName = "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1750 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_150")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_1550")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5ttttZg_1800_1000")){lName = "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)";}

  else if(lName.Contains("T6ttZg_1000_100")){lName = "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)";}
  else if(lName.Contains("T6ttZg_1000_900")){lName = "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 900 GeV)";}

  // else if(lName.Contains("T5bbbbZg_1600_150")){lName = "T5bbbbZG (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 150 GeV)";}
  // else if(lName.Contains("T5bbbbZg_1600_1550")){lName = "T5bbbbZG (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 1550 GeV)";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName="T5bbbbZg 1600, 1550";}
  else if(lName.Contains("T5qqqqHg_1600_1550")){lName="T5qqqqHg_1550";}
  else if(lName.Contains("T5qqqqHg_1600_150")){lName="T5qqqqHg_150";}
  else if(lName.Contains("TChiWg_0_400")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 400 GeV";}
  else if(lName.Contains("TChiWg_0_800")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 800 GeV";}
  else if(lName.Contains("TChiWg_0_1200")){lName="TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 1200 GeV";}
  else if(lName.Contains("T5bbbbZG_2350_10")){lName="#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2350 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)";}
  else if(lName.Contains("T5bbbbZG_2350_50")){lName="#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2350 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)";}
  else if(lName.Contains("T5bbbbZG_2350_200")){lName="#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2350 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)";}
  else if(lName.Contains("T5bbbbZG_2350_1500")){lName="#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2350 GeV, m_{#tilde{#chi}_{1}^{0}} = 1500 GeV)";}


  else if(lName.Contains("pMSSM_MCMC_106_19786")){ lName = "pMSSM_MCMC_106_19786";}
    else if(lName.Contains("pMSSM_MCMC_473_54451")){ lName = "pMSSM_MCMC_473_54451";}
   else if(lName.Contains("pMSSM_MCMC_70_90438")){ lName = "pMSSM_MCMC_70_90438";}
  else if(lName.Contains("pMSSM_MCMC_86_7257_Mg1400")){ lName = "pMSSM_86_7257_mg1400";}
  else if(lName.Contains("pMSSM_MCMC_86_7257_Mg1500")){ lName = "pMSSM_86_7257_mg1500";}
else if(lName.Contains("pMSSM_MCMC_86_7257_Mg1600")){ lName = "pMSSM_86_7257_mg1600";}
  else if(lName.Contains("pMSSM_MCMC_86_7257_Mg1700")){ lName = "pMSSM_86_7257_mg1700";}


  // const char *l_name=lName.c_str();
  if(i<nBG)legend1->AddEntry(hist,lName,"f");
  else legend2->AddEntry(hist,lName,"l");
  // legend1->SetTextSize(0.04);
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

TH1D* setMyRange(TH1D *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;
  // h1->Print("all");
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);
  //  cout<<nMax<<endl;
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
    //    cout<<":";
  }
  //  h1->Print("all");
  //  cout<<endl;
  return h1;
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);                                                                                                                      
}

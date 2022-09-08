const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[4] = {"(0#mu,0e) SR","(0#mu,1e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[4] = {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
int line_color[9] = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
TH1F* setLastBinAsOverFlow(TH1F*, int);
TH1F* setMyRange(TH1F*,double,double);
TH1F* DrawOverflow(TH1F*);
TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){
    //function to paint the histogram h with an extra bin for overflows
       // This function paint the histogram h with an extra bin for overflows
   UInt_t nx    = h->GetNbinsX()+1;
   Double_t *xbins= new Double_t[nx+1];
   for (UInt_t i=0;i<nx;i++)
     xbins[i]=h->GetBinLowEdge(i+1);
   xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
   char *tempName= new char[strlen(h->GetName())+10];
   sprintf(tempName,"%swtOverFlow",h->GetName());
   h->GetXaxis()->SetLimits(xmin,xrange);
   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
   htmp->GetXaxis()->SetRange(xmin,xrange);
   // Reset the axis labels
   htmp->SetXTitle(h->GetXaxis()->GetTitle());
   htmp->SetYTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++)
     htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   // Fill the underflows
   htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   // Restore the number of entries
   htmp->SetEntries(h->GetEntries());
   // FillStyle and color
   // htmp->SetFillStyle(h->GetFillStyle());
   // htmp->SetFillColor(h->GetFillColor());
   return htmp;
}
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
TH1F* setLastBinAsOverFlow(TH1F* h_hist, int xrange){
  //     h_hist = setMyRange(h_hist,0,xrange);
  //  h_hist->GetXaxis()->SetRangeUser(0,xrange);
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX());
  //  cout<<h_hist->GetNbinsX()<<"\t"<<lastBinCt<<"\t"<<overflCt<<endl;

  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;
  //h_temp->GetXaxis()->SetRangeUser(0,xrange);

  lastBinCt = lastBinCt+overflCt;
  //  cout<<lastBinCt<<endl;
  TH1F* h_temp = (TH1F*)h_hist->Clone();
  h_temp->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_temp->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  //  h_temp->GetXaxis()->SetRangeUser(0,xrange);

  // h_hist = setMyRange(h_hist,0,xrange);
  //
  return h_temp;
}


TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow                                                                                                    
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;

  // h1->Print("all");
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);  
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);

  //  cout<<nMax<<endl;
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
    //    cout<<":";
    //h1->GetXaxis()->SetRangeUser(xLow,xHigh); 
  }
  return h1;
}
void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", float energy=-1, int xrange=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title=""){  
  

  TCanvas *c = new TCanvas(tag_name, tag_name, 700, 600);
  c->SetRightMargin(0.009);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.025);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.06);
   c->SetBottomMargin(0.12);

    gStyle->SetOptStat(1111111);
  //   gStyle->SetOptStat(0);
  //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  gStyle->SetOptStat(0);

  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.65,0.65,0.9,0.9);  
  legend->SetTextSize(0.045);
  legend->SetLineColor(kWhite);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
 legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;


  for(int i = 0; i < (int)hist.size(); i++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    //    normalize = true;
    if(normalize) {
      	hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
    //    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange);
    hist.at(i)->SetLineWidth(line_width[i]);
    
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    //    hist.at(i)->GetXaxis()->SetRange(0,1500);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.2);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);
    //   //hist.at(i)->Rebin(1);
    // }

  //     double lastBinCt =hist.at(i)->GetBinContent(hist.at(i)->GetNbinsX()),overflCt =hist.at(i)->GetBinContent(hist.at(i)->GetNbinsX()+1);
  // double lastBinErr=hist.at(i)->GetBinError(hist.at(i)->GetNbinsX()),  overflErr=hist.at(i)->GetBinError(hist.at(i)->GetNbinsX()+1);
  // if(lastBinCt!=0 && overflCt!=0)
  //   lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  // else if(lastBinCt==0 && overflCt!=0)
  //   lastBinErr = overflErr;
  // else if(lastBinCt!=0 && overflCt==0)
  //   lastBinErr = lastBinErr;
  // else lastBinErr=0;

  // lastBinCt = lastBinCt+overflCt;
  // hist.at(i)->SetBinContent(hist.at(i)->GetNbinsX(),lastBinCt);
  // hist.at(i)->SetBinError(hist.at(i)->GetNbinsX(),lastBinErr);
  //  hist.at(i)->GetXaxis()->SetRange(1, hist.at(i)->GetNbinsX() + 1);
    /* hist.at(i)->GetXaxis()->SetRangeUser(x_min[energy],x_max[energy]); */
    //    hist.at(i)= DrawOverflow(hist.at(i));
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //setLastBinAsOverFlow(hist.at(i),0);
    
    

  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.001,ymax*10);
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.00001,ymax*10.0);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
	
      }
  	if(!i) hist.at(i)->Draw("hist");
    else   hist.at(i)->Draw("hist sames");

  }
  
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
    sprintf(canvas_name,"./Results/e_overlayPlots/%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
     c->SaveAs(canvas_name);   
     sprintf(canvas_name,"./Results/e_overlayPlots/%s.pdf",tag_name);
    c->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=10,nBG=6;                                                                                                                                                              
TFile *f[nfiles];


void overlayPlots_v1()
{
  char* hname = new char[200];
  char* hist_name  = new char[200];
  char* hist_name1 = new char[200];
  char* hist_name2 = new char[200];
  char* hist_name3 = new char[200];
  char* hist_name4 = new char[200];
  char* hist_name5 = new char[200];
  char* hist_name6 = new char[200];
  char* hist_name7 = new char[200];
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
  sprintf(path2,"Results/dPhi_METHadJets");
  // TFile *f[nfiles];                                                                                                                                                                                    
  /* f[0] = new TFile("Out_TTGJets_v18_Dphi_MetHadJets_corr_v2.root"); */
  /* f[1] = new TFile("Out_TTjets_v18_Dphi_MetHadJets_corr_v2.root");//ZG                                                                                                                                         */
  /* f[2]= new TFile("Out_ZJets_Gamma_v18_Dphi_MetHadJets_corr_v2.root"); */
  /* f[5] = new TFile("Out_GJets_DR_combined_v18_Dphi_MetHadJets_corr_v2.root"); */
  f[0] = new TFile("Out_Autumn18_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[1] = new TFile("Out_Fall17_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[2] = new TFile("Out_Summer16v3_TTGJets_LL_estimation_v18_Electron_Aug26.root");
  f[3]= new TFile("Out_FullRun2_TTGJets_LL_estimation_v18_Electron_Aug26.root");

  f[4] = new TFile("Out_Autumn18_WGJets_LL_estimation_v18_Electron_Aug26.root");
  f[5] = new TFile("Out_Fall17_WGJets_LL_estimation_v18_Electron_Aug26.root");
  f[6] = new TFile("Out_Summer16v3_WGJets_LL_estimation_v18_Electron_Aug26.root");

  f[7]= new TFile("Out_FullRun2_WGJets_LL_estimation_v18_Electron_Aug26.root");

  //f[0] = new TFile("Out_Autumn18_TTGJets_LL_estimation_v18_wBDTcut_STcorr.root");
  //  f[0] = new TFile("Out_Autumn18_WGJets_LL_estimation_v18_Jul22.root");//Autumn18_lowPhotpT_WGJets_MonoPhoton_v18.root");
  /* f[4] = new TFile("Out_WJets_v18_Dphi_MetHadJets_corr_v2.root"); */
  /* f[6] = new TFile("OutT5bbbbZg_2350_10_Dphi_MetHadJets_corr_v2.root"); */
  /* f[7] = new TFile("OutT5bbbbZg_2350_50_Dphi_MetHadJets_corr_v2.root"); */
  /* f[8] = new TFile("OutT5bbbbZg_2350_200_Dphi_MetHadJets_corr_v2.root"); */
  /* f[9] = new TFile("OutT5bbbbZg_2350_1500_Dphi_MetHadJets_corr_v2.root"); */
  // f[0] = new TFile("Out_pMSSM_MCMC_86_7257.root");
  // f[1] = new TFile("Out_pMSSM_MCMC_70_90438.root");
  // f[2] = new TFile("Out_pMSSM_MCMC_106_19786.root");
  // f[3] = new TFile("Out_pMSSM_MCMC_473_54451.root");

  const char* baseline[11] = {"nocut","bkg_comp","veto-lep","veto-iso","Phot-pt_20","NhadJetsCut","dPhi_Met","Met_100","Met_250","st_300_Met250","pt_st_Met_250"};
  //  const char* baseline[9]= {"Nocut","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250","nocut"};
  //const char* baseline1[11]= {"Nocut","nocut_sam","nocut","Met_100","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250"};
  //const char *baseline1[25]={"Nocut","photon_selec","Phot_pT_20","nHadJets_2","MET_100","ST_300","bkg_comp","Met_cleaning","lept_veto","veto_chargedTracks","dPhi_MET","jet_pT_Pho_pT","MET_250","pho_pt_100","Final","Pho_pT_30","HT_1TeV_Met250","HT_1TeV_Met250_pt_100","HT_15TeV_Met100","HT_15TeV_Met250","HT_15TeV_Met250_pt_100","HT_175TeV_Met100","nocut_sam","basic_sam"};
    const char *baseline1[3]={"Nocut","Elec_SR","Elec_CR"};//
  //   const char *baseline1[3]={"Nocut","Mu_SR","Mu_CR"};
      //   const char *baseline1[9]={"Nocut","SignalRegion","lostElec_SR","lostMu_SR","lostTau_SR","lostElec_SR_iso","lostElec_SR_Accept","lostElec_SR_ident"};
  const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
  float energyy[8]={59.74,41.529,35.922,137.19,59.74,41.529,35.922,137.19};
  bool flag=false;

  //  const char* filetag[10]={"TTGJets","pMSSM_MCMC_70_90438","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451","WJets","GJets","T5bbbbZG_10","T5bbbbZG_50","T5bbbbZG_200","T5bbbbZG_1500"};
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
      vector<TH1F*> hist_BDTresp;
      for(int i_cut=1; i_cut<3;i_cut++)
	{
	  //if(i_cut>0 && i_cut!=6)continue;
	  //	  if(i_file<6 )
	     sprintf(hist_name,"h_NhadJets_%s",baseline1[i_cut]);
	      sprintf(hist_name1,"h_NBJets_%s",baseline1[i_cut]);
	      sprintf(hist_name2,"h_MET_%s",baseline1[i_cut]);
	      sprintf(hist_name3,"h_PhoPt_%s",baseline1[i_cut]);
	      //sprintf(hist_name4,"h_Mt_phoMET_%s",baseline1[i_cut]);
	      sprintf(hist_name5,"h_HT_%s",baseline1[i_cut]);
	      sprintf(hist_name6,"h_St_%s",baseline1[i_cut]);
	      sprintf(hist_name7,"h_BDT_response_%s",baseline1[i_cut]);

	    // }
	  // else
	    
	      // sprintf(hist_name,"h_NhadJets_%s",baseline[i_cut]);
	      // sprintf(hist_name1,"h_NBJets_%s",baseline[i_cut]);
	      // sprintf(hist_name2,"h_MET_%s",baseline[i_cut]);
	      // sprintf(hist_name3,"h_PhoPt_%s",baseline[i_cut]);
	      // sprintf(hist_name4,"h_Mt_phoMET_%s",baseline[i_cut]);
	      // sprintf(hist_name5,"h_dPhi_phoMet_%s",baseline[i_cut]);
	      // sprintf(hist_name6,"h_St_%s",baseline[i_cut]);
	      // sprintf(hist_name7,"h_BDT_response_%s",baseline[i_cut]);
	    
	  TH1F* hist_Njets_temp = (TH1F*)f[i_file]->Get(hist_name);
	  TH1F* hist_Bjets_temp = (TH1F*)f[i_file]->Get(hist_name1);
	  TH1F* hist_MET_temp = (TH1F*)f[i_file]->Get(hist_name2);
	  TH1F* hist_phopt_temp = (TH1F*)f[i_file]->Get(hist_name3);
	  //	  TH1F* hist_Mt_temp = (TH1F*)f[i_file]->Get(hist_name4);
	  TH1F* hist_Ht_temp = (TH1F*)f[i_file]->Get(hist_name5);
	  TH1F* hist_ST_temp = (TH1F*)f[i_file]->Get(hist_name6);
	  TH1F* hist_BDTresp_temp = (TH1F*)f[i_file]->Get(hist_name7);
	  hist_BDTresp_temp->Rebin(4);
	  //	  hist_ST_temp->GetXaxis()->SetRangeUser(0,xrange);
	  // for (int ib=0; ib<hist_ST_temp->GetNbinsX();ib++)
	  //   {
	  //     double lastBinCt =hist_ST_temp->GetBinContent(ib);
	  //     //    hist_ST_temp->GetNbinsX()),
	  //     double overflCt =hist_ST_temp->GetBinContent(hist_ST_temp->GetNbinsX());
	  //     cout<<"Alps-0"<<"\t"<<ib<<"\t"<<i_file<<"\t"<<lastBinCt<<"\t"<<hist_ST_temp->Integral()<<endl;
	  //   }
	  hist_BDTresp_temp->GetXaxis()->SetTitle("BDT response");
	  hist_MET_temp->Rebin(2);
	  hist_phopt_temp->Rebin(2);
	  hist_Ht_temp->Rebin(2);
	  hist_ST_temp->Rebin(2);
	  hist_Njets_temp->GetXaxis()->SetTitle("N_{jets}");
	  //sprintf(title,"%s",f[i_file]);
	  //hist_Njets_temp->SetTitle("");
	  //	   hist_Bjets_temp->
	  hist_Bjets_temp->GetXaxis()->SetTitle("N_{B-jets}");
	  hist_MET_temp->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
	  hist_phopt_temp->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
	   // for (int ib=0; ib<hist_ST_temp->GetNbinsX();ib++)
           //  {
           //    double lastBinCt =hist_ST_temp->GetBinContent(ib);
           //    //    hist_ST_temp->GetNbinsX()),                                                                                                   
           //    double overflCt =hist_ST_temp->GetBinContent(hist_ST_temp->GetNbinsX());
           //    cout<<"Alps-0"<<"\t"<<ib<<"\t"<<i_file<<"\t"<<lastBinCt<<"\t"<<hist_ST_temp->Integral()<<endl;
           //  }

	  //	  hist_Bjets_temp = setMyRange(hist_Bjets_temp,0,12);
	  
	  // setLastBinAsOverFlow(hist_Bjets_temp,12);
	  // setLastBinAsOverFlow(hist_Njets_temp,16);
	  // setLastBinAsOverFlow(hist_MET_temp,1500);
	  // setLastBinAsOverFlow(hist_phopt_temp,1000);
	  // setLastBinAsOverFlow(hist_Ht_temp,3000);
	  // setLastBinAsOverFlow(hist_ST_temp,3000);

	  //	  hist_Mt_temp->GetXaxis()->SetTitle("M_{T}^{#gamma & MET}");
	  hist_Ht_temp->GetXaxis()->SetTitle("HT");
	  hist_Ht_temp->GetXaxis()->SetRangeUser(0,3000);
	  hist_ST_temp->GetXaxis()->SetTitle("EnSum for jets");
	  hist_ST_temp->GetXaxis()->SetRangeUser(0,3000);
	  //	  lastBinCt =float(hist_ST_temp->GetBinContent(hist_ST_temp->GetNbinsX())),overflCt =float(hist_ST_temp->GetBinContent(hist_ST_temp->GetNbinsX()+1));
	  //cout<<"Alps"<<"\t"<<hist_ST_temp->GetNbinsX()<<"\t"<<lastBinCt<<"\t"<<overflCt<<endl;
	  // // if(i_file==0)
	  // //   {
	  hist_Njets_temp->GetXaxis()->SetRangeUser(0,16);
	      hist_Bjets_temp->GetXaxis()->SetRangeUser(0,10);
	  //   }
	  // else if(i_file==1)
          //   {
	  //  hist_Njets_temp->GetXaxis()->SetRangeUser(0,16);
	  //  hist_Bjets_temp->GetXaxis()->SetRangeUser(0,10);
          //   }

	  // else
	  //   {
	  //     hist_Njets_temp->GetXaxis()->SetRangeUser(0,10);
          //     hist_Bjets_temp->GetXaxis()->SetRangeUser(0,6);

	  //   }
	  hist_MET_temp->GetXaxis()->SetRangeUser(0,2000);
	  hist_phopt_temp->GetXaxis()->SetRangeUser(0,1500);
	  //  setLastBinAsOverFlow(hist_Bjets_temp,12);
          // setLastBinAsOverFlow(hist_Njets_temp,16);
          // setLastBinAsOverFlow(hist_MET_temp,1500);
          // setLastBinAsOverFlow(hist_phopt_temp,1000);
          // setLastBinAsOverFlow(hist_Ht_temp,3000);
          // setLastBinAsOverFlow(hist_ST_temp,3000);
	  // TH1F* hist_ST_temp1 =  DrawOverflow(hist_ST_temp);
	  // hist_ST_temp1->GetXaxis()->SetTitle("EnSum for jets");
	  
	  // hist_ST_temp1->GetXaxis()->SetRangeUser(0,3000);

	  hist_list_Njets.push_back(hist_Njets_temp);
	  hist_list_Bjets.push_back(hist_Bjets_temp);
	  hist_list_MET.push_back(hist_MET_temp);
	  hist_list_PhoPt.push_back(hist_phopt_temp);
	  hist_list_HT.push_back(hist_Ht_temp);
	  hist_list_ST.push_back(hist_ST_temp);
	  hist_BDTresp.push_back(hist_BDTresp_temp);
	  // setLastBinAsOverFlow(hist_Bjets_temp,12);
          // setLastBinAsOverFlow(hist_Njets_temp,16);
          // setLastBinAsOverFlow(hist_MET_temp,1500);
          // setLastBinAsOverFlow(hist_phopt_temp,1000);
          // setLastBinAsOverFlow(hist_Ht_temp,3000);
          // setLastBinAsOverFlow(hist_ST_temp,3000);


	  //hist_list_dPhi.push_back(hist_dPhi_temp);
	}
      // }
      //path to save the png file
      float energy=energyy[i_file];
      int xrange=0.0;
      //path to save the file
      /* sprintf(full_path6,"Overlay_dPhi_photon_MET_%s",filetag[i_file]); */
      /* generate_1Dplot(hist_list_dPhi,full_path6,energy,leg_head,false,true,false,true); */
      sprintf(full_path,"Overlay_NJets_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_Njets,full_path,energy,16,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path1,"Overlay_N-BJets_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_Bjets,full_path1,energy,10,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path2,"Overlay_MET_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_MET,full_path2,energy,1500,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path3,"Overlay_PhoPt_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_PhoPt,full_path3,energy,1000,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path4,"Overlay_Ht_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_HT,full_path4,energy,3000,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path5,"Overlay_ST_%s_",filetag[i_file]);
      generate_1Dplot(hist_list_ST,full_path5,energy,3000,0,leg_head,false,true,false,true,filetag[i_file]);
      sprintf(full_path5,"Overlay_BDT_resp_%s_",filetag[i_file]);
      generate_1Dplot(hist_BDTresp,full_path5,energy,1,-1,leg_head,false,true,true,true,filetag[i_file]);


      sprintf(full_path,"Overlay_NJets_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_Njets,full_path,energy,16,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path1,"Overlay_N-BJets_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_Bjets,full_path1,energy,10,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path2,"Overlay_MET_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_MET,full_path2,energy,1500,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path3,"Overlay_PhoPt_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_PhoPt,full_path3,energy,1000,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path4,"Overlay_Ht_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_HT,full_path4,energy,3000,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path5,"Overlay_ST_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_list_ST,full_path5,energy,3000,0,leg_head,true,true,false,true,filetag[i_file]);
      sprintf(full_path5,"Overlay_BDT_resp_%s_normalized",filetag[i_file]);
      generate_1Dplot(hist_BDTresp,full_path5,energy,1,-1,leg_head,true,true,true,true,filetag[i_file]);

      
      }
      

}






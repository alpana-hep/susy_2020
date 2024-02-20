const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[5] ={"t #bar{t} + #gamma","t #bar{t}","W(l#nu) + jets","W(l#nu) + #gamma","Data"};//(1#mu,1#gamma) CR","(0#mu,1#gamma) SR","Failed Id","Failed Iso"};//,"Failed acceptance","1e CR"};// {"(0e,0e) SR","(1e,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//  if(lName.Contains("ZGZJ")){lName="Z(#nu#bar{#nu}) + #gamma";}
  //  else if(lName.Contains("ZJets")){lName="Z(#nu#bar{#nu}) + jets";}                                                                           
  // else if(lName.Contains("DYJetsToLL")){lName="DY(l^{+}l^{-})";}
  // else if(lName.Contains("WJets")){lName="W(l#nu) + jets";}
  // else if(lName.Contains("RareProcess")){}
  // else if(lName.Contains("TTJets")){lName="t #bar{t}";}
  // else if(lName.Contains("WGJets")){lName="W(l#nu) + #gamma";}
  // else if(lName.Contains("ZNuNuGJets")){lName="Z(#nu#bar{#nu}) + #gamma";}
  // else if(lName.Contains("TTGJets")){lName="t #bar{t} + #gamma";}

//TString legend_text[4] = {"(0e,0e) SR","(1e,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,3,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};                                                                               
int line_color[9] = {9,kCyan+2,45,kMagenta,kGray+1,kRed,kBlue+2,kMagenta,kCyan};
//vector <int>col;
//col.resize(0);
vector<int> col={kGray,kTeal+9,kOrange,kCyan-1,kBlack,kBlue,kMagenta,kBlack,kPink+4,kGreen+2,kOrange+1,kBlue+2};

//int line_color[9] = {kMagenta+2, kGray+2, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
//vector<int> col={kBlue,kBlack,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={1001,3021,3144,3244};
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
void decorate(TH1D*,int );
 
void decorate(TH1D* hist,int i){
  hist->SetLineColor(col[i]);
  //  if(i<nBG) {                                                                                                                                 
  hist->SetFillColor(col[i]);
    // if(i!=0)
    //   {
    // 	//	hist->SetFillColor(col[i]);

    // 	hist->SetFillStyle(Style[i]);
    //   }
    // else
    //   {
    // 	hist->SetFillColor(kGray+2);

    // 	hist->SetFillStyle(1001);
    // 	}
    hist->SetLineWidth(2);
    //}                                                                                                                                           
}

void setLastBinAsOverFlow(TH1D*);
TH1D* setMyRange(TH1D*,double,double);
TH1D* DrawOverflow(TH1D*);
TH1D* DrawOverflow(TH1D* h,int xmin, int xrange){
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
   TH1D *htmp = new TH1D(tempName, h->GetTitle(), nx, xbins);
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

// TH1D* setLastBinAsOverFlow(TH1D* h_hist, int xrange){
//   //     h_hist = setMyRange(h_hist,0,xrange);setLastBinAsOverFlow
//   //  h_hist->GetXaxis()->SetRangeUser(0,xrange);
//   double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX());
//   //  cout<<h_hist->GetNbinsX()<<"\t"<<lastBinCt<<"\t"<<overflCt<<endl;

//   double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
//   if(lastBinCt!=0 && overflCt!=0)
//     lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

//   else if(lastBinCt==0 && overflCt!=0)
//     lastBinErr = overflErr;
//   else if(lastBinCt!=0 && overflCt==0)
//     lastBinErr = lastBinErr;
//   else lastBinErr=0;
//   //h_temp->GetXaxis()->SetRangeUser(0,xrange);

//   lastBinCt = lastBinCt+overflCt;
//   //  cout<<lastBinCt<<endl;
//   TH1D* h_temp = (TH1D*)h_hist->Clone();
//   h_temp->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
//   h_temp->SetBinError(h_hist->GetNbinsX(),lastBinErr);
//   //  h_temp->GetXaxis()->SetRangeUser(0,xrange);

//   // h_hist = setMyRange(h_hist,0,xrange);
//   //
//   return h_temp;
// }


TH1D* setMyRange(TH1D *h1,double xLow,double xHigh){
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
void generate_1Dplot(vector<TH1D*> hist, TH1D* hist_ratio, char const *tag_name="", char const *xlabel="",char const *ylabel="",float energyy=0.1,int rebin=-1,double ymin=0,double ymax=0, double xmin=0,double xmax=0,char const *leg_head="",bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title=""){    
   TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
   canvas_n1->Range(-60.25,-0.625,562.25,0.625);
   canvas_n1->SetFillColor(0);
   canvas_n1->SetBorderMode(0);
   canvas_n1->SetBorderSize(2);
   auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
   pad_1->SetTopMargin(0.013);
   pad_1->SetBottomMargin(0.3);
   pad_1->SetRightMargin(0.018);
   pad_1->SetLeftMargin(0.12);
   
   auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
   p1->SetBottomMargin(0.01);
   p1->SetRightMargin(0.018);
   p1->SetLeftMargin(0.12);
   p1->SetTopMargin(0.07);
   p1->cd();
   // p1->SetGrid();
   
//   TCanvas *c = new TCanvas(tag_name, tag_name, 700, 600);
//   c->SetRightMargin(0.009);
//   c->SetLeftMargin(0.11);
//   c->SetTopMargin(0.05);
// c->SetLeftMargin(0.11);
//   c->SetTopMargin(0.05);
//   c->SetRightMargin(0.025);
//   c->SetLeftMargin(0.13);
//   c->SetTopMargin(0.06);
//    c->SetBottomMargin(0.12);
     THStack *hs_var=new THStack("var_Stack","");
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
  legend = new TLegend(0.6,0.6,0.96,0.91);  
  legend->SetTextSize(0.05);
  //  legend->SetLineColor(kWhite);
  legend->SetNColumns(2);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);
   legend->SetFillStyle(0);
  TLegendEntry* leg_entry[11];
  float x_label_size = 0.045;
  // double ymin = 100000.0;
  // double ymax = 0.0;
  double xrange = xmax;
  float energy = energyy;
  cout<<" hist.size() = "<<hist.size()<<endl;
  // if(hist.size()==2)
  //   {  col[1]=col[1]-41;
  //     legend_text[1] = "(0e,0e) SR";
  //   }
  // else
  //   {
  //     col[1]=45;
  //     legend_text[1] = "Failed acceptance";
  //   }
  for(int i = 0; i < (int)hist.size(); i++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(rebin);

    // }
    //        hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
    //    setLastBinAsOverFlow(hist.at(i));

    //    normalize = true;
    if(normalize) {
      	hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
     hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);
    hist.at(i)->SetLineWidth(line_width[i]);
    //    setLastBinAsOverFlow(hist.at(i));
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    //    hist.at(i)->GetXaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    decorate(hist.at(i),i);
    if(i<5){
    legName.push_back(hist.at(i)->GetName());
    //leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"f");
    if(i==4)
      leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"e2p");
    else
      leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"f");
    leg_entry[i]->SetTextColor(kBlack);}//hist.at(i)->GetLineColor());}
    // if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    // if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //setLastBinAsOverFlow(hist.at(i),0);
    
    

  
  // if(ymin == 0.0) ymin = 1e-3;
  // if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  // for(int i = 0; i < (int)hist.size(); i++) {
  //   if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.0001,100.0*ymax);
  //   else
  //     {  hist.at(i)->GetYaxis()->SetRangeUser(0.000001,ymax*60.0);
  // 	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
	
  //     }
    //    p1->SetGrid();
    hist.at(i)= setMyRange(hist.at(i),xmin,xmax+4);
    setLastBinAsOverFlow(hist.at(i));
    if(i<4)
      hs_var->Add(hist.at(i));
    hs_var->SetMinimum(ymin);
    hs_var->SetMaximum(ymax);
  }
    // 	if(!i) hist.at(i)->Draw("");
    // else   hist.at(i)->Draw(" sames");
	

  // hs_var->SetMinimum(ymin);
  // hs_var->SetMaximum(ymax);

  //hs_var->GetXaxis()->SetTitle(xlabel);
  //hs_var->GetYaxis()->SetTitle(ylabel);
  // hs_var->Draw("BAR HIST");
  //  hs_var->GetYaxis()->SetTitle("Entries");

  hs_var->Draw("HIST");
  hs_var->GetYaxis()->SetTitle("Entries");

  hs_var->GetXaxis()->SetRangeUser(xmin,xmax+4);
  hs_var->GetXaxis()->SetTitle(xlabel);
  hs_var->GetXaxis()->SetTitleOffset(1.2);
  //  hist.at(i)->GetXaxis()->SetLabelOffset(1.2);
  gPad->Modified(); gPad->Update();
  hs_var->GetXaxis()->SetTitle(xlabel);
  hs_var->GetXaxis()->SetTitle(ylabel);
  
  //  hs_var->GetYaxis()->SetTitleSize(
  //hs_var->GetYaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");//hs_var->GetYaxis()->SetTitle("Events");
  hs_var->SetTitle(0);
  hs_var->GetYaxis()->SetTitleOffset(1.2);
  hs_var->GetXaxis()->SetTitleSize(00.05);
  hs_var->GetXaxis()->SetLabelSize(0.04);
  hs_var->GetYaxis()->SetLabelSize(0.04);
  hs_var->GetYaxis()->SetTitleSize(00.065);
  hs_var->GetYaxis()->SetTitleOffset(0.8);
  

  legend->Draw();
  hist.at(4)->SetLineWidth(1);
  hist.at(4)->SetMarkerStyle(20);
  hist.at(4)->GetYaxis()->SetTitle("Entries");
  //  hist.at(4)->SetLineStyle(2);
  hist.at(4)->Draw("E1same");

  if(log_flag) {
      gPad->SetLogy();
    }
  // if(logx)
  //   gPad->SetLogx();

  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.045);
  textOnTop->DrawLatexNDC(0.12,0.94,"CMS #it{#bf{Preliminary}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.045);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.78,0.94,en_lat);
   TArrow *arrow1 = new TArrow( 1.0,0.10, 2.0,0.1,0.01,"<|>");
    TArrow *arrow2 = new TArrow( 2.0,0.10,3.0,0.1,0.01,"<|>");
    TArrow *arrow3 = new TArrow(3.0,0.10,4.0,0.1,0.01,"<|>");
    TArrow *arrow4 = new TArrow(4.0,0.10, 5.0,0.1,0.01,"<|>");
    TArrow *arrow5 = new TArrow(5.0,0.10, 7.0,0.1,0.01,"<|>");
    TArrow *arrow6 = new TArrow(7.0,0.10, 8.0,0.1,0.01,"<|>");
TArrow *arrow7 = new TArrow(8.0,0.1, 9.0,0.1,0.01,"<|>");
 TArrow *arrow8 = new TArrow(15.0,0.1, 17.0,0.1,0.01,"<|>");

    // arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
    // arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
    // arrow7->Draw(); //arrow8->Draw();
    // TLatex Tl;
    // Tl.SetTextSize(0.055);
    // Tl.SetTextColor(kRed);
    // Tl.DrawLatex(1.0,0.2,"N^{ 0}_{ 2}");
    // Tl.DrawLatex(2.0,0.2,"N^{ 0}_{ 3}");
    // Tl.DrawLatex(3.0,0.2,"N^{ 0}_{ 4}");
    // Tl.DrawLatex(4.0,0.2,"N^{ 0}_{ 5,6}");
    // Tl.DrawLatex(5.0,0.2,"N^{ 0}_{>=7}");
    // Tl.DrawLatex(6.0,0.2,"N^{ #geq1}_{ #geq2}");
    // Tl.DrawLatex(7.0,0.2,"N^{ #geq1}_{ 5,6}");
    // Tl.DrawLatex(8.0,0.2,"N^{ #geq1}_{ #geq7}");

//     TArrow *arrow1 = new TArrow( 0.0,300, 7.0,300,0.01,"<|>");
//     TArrow *arrow2 = new TArrow( 7.0,300,13.0,300,0.01,"<|>");
//     TArrow *arrow3 = new TArrow(13.0,300,19.0,300,0.01,"<|>");
//     TArrow *arrow4 = new TArrow(19.0,300, 25.0,300,0.01,"<|>");
//     TArrow *arrow5 = new TArrow(25.0,300, 31.0,300,0.01,"<|>");
//     TArrow *arrow6 = new TArrow(31.0,300, 38.0,300,0.01,"<|>");

// //     TArrow *arrow1 = new TArrow( 1.0,0.10, 3.0,0.1,0.01,"<|>");
// //     TArrow *arrow2 = new TArrow( 3.0,0.10,5.0,0.1,0.01,"<|>");
// //     TArrow *arrow3 = new TArrow(5.0,0.10,7.0,0.1,0.01,"<|>");
// //     TArrow *arrow4 = new TArrow(7.0,0.10, 9.0,0.1,0.01,"<|>");
// //     TArrow *arrow5 = new TArrow(9.0,0.10, 11.0,0.1,0.01,"<|>");
// //     TArrow *arrow6 = new TArrow(11.0,0.10, 13.0,0.1,0.01,"<|>");
// // TArrow *arrow7 = new TArrow(13.0,0.1, 15.0,0.1,0.01,"<|>");
// //  TArrow *arrow8 = new TArrow(15.0,0.1, 17.0,0.1,0.01,"<|>");

//     arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
//     arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
//     // arrow7->Draw(); arrow8->Draw();
//     TLatex Tl;
//     Tl.SetTextSize(0.055);
    //    Tl.SetTextColor(kWhite);
    // Tl.DrawLatex(3.5,100,"N^{ 0}_{ 2-4}");
    // Tl.DrawLatex(10.5,100,"N^{ 0}_{ 5-6}");
    // Tl.DrawLatex(15.5,100,"N^{ 0}_{ #geq7}");
    // Tl.DrawLatex(21.5,100,"N^{ #geq1}_{ 2-4}");
    // Tl.DrawLatex(26.5,100,"N^{ #geq1}_{ 5-6}");
    // Tl.DrawLatex(33.5,100,"N^{ #geq1}_{ #geq7}");
    // Tl.DrawLatex(1.5,0.1,"N^{ 0}_{ 2}");
    // Tl.DrawLatex(3.5,0.1,"N^{ 0}_{ 3}");
    // Tl.DrawLatex(5.5,0.1,"N^{ 0}_{ 4}");
    // Tl.DrawLatex(7.5,0.1,"N^{ 0}_{ 5,6}");
    // Tl.DrawLatex(9.5,0.1,"N^{ 0}_{>=7}");
    // Tl.DrawLatex(11.5,0.1,"N^{ #geq1}_{ #geq2}");
    // Tl.DrawLatex(13.5,0.1,"N^{ #geq1}_{ 5,6}");
    // Tl.DrawLatex(15.5,0.1,"N^{ #geq1}_{ #geq7}");


  gPad->Modified();
                                                                                       
  hist_ratio->SetLineWidth(1);//line_width[0]+1);
    hist_ratio->SetLineStyle(1);
    hist_ratio->SetLineColor(kBlack);
    hist_ratio->SetTitle(" ");
    hist_ratio->GetXaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetTitle("Data/MC");//TF = #frac{(0#mu,1#gamma)}{(1#mu,1#gamma)}");
    hist_ratio->GetXaxis()->SetLabelSize(0.1);
    hist_ratio->GetYaxis()->SetRangeUser(-0.5,3);
    hist_ratio->GetXaxis()->SetTitle(xlabel);
    // if(DoRebin)
    //   hist_ratio->Rebin(rebin);

    //    hist_ratio->GetXaxis()->SetRangeUser(xmin,xmax+4);
    //hist_ratio =setMyRange(hist_ratio,xmin,xmax+4);
    setLastBinAsOverFlow(hist_ratio);
    //hist_ratio->GetXaxis()->SetRangeUser(xmin,xmax+4);

    //    hist_ratio->GetXaxis()->SetLabelSize(0.0450);
    hist_ratio->GetYaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetLabelSize(0.08);
    hist_ratio->GetYaxis()->SetTitleOffset(.3);
    hist_ratio->GetYaxis()->SetNdivisions(505);
    //    hist_ratio->GetYaxis()->SetLabelSize(x_label_size);
   pad_1->cd();
   pad_1->SetGrid();
   TLine *l =new TLine(xmin,1.0,xmax+4,1.0);
   hist_ratio->SetMarkerStyle(20);
   hist_ratio->Draw("E1");
   l->Draw("sames");
   TLine *l1 =new TLine(xmin,1.5,xmax+4,1.5);
   l1->SetLineStyle(7);
   l1->Draw("sames");
   TLine *l2 =new TLine(xmin,0.5,xmax+4,0.5);
   l2->SetLineStyle(7);

   l2->Draw("sames");

  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
     canvas_n1->SaveAs(canvas_name);   
     sprintf(canvas_name,"%s.pdf",tag_name);
    canvas_n1->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
TFile *f[nfiles];
TFile *f1[nfiles];


void plot_DataValida_kinem(string pathname, int which_Lept, int which_year, int which_TFBins)
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
  char *dataset=new char[200];
  char *dataset1= new char[200];
  char *year =new char[200];
  float energyy[2]={};
  char *string_png = new char[200];
 char *TFbins_str= new char[2000];
  if(which_TFBins==1)
    sprintf(TFbins_str,"TFbins_v1_nJetsBjets");
  else if (which_TFBins==2)
    sprintf(TFbins_str,"TFbins_v2_nJetsBjets_PhoPt");
  else if(which_TFBins==3)
    sprintf(TFbins_str,"TFbins_v3_nJetsBjets_MET");

  cout<<which_TFBins<<"\t"<<which_Lept<<"\t"<<which_year<<endl;
  if(which_Lept==1){
    sprintf(string_png,"Electron_LL_%s",TFbins_str);
    if(which_year==0) {sprintf(dataset,"2016: (0e, 1#gamma)"); sprintf(year,"SummerUL16");energyy[0]=35.922;}    
    if(which_year==1) {sprintf(dataset,"2017: (0e, 1#gamma)");sprintf(year,"Summer20UL17");energyy[0]=41.529;}
    if(which_year==2) {sprintf(dataset,"2018: (0e, 1#gamma)");sprintf(year,"Summer20UL18");energyy[0]=59.74;}
    if(which_year==3){sprintf(dataset,"2016preVFP: (0e, 1#gamma)"); sprintf(year,"SummerUL16preVFP");energyy[0]=19.5;}
    if(which_year==4){sprintf(dataset,"2016postVFP: (0e, 1#gamma)"); sprintf(year,"SummerUL16postVFP");energyy[0]=16.5;}
    if(which_year==5){sprintf(dataset,"Run2: (0e, 1#gamma)"); sprintf(year,"FullRun2");energyy[0]=137.19;}
  }
  else if(which_Lept==3){
    sprintf(string_png,"Lepton_LL_%s",TFbins_str);
    if(which_year==0) {sprintf(dataset,"2016: (0l, 1#gamma)"); sprintf(year,"SummerUL16");energyy[0]=35.922;}
    if(which_year==1) {sprintf(dataset,"2017: (0l, 1#gamma)");sprintf(year,"Summer20UL17");energyy[0]=41.529;}
    if(which_year==2) {sprintf(dataset,"2018: (0l, 1#gamma)");sprintf(year,"Summer20UL18");energyy[0]=59.74;}
    if(which_year==3){sprintf(dataset,"2016preVFP: (0l, 1#gamma)"); sprintf(year,"SummerUL16preVFP");energyy[0]=19.5;}
    if(which_year==4){sprintf(dataset,"2016postVFP: (0l, 1#gamma)"); sprintf(year,"SummerUL16postVFP");energyy[0]=16.5;}
    if(which_year==5){sprintf(dataset,"Run2: (0l, 1#gamma)"); sprintf(year,"FullRun2");energyy[0]=137.19;}

  }
  else if(which_Lept==2)
    {
      sprintf(string_png,"Muon_LL_%s",TFbins_str);
      if(which_year==0) {sprintf(dataset,"2016: (0#mu, 1#gamma)"); sprintf(year,"Summer20UL16");energyy[0]=35.922;}
      if(which_year==1) {sprintf(dataset,"2017: (0#mu, 1#gamma)");sprintf(year,"Summer20UL17");energyy[0]=41.529;}
      if(which_year==2) {sprintf(dataset,"2018: (0#mu, 1#gamma)");sprintf(year,"Summer20UL18");energyy[0]=59.74;}
      if(which_year==3){sprintf(dataset,"2016preVFP: (0#mu, 1#gamma)"); sprintf(year,"SummerUL16preVFP");energyy[0]=19.5;}
      if(which_year==4){sprintf(dataset,"2016postVFP: (0#mu, 1#gamma)"); sprintf(year,"SummerUL16postVFP");energyy[0]=16.5;}
       if(which_year==5){sprintf(dataset,"Run2: (0#mu, 1#gamma)"); sprintf(year,"FullRun2");energyy[0]=137.19;}
      
    }

  if(which_Lept==1)
    {
      if(which_year==2){
       f[0] = new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[1] = new TFile("./Summer20UL18_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[3] = new TFile("./Summer20UL18_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[2] = new TFile("./Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[4]= new TFile("./out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
       f[5]= new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");}
       if(which_year==1)
	 {
	   f[0] = new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[1] = new TFile("./Summer20UL17_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[3] = new TFile("./Summer20UL17_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[2] = new TFile("./Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[4]= new TFile("./out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
       f[5]= new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
	 }

       if(which_year==3){
       f[0] = new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[1] = new TFile("./Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[3] = new TFile("./Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[2] = new TFile("./Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[4]= new TFile("./out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
       f[5]= new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");}
       if(which_year==4)
         {
           f[0] = new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[1] = new TFile("./Summer20UL16_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[3] = new TFile("./Summer20UL16_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[2] = new TFile("./Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
       f[4]= new TFile("./out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
       f[5]= new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
         }

       if(which_year==0)
         {
	   f[0] = new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
	   f[1] = new TFile("./Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
	   f[3] = new TFile("./Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
	   f[2] = new TFile("./Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
	   f[4]= new TFile("./out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
	   f[5]= new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
	 }
       if(which_year==5)
        {
          f[0] = new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");
          f[1] = new TFile("./FullRun2_TTJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
          f[3] = new TFile("./FullRun2_WGJets_PhoIdloose_phopt40_MET200_ElectronLL.root");
          f[2] = new TFile("./FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_ElectronLL.root");
          f[4]= new TFile("./out_Data_FullRun2_Allruns_MET_phoID_loose_pt40_MET200_Electron.root");
          f[5]= new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_ElectronLL.root");

         }


    }
  else if(which_Lept==3)
    {
      if(which_year==0)
         {
	   f[0] = new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
	   f[1] = new TFile("./Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
	   f[3] = new TFile("./Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
	   f[2] = new TFile("./Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
	   f[4]= new TFile("./out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
	   f[5]= new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
	 }
       if(which_year==2)
        {
            f[0] = new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
           f[1] = new TFile("./Summer20UL18_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
           f[3] = new TFile("./Summer20UL18_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[2] = new TFile("./Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[4]= new TFile("./out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
       f[5]= new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");

         }
         if(which_year==1)
        {
	  f[0] = new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
	  f[1] = new TFile("./Summer20UL17_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
	  f[3] = new TFile("./Summer20UL17_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
	  f[2] = new TFile("./Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
	  f[4]= new TFile("./out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
	  f[5]= new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");

         }
	    if(which_year==5)
        {
          f[0] = new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
          f[1] = new TFile("./FullRun2_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
          f[3] = new TFile("./FullRun2_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
          f[2] = new TFile("./FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
          f[4]= new TFile("./out_Data_FullRun2_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
          f[5]= new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");

         }
	    if(which_year==3){
       f[0] = new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[1] = new TFile("./Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[3] = new TFile("./Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[2] = new TFile("./Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[4]= new TFile("./out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
       f[5]= new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");}
       if(which_year==4)
         {
           f[0] = new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[1] = new TFile("./Summer20UL16_TTJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[3] = new TFile("./Summer20UL16_WGJets_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[2] = new TFile("./Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_LeptonLL.root");
       f[4]= new TFile("./out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root");
       f[5]= new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_LeptonLL.root");
         }




    }
  else if (which_Lept==2)
    {
      if(which_year==2)
	{
	    f[0] = new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
           f[1] = new TFile("./Summer20UL18_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
           f[3] = new TFile("./Summer20UL18_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[2] = new TFile("./Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[4]= new TFile("./out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
       f[5]= new TFile("./Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");

	 }
      if(which_year==5)
        {
          f[0] = new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
          f[1] = new TFile("./FullRun2_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
          f[3] = new TFile("./FullRun2_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
          f[2] = new TFile("./FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
          f[4]= new TFile("./out_Data_FullRun2_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
          f[5]= new TFile("./FullRun2_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");

         }

      if(which_year==1)
	{
	  f[0] = new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[1] = new TFile("./Summer20UL17_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[3] = new TFile("./Summer20UL17_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[2] = new TFile("./Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[4]= new TFile("./out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
       f[5]= new TFile("./Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
	}
      if(which_year==0)
	{
	  f[0] = new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[1] = new TFile("./Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[3] = new TFile("./Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[2] = new TFile("./Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[4]= new TFile("./out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
       f[5]= new TFile("./Summer20UL_total2016_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");

	}
      if(which_year==3){
       f[0] = new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[1] = new TFile("./Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[3] = new TFile("./Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[2] = new TFile("./Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[4]= new TFile("./out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
       f[5]= new TFile("./Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");}
       if(which_year==4)
         {
           f[0] = new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[1] = new TFile("./Summer20UL16_TTJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[3] = new TFile("./Summer20UL16_WGJets_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[2] = new TFile("./Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200_MuonLL.root");
       f[4]= new TFile("./out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root");
       f[5]= new TFile("./Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_MET200_MuonLL.root");
         }

    
    }

  vector<string>varName;
  vector<string>varName1;
   vector<string>varName2;
  vector<string>varName3;
  vector<string>varName4;
  vector<string>varName5;

  if(which_TFBins==1){
    varName ={"h_St_validation_Elec_CR","h_HT_validation_Elec_CR","h_NhadJets_validation_Elec_CR","h_NBJets_validation_Elec_CR","h_MET_validation_Elec_CR","h_PhoPt_validation_Elec_CR"};
    varName1 ={"h_St_validation_Mu_CR","h_HT_validation_Mu_CR","h_NhadJets_validation_Mu_CR","h_NBJets_validation_Mu_CR","h_MET_validation_Mu_CR","h_PhoPt_validation_Mu_CR"};
    varName2={"h_St_Elec_SR","h_HT_Elec_SR","h_NhadJets_Elec_SR","h_NBJets_Elec_SR","h_MET_Elec_SR","h_PhoPt_Elec_SR"};
    varName3={"h_St_Mu_SR","h_HT_Mu_SR","h_NhadJets_Mu_SR","h_NBJets_Mu_SR","h_MET_Mu_SR","h_PhoPt_Mu_SR"};
    varName4={"h_St_TauHad_SR","h_HT_TauHad_SR","h_NhadJets_TauHad_SR","h_NBJets_TauHad_SR","h_MET_TauHad_SR","h_PhoPt_TauHad_SR"};
  }

  else if(which_TFBins==2)
    {

      varName ={"h_St_validation_TFbins_v2_Elec_CR","h_HT_validation_TFbins_v2_Elec_CR","h_NhadJets_validation_TFbins_v2_Elec_CR","h_NBJets_validation_TFbins_v2_Elec_CR","h_MET_validation_TFbins_v2_Elec_CR","h_PhoPt_validation_TFbins_v2_Elec_CR"};
      varName1 ={"h_St_validation_TFbins_v2_Mu_CR","h_HT_validation_TFbins_v2_Mu_CR","h_NhadJets_validation_TFbins_v2_Mu_CR","h_NBJets_validation_TFbins_v2_Mu_CR","h_MET_validation_TFbins_v2_Mu_CR","h_PhoPt_validation_TFbins_v2_Mu_CR"};
      varName2={"h_St_Elec_SR","h_HT_Elec_SR","h_NhadJets_Elec_SR","h_NBJets_Elec_SR","h_MET_Elec_SR","h_PhoPt_Elec_SR"};
      varName3={"h_St_Mu_SR","h_HT_Mu_SR","h_NhadJets_Mu_SR","h_NBJets_Mu_SR","h_MET_Mu_SR","h_PhoPt_Mu_SR"};
      varName4={"h_St_TauHad_SR","h_HT_TauHad_SR","h_NhadJets_TauHad_SR","h_NBJets_TauHad_SR","h_MET_TauHad_SR","h_PhoPt_TauHad_SR"};

    }
  else if(which_TFBins==3){
    varName ={"h_St_validation_TFbins_v3_Elec_CR","h_HT_validation_TFbins_v3_Elec_CR","h_NhadJets_validation_TFbins_v3_Elec_CR","h_NBJets_validation_TFbins_v3_Elec_CR","h_MET_validation_TFbins_v3_Elec_CR","h_PhoPt_validation_TFbins_v3_Elec_CR"};
    varName1 ={"h_St_validation_TFbins_v3_Mu_CR","h_HT_validation_TFbins_v3_Mu_CR","h_NhadJets_validation_TFbins_v3_Mu_CR","h_NBJets_validation_TFbins_v3_Mu_CR","h_MET_validation_TFbins_v3_Mu_CR","h_PhoPt_validation_TFbins_v3_Mu_CR"};
    varName2={"h_St_Elec_SR","h_HT_Elec_SR","h_NhadJets_Elec_SR","h_NBJets_Elec_SR","h_MET_Elec_SR","h_PhoPt_Elec_SR"};
    varName3={"h_St_Mu_SR","h_HT_Mu_SR","h_NhadJets_Mu_SR","h_NBJets_Mu_SR","h_MET_Mu_SR","h_PhoPt_Mu_SR"};
    varName4={"h_St_TauHad_SR","h_HT_TauHad_SR","h_NhadJets_TauHad_SR","h_NBJets_TauHad_SR","h_MET_TauHad_SR","h_PhoPt_TauHad_SR"};

  }
  vector <string>  xLabel;
  xLabel={"Sum of P_{T}^{Jets} & P_{T}^{#gamma} [GeV]","HT[GeV]","N_{jets}","N_{ b-jets}","p_{T}^{miss} [GeV]","p_{T}^{#gamma} [GeV]"};//,"BDT response","M_{T}^{miss & #gamma} [GeV]","dPhi(#gamma,MET)"};
  vector <int> rebin;
  rebin={4,4,1,1,8,8,8,4,4};
  vector<double> ymin ={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  vector<double> ymax={100000,100000,100000,100000,100000,100000,100000,100000,100000};
  vector<double> xmin ={300,300,2,0,100,20,0,0,0};
  vector<double> xmax={1500,1500,16,12,1000,400,1,700,4};
     
  vector<string>baseline1;
  vector<string> baseline = {"Nocut", "PreSmuontion","Electron_CR","Electron_SR","FailAcep_ElectronSR","FailId_ElectronSR","FailIso_ElectronSR","Mu_CR","Mu_SR","FailAcep_MuSR","FailId_MuSR","FailIso_MuSR","Electron_CR_BDTcut1","Electron_SR_BDTcut1","FailAcep_ElectronSR_BDTcut1","FailId_ElectronSR_BDTcut1","FailIso_ElectronSR_BDTcut1","Mu_CR_BDTcut1","Mu_SR_BDTcut1","FailAcep_MuSR_BDTcut1","FailId_MuSR_BDTcut1","FailIso_MuSR_BDTcut1","Electron_CR_BDTcut2","Electron_SR_BDTcut2","FailAcep_ElectronSR_BDTcut2","FailId_ElectronSR_BDTcut2","FailIso_ElectronSR_BDTcut2","Mu_CR_BDTcut2","Mu_SR_BDTcut2","FailAcep_MuSR_BDTcut2","FailId_MuSR_BDTcut2","FailIso_MuSR_BDTcut2"};

  //     vector<string> baseline1 = {"FailId_ElecSR","FailIso_ElecSR","FailAcep_ElecSR"};
  // if(which_Lept)
  //   baseline1={"Elec_CR","Mu_CR","Elec_SR","Elec_SR","Mu_SR","FailAcep_ElecSR","FailId_ElecSR","FailIso_ElecSR","Elec_SR","Elec_"};//
  // else
  //   baseline1 = {"Mu_CR","Mu_SR","FailAcep_MuSR","FailId_MuSR","FailIso_MuSR","FailAcep_MuSR","Mu_SR","Mu_"};
  //   const char *baseline1[3]={"Nocut","Mu_SR","Mu_CR"};
      //   const char *baseline1[9]={"Nocut","SignalRegion","lostElec_SR","lostMu_SR","lostTau_SR","lostElec_SR_iso","lostElec_SR_Accept","lostElec_SR_ident"};
  //  const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
  //  float energyy[8]={59.74,41.529,35.922,137.19,59.74,41.529,35.922,137.19};
  bool flag=false;
  const char* filetag[4]={"2018","2017","2016","Run2"};
  //  const char* filetag[10]={"TTGJets","pMSSM_MCMC_70_90438","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451","WJets","GJets","T5bbbbZG_10","T5bbbbZG_50","T5bbbbZG_200","T5bbbbZG_1500"};

  //  vector<TH1D*> hist_list_Njets;
  // vector<TH1D*> hist_list_Bjets;
  // vector<TH1D*> hist_list_MET;
  // vector<TH1D*> hist_list_PhoPt;
  // vector<TH1D*> hist_list_Mt;
  // vector<TH1D*> hist_list_ST;
  // vector<TH1D*> hist_list_HT;
  // vector<TH1D*> hist_list_dPhiPhoMET;
  // vector<TH1D*> hist_list_BDTresponse;
  for(int i_cut=0; i_cut<varName.size();i_cut++)
    {
      if(i_cut==1) continue;
      vector<TH1D*> hist_list_Njets;
      sprintf(hist_name,"%s",varName[i_cut].c_str());
      sprintf(hist_name1,"%s",varName1[i_cut].c_str());
      sprintf(hist_name2,"%s",varName2[i_cut].c_str());
      sprintf(hist_name3,"%s",varName3[i_cut].c_str());
      sprintf(hist_name4,"%s",varName4[i_cut].c_str());
      //vector<TH1D*> hist_list_Njets;
      for(int i_file=0; i_file<6;i_file++)
	{
	  cout<<hist_name<<"\t"<<i_cut<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
	  TH1D* h_resp = (TH1D*)f[i_file]->Get(hist_name);
          TH1D* h_resp1 = (TH1D*)f[i_file]->Get(hist_name1);
          TH1D* h_resp2 = (TH1D*)f[i_file]->Get(hist_name2);
          TH1D* h_resp3 = (TH1D*)f[i_file]->Get(hist_name3);
          TH1D* h_resp4 = (TH1D*)f[i_file]->Get(hist_name4);
          cout<<"resp "<<h_resp->Integral()<<"\t"<<"resp1 "<<h_resp1->Integral()<<"\t"<<"resp2 "<<h_resp2->Integral()<<"\t"<<"resp3 "<<h_resp3->Integral()<<"\t"<<"resp4 "<<h_resp4->Integral()<<"\t"<<endl;
	  if(i_file!=4){
            if(which_Lept==1){
              h_resp2->Rebin(rebin[i_cut]);
              hist_list_Njets.push_back(h_resp2);
            }
            else if(which_Lept==2){
              h_resp4->Add(h_resp3);
              h_resp4->Rebin(rebin[i_cut]);
              hist_list_Njets.push_back(h_resp4);
            }
            else if(which_Lept==3){
              h_resp2->Add(h_resp3);
              h_resp2->Add(h_resp4);
              h_resp2->Rebin(rebin[i_cut]);
              hist_list_Njets.push_back(h_resp2);
            }
          }
          else
            {
              if(which_Lept==1){
                h_resp->Rebin(rebin[i_cut]);
                hist_list_Njets.push_back(h_resp);
              }
              else if(which_Lept==2){
                h_resp1->Rebin(rebin[i_cut]);
                hist_list_Njets.push_back(h_resp1);
              }
              else if(which_Lept==3){
                h_resp->Add(h_resp1);
                h_resp->Rebin(rebin[i_cut]);
                hist_list_Njets.push_back(h_resp);
              }
            }
	}
      //path to save the png file
      float energy=energyy[0];
      int xrange=0.0;
      TH1D* hNjets_total =(TH1D*)hist_list_Njets.at(0)->Clone();
      hNjets_total->Add(hist_list_Njets.at(1));
      hNjets_total->Add(hist_list_Njets.at(2));
      hNjets_total->Add(hist_list_Njets.at(3));
      TH1D* hNjets_ratio = (TH1D*)hist_list_Njets.at(4)->Clone();
      hNjets_ratio->Divide(hNjets_total);
      //hNjets_ratio->Rebin(rebin[i_cut]);
      sprintf(full_path,"%s/%s_%s_%s_DataMC_PredVsExp_compare",pathname.c_str(),string_png,year,varName[i_cut].c_str());
      generate_1Dplot(hist_list_Njets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Entries",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,false,true,false,true,dataset);
      
      
    }
      

}






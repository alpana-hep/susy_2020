const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[5] ={"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR","Failed Iso"};//,"Failed acceptance","1e CR"};// {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[4] = {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};                                                                               
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};//{9,kCyan+2,45,kMagenta,kGray+1,kRed,kBlue+2,kMagenta,kCyan};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
//int line_color[9] = {kMagenta+2, kGray+2, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
void decorate(TH1D*,int,int );
 
void decorate(TH1D* hist,int i, int j){
  //  hist->SetLineColor(col[i]);
  // hist->SetFillColor(col[i]);
  // vector<int> col;
  // vector<int> Style;
  // if(j==1){
  //   col={kViolet,kGray+1,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
  //   Style= {3008,1001,3019,3244};
    
  //  }
  // else if(j==2){
  //    col={kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
  //    Style={1001,3008,1001,3244};
  // }
  // else if(j==3){
  //   col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
  //   Style={3008,1001,3008,1001};
  // }
  //   if(i!=4)
  //     {
  //       hist->SetFillColor(col[i]);

  //       hist->SetFillStyle(Style[i]);
  //     }
  //   else 
  //     {
  //       hist->SetFillColor(kGray+1);

  //       hist->SetFillStyle(1001);
  //       }
    hist->SetLineWidth(3);

  //  if(i<nBG) {                                                                                                                                 
  //hist->SetFillColor(col[i]);
    // if(i!=0)
    //   {
    // 	//	hist->SetFillColor(col[i]);

    //  	hist->SetFillStyle(Style[i]);
    // //   }
    // // else
    // //   {
    // // 	hist->SetFillColor(kGray+2);

    // // 	hist->SetFillStyle(1001);
    // // 	}
    // hist->SetLineWidth(2);
    //}                                                                                                                                           
}

void setLastBinAsOverFlow(TH1D*);
TH1D* setMyRange(TH1D*,double,double);
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
  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
   return h1;

}

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
  cout<<lastBinCt<<"\t"<<"Last bin values"<<endl;

}

// TH1D* setLastBinAsOverFlow(TH1D* h_hist, int xrange){
//   //     h_hist = setMyRange(h_hist,0,xrange);
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


// TH1D* setMyRange(TH1D *h1,double xLow,double xHigh){
//   //call it after setting last bin as overflow                                                                                                    
//   double err=0;
//   if(xHigh > 13000) return h1;
//   if(xLow < -13000) return h1;

//   // h1->Print("all");
//   //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);  
//   int nMax=h1->FindBin(xHigh);
//   h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
//   h1->SetBinError(nMax,err);

//   //  cout<<nMax<<endl;
//   //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
//   for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
//     h1->SetBinContent(i,0);
//     h1->SetBinError(i,0);
//     //    cout<<":";
//     //h1->GetXaxis()->SetRangeUser(xLow,xHigh); 
//   }
//   return h1;
// }

void generate_1Dplot(vector<TH1D*> hist, TH1D* hist_ratio, char const *tag_name="",char const *xlabel="",char const *ylabel="", float energy=-1, int rebin=-1,double ymin=0,double ymax=0,int xmin=-1,int xmax=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"}, int which_TFbins=-1, int which_Lept=-1){  
     TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,950,850);
       canvas_n1->Range(-60.25,-0.625,562.25,0.625);
       canvas_n1->SetFillColor(0);
       canvas_n1->SetBorderMode(0);
       canvas_n1->SetBorderSize(2);
       auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
       pad_1->SetTopMargin(0.04);
       pad_1->SetBottomMargin(0.33);
       pad_1->SetRightMargin(0.035);
       pad_1->SetLeftMargin(0.13);
       auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
       p1->SetBottomMargin(0.026);
       p1->SetRightMargin(0.035);
       p1->SetLeftMargin(0.13);
       p1->SetTopMargin(0.1);
       p1->cd();
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
  gROOT->ForceStyle();
  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.2,0.75,0.65,0.88);  
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(4);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.045;
  //  double ymin = 100000.0;
  //double ymax = 0.0;
  double xrange = xmax;
  // float energy = energyy;
  vector<TH1D*> hist_list_temp;
  cout<<" hist.size() = "<<hist.size()<<endl;
  // if(hist.size()==2)
  //   {  col[1]=col[1]-41;
  //     legend_text[1] = "(0e,0#mu) SR";
  //   }
  // else
  //   {
  //     col[1]=45;
  //     legend_text[1] = "Failed acceptance";
  //   }
  for(int i =0;i<(int)hist.size(); i ++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    //    hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    //    normalize = true;
    if(normalize) {
      //hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
    //   hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);
    hist.at(i)->SetLineWidth(line_width[i]);
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    //    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleSize(0.06);
    hist.at(i)->GetYaxis()->SetLabelSize(0.06);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.);
     decorate(hist.at(i),i, which_Lept);
    hist.at(i)->SetMarkerSize(0.8);
    hist.at(i)->SetMarkerStyle(20);
    hist.at(i)->SetMarkerColor(line_color[i]);
    //new ones
    // hist.at(i)->GetXaxis()->SetTitleSize(0.08);
    // hist.at(i)->GetXaxis()->SetLabelSize(0.06);

    // hist.at(i)->GetYaxis()->SetTitleSize(0.07);
    //hist.at(i)->GetYaxis()->SetLabelSize(0.06);

    // hist.at(i)->GetXaxis()->SetTitleOffset(3);
 hist.at(i)->GetXaxis()->SetLabelOffset(1.6);

    // hist.at(i)->GetYaxis()->SetTitleOffset(0.9);

    // //decorate(hist.at(i),i, which_Lept);
    // hist.at(i)->GetYaxis()->SetNdivisions(506);
    // hist.at(i)->GetXaxis()->SetTitle(title);
    
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
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"e2p");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    // hist.at(i)= setMyRange(hist.at(i),xmin,xmax+4);
    //setLastBinAsOverFlow(hist.at(i));
    //    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);

    

  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0;i<(int)hist.size(); i++) {
    if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.0001,10.0*ymax);
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.00001,5.0);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
      }
    //    p1->SetGrid();
    //hs_var->Add(hist.at(i));
    // hs_var->SetMinimum(0.00001);
    // hs_var->SetMaximum(ymax*60);
    //gPad->SetLogu
    cout<<"i Alps "<<i<<endl;
    if(i==0) hist.at(i)->Draw("hist ");
    else hist.at(i)->Draw("hist sames");
	
  }
//  hs_var->SetMinimum(0.0);
//  hs_var->SetMaximum(ymax+0.5);


//   hs_var->Draw("BAR HIST");
//   hs_var->Draw("HIST");
//   hs_var->GetXaxis()->SetRangeUser(0,xrange);
//  hs_var->GetXaxis()->SetTitle(title);
// hs_var->GetXaxis()->SetTitleOffset(1.0);
//   gPad->Modified(); gPad->Update();
//   hs_var->GetXaxis()->SetTitle(title);
//   //  hs_var->GetYaxis()->SetTitleSize(
//   hs_var->GetYaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");//hs_var->GetYaxis()->SetTitle("Events");
//   hs_var->SetTitle(0);
//    hs_var->GetYaxis()->SetTitleOffset(1.2);
// hs_var->GetXaxis()->SetTitleSize(00.05);
//     hs_var->GetXaxis()->SetLabelSize(0.04);
//     hs_var->GetYaxis()->SetLabelSize(0.04);
//     hs_var->GetYaxis()->SetTitleSize(00.055);
//     hs_var->GetYaxis()->SetTitleOffset(1.0);

  // hs_var->SetMinimum(0.0);
  // hs_var->SetMaximum(1.5);
  
  
  // hs_var->Draw("BAR HIST");
  // hs_var->Draw("HIST");
  // if(which_TFbins==1) //default 8 bins                                                                                                                                     
  //   hs_var->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                    
  // else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100                                                                                                
  //   hs_var->GetXaxis()->SetRangeUser(0,18);
  // else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300                                                                                                     
  //   hs_var->GetXaxis()->SetRangeUser(0,18);
  // hs_var->GetXaxis()->SetTitle(title);
  // hs_var->GetXaxis()->SetTitleOffset(1.0);
  // gPad->Modified(); gPad->Update();
  // hs_var->GetXaxis()->SetTitle(title);
  // //  hs_var->GetYaxis()->SetTitleSize(
  // hs_var->GetYaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");//hs_var->GetYaxis()->SetTitle("Events");
  // hs_var->SetTitle(0);
  // hs_var->GetYaxis()->SetTitleOffset(1.2);
  // hs_var->GetXaxis()->SetTitleSize(00.05);
  // hs_var->GetXaxis()->SetLabelSize(0.04);
  // hs_var->GetYaxis()->SetLabelSize(0.04);
  // hs_var->GetYaxis()->SetTitleSize(00.055);
  // hs_var->GetYaxis()->SetTitleOffset(1.0);


  legend->Draw();

    //    legend->Draw();


  if(log_flag) {
      gPad->SetLogy();
    }
  // if(logx)
  //   gPad->SetLogx();
  gPad->Update(); 
  TLatex* textOnTop = new TLatex();
  //new
    textOnTop->SetTextSize(0.054);
  textOnTop->DrawLatexNDC(0.146,0.925,"CMS #it{#bf{Simulation Preliminary}}");

  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.054);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.72,0.925,en_lat);


  gPad->Modified();
                                                                                       
    hist_ratio->SetLineWidth(2);
    hist_ratio->SetLineStyle(1);
    hist_ratio->SetMarkerSize(0.2);
    hist_ratio->SetLineColor(kBlack);
    hist_ratio->SetTitle(" ");
    hist_ratio->GetXaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetTitle("Z(#nu#nu)/Z(ll)");//TF = #frac{N_{SR}}{N_{CR}}");//(0#mu,1#gamma)}{(1#mu,1#gamma)}");
    //    hist_ratio->GetXaxis()->SetLabelSize(0.1);
    hist_ratio->GetYaxis()->SetRangeUser(-0.001,0.03);
    //hist_ratio->GetXaxis()->SetRangeUser(xmin,xmax+4);
    // hist_ratio= setMyRange(hist_ratio,xmin,xmax+6);
    //setLastBinAsOverFlow(hist_ratio);

    // if(which_TFbins==1) //default 8 bins                                                                                                    
    //   hist_ratio->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                 
    // else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100
    //   hist_ratio->GetXaxis()->SetRangeUser(0,18);
    // else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300                                                                                        
    //    hist_ratio->GetXaxis()->SetRangeUser(0,39);
    
    hist_ratio->GetXaxis()->SetLabelSize(0.0450);
    hist_ratio->GetYaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetLabelSize(0.08);
    hist_ratio->GetYaxis()->SetTitleOffset(.4);
    hist_ratio->SetMarkerSize(1.0);
    hist_ratio->SetMarkerStyle(20);
    hist_ratio->SetMarkerColor(kBlue);
    hist_ratio->GetXaxis()->SetTitle(xlabel);
    hist_ratio->GetYaxis()->SetNdivisions(505);
    //new
    hist_ratio->GetXaxis()->SetTitleSize(0.05);
    hist_ratio->GetXaxis()->SetLabelSize(0.11);
    hist_ratio->GetYaxis()->SetTitleSize(0.125);
    hist_ratio->GetYaxis()->SetNdivisions(505);

    hist_ratio->GetXaxis()->SetTitleOffset(1);
    hist_ratio->GetYaxis()->SetTitleOffset(0.41);
    hist_ratio->GetXaxis()->SetTitleSize(0.14);

    hist_ratio->GetYaxis()->SetLabelSize(0.11);
    
    if(normalize){
      hist_ratio->GetYaxis()->SetRangeUser(0,2);
      hist_ratio->GetXaxis()->SetLabelOffset(0);
      //      hist_ratio->GetXaxis()->SetLabelSize(0.4);
    }
    else {
      hist_ratio->GetYaxis()->SetRangeUser(0,10);
    }
    // gStyle->SetLabelOffset(1.2);
    // gStyle->SetLabelSize(1.2);
     //     hist_ratio->GetXaxis()->SetLabelOffset();
     //hist_ratio->GetYaxis()->SetLabelSize(x_label_size);
   pad_1->cd();
   //   pad_1->SetGrid();
   // if(which_TFbins==1){
   // TLine *l =new TLine(0,1.0,10,1.0);   
   // hist_ratio->Draw("");
   // l->Draw("sames");
   // TLine *l1 =new TLine(0,1.5,10,1.5);
   // l1->SetLineStyle(7);
   // l1->Draw("sames");
   // TLine *l2 =new TLine(0,0.5,10,0.5);
   // l2->SetLineStyle(7);

   // l2->Draw("sames");
   // }

   // else{
     
   //    TLine *l =new TLine(0,1.0,18,1.0);
   // hist_ratio->Draw("");
   // l->Draw("sames");
   // TLine *l1 =new TLine(0,1.5,18,1.5);
   // l1->SetLineStyle(7);
   // l1->Draw("sames");
   // TLine *l2 =new TLine(0,0.5,18,0.5);
   // l2->SetLineStyle(7);

   // l2->Draw("sames");
   // }
   // TLine *l =new TLine(xmin,.01,xrange,.01);
   // hist_ratio->Draw("");
   // l->Draw("sames");
   // TLine *l1 =new TLine(xmin,0.02,xrange,0.02);
   // l1->SetLineStyle(7);
   // l1->Draw("sames");
   // TLine *l2 =new TLine(xmin,0.03,xrange,0.03);
   // l2->SetLineStyle(7);

   // l2->Draw("sames");

     if(normalize){
       //       hist_ratio->Scale(1.0/hist_ratio->Integral());
       TLine *l =new TLine(xmin,1,xrange,1);
       hist_ratio->Draw("");
       l->Draw("sames");
       TLine *l1 =new TLine(xmin,1.5,xrange,1.5);
       l1->SetLineStyle(7);
       l1->Draw("sames");
       TLine *l2 =new TLine(xmin,0.5,xrange,0.5);
       l2->SetLineStyle(7);
       l2->Draw("sames");
   }

   else{

      TLine *l =new TLine(xmin,5,xrange,5);
      hist_ratio->Draw("");
       l->Draw("sames");
       TLine *l1 =new TLine(xmin,10,xrange,10);
       l1->SetLineStyle(7);
       l1->Draw("sames");
       TLine *l2 =new TLine(xmin,15,xrange,15);
       l2->SetLineStyle(7);
      l2->Draw("sames");

   //   TLine *l =new TLine(0,10,50,10);
   // hist_ratio->Draw("");
   // l->Draw("sames");
   // TLine *l1 =new TLine(0,12,50,12);
   // l1->SetLineStyle(7);
   // l1->Draw("sames");
   // TLine *l2 =new TLine(0,5,50,5);
   // l2->SetLineStyle(7);

   // l2->Draw("sames");
   }
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
     canvas_n1->SaveAs(canvas_name);   
     sprintf(canvas_name,"%s.pdf",tag_name);
    canvas_n1->SaveAs(canvas_name);
  // sprintf(canvas_name,"%s.root",tag_name);
  //   canvas_n1->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
TFile *f[nfiles];
TFile *f1[nfiles];


void TF_varKinem_ratio(string pathname, int which_plots)
{
  char* hname = new char[200];
  char* hname1 = new char[200];
  char* hname2 = new char[200];
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
  int n=0;
  char *dataset=new char[200];
  char *year =new char[200];
  //  float energyy[2]={};
  int n_files=42;
  char *string_png = new char[200];
  vector<string>baseline1;
  vector<string> legend_texts;
  char *TFbins_str= new char[2000];
  int which_Lept=1;
  int which_TFBins=1;

  // if(which_TFBins==1)
  //   sprintf(TFbins_str,"TFbins_v1_nJetsBjets");
  // else if (which_TFBins==2)
  //   sprintf(TFbins_str,"TFbins_v2_nJetsBjets_PhoPt");
  // else if(which_TFBins==3)
  sprintf(TFbins_str,"Default");
  //if(which_Lept==1){
    //    sprintf(string_png,"Electron_LL");
    baseline1={"Pho_SR","Elec_CR"};//,"TauHad_SR","Mu_SR","Elec_SR","FailAcep_ElecSR","FailId_ElecSR","FailIso_ElecSR","Elec_SR","Elec_"};//
    
    //legend_texts ={"(0e,1#gamma) SR","(1e,0#gamma) CR"};//,"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR","Failed Iso"};
    legend_texts ={"Z -> #nu#nu + #gamma","Z -> ll + #gamma"};
    //sprintf(hname,"");
    sprintf(string_png,"Zinv_CRvsSR_%s",TFbins_str);
    sprintf(hname,"%s_phoID_loose_29Jan24",string_png);
    f[0] = new TFile("Summer20UL18_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");
    f[1] = new TFile("Summer20UL17_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");
    f[2] = new TFile("Summer20UL16_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");
    f[3] = new TFile("Summer20UL16APV_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");
    f[4] = new TFile("Summer20UL_total2016_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");
    f[5] = new TFile("FullRun2_totalZtoNuNu_PhoIdloose_phopt40_MET200.root");

    f1[0] = new TFile("Summer20UL18_totalZtoLL_PhoIdloose_phopt40_MET200.root");
    f1[1] = new TFile("Summer20UL17_totalZtoLL_PhoIdloose_phopt40_MET200.root");
    f1[2] = new TFile("Summer20UL16_totalZtoLL_PhoIdloose_phopt40_MET200.root");
    f1[3] = new TFile("Summer20UL16APV_totalZtoLL_PhoIdloose_phopt40_MET200.root");
    f1[4] = new TFile("Summer20UL_total2016_totalZtoLL_PhoIdloose_phopt40_MET200.root");
    f1[5] = new TFile("FullRun2_totalZtoLL_PhoIdloose_phopt40_MET200.root");

    vector<string> filetag=  {"2018","2017","2016postVFP","2016preVFP","2016","FullRun2"};
    vector<float> energyy={59.74,41.53,16.5,19.5,36,137.19};

    // }
  // if(which_Lept==2){
  //   sprintf(string_png,"Muon_LL_Sbins_Valid_%s",TFbins_str);
  //   sprintf(hname,"%s_phoID_loose_09Jan24",string_png);
  //   baseline1 = {"Mu_SR","TauHad_SR","Mu_CR"};
  //   //  legend_texts ={"(0#mu,1#gamma) SR","(#tau-had, 1#gamma) SR","(1#mu,1#gamma) CR"};//,"(#tau-had, 1#gamma) SR"};
  //    legend_texts ={"(0#mu,1#gamma) SR","(1#mu,1#gamma) CR"};
  // }
  // if(which_Lept==3){
  //   sprintf(string_png,"Lepton_LL_Sbins_Valid_%s",TFbins_str);
  //   sprintf(hname,"%s_phoID_loose_09Jan24",string_png);
  //   if(which_TFBins==1)
  //     baseline1 = {"Validation_Elec_CR","Validation_Mu_CR","TauHad_SR","Mu_SR","Elec_SR"};
  //   else if(which_TFBins==2)
  //     baseline1 = {"Validation_TFbins_V2_Elec_CR","Validation_TFbins_V2_Mu_CR","TauHad_SR","Mu_SR","Elec_SR"};
  //   else if(which_TFBins==3)
  //     baseline1 = {"Validation_TFbins_V3_Elec_CR","Validation_TFbins_V3_Mu_CR","TauHad_SR","Mu_SR","Elec_SR"};
  //   legend_texts = {"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR"};
  //   legend_texts ={"(0l,1#gamma) SR","(1l,1#gamma) CR"};
  //   n_files = 18;
  // }
  
  // cout<<string_png<<"\t"<<TFbins_str<<"\t"<<which_TFBins<<endl;
  // if(which_Lept==1)
  //   {

  vector<string>varName;
  vector<string>varName1;
  vector<string>varName2;
  vector<string>varName3;
  vector<string>varName4;
  vector<string>varName5;

    if(which_plots==0){
      varName ={"h_St_Elec_CR","h_HT_Elec_CR","h_NhadJets_Elec_CR","h_NBJets_Elec_CR","h_MET_Elec_CR","h_PhoPt_Elec_CR","h_qmulti_Elec_CR","h_Photon_Eta_Elec_CR","h_Photon_Phi_Elec_CR","h_MET_Phi_Elec_CR","h_leadJets_qmulti_Elec_CR","h_leadJet_Pt_Elec_CR","h_leadbjet_tag_Elec_CR","h_nvrtx_Elec_CR","h_minDR_Jets_EMObject_Elec_CR","h_Phi_leadJet1_Elec_CR","h_Eta_leadJet1_Elec_CR","h_Pt_leadJet1_Elec_CR","h_dPhi_METJet1_Elec_CR","h_Phi_leadJet2_Elec_CR","h_Eta_leadJet2_Elec_CR","h_Pt_leadJet2_Elec_CR","h_dPhi_METJet2_Elec_CR","h_Phi_leadJet3_Elec_CR","h_Eta_leadJet3_Elec_CR","h_Pt_leadJet3_Elec_CR","h_dPhi_METJet3_Elec_CR","h_Phi_leadJet4_Elec_CR","h_Eta_leadJet4_Elec_CR","h_Pt_leadJet4_Elec_CR","h_dPhi_METJet4_Elec_CR","h_Phi_matchedJet_Elec_CR","h_Eta_matchedJet_Elec_CR","h_Pt_matchedJet_Elec_CR","h_HT5HT_Elec_CR","h_invariantMass_Elec_CR","h_ZpT_Elec_CR"};                        

      varName2={"h_St_Pho_SR","h_HT_Pho_SR","h_NhadJets_Pho_SR","h_NBJets_Pho_SR","h_MET_Pho_SR","h_PhoPt_Pho_SR","h_qmulti_Pho_SR","h_Photon_Eta_Pho_SR","h_Photon_Phi_Pho_SR","h_MET_Phi_Pho_SR","h_leadJets_qmulti_Pho_SR","h_leadJet_Pt_Pho_SR","h_leadbjet_tag_Pho_SR","h_nvrtx_Pho_SR","h_minDR_Jets_EMObject_Pho_SR","h_Phi_leadJet1_Pho_SR","h_Eta_leadJet1_Pho_SR","h_Pt_leadJet1_Pho_SR","h_dPhi_METJet1_Pho_SR","h_Phi_leadJet2_Pho_SR","h_Eta_leadJet2_Pho_SR","h_Pt_leadJet2_Pho_SR","h_dPhi_METJet2_Pho_SR","h_Phi_leadJet3_Pho_SR","h_Eta_leadJet3_Pho_SR","h_Pt_leadJet3_Pho_SR","h_dPhi_METJet3_Pho_SR","h_Phi_leadJet4_Pho_SR","h_Eta_leadJet4_Pho_SR","h_Pt_leadJet4_Pho_SR","h_dPhi_METJet4_Pho_SR","h_Phi_matchedJet_Pho_SR","h_Eta_matchedJet_Pho_SR","h_Pt_matchedJet_Pho_SR","h_HT5HT_Pho_SR", "h_invariantMass_Pho_SR","h_ZpT_Pho_SR"};

    }
    else if (which_plots==1){
      varName ={"h_St_L1Trig_Elec_CR","h_HT_L1Trig_Elec_CR","h_NhadJets_L1Trig_Elec_CR","h_NBJets_L1Trig_Elec_CR","h_MET_L1Trig_Elec_CR","h_PhoPt_L1Trig_Elec_CR","h_qmulti_L1Trig_Elec_CR","h_Photon_Eta_L1Trig_Elec_CR","h_Photon_Phi_L1Trig_Elec_CR","h_MET_Phi_L1Trig_Elec_CR","h_leadJets_qmulti_L1Trig_Elec_CR","h_leadJet_Pt_L1Trig_Elec_CR","h_leadbjet_tag_L1Trig_Elec_CR","h_nvrtx_L1Trig_Elec_CR","h_minDR_Jets_EMObject_L1Trig_Elec_CR","h_Phi_leadJet1_L1Trig_Elec_CR","h_Eta_leadJet1_L1Trig_Elec_CR","h_Pt_leadJet1_L1Trig_Elec_CR","h_dPhi_METJet1_L1Trig_Elec_CR","h_Phi_leadJet2_L1Trig_Elec_CR","h_Eta_leadJet2_L1Trig_Elec_CR","h_Pt_leadJet2_L1Trig_Elec_CR","h_dPhi_METJet2_L1Trig_Elec_CR","h_Phi_leadJet3_L1Trig_Elec_CR","h_Eta_leadJet3_L1Trig_Elec_CR","h_Pt_leadJet3_L1Trig_Elec_CR","h_dPhi_METJet3_L1Trig_Elec_CR","h_Phi_leadJet4_L1Trig_Elec_CR","h_Eta_leadJet4_L1Trig_Elec_CR","h_Pt_leadJet4_L1Trig_Elec_CR","h_dPhi_METJet4_L1Trig_Elec_CR","h_Phi_matchedJet_L1Trig_Elec_CR","h_Eta_matchedJet_L1Trig_Elec_CR","h_Pt_matchedJet_L1Trig_Elec_CR","h_HT5HT_L1Trig_Elec_CR","h_Mt_phoMET_L1Trig_Elec_CR","h_dPhi_phoMet_L1Trig_Elec_CR"};

  varName2={"h_St_L1Trig_Pho_SR","h_HT_L1Trig_Pho_SR","h_NhadJets_L1Trig_Pho_SR","h_NBJets_L1Trig_Pho_SR","h_MET_L1Trig_Pho_SR","h_PhoPt_L1Trig_Pho_SR","h_qmulti_L1Trig_Pho_SR","h_Photon_Eta_L1Trig_Pho_SR","h_Photon_Phi_L1Trig_Pho_SR","h_MET_Phi_L1Trig_Pho_SR","h_leadJets_qmulti_L1Trig_Pho_SR","h_leadJet_Pt_L1Trig_Pho_SR","h_leadbjet_tag_L1Trig_Pho_SR","h_nvrtx_L1Trig_Pho_SR","h_minDR_Jets_EMObject_L1Trig_Pho_SR","h_Phi_leadJet1_L1Trig_Pho_SR","h_Eta_leadJet1_L1Trig_Pho_SR","h_Pt_leadJet1_L1Trig_Pho_SR","h_dPhi_METJet1_L1Trig_Pho_SR","h_Phi_leadJet2_L1Trig_Pho_SR","h_Eta_leadJet2_L1Trig_Pho_SR","h_Pt_leadJet2_L1Trig_Pho_SR","h_dPhi_METJet2_L1Trig_Pho_SR","h_Phi_leadJet3_L1Trig_Pho_SR","h_Eta_leadJet3_L1Trig_Pho_SR","h_Pt_leadJet3_L1Trig_Pho_SR","h_dPhi_METJet3_L1Trig_Pho_SR","h_Phi_leadJet4_L1Trig_Pho_SR","h_Eta_leadJet4_L1Trig_Pho_SR","h_Pt_leadJet4_L1Trig_Pho_SR","h_dPhi_METJet4_L1Trig_Pho_SR","h_Phi_matchedJet_L1Trig_Pho_SR","h_Eta_matchedJet_L1Trig_Pho_SR","h_Pt_matchedJet_L1Trig_Pho_SR","h_HT5HT_L1Trig_Pho_SR","h_Mt_phoMET_L1Trig_Pho_SR","h_dPhi_phoMet_L1Trig_Pho_SR"};
    }

    else if (which_plots==2){
      varName ={"h_St_ProbL1Trig_Elec_CR","h_HT_ProbL1Trig_Elec_CR","h_NhadJets_ProbL1Trig_Elec_CR","h_NBJets_ProbL1Trig_Elec_CR","h_MET_ProbL1Trig_Elec_CR","h_PhoPt_ProbL1Trig_Elec_CR","h_qmulti_ProbL1Trig_Elec_CR","h_Photon_Eta_ProbL1Trig_Elec_CR","h_Photon_Phi_ProbL1Trig_Elec_CR","h_MET_Phi_ProbL1Trig_Elec_CR","h_leadJets_qmulti_ProbL1Trig_Elec_CR","h_leadJet_Pt_ProbL1Trig_Elec_CR","h_leadbjet_tag_ProbL1Trig_Elec_CR","h_nvrtx_ProbL1Trig_Elec_CR","h_minDR_Jets_EMObject_ProbL1Trig_Elec_CR","h_Phi_leadJet1_ProbL1Trig_Elec_CR","h_Eta_leadJet1_ProbL1Trig_Elec_CR","h_Pt_leadJet1_ProbL1Trig_Elec_CR","h_dPhi_METJet1_ProbL1Trig_Elec_CR","h_Phi_leadJet2_ProbL1Trig_Elec_CR","h_Eta_leadJet2_ProbL1Trig_Elec_CR","h_Pt_leadJet2_ProbL1Trig_Elec_CR","h_dPhi_METJet2_ProbL1Trig_Elec_CR","h_Phi_leadJet3_ProbL1Trig_Elec_CR","h_Eta_leadJet3_ProbL1Trig_Elec_CR","h_Pt_leadJet3_ProbL1Trig_Elec_CR","h_dPhi_METJet3_ProbL1Trig_Elec_CR","h_Phi_leadJet4_ProbL1Trig_Elec_CR","h_Eta_leadJet4_ProbL1Trig_Elec_CR","h_Pt_leadJet4_ProbL1Trig_Elec_CR","h_dPhi_METJet4_ProbL1Trig_Elec_CR","h_Phi_matchedJet_ProbL1Trig_Elec_CR","h_Eta_matchedJet_ProbL1Trig_Elec_CR","h_Pt_matchedJet_ProbL1Trig_Elec_CR","h_HT5HT_ProbL1Trig_Elec_CR","h_Mt_phoMET_ProbL1Trig_Elec_CR","h_dPhi_phoMet_ProbL1Trig_Elec_CR"};

  varName2={"h_St_ProbL1Trig_Pho_SR","h_HT_ProbL1Trig_Pho_SR","h_NhadJets_ProbL1Trig_Pho_SR","h_NBJets_ProbL1Trig_Pho_SR","h_MET_ProbL1Trig_Pho_SR","h_PhoPt_ProbL1Trig_Pho_SR","h_qmulti_ProbL1Trig_Pho_SR","h_Photon_Eta_ProbL1Trig_Pho_SR","h_Photon_Phi_ProbL1Trig_Pho_SR","h_MET_Phi_ProbL1Trig_Pho_SR","h_leadJets_qmulti_ProbL1Trig_Pho_SR","h_leadJet_Pt_ProbL1Trig_Pho_SR","h_leadbjet_tag_ProbL1Trig_Pho_SR","h_nvrtx_ProbL1Trig_Pho_SR","h_minDR_Jets_EMObject_ProbL1Trig_Pho_SR","h_Phi_leadJet1_ProbL1Trig_Pho_SR","h_Eta_leadJet1_ProbL1Trig_Pho_SR","h_Pt_leadJet1_ProbL1Trig_Pho_SR","h_dPhi_METJet1_ProbL1Trig_Pho_SR","h_Phi_leadJet2_ProbL1Trig_Pho_SR","h_Eta_leadJet2_ProbL1Trig_Pho_SR","h_Pt_leadJet2_ProbL1Trig_Pho_SR","h_dPhi_METJet2_ProbL1Trig_Pho_SR","h_Phi_leadJet3_ProbL1Trig_Pho_SR","h_Eta_leadJet3_ProbL1Trig_Pho_SR","h_Pt_leadJet3_ProbL1Trig_Pho_SR","h_dPhi_METJet3_ProbL1Trig_Pho_SR","h_Phi_leadJet4_ProbL1Trig_Pho_SR","h_Eta_leadJet4_ProbL1Trig_Pho_SR","h_Pt_leadJet4_ProbL1Trig_Pho_SR","h_dPhi_METJet4_ProbL1Trig_Pho_SR","h_Phi_matchedJet_ProbL1Trig_Pho_SR","h_Eta_matchedJet_ProbL1Trig_Pho_SR","h_Pt_matchedJet_ProbL1Trig_Pho_SR","h_HT5HT_ProbL1Trig_Pho_SR","h_Mt_phoMET_ProbL1Trig_Pho_SR","h_dPhi_phoMet_ProbL1Trig_Pho_SR"};
    }

   vector <string>  xLabel;

   xLabel={"Sum of P_{T}^{Jets} & P_{T}^{#gamma} [GeV]","HT[GeV]","N_{jets}","N_{ b-jets}","p_{T}^{miss} [GeV]","p_{T}^{#gamma} [GeV]","q-multi","Eta coordinate of EM-obj","Phi coordinate of EM-obj","MET phi","q multi of leading jets","pT of leading jets","b-tagger deep csv value","number of vertices","mindR(matched Jet, EM obj)","Phi coordinate of leading Jet1","Eta coordinate of leading Jet1","P_{T} of leading Jet1","d#phi(P_{T}^{miss},lead Jet1)","Phi coordinate of leading Jet2","Eta coordinate of leading Jet2","P_{T} of leading Jet2","d#phi(P_{T}^{miss},lead Jet2)","Phi coordinate of leading Jet3","Eta coordinate of leading Jet3","P_{T} of leading Jet3","d#phi(P_{T}^{miss},lead Jet3)","Phi coordinate of leading Jet4","Eta coordinate of leading Jet4","P_{T} of leading Jet4","d#phi(P_{T}^{miss},lead Jet4)","Phi coordinate of matched Jet","Eta coordinate of matched Jet","P_{T} of matched Jet","HT5/HT","M_{T}^{miss & elec} [GeV]","dPhi(elec,MET)","M^{inv}_{lep1 & lep2}","P_{T}^{Z}"};                                                                                                                                              
  vector <int> rebin;
  rebin={5,5,1,1,5,10,1,4,4,4,1,4,4,1,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,2,4,4,4};
  vector<double> ymin ={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  vector<double> ymax={100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000};
  vector<double> xmin ={200,300,2,0,200,20,0,-10,-10,-10,0,0,0,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0.9,0,0,80,0};
  vector<double> xmax={1500,1500,11,8,1000,800,100,10,10,10,100,1000,1,50,0.5,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,2,400,3,100,800};

  cout<<"different vector sizes "<<endl;
  cout<<varName.size()<<"\t"<<varName2.size()<<"\t"<<xLabel.size()<<"\t"<<rebin.size()<<"\t"<<xmax.size()<<"\t"<<xmin.size()<<endl;
  bool flag=false;
  //  const char* filetag[10]={"TTGJets","pMSSM_MCMC_70_90438","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451","WJets","GJets","T5bbbbZG_10","T5bbbbZG_50","T5bbbbZG_200","T5bbbbZG_1500"};
   /* vector<TH1D*> hist_list_Njets; */
   /*    vector<TH1D*> hist_list_Bjets; */
   /*    vector<TH1D*> hist_list_MET; */
   /*    vector<TH1D*> hist_list_PhoPt; */
   /*    //vector<TH1D*> hist_list_Mt;                                                                                                                                    */
   /*    vector<TH1D*> hist_list_ST; */
   /*    vector<TH1D*> hist_list_HT; */

  //  TFile* fout = new TFile("TF_allin1_LLEstimation_electron.root","RECREATE");
  sprintf(hname1,"temp.root");//,hname); 
  TFile* fout = new TFile(hname1,"RECREATE");
  // sprintf(hname,"EventYields_TF_LL_muon_allProcess_binsV3_phoID_loose_09Jan24.txt");
  // std::ofstream file_;
  // file_.open(hname,ios::out);
    n_files=6;  
  for(int i_file=0; i_file<n_files;i_file++)
    {      
      //      vector<TH1D*> hist_list_Njets;
      //      vector<TH1D*> hist_list_Bjets;
      vector<TH1D*> hist_list_MET;
      vector<TH1D*> hist_list_PhoPt;
      vector<TH1D*> hist_list_ST;
      for(int i_cut=0; i_cut<varName.size();i_cut++)
	{
	  //	  	  if(i_cut!=4) continue;
	  vector<TH1D*> hist_list_Njets;
	  vector<TH1D*> hist_list_Bjets;

	  sprintf(hist_name,"%s",varName[i_cut].c_str());
	  sprintf(hist_name2,"%s",varName2[i_cut].c_str());
	  cout<<hist_name<<"\t"<<hist_name2<<"\t"<<i_cut<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
          TH1D* h_resp = (TH1D*)f1[i_file]->Get(hist_name);//CR
	  TH1D* h_resp2 = (TH1D*)f[i_file]->Get(hist_name2); // SR
	  
	  cout<<"resp "<<h_resp->Integral()<<"\t"<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin[i_cut]<<"\t"<<xmin[i_cut]<<"\t"<<xmax[i_cut]<<endl;
	  h_resp->Rebin(rebin[i_cut]);
	  h_resp2->Rebin(rebin[i_cut]);
	  if(rebin[i_cut]!=1){
	  h_resp->Rebin(2);
          h_resp2->Rebin(2);
	  }
	  
	  h_resp2= setMyRange(h_resp2,xmin[i_cut],xmax[i_cut]+0.01*xmax[i_cut]);//xmin,xmax+6);
	  setLastBinAsOverFlow(h_resp2);
	  h_resp= setMyRange(h_resp,xmin[i_cut],xmax[i_cut]+0.01*xmax[i_cut]);
	  setLastBinAsOverFlow(h_resp);
	  
	  hist_list_Njets.push_back(h_resp2); // 0th is SR
	  hist_list_Njets.push_back(h_resp); //1st is CR
	  double factor=1.0;
          // h_resp->Scale(factor/h_resp->Integral());
          // h_resp2->Scale(factor/h_resp2->Integral());
	  hist_list_Bjets.push_back(h_resp2);
	  hist_list_Bjets.push_back(h_resp); 
	  
	  cout<<h_resp->Integral()<<"\t"<<f[i_file]->GetName()<<endl;	 
	  cout<<" hist_list_Njets.size() "<<hist_list_Njets.size()<<"\t "<<endl;//hist_list_Bjets.size() "<<hist_list_Bjets.size()<<endl;
	  //path to save the png file
	  float energy=energyy[i_file];
	  int xrange=0.0;
	  TH1D* hNjets_total =(TH1D*)hist_list_Njets.at(1)->Clone();
	  TH1D* hNjets_ratio = (TH1D*)hist_list_Njets.at(0)->Clone();
	  hNjets_ratio = setMyRange(hNjets_ratio,xmin[i_cut],xmax[i_cut]+0.01*xmax[i_cut]);
	  hNjets_ratio->Divide(hNjets_total);
	  gROOT->ForceStyle();
	  setLastBinAsOverFlow(hNjets_ratio);
	  sprintf(full_path,"%s/%s_%s_%s",pathname.c_str(),string_png,varName[i_cut].c_str(),filetag[i_file].c_str());
	  if(i_cut==2 || i_cut==3)
	    generate_1Dplot(hist_list_Njets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Entries",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);
	  else
	    generate_1Dplot(hist_list_Njets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Entries",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);

	  sprintf(full_path,"%s/%s_%s_%s_normalize",pathname.c_str(),string_png,varName[i_cut].c_str(),filetag[i_file].c_str());

	  
	  // hist_list_Njets.at(0)->GetXaxis()->SetLabelSize(0.045);
          // hist_list_Njets.at(1)->GetXaxis()->SetLabelSize(0.045);
	  // double factor=1.0;
	  // hist_list_Njets.at(0)->Scale(factor/hist_list_Njets.at(0)->Integral());
	  // hist_list_Njets.at(1)->Scale(factor/hist_list_Njets.at(1)->Integral());
	  // hist_list_Njets.at(0)->GetXaxis()->SetLabelSize(0.045);
	  // hist_list_Njets.at(1)->GetXaxis()->SetLabelSize(0.045);
	  //	  TH1D* h_temp = (TH1D*)hist_list_Bjets.at(0)->Clone();	  
	  hist_list_Bjets.at(0)->Scale(factor/hist_list_Bjets.at(0)->Integral());//(factor/hist_list_Njets.at(0)->Integral())*hist_list_Njets.at(0);//->Scale(factor/hist_list_Njets.at(0)->Integral());
	  //	  hist_list_Bjets.push_back(h_temp);
	  //	  TH1D*	h_temp1 = (TH1D*)hist_list_Bjets.at(1)->Clone();
	  hist_list_Bjets.at(1)->Scale(factor/hist_list_Bjets.at(1)->Integral());//(factor/hist_list_Njets.at(1)->Integral())*hist_list_Njets.at(1);
          //h_temp1->Scale(factor/hist_list_Njets.at(1)->Integral());
	  //          hist_list_Bjets.push_back(h_temp1);

	  hNjets_ratio =(TH1D*)hist_list_Njets.at(1)->Clone();//(TH1D*)hist_list_Bjets.at(1)->Clone();
	  //hNjets_ratio = hist_list_Bjets.at(1);
	  hNjets_ratio->Divide(hist_list_Bjets.at(0));
	  // hist_list_Bjets.at(1)->Draw("Hist");
	  // hist_list_Bjets.at(1)->SetLineColor(kRed);
	  // hist_list_Bjets.at(0)->SetLineColor(kGreen);
	  // hist_list_Bjets.at(1)->SetLineStyle(3);
	  // hist_list_Bjets.at(0)->Draw("Hist sames");
	  //	  hNjets_ratio->Draw("Hist");
	  // setLastBinAsOverFlow(hNjets_ratio);
	  if(i_cut==2 || i_cut==3)
            generate_1Dplot(hist_list_Bjets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Normalize",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);
          else
            generate_1Dplot(hist_list_Bjets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Normalize",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);

	}
      //fout->Close();
      
    }
}







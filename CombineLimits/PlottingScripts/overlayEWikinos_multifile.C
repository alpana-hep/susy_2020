const int n_pl = 4;
bool logx = false;
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[5] ={"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR","Failed Iso"};//,"Failed acceptance","1e CR"};// {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[4] = {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {4,4,4,4,4,2,4,2,4,2,4,2};
int line_style[12] = {1,1,1,1,2,1,2,1,2,1,2,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};                                                                               
int line_color[9] = {kBlue,kBlack,kYellow+2,kYellow+2,kMagenta,kMagenta, kRed+2,9,kCyan+2};//,45,kMagenta,kGray+1,kRed,kBlue+2,kMagenta,kCyan};
int line_color1[9]= {kBlue,kBlack,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
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

void generate_1Dplot(vector<TGraph*> hist, char const *tag_name="",char const *xlabel="",char const *ylabel="", float energy=-1, int rebin=-1,double ymin=0,double ymax=0,int xmin=-1,int xmax=-1, char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true,  vector<string> legend_texts={"nil"}, char const *legend_title="", TString model="", vector<int> linecolorlist={}){  

  cout<<" inside generate 1D plot "<<"\t"<<legend_title<<"\t"<<endl;
  TCanvas *canvas_n1 =      new TCanvas(tag_name, tag_name,1200,900);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  canvas_n1->SetTopMargin(0.05);
  canvas_n1->SetRightMargin(0.035);
  canvas_n1->SetLeftMargin(0.13);
canvas_n1->SetBottomMargin(0.13);
  // auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
  // pad_1->SetTopMargin(0.04);
  // pad_1->SetBottomMargin(0.33);
  // pad_1->SetRightMargin(0.035);
  // pad_1->SetLeftMargin(0.13);
  // auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
  // p1->SetBottomMargin(0.026);
  // p1->SetRightMargin(0.035);
  // p1->SetLeftMargin(0.13);
  // p1->SetTopMargin(0.1);
  // //  p1->cd();
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
  //  std::string leg_head_str = ;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.15,0.77,0.9,0.945);  
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
    legend->SetNColumns(2);
  char* lhead = new char[100];
  //  cout<<"before legend fixing "<<"\t"<<leg_head<<"\t"<<legend_title<<endl; 
  //  sprintf(lhead,"%s ",leg_head);
 auto  legend1 = new TLegend(0.15,0.7,0.3,0.65);
 legend1->SetTextSize(0.035);
 legend1->SetLineColor(kWhite);
 legend1->SetHeader(legend_title);
 //   cout<<"before legend fixing "<<"\t"<<leg_head<<"\t"<<legend_title<<endl;

 legend->SetHeader(legend_title);
 legend->SetLineColor(kWhite);
  cout<<"after legend fixing "<<endl;
  TLegendEntry* leg_entry[11];
  float x_label_size = 0.045;
  //  double ymin = 100000.0;
  //double ymax = 0.0;
  double xrange = xmax;
  
  // float energy = energyy;
  vector<TGraph*> hist_list_temp;
  TMultiGraph* mg = new TMultiGraph();
  cout<<" hist.size() = "<<hist.size()<<endl;
  for(int i =0;i<(int)hist.size(); i ++) {
    
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    //    hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    //    normalize = true;
     hist.at(i)->GetXaxis()->SetTitle(xlabel);
     hist.at(i)->GetYaxis()->SetTitle(ylabel);
     //     cout<<" i"<<i<<endl;
     //   hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);
    hist.at(i)->SetLineWidth(line_width[i]);
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(linecolorlist[i]);
    hist.at(i)->SetTitle(" ");
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->SetLineColor(linecolorlist[i]);
    hist.at(i)->SetTitle(" ");
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleSize(0.06);
    hist.at(i)->GetYaxis()->SetLabelSize(0.06);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.);
    //decorate(hist.at(i),i, 0);
    hist.at(i)->SetMarkerSize(1.1);
    hist.at(i)->SetMarkerStyle(20);
    hist.at(i)->SetMarkerColor(linecolorlist[i]);
    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);
    hist.at(i)->GetYaxis()->SetRangeUser(ymin,ymax);
    //new ones
    // hist.at(i)->GetXaxis()->SetTitleSize(0.08);
    // hist.at(i)->GetXaxis()->SetLabelSize(0.06);

    // hist.at(i)->GetYaxis()->SetTitleSize(0.07);
    //hist.at(i)->GetYaxis()->SetLabelSize(0.06);

    // hist.at(i)->GetXaxis()->SetTitleOffset(3);
    // hist.at(i)->GetXaxis()->SetLabelOffset(1.6);

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
    //cout<<"before setting legends  "<<endl;
    
    legName.push_back(hist.at(i)->GetName());
    if(i!=0 && !model.Contains("T5gg")){
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    }
    else if(model.Contains("T5gg")){
      leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"l");
      leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());

    }
    // if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    // if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    // hist.at(i)= setMyRange(hist.at(i),xmin,xmax+4);
    //setLastBinAsOverFlow(hist.at(i));
    //    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);
    if(i!=0 && !model.Contains("T5gg"))
      mg->Add(hist.at(i));
    else if(model.Contains("T5gg"))
      mg->Add(hist.at(i));
  }
  cout<<"outside  the hist loop "<<endl;
  // if(ymin == 0.0) ymin = 1e-3;
  // if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  // for(int i = 0;i<(int)hist.size(); i++) {
  //   if(!normalize) {
  //   //   if(model.Contains("TChiWG") || model.Contains("T6ttZg") || model.Contains("TChiNG") || model.Contains("WGJets") || model.Contains("WlnuJets") || model.Contains("ttbarG")|| model.Contains("ttbarJets") || model.Contains("GJets") || model.Contains("ZnunuGJets") ) hist.at(i)->GetYaxis()->SetRangeUser(0.01,1000*ymax);
  //   //   else       hist.at(i)->GetYaxis()->SetRangeUser(0.01,0.01*ymax);
  //   // }
  //   // else
  //   //   {  hist.at(i)->GetYaxis()->SetRangeUser(0.00001,5.0);
  //   // 	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
  //   //   }
  //   //    p1->SetGrid();
  //   //hs_var->Add(hist.at(i));
  //   // hs_var->SetMinimum(0.00001);
  //   // hs_var->SetMaximum(ymax*60);
  //   //gPad->SetLogu
  //   cout<<"i Alps "<<i<<endl;
  //   if(i==0) hist.at(i)->Draw("hist ");
  //   else hist.at(i)->Draw("hist sames");
	
  // }
 // hs_var->SetMinimum(0.0);
 // hs_var->SetMaximum(ymax+0.5);
 // TAxis *axis51=  mg->GetXaxis();
 // axis51->SetRangeUser(0,2000);
 
  mg->GetXaxis()->SetLimits(xmin,xmax);
  mg->GetYaxis()->SetRangeUser(ymin,ymax);
  mg->SetMaximum(ymax);
  mg->SetMinimum(ymin);
  gPad->SetLogy();
  mg->Draw("ALP");
   if(!model.Contains("T5gg")){

 hist.at(0)->SetMarkerSize(1.2);
 hist.at(0)->SetMarkerStyle(20);
 hist.at(0)->SetMarkerColor(kRed+1);
 hist.at(0)->SetLineColor(kRed+1);
 hist.at(0)->SetFillColor(kRed-9);
 hist.at(0)->Draw("3same");
 hist.at(0)->Draw("LXsame");
 }
 //mg->Draw("L same");
 //hist.at(0)->Draw("3same");

  // TAxis *axis51=  mg->GetXaxis();
  // axis51->SetRangeUser(0,2000);
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(xmin,2000);
  
  mg->GetXaxis()->SetTitleOffset(1.0);
  cout<<"after drawing "<<endl;
  gPad->Modified(); gPad->Update();
  mg->GetXaxis()->SetTitle(xlabel);
  mg->GetYaxis()->SetTitle(ylabel);
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->SetTitleSize(00.05);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.0);

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

  if(!model.Contains("T5gg")){
  leg_entry[7] = legend->AddEntry(hist.at(0),legend_texts[0].c_str(),"ep");
  leg_entry[7]->SetTextColor(hist.at(0)->GetLineColor());
  }
 TLine *l =new TLine(xmin,1.0,xmax,1.0);
 if(model.Contains("T5gg"))
   l->Draw("sames");             
  legend->Draw();

  //  legend1->Draw();


  // if(log_flag) {
  //     gPad->SetLogy();
  //   }
  // if(logx)
  //   gPad->SetLogx();
  gPad->Update(); 
  TLatex* textOnTop = new TLatex();
  //new
  textOnTop->SetTextSize(0.04);
  //  textOnTop->DrawLatexNDC(0.13,0.965,"CMS #it{#bf{Simulation Preliminary}}");

  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.73,0.965,en_lat);


  gPad->Modified();
                                                                                       
   //  hist_ratio->SetLineWidth(2);
   //  hist_ratio->SetLineStyle(1);
   //  hist_ratio->SetMarkerSize(0.2);
   //  hist_ratio->SetLineColor(kBlack);
   //  hist_ratio->SetTitle(" ");
   //  hist_ratio->GetXaxis()->SetTitleSize(0.13);
   //  hist_ratio->GetYaxis()->SetTitle("Z(#nu#nu)/Z(ll)");//TF = #frac{N_{SR}}{N_{CR}}");//(0#mu,1#gamma)}{(1#mu,1#gamma)}");
   //  //    hist_ratio->GetXaxis()->SetLabelSize(0.1);
   //  hist_ratio->GetYaxis()->SetRangeUser(-0.001,0.03);
   //  //hist_ratio->GetXaxis()->SetRangeUser(xmin,xmax+4);
   //  // hist_ratio= setMyRange(hist_ratio,xmin,xmax+6);
   //  //setLastBinAsOverFlow(hist_ratio);

   //  // if(which_TFbins==1) //default 8 bins                                                                                                    
   //  //   hist_ratio->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                 
   //  // else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100
   //  //   hist_ratio->GetXaxis()->SetRangeUser(0,18);
   //  // else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300                                                                                        
   //  //    hist_ratio->GetXaxis()->SetRangeUser(0,39);
    
   //  hist_ratio->GetXaxis()->SetLabelSize(0.0450);
   //  hist_ratio->GetYaxis()->SetTitleSize(0.13);
   //  hist_ratio->GetYaxis()->SetLabelSize(0.08);
   //  hist_ratio->GetYaxis()->SetTitleOffset(.4);
   //  hist_ratio->SetMarkerSize(1.0);
   //  hist_ratio->SetMarkerStyle(20);
   //  hist_ratio->SetMarkerColor(kBlue);
   //  hist_ratio->GetXaxis()->SetTitle(xlabel);
   //  hist_ratio->GetYaxis()->SetNdivisions(505);
   //  //new
   //  hist_ratio->GetXaxis()->SetTitleSize(0.05);
   //  hist_ratio->GetXaxis()->SetLabelSize(0.11);
   //  hist_ratio->GetYaxis()->SetTitleSize(0.125);
   //  hist_ratio->GetYaxis()->SetNdivisions(505);

   //  hist_ratio->GetXaxis()->SetTitleOffset(1);
   //  hist_ratio->GetYaxis()->SetTitleOffset(0.41);
   //  hist_ratio->GetXaxis()->SetTitleSize(0.14);

   //  hist_ratio->GetYaxis()->SetLabelSize(0.11);
    
   //  if(normalize){
   //    hist_ratio->GetYaxis()->SetRangeUser(0,2);
   //    hist_ratio->GetXaxis()->SetLabelOffset(0);
   //    //      hist_ratio->GetXaxis()->SetLabelSize(0.4);
   //  }
   //  else {
   //    hist_ratio->GetYaxis()->SetRangeUser(0,10);
   //  }
   //  // gStyle->SetLabelOffset(1.2);
   //  // gStyle->SetLabelSize(1.2);
   //   //     hist_ratio->GetXaxis()->SetLabelOffset();
   //   //hist_ratio->GetYaxis()->SetLabelSize(x_label_size);
   // pad_1->cd();
   // //   pad_1->SetGrid();
   // // if(which_TFbins==1){
   // // TLine *l =new TLine(0,1.0,10,1.0);   
   // // hist_ratio->Draw("");
   // // l->Draw("sames");
   // // TLine *l1 =new TLine(0,1.5,10,1.5);
   // // l1->SetLineStyle(7);
   // // l1->Draw("sames");
   // // TLine *l2 =new TLine(0,0.5,10,0.5);
   // // l2->SetLineStyle(7);

   // // l2->Draw("sames");
   // // }

   // // else{
     
   // //    TLine *l =new TLine(0,1.0,18,1.0);
   // // hist_ratio->Draw("");
   // // l->Draw("sames");
   // // TLine *l1 =new TLine(0,1.5,18,1.5);
   // // l1->SetLineStyle(7);
   // // l1->Draw("sames");
   // // TLine *l2 =new TLine(0,0.5,18,0.5);
   // // l2->SetLineStyle(7);

   // // l2->Draw("sames");
   // // }
   // // TLine *l =new TLine(xmin,.01,xrange,.01);
   // // hist_ratio->Draw("");
   // // l->Draw("sames");
   // // TLine *l1 =new TLine(xmin,0.02,xrange,0.02);
   // // l1->SetLineStyle(7);
   // // l1->Draw("sames");
   // // TLine *l2 =new TLine(xmin,0.03,xrange,0.03);
   // // l2->SetLineStyle(7);

   // // l2->Draw("sames");

   //   if(normalize){
   //     //       hist_ratio->Scale(1.0/hist_ratio->Integral());
   //     TLine *l =new TLine(xmin,1,xrange,1);
   //     hist_ratio->Draw("");
   //     l->Draw("sames");
   //     TLine *l1 =new TLine(xmin,1.5,xrange,1.5);
   //     l1->SetLineStyle(7);
   //     l1->Draw("sames");
   //     TLine *l2 =new TLine(xmin,0.5,xrange,0.5);
   //     l2->SetLineStyle(7);
   //     l2->Draw("sames");
   // }

   // else{

   //    TLine *l =new TLine(xmin,5,xrange,5);
   //    hist_ratio->Draw("");
   //     l->Draw("sames");
   //     TLine *l1 =new TLine(xmin,10,xrange,10);
   //     l1->SetLineStyle(7);
   //     l1->Draw("sames");
   //     TLine *l2 =new TLine(xmin,15,xrange,15);
   //     l2->SetLineStyle(7);
   //    l2->Draw("sames");

   // //   TLine *l =new TLine(0,10,50,10);
   // // hist_ratio->Draw("");
   // // l->Draw("sames");
   // // TLine *l1 =new TLine(0,12,50,12);
   // // l1->SetLineStyle(7);
   // // l1->Draw("sames");
   // // TLine *l2 =new TLine(0,5,50,5);
   // // l2->SetLineStyle(7);

   // // l2->Draw("sames");
   // }
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


void overlayEWikinos_multifile(string pathname, string model, int which_plot)
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
  int n=0;
  char *dataset=new char[200];
  char *year =new char[200];
  //  float energyy[2]={};
  int n_files=42;
  char *string_png = new char[200];
  vector<string>baseline1, baseline;
  vector<string> legend_texts, legend_texts_v1;
  //  char *legend_title= new char[2000];
  int which_Lept=1; 
 int which_TFBins=1;
  vector<string>  Nlsp_m;
  string png_string = "";//gluino_m.c_str();
  TString Model = model.c_str();
  baseline = {"gr1d_xseclimit"};//No_selection","veto_chargTracks","photon_pT40","MET_100","nJets_2","ST_300","dPhi_cut","overlap_removal","MET_200","MET_300"}; //lepton_veto
  vector <string > varName ={"gr1d_xseclimit"};//h_St","h_HT","h_NhadJets","h_NBJets","h_MET","h_PhoPt","h_Photon_Eta","h_Photon_Phi","h_MET_Phi","FR_nbtagBins","h_Mt_phoMET","h_dPhi_phoMet","h_Nphotons","h_dPhi_METLeadJet"};
  vector <string> xlabel = {"m_{#tilde{g}} [GeV]"};//Sum of p_{T}^{Jets} & p_{T}^{#gamma} [GeV]","HT[GeV]","N_{jets}","N_{ b-jets}","p_{T}^{miss} [GeV]","p_{T}^{#gamma} [GeV]","#eta^{#gamma}","#phi^{#gamma}","#phi^{MET}","Bin No.","M_{T}^{miss & #gamma} [GeV]","dPhi(#gamma,MET)","N_{#gamma}","d#phi(p_{T}^{miss},lead Jet1)"};
  legend_texts ={""};
  vector<int>  linecolorlist = {kBlue,kBlack,kBlue,kBlack,kRed,kMagenta, kRed+2,9,kCyan+2};
  vector <string> ylabel = {"m_{#tilde{#chi}_{1}^{0}} [GeV]"};
  if(Model.Contains("T5bbbbZg")  )
    legend_texts ={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};//,"veto lept+isoTracks","p_{T}^{#gamma} > 40","p_{T}^{miss} [GeV] > 100","N_{jets} #geq 2","Sum of p_{T}^{Jets} & p_{T}^{#gamma} > 300 ","d#phi(p_{T}^{miss},lead Jets) > 0.3 ","Preselection","Preselection + p_{T}^{miss} [GeV] > 200","Preselection + p_{T}^{miss} [GeV] > 300"};

  if(Model.Contains("T5qqqqHg"))
      legend_texts ={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}"};
  if(Model.Contains("T5ttttZg")  )
    legend_texts ={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
  if(Model.Contains("T6ttZg")  )
    legend_texts ={"#tilde{t} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
  vector <string>legend_title;
  vector<string> filetag= {"FullRun2"};// {"2018","2017","2016postVFP","2016preVFP","2016","FullRun2"};
  vector<float> energyy={137.19};//{ 59.74,41.53,16.5,19.5,36,137.19};
  vector <int> rebin;
  rebin={1};//,4,1,1,5,5,4,4,4,1,4,4,1,4};//,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,2,4,4,4};
  vector<double> ymin ={0.0};//,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};//,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  vector<double> ymax={3000};//,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000};//,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000};
  vector<double> xmin ={0};//,0,0,0,0,0,-5,-5,-5,0,0,0,0,0,0};//-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0.9,0,0,80,0};
  vector<double> xmax={3000};//,2500,20,16,1000,900,5,5,5,3,800,5,5,5};//,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,2,400,3,100,800};

  cout<<"different vector sizes "<<endl;
  cout<<varName.size()<<"\t"<<baseline.size()<<"\t"<<xlabel.size()<<"\t"<<rebin.size()<<"\t"<<xmax.size()<<"\t"<<xmin.size()<<"\t"<<legend_texts.size()<<endl;
  bool flag=false;
  n_files=2;
    if(Model.Contains("T5bbbbZg") && which_plot==1 ){//&& Gluino_m.Contains("2200")){    
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};//, m_{#tilde{g}} = 2200 "};
      legend_texts_v1 = {"SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) : p_{T}^{#gamma}>40", "SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) : p_{T}^{#gamma}>100" };//m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};
    f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");//out_T5bbbbZg_2200_10.root");
    f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      ymax={3000};
      xmax={2800};
      ymin={0};
      xmin={1500};

    n_files=2;
    png_string = "pt40vspt100";
    
  }

    if(Model.Contains("T5bbbbZg") && which_plot==2 ){
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
      legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0" };
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");   
    f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
    f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_withMvaCut.root");
    f[3] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
    n_files=4;
    png_string = "pt40vspt100_withMvaCut";
  }

    if(Model.Contains("T6ttZg") && which_plot==1 ){
      legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
      legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0" };
      f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100_withMvaCut.root");
      n_files=4;
      png_string = "pt40vspt100_withMvaCut_T6ttZg";
      ymax={2000};
      xmax={1600};
      ymin={0};
      xmin={800};

  }
    if(Model.Contains("T6ttZg") && which_plot==2 ){
      legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
      //f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");      
      f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      linecolorlist = {kBlue, kBlue+1, kYellow+2, kCyan, kViolet,kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=4;
      png_string = "STNbjets_METBins_withMvaCut_T6ttZg";
        ymax={2000};
      xmax={1600};
      ymin={0};
      xmin={800};

  }
    if(Model.Contains("T6ttZg") && which_plot==3 ){
      legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");
      f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
linecolorlist = {kBlue, kBlue+1, kYellow+2, kCyan, kViolet,kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=4;
       ymax={2000};
      xmax={1600};
      ymin={0};
      xmin={800};

      png_string = "STMETNbjet_phoptBins_withMvaCut_T6ttZg";
  }
    if(Model.Contains("T6ttZg") && which_plot==4 ){
      legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      //f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");
      f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");

      f[4] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=6;
      ymax={2000};
      xmax={1600};
      ymin={0};
      xmin={800};
       xlabel = {"m_{#tilde{t}} [GeV]"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut_T6ttZg";
  }

  if(Model.Contains("T6ttZg") && which_plot==1 ){
    legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
    legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamm\
a}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
    f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
    f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
    // f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");
    // f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
    // f[4] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
    //    varName ={"sigmatheory","exp"};
    ymax={2000};
    xmax={1600};
    ymin={0};
    xmin={800};
    n_files=5;
    png_string = "STMETNbjet_phoptBins_withMvaCut_T6ttZg";
  }

    
    if(Model.Contains("T6ttZg") && which_plot==5 ){
      legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
      legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[4] = new TFile("Excl_out_T6ttZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");

      ymax={2000};
      xmax={1600};
      ymin={0};
      xmin={800};
      n_files=5;
      png_string = "STMETNbjet_phoptBins_withMvaCut_T6ttZg";
  }

    if(Model.Contains("T5gg") && which_plot==11 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=10GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v3_MET_200.root");
      //      f[3] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v7_MET_200.root");
      // f[5] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM10_final";
      linecolorlist = {kBlue, kBlack, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }

    if(Model.Contains("T5gg") && which_plot==111 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=10GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#g\
amma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM10_final_update";
      linecolorlist = {kBlue, kBlack, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }


     if(Model.Contains("T5gg") && which_plot==1 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=10GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root");                                                                      
      f[4] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
       f[5] = new TFile("Excl_T5gg_2200_delM10_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut_v1.root");                                                                           
      varName ={"exp"};
      n_files=6;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut_delM10";

      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }


    if(Model.Contains("T5gg") && which_plot==2 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=30GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root");
      f[4] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      f[5] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut_v1.root");
      varName ={"exp"};
      n_files=6;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut_delM30";

      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }

    
    if(Model.Contains("T5gg") && which_plot==22 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=30GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200.root");
      //      f[3] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200.root");
      //      f[5] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM30_final";

      linecolorlist = {kBlue, kBlack, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }

    if(Model.Contains("T5gg") && which_plot==222 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=30GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      //      f[3] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");                                                                                                    
      f[3] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      //      f[5] = new TFile("Excl_T5gg_2200_delM30_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");                                                                                                    
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM30_final_update";

      linecolorlist = {kBlue, kBlack, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }


    if(Model.Contains("T5gg") && which_plot==3 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=100GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root");
      f[4] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      f[5] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut_v1.root");
      varName ={"exp"};
      n_files=6;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut_delM100";

      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }
    if(Model.Contains("T5gg") && which_plot==33 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=100GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamm\
a}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200.root");
      //      f[3] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200.root");
      //      f[5] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM100_final";

      linecolorlist = {kBlue, kBlack,  kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }

    if(Model.Contains("T5gg") && which_plot==333 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{2}^{0}, #tilde{#chi}_{2}^{0} #rightarrow #gamma #tilde{#chi}_{1}^{0}, #delta_{m}=100GeV"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      //      f[3] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");                                                             
      f[3] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      //      f[5] = new TFile("Excl_T5gg_2200_delM100_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");                                                             
      varName ={"exp"};
      n_files=4;
      ymax={1000};
      xmax={2000};
      ymin={0.01};
      xmin={0};
      xlabel = {"m_{#tilde{#chi}_{2}^{0}} [GeV]"};
      ylabel={"r value"};
      png_string = "overlay_v3vsV7SRbins_delM100_final_update";

      linecolorlist = {kBlue, kBlack,  kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
  }





    if(Model.Contains("T5qqqqHg") && which_plot==4 ){
      legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_out_T5qqqqHg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=6;
      ymax={3000};
      xmax={2800};
      ymin={0};
      xmin={1500};
       xlabel = {"m_{#tilde{g}} [GeV]"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut_T5qqqqHg";
  }

    if(Model.Contains("T5ttttZg") && which_plot==4 ){
      legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_out_T5ttttZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=6;
      ymax={3000};
      xmax={2800};
      ymin={0};
      xmin={1400};
      xlabel = {"m_{#tilde{g}} [GeV]"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut";
  }
    
    if(Model.Contains("T5bbbbZg") && which_plot==41 ){
      legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      linecolorlist = {kBlue, kBlack, kYellow+2, kYellow+2, kViolet,kViolet,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      n_files=6;
      ymax={3000};
      xmax={2800};
      ymin={0};
      xmin={1500};
      xlabel = {"m_{#tilde{g}} [GeV]"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut";
  }

    if(Model.Contains("TChiWG") && which_plot==4 ){
      legend_title={"TChiWG , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-}"};//#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"Theory cross section","SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_TChiWG_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_TChiWG_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root");
      f[4] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      f[5] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200_v1.root");
      linecolorlist = {kGreen, kBlue, kBlack,kYellow+2, kYellow+2, kMagenta+1,kMagenta+1,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      varName ={"exp"}; 
      n_files=4;
      ymax={10};
      xmax={1400};
      ymin={0.0001};
      xmin={300};
      xlabel = {"m_{#tilde{#chi}_{1}^{0}} [GeV]"};
      ylabel={"Cross section (pb)"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut";
  }

    if(Model.Contains("TChiWG") && which_plot==411 ){
      legend_title={"TChiWG , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-}"};//#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};                                                                                                                                             
      legend_texts_v1 = {"Theory cross section","SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma})-Merge","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma})-Merge with BDT score>"};
      f[0] = new TFile("Excl_TChiWG_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_TChiWG_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[4] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      f[6] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[7] = new TFile("Excl_TChiWG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");

      linecolorlist = {kGreen, kBlue, kBlack,kYellow+2, kYellow+2, kMagenta+1,kMagenta+1,kRed, kRed+2,9, kGreen+2, kYellow+2};
      varName ={"exp"};
      n_files=8;
      ymax={10};
      xmax={1400};
      ymin={0.0001};
      xmin={300};
      xlabel = {"m_{#tilde{#chi}_{1}^{0}} [GeV]"};
      ylabel={"Cross section (pb)"};
      png_string = "overlay_v3vsV7merge_STMETNbjet_phoptBins_withMvaCut";
  }

    if(Model.Contains("TChiNG") && which_plot==4 ){
      legend_title={"TChiNG , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-}"}; 
      legend_texts_v1 = {"Theory cross section","SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
      f[0] = new TFile("Excl_TChiNG_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_TChiNG_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v3_MET_200_v1.root");
      f[3] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut_v1.root");
      f[4] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      linecolorlist = {kGreen, kBlue, kBlack,kYellow+2, kYellow+2, kMagenta+1,kMagenta+1,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
      varName ={"exp"};
      n_files=4;
      ymax={10};
      xmax={1400};
      ymin={0.0001};
      xmin={300};
      xlabel = {"m_{#tilde{#chi}_{1}^{0}} [GeV]"};
      ylabel={"Cross section (pb)"};
      png_string = "overlay_v3vsV7STMETNbjet_phoptBins_withMvaCut";
  }

     if(Model.Contains("TChiNG") && which_plot==411){
      legend_title={"TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-}"};
      legend_texts_v1 = {"Theory cross section","SP bins : p_{T}^{#gamma}>40","SP bins : p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0", "Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma})-Merge","Bins(p_{T}^{miss},ST,N_{jets}^{b},p_{T}^{#gamma})-Merge with BDT score>"};
      f[0] = new TFile("Excl_TChiNG_h_Sbins_LL_MET_200.root");
      f[1] = new TFile("Excl_TChiNG_h_Sbins_LL_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v3_MET_200.root");
      f[3] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
      f[4] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[5] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
      f[6] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200.root");
      f[7] = new TFile("Excl_TChiNG_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");

      linecolorlist = {kGreen, kBlue, kBlack,kYellow+2, kYellow+2, kMagenta+1,kMagenta+1,kRed, kRed+2,9, kGreen+2, kYellow+2};
      varName ={"exp"};
      n_files=8;
      ymax={10};
      xmax={1400};
      ymin={0.0001};
      xmin={300};
      xlabel = {"m_{#tilde{#chi}_{1}^{0}} [GeV]"};
      ylabel={"Cross section (pb)"};
      png_string = "overlay_v3vsV7_STMETNbjet_phoptBins_withMvaCut";
  }




    if(Model.Contains("T5bbbbZg") && which_plot==3 ){
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
      legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0","Add bins in p_{T}^{#gamma}","Add bins in p_{T}^{#gamma} & BDT score>0" };
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5bbbbZg_combine_v3sbins_MET_200.root");
      f[5] = new TFile("Excl_out_T5bbbbZg_combine_v3sbins_MET_200_withMvaCut.root");
      n_files=6;
      png_string = "withPtbins_withMvaCut";
  }

      if(Model.Contains("T5bbbbZg") && which_plot==4 ){
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
      legend_texts_v1 = {"p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0","Add bins in ST","Add bins in ST & BDT score>0" };
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200.root");
      f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_withMvaCut.root");
      n_files=6;
      png_string = "withSTbins_withMvaCut";
  }

      if(Model.Contains("T5bbbbZg") && which_plot==5 ){
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) + ST bins "};
      legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40","SP bins: p_{T}^{#gamma}>100","SP bins: p_{T}^{#gamma}>100 & BDT score>0","p_{T}^{#gamma}>40", "p_{T}^{#gamma}>100","p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0"};//,"Add bins in ST","Add bins in ST & BDT score>0" };
      // f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      // f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      // f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_withMvaCut.root");
      // f[3] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
      f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200.root");
      f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100.root");
      f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_withMvaCut.root");
      f[6] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100_withMvaCut.root");

      n_files=7;
      png_string = "withSTbins_highlowpT_withMvaCut";
  }
      if(Model.Contains("T5bbbbZg") && which_plot==10 ){
	legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) + ST bins "};
	legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40","SP bins: p_{T}^{#gamma}>100","SP bins: p_{T}^{#gamma}>40 & BDT score>0","SP bins: p_{T}^{#gamma}>100 & BDT score>0", "p_{T}^{#gamma}>100", "p_{T}^{#gamma}>100 & BDT score>0"};
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_withMvaCut.root");
      f[3] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100_withMvaCut.root");
      f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100.root");
      f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100_withMvaCut.root");

      n_files=6;
      png_string = "withSTbins_OnlyhighlowpT_withMvaCut";
  }


      if(Model.Contains("T5bbbbZg") && which_plot==6 ){
	legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};// : SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss})"};
	legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>100 & BDT score>0"};//,"p_{T}^{#gamma}>40 & BDT score>0", "p_{T}^{#gamma}>100 & BDT score>0","Add bins in ST","Add bins in ST & BDT score>0" };
	f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
	f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
	f[2] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200.root");
	f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200_phopt_100.root");
	f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200_withMvaCut.root");
	f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200_phopt_100_withMvaCut.root");
      
	n_files=6;
	png_string = "withSTandPtbins_withMvaCut";
      }
      if(Model.Contains("T5bbbbZg") && which_plot==8 ){
        legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : p_{T}^{miss}>300"};
        legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>100 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300_phopt_100.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300_withMvaCut.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300_phopt_100_withMvaCut.root");

        n_files=6;
        png_string = "withSTandPtbins_withMvaCut_highMET";
      }

      if(Model.Contains("T5bbbbZg") && which_plot==7 ){
      legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
      legend_texts_v1 = {"SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) : p_{T}^{#gamma}>40", "SP bins(N_{jets},N_{jets}^{b},p_{T}^{miss}) : p_{T}^{#gamma}>100", "Bin in BDT score"};
      f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
      f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
      f[2] = new TFile("Excl_out_T5bbbbZg_combine_v4sbins_bdt_MET_200.root");
    n_files=3;
    png_string = "pt40vspt100_andBDTbins";

  }
      if(Model.Contains("T5bbbbZg") && which_plot==9 ){
	legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : p_{T}^{miss}>300"};
        legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300 & BDT score>0","SP+STbins:p_{T}^{#gamma}>40","SP+STbins:p_{T}^{#gamma}>40 & BDT score>0"};//,"Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>100 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300_withMvaCut.root");
	f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200.root");
	f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_withMvaCut.root");

        n_files=6;
        png_string = "withSTandPtbins_vsSt+SPbins_withMvaCut";
      }

        if(Model.Contains("T5bbbbZg") && which_plot==11){
	  legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} : p_{T}^{miss}>300"};
	  legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300 & BDT score>0","SP+STbins:p_{T}^{#gamma}>100","SP+STbins:p_{T}^{#gamma}>100 & BDT score>0"};//,"Bins(p_{T}^{#gamma},ST):
	  f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
	f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_300_withMvaCut.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v4_MET_200_phopt_100_withMvaCut.root");

        n_files=6;
        png_string = "withSTandPtbins_vsSt+SPbins_withMvaCut_highpT";
      }

	if(Model.Contains("T5qqqqHg") && which_plot==1){
	   legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}"};
	   legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300", "Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300 & BDT score>0","bins in BDT score"};
	   f[0] = new TFile("Excl_out_T5qqqqHg_combine_MET_200_phopt_40.root");
	   f[1] = new TFile("Excl_out_T5qqqqHg_combine_MET_200_phopt_100.root");
	   f[2] = new TFile("Excl_out_T5qqqqHg_combine_T5qqqqHg_v5_MET_300.root");
	   f[3] = new TFile("Excl_out_T5qqqqHg_combine_T5qqqqHg_v5_MET_300_withMvaCut.root");
	   f[4]= new TFile("Excl_out_T5qqqqHg_combine_T5qqqqHg_v4sbins_bdt_MET_200.root");
	   n_files=5;
	   png_string = "DifferentSRbins_firstTest";

	}
	
	if(Model.Contains("T5ttttZg") && which_plot==1){
           legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
           legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300", "Bins(p_{T}{#gamma},ST):p_{T}^{#gamma}>40 & p_{T}^{miss}>300 & BDT score>0","bins in BDT score"};
           f[0] = new TFile("Excl_out_T5ttttZg_combine_MET_200_phopt_40.root");
           f[1] = new TFile("Excl_out_T5ttttZg_combine_MET_200_phopt_100.root");
           f[2] = new TFile("Excl_out_T5ttttZg_combine_T5ttttZg_v5_MET_300.root");
           f[3] = new TFile("Excl_out_T5ttttZg_combine_T5ttttZg_v5_MET_300_withMvaCut.root");
           //f[4]= new TFile("Excl_out_T5ttttZg_combine_T5ttttZg_v4sbins_bdt_MET_200.root");
           n_files=4;
	   png_string = "DifferentSRbins_firstTest";

        }

	if(Model.Contains("T5bbbbZg") && which_plot==12){
	  legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
	  legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_T5bbbbZg_v5_MET_200_withMvaCut.root");
	f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200.root");
	f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200_withMvaCut.root");
        n_files=6;
        png_string = "withSTandPt_bjetsbins_withMvaCut";
	}
	if(Model.Contains("T5bbbbZg") && which_plot==13){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
	f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200_withMvaCut.root");
	f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200_withMvaCut.root");
	
        n_files=6;
        png_string = "withSTandMETbins_withMvaCut";
	}
	if(Model.Contains("T5bbbbZg") && which_plot==13){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
	f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200_withMvaCut.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200_withMvaCut.root");
	linecolorlist = {kBlue,kBlack,kRed,kMagenta, kRed+2,9,kCyan+2};

        n_files=6;
        png_string = "withSTandMETbins_withMvaCut";
	}

	if(Model.Contains("T5bbbbZg") && which_plot==14){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{#gamma},ST, N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST):p_{T}^{#gamma}>40 & BDT score>0","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
	f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");      
	f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v1_MET_200_withMvaCut.root");
	f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v2_MET_200_withMvaCut.root");
	f[6] = new  TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
        f[7] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");

        linecolorlist = {kBlue,kBlack,kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2};
	
        n_files=8;
        png_string = "withSTandMET_nbjetsbins_withMvaCut";
        }

	if(Model.Contains("T5bbbbZg") && which_plot==15){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},N_{jets}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v4_MET_200.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v4_MET_200_withMvaCut.root");
        linecolorlist = {kBlue,kBlack, kGreen+2, kYellow+2, kCyan, kViolet,kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2};

        n_files=6;
        png_string = "withSTNjets_nbjetsbins_withMvaCut";
        }
	
	if(Model.Contains("T5bbbbZg") && which_plot==16){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},N_{jets}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets}, p_{T}^{miss}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
	f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
	f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v4_MET_200.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v4_MET_200_withMvaCut.root");
	f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v5_MET_200.root");
        f[6] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v5_MET_200_withMvaCut.root");

	linecolorlist = {kBlue,kBlack, kYellow+2, kCyan, kViolet,kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2};

        n_files=7;
        png_string = "withSTNjets_nbjet_METbins_withMvaCut";
	}

	if(Model.Contains("T5bbbbZg") && which_plot==17){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets}, p_{T}^{miss}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{#gamma}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{#gamma}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v5_MET_200.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v5_MET_200_withMvaCut.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v6_MET_200.root");
        f[6] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v6_MET_200_withMvaCut.root");

        linecolorlist = {kBlue,kBlack, kYellow+2, kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2,kCyan, kViolet,};

        n_files=7;
        png_string = "withSTNjets_nbjet_phoptbins_withMvaCut";
        }

	   if(Model.Contains("T5bbbbZg") && which_plot==18){
          legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
          legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},p_{T}^{#gamma},p_{T}^{miss}):p_{T}^{#gamma}>40","Bins(ST,N_{jets}^{b},p_{T}^{#gamma},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0"};
        f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
        f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
        f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200.root");
        f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
        f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v7_MET_200.root");
        f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");

        linecolorlist = {kBlue,kBlack, kYellow+2, kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2,kCyan, kViolet,};

        n_files=6;
        png_string = "withSTnbjet_phoptandMETbins_withMvaCut";
        }

	   if(Model.Contains("T5bbbbZg") && which_plot==19){
	     legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}"};
	     legend_texts_v1 = {"SP bins: p_{T}^{#gamma}>40", "SP bins: p_{T}^{#gamma}>100","Bins(p_{T}^{miss},ST,N_{jets}^{b}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},N_{jets},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0","Bins(ST,N_{jets}^{b},p_{T}^{#gamma},p_{T}^{miss}):p_{T}^{#gamma}>40 & BDT score>0","SP bins (N_{jets}^{b},N_{jets},p_{T}^{miss},ST,p_{T}^{#gamma})","SP bins (N_{jets}^{b},N_{jets},p_{T}^{miss},ST,p_{T}^{#gamma}): BDT score>0" };
	     f[0] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_40.root");
	     f[1] = new TFile("Excl_out_T5bbbbZg_combine_MET_200_phopt_100.root");
	     f[2] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v3_MET_200_withMvaCut.root");
	     f[3] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v5_MET_200_withMvaCut.root");
	     f[4] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v7_MET_200_withMvaCut.root");
	     f[5] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v8_MET_200.root");
	     f[6] = new TFile("Excl_out_T5bbbbZg_combine_h_Sbins_LL_newSbins_v8_MET_200_withMvaCut.root");
	     
	     linecolorlist = {kBlue,kBlack, kYellow+2, kRed,kMagenta, kRed+2,9, kGreen+2, kYellow+2,kCyan, kViolet,};

	     n_files=7;
        png_string = "withMETNjets_nbjet_phoptandSTbins_withMvaCut";
        }



	
  int  n_var = varName.size();
  int   n_cut =baseline.size();
 // for(int i_file=0; i_file<n_files;i_file++)
 //   {      
 //     for(int i_var=0; i_var<n_var;i_var++)
 // 	{
 // 	  vector<TGraph*> hist_list_Njets;
 // 	  vector<TGraph*> hist_list_Bjets;
 for(int i_cut=0; i_cut<n_cut;i_cut++){
   // vector<TGraph*> hist_list_Njets;
   // vector<TGraph*> hist_list_Bjets;
   for(int i_var=0; i_var<n_var;i_var++)
     {
       vector<TGraph*> hist_list_Njets;
       vector<TGraph*> hist_list_Bjets; 
       for(int i_file=0; i_file<n_files;i_file++)
	 {
	   if(i_file==0 && !Model.Contains("T5gg"))
             {
	       sprintf(hist_name,"sigmatheory");//,varName[0].c_str());
	        TGraph* h_resp = (TGraph*)f[i_file]->Get(hist_name);
		h_resp->GetXaxis()->SetRangeUser(xmin[i_var],xmax[i_var]);
		h_resp->GetYaxis()->SetRangeUser(ymin[i_var],ymax[i_var]);
		hist_list_Njets.push_back(h_resp);
             }
	   
	   sprintf(hist_name,"%s",varName[i_var].c_str());//,baseline[i_cut].c_str());
	   cout<<hist_name<<"\t"<<i_cut<<"\t"<<i_var<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
	   TGraph* h_resp = (TGraph*)f[i_file]->Get(hist_name);
	   cout<<"resp "<<h_resp->Integral()<<"\t"<<rebin[i_var]<<"\t"<<xmin[i_var]<<"\t"<<xmax[i_var]<<endl;
	   h_resp->GetXaxis()->SetRangeUser(xmin[i_var],xmax[i_var]);
	   h_resp->GetYaxis()->SetRangeUser(ymin[i_var],ymax[i_var]);
	   hist_list_Njets.push_back(h_resp); 
	   double factor=1.0;
	   
	 }	
	   cout<<" hist_list_Njets.size() "<<hist_list_Njets.size()<<"\t "<<"baseline.size()  "<<baseline.size()<<endl;//hist_list_Bjets.size() "<<hist_list_Bjets.size()<<endl;
	  float energy=energyy[0];
	  int xrange=0.0;
	  sprintf(full_path,"%s/overlayLimits_%s_%s_%s",pathname.c_str(),model.c_str(),varName[i_var].c_str(),png_string.c_str());
	  cout<<"varName "<< varName[i_var].c_str() <<"\t"<<i_var<<"\t"<<xlabel[i_var].c_str()<<"\t"<<rebin[i_var]<<"\t"<<ymin[i_var]<<"\t"<<ymax[i_var]<<"\t"<<xmin[i_var]<<"\t"<<xmax[i_var]<<"\t"<<legend_title[0].c_str()<<"\t"<<legend_texts_v1.size()<<endl;
	  // if(i_var==2 || i_var==3)
	  //   generate_1Dplot(hist_list_Njets,full_path,xlabel[i_var].c_str(),"Entries",energy,rebin[i_var],ymin[i_var],ymax[i_var],xmin[i_var],xmax[i_var],legend_texts[i_cut].c_str(),false,true,false,true,legend_texts_v1,legend_title[0].c_str(), Model);
	  // else
	  cout<<"temp"<<n_var<<"\t"<<n_cut<<"\t"<<xlabel[i_var].c_str()<<"\t"<<ylabel[i_var].c_str()<<"\t"<<energy<<"\t"<<legend_texts[i_cut].c_str()<<endl;
	  generate_1Dplot(hist_list_Njets,full_path,xlabel[i_var].c_str(),ylabel[i_var].c_str(),energy,rebin[i_var],ymin[i_var],ymax[i_var],xmin[i_var],xmax[i_var],legend_texts[i_cut].c_str(),false,true,false,true,legend_texts_v1,legend_title[0].c_str(), Model, linecolorlist);
	  cout<<"\temp2"<<"what the fuck "<<endl;
	}

    }
}

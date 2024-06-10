const int n_pl = 4;
bool logx = false;
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[6] ={"t #bar{t} + #gamma","t #bar{t} + jets", "Z(#nu#nu)+ #gamma", "W(l#nu) + #gamma","W(l#nu) + jets","#gamma + jets"};//,"Failed acceptance","1e CR"};// {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[4] = {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,3,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,2,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};                                                                               
int line_color[9] = {kBlue,9,kCyan+3,kGreen,kRed,kMagenta, kRed+2,9,kCyan+2};//,45,kMagenta,kGray+1,kRed,kBlue+2,kMagenta,kCyan};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
//int line_color[9] = {kMagenta+2, kGray+2, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
vector<int> col={kGray,kTeal+9,kRed,kOrange,kCyan-1,kCyan,kBlue};//kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={1001,1001,1001,1001,1001,1001};
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

void generate_1Dplot(vector<TH1D*> hist, char const *tag_name="",char const *xlabel="",char const *ylabel="", float energy=-1, int rebin=-1,double ymin=0,double ymax=0,int xmin=-1,int xmax=-1, char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true,  vector<string> legend_texts={"nil"}, char const *legend_title="", TString model=""){  

  cout<<" inside generate 1D plot "<<"\t"<<legend_title<<"\t"<<endl;
  TCanvas *canvas_n1 =      new TCanvas(tag_name, tag_name,1500,1000);
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
  legend = new TLegend(0.15,0.745,0.5,0.945);  
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(2);
  char* lhead = new char[100];
 cout<<"before legend fixing "<<endl; 
  sprintf(lhead,"%s ",leg_head);
 auto  legend1 = new TLegend(0.7,0.7,0.95,0.945);
  legend1->SetTextSize(0.03);
  legend1->SetLineColor(kWhite);
  legend1->SetNColumns(2);

  legend1->SetHeader(lhead);

  legend->SetHeader(legend_title);
  legend->SetLineColor(kWhite);
  cout<<"after legend fixing "<<endl;
  TLegendEntry* leg_entry[20];
  float x_label_size = 0.045;
  //  double ymin = 100000.0;
  //double ymax = 0.0;
  double xrange = xmax;
  // float energy = energyy;
  vector<TH1D*> hist_list_temp;
  cout<<" hist.size() = "<<hist.size()<<endl;
  for(int i =0;i<(int)hist.size(); i ++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    //    hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    //    normalize = true;
    if(normalize) {
      hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
     hist.at(i)->GetXaxis()->SetTitle(xlabel);
    //   hist.at(i)->GetXaxis()->SetRangeUser(xmin,xrange+4);
     //     cout<<i<<"\t"<<"oinside loop "<<endl;
    hist.at(i)->SetLineWidth(line_width[i]);
    if(i>5){
    hist.at(i)->SetLineStyle(line_style[i-6]);
  
    hist.at(i)->SetLineColor(line_color[i-6]);
    }
    else
      {
	hist.at(i)->SetLineColor(col[i]);
	hist.at(i)->SetFillColor(col[i]);
      }
    //    cout<<i<<"\t"<<"oinside loop "<<endl;

    hist.at(i)->SetTitle(" ");
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    //    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleSize(0.06);
    hist.at(i)->GetYaxis()->SetLabelSize(0.06);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.);
     decorate(hist.at(i),i, 0);
    hist.at(i)->SetMarkerSize(0.8);
    hist.at(i)->SetMarkerStyle(20);
    //    hist.at(i)->SetMarkerColor(line_color[i]);
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
    legName.push_back(hist.at(i)->GetName());
    if(i>5){
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i-6].c_str(),"ep");
    }
    else
      leg_entry[i] = legend1->AddEntry(hist.at(i),legend_text[i],"f"); 
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(i==0 || i==5)
      leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor()+1);
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
    if(!normalize) {
      if(model.Contains("TChiWG") || model.Contains("T6ttZg") || model.Contains("TChiNG") || model.Contains("WGJets") || model.Contains("WlnuJets") || model.Contains("ttbarG")|| model.Contains("ttbarJets") || model.Contains("GJets") || model.Contains("ZnunuGJets") ) hist.at(i)->GetYaxis()->SetRangeUser(0.01,1000*ymax);
      else       hist.at(i)->GetYaxis()->SetRangeUser(0.001,100*ymax);    }
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.00001,5.0);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
      }
    //    p1->SetGrid();
    // if(i==6)
    //   hist.at(i)->Draw("hist ");
    // else if (i>6)
    //    hist.at(i)->Draw("hist sames ");
    // else
      if(i<=5)
	hs_var->Add(hist.at(i));
    hs_var->SetMinimum(0.01);
    hs_var->SetMaximum(ymax*1000);
    //gPad->SetLogu
    //    cout<<"i Alps "<<i<<endl;
    // if(i>=0) hist.at(i)->Draw("hist ");
    // else hist.at(i)->Draw("hist sames");
	
  }
//  hs_var->SetMinimum(0.0);
//  hs_var->SetMaximum(ymax+0.5);


  hs_var->Draw("BAR HIST");
//   hs_var->Draw("HIST");
  hs_var->GetXaxis()->SetTitleOffset(1.0);
  gPad->Modified(); gPad->Update();
  hs_var->GetXaxis()->SetTitle(xlabel);
  hs_var->GetXaxis()->SetRangeUser(xmin,xmax+0.01*xmax);
  //  hs_var->GetYaxis()->SetTitleSize(
  hs_var->GetYaxis()->SetTitle("Entries");//hs_var->GetYaxis()->SetTitle("Events");
  hs_var->SetTitle(0);
  hs_var->GetYaxis()->SetTitleOffset(1.2);
  hs_var->GetXaxis()->SetTitleSize(00.05);
  hs_var->GetXaxis()->SetLabelSize(0.04);
  hs_var->GetYaxis()->SetLabelSize(0.04);
  hs_var->GetYaxis()->SetTitleSize(00.055);
  hs_var->GetYaxis()->SetTitleOffset(1.0);

  for(int j=6;j<hist.size();j++){
    hist.at(j)->Draw("Hist sames");
  }
  // // hs_var->SetMinimum(0.0);
  // // hs_var->SetMaximum(1.5);
  
  
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

  legend1->Draw();


  if(log_flag) {
      gPad->SetLogy();
    }
  // if(logx)
  //   gPad->SetLogx();
  gPad->Update(); 
  TLatex* textOnTop = new TLatex();
  //new
  textOnTop->SetTextSize(0.04);
  textOnTop->DrawLatexNDC(0.13,0.965,"CMS #it{#bf{Simulation Preliminary}}");

  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.75,0.965,en_lat);


  gPad->Modified();
  TString Tag = tag_name;
  if(Tag.Contains("Sbins")){
      
    TLine *line1V7=new TLine( 8.0,0.01,  8.0,10000);
    TLine *line2V7=new TLine(14.0,0.01, 14.0,10000);
    TLine *line3V7=new TLine(19.0,0.01, 19.0,10000);
    TLine *line4V7=new TLine(24.0,0.01, 24.0,10000);
  TLine *line5V7=new TLine(29.0,0.1, 29.0,10000);
  line1V7->Draw();      line2V7->Draw();  line3V7->Draw();
    line4V7->Draw();      line5V7->Draw();
  TArrow *arrow1 = new TArrow( 1.0,600, 8.0,600,0.01,"<|>");
    TArrow *arrow2 = new TArrow( 8.0,600,14.0,600,0.01,"<|>");
    TArrow *arrow3 = new TArrow(14.0,600,19.0,600,0.01,"<|>");
    TArrow *arrow4 = new TArrow(19.0,600, 24.0,600,0.01,"<|>");
    TArrow *arrow5 = new TArrow(24.0,600, 29.0,600,0.01,"<|>");
    TArrow *arrow6 = new TArrow(29.0,600, 34.0,600,0.01,"<|>");

    arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
    arrow4->Draw(); arrow5->Draw(); arrow6->Draw();

    TLatex Tl;
    Tl.SetTextSize(0.04);
    Tl.DrawLatex(3.5,10000,"N^{ 0}_{ 2-4}");
    Tl.DrawLatex(9.5,10000,"N^{ 0}_{ 5-6}");
    Tl.DrawLatex(15.5,10000,"N^{ 0}_{ #geq7}");
    Tl.DrawLatex(19.5,10000,"N^{ #geq1}_{ 2-4}");
    Tl.DrawLatex(25.5,10000,"N^{ #geq1}_{ 5-6}");
    Tl.DrawLatex(30.5,10000,"N^{ #geq1}_{ #geq7}");

    }
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


void StackPlots_multifile(string pathname, string model, string gluino_m)
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
  TString Gluino_m = gluino_m.c_str();
  TString Model = model.c_str();
  baseline = {"overlap_removal","MET_200","MET_300","MET_300_phopt_100","MET_200_phopt_100"}; //lepton_veto
  vector <string > varName ={"h_St","h_HT","h_NhadJets","h_NBJets","h_MET","h_PhoPt","h_Photon_Eta","h_Photon_Phi","h_MET_Phi","FR_nbtagBins","h_Mt_phoMET","h_dPhi_phoMet","h_Nphotons","h_dPhi_METLeadJet","h_Sbins_LL"};
  vector <string> xlabel = {"Sum of p_{T}^{Jets} & p_{T}^{#gamma} [GeV]","HT[GeV]","N_{jets}","N_{ b-jets}","p_{T}^{miss} [GeV]","p_{T}^{#gamma} [GeV]","#eta^{#gamma}","#phi^{#gamma}","#phi^{MET}","Bin No.","M_{T}^{miss & #gamma} [GeV]","dPhi(#gamma,MET)","N_{#gamma}","d#phi(p_{T}^{miss},lead Jet1)","Bin No."};

  legend_texts ={"p_{T}^{miss} [GeV] > 100","p_{T}^{miss} [GeV] > 200","p_{T}^{miss} [GeV] > 300","p_{T}^{miss} [GeV] > 300 & p_{T}^{#gamma} >100","p_{T}^{miss} [GeV] > 200 & p_{T}^{#gamma} >100"};

  vector <string>legend_title;
  vector<string> filetag= {"FullRun2"};// {"2018","2017","2016postVFP","2016preVFP","2016","FullRun2"};
  vector<float> energyy={137.19};//{ 59.74,41.53,16.5,19.5,36,137.19};
  vector <int> rebin;
  rebin={5,5,1,1,5,5,4,4,4,1,4,4,1,4,1};//,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,2,4,4,4};
  vector<double> ymin ={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};//,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  vector<double> ymax={100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000};//,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000,100000};
  vector<double> xmin ={0,0,0,0,0,0,-5,-5,-5,0,0,0,0,0,0,0};//-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0,-5,-5,0,0.9,0,0,80,0};
  vector<double> xmax={2500,2500,20,16,1000,900,5,5,5,3,800,5,5,5,40};//,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,5,5,5,1000,2,400,3,100,800};

  cout<<"different vector sizes "<<endl;
  cout<<varName.size()<<"\t"<<baseline.size()<<"\t"<<xlabel.size()<<"\t"<<rebin.size()<<"\t"<<xmax.size()<<"\t"<<xmin.size()<<"\t"<<legend_texts.size()<<endl;
  bool flag=false;
  n_files=7;
  if(Model.Contains("T5bbbbZg") && Gluino_m.Contains("2200")){
    Nlsp_m = {"10","50","100","200","400","600","1000"};
    legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)"};

    legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2200 "};
    legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};
    f[0] = new TFile("out_T5bbbbZg_2200_10.root");
    f[1] = new TFile("out_T5bbbbZg_2200_50.root");
    f[2] = new TFile("out_T5bbbbZg_2200_100.root");
    f[3] = new TFile("out_T5bbbbZg_2200_200.root");
    f[4] = new TFile("out_T5bbbbZg_2200_400.root");
    f[5] = new TFile("out_T5bbbbZg_2200_600.root");
    f[6] = new TFile("out_T5bbbbZg_2200_1000.root");
  }
  else if(Model.Contains("T5bbbbZg") && Gluino_m.Contains("2300")){
    Nlsp_m = {"10","50","100","200","400","600","1000"};
   legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)"};
   legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2300 "};
   legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};


    f[0] = new TFile("out_T5bbbbZg_2300_10.root");
    f[1] = new TFile("out_T5bbbbZg_2300_50.root");
    f[2] = new TFile("out_T5bbbbZg_2300_100.root");
    f[3] = new TFile("out_T5bbbbZg_2300_200.root");
    f[4] = new TFile("out_T5bbbbZg_2300_400.root");
    f[5] = new TFile("out_T5bbbbZg_2300_600.root");
    f[6] = new TFile("out_T5bbbbZg_2300_1000.root");
  }
  else if(Model.Contains("T5bbbbZg") && Gluino_m.Contains( "2400")){
    Nlsp_m = {"10","50","100","200","400","600","1000"};
    legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)"};
 legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2400 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};


    f[0] = new TFile("out_T5bbbbZg_2400_10.root");
    f[1] = new TFile("out_T5bbbbZg_2400_50.root");
    f[2] = new TFile("out_T5bbbbZg_2400_100.root");
    f[3] = new TFile("out_T5bbbbZg_2400_200.root");
    f[4] = new TFile("out_T5bbbbZg_2400_400.root");
    f[5] = new TFile("out_T5bbbbZg_2400_600.root");
    f[6] = new TFile("out_T5bbbbZg_2400_1000.root");
  }
  else if(Model.Contains("T5bbbbZg") && Gluino_m.Contains("2500")){
    Nlsp_m = {"10","50","100","200","400","600","1000"};
        legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 1000 GeV)"};

 legend_title={"#tilde{g} #rightarrow b #bar{b} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2500 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5bbbbZg_2500_10.root");
    f[1] = new TFile("out_T5bbbbZg_2500_50.root");
    f[2] = new TFile("out_T5bbbbZg_2500_100.root");
    f[3] = new TFile("out_T5bbbbZg_2500_200.root");
    f[4] = new TFile("out_T5bbbbZg_2500_400.root");
    f[5] = new TFile("out_T5bbbbZg_2500_600.root");
    f[6] = new TFile("out_T5bbbbZg_2500_1000.root");
  }

  else if(Model.Contains("T5qqqqHg") && Gluino_m.Contains("2500")){
     Nlsp_m = {"127","150","200","400","600","1000"};
     n_files=6;
     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 127 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 150 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

 legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}, m_{#tilde{g}} = 2500 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 127 GeV","m_{#tilde{#chi}_{1}^{0}} = 150 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};
     
    f[0] = new TFile("out_T5qqqqHg_2500_127.root");
    f[1] = new TFile("out_T5qqqqHg_2500_150.root");
    f[2] = new TFile("out_T5qqqqHg_2500_200.root");
    f[3] = new TFile("out_T5qqqqHg_2500_400.root");
    f[4] = new TFile("out_T5qqqqHg_2500_600.root");
    f[5] = new TFile("out_T5qqqqHg_2500_1000.root");

  }

   else if(Model.Contains("T5qqqqHg") && Gluino_m.Contains("2400")){
     Nlsp_m = {"127","150","200","400","600","1000"};
     n_files=6;
     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 127 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 150 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}, m_{#tilde{g}} = 2400 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 127 GeV","m_{#tilde{#chi}_{1}^{0}} = 150 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5qqqqHg_2400_127.root");
    f[1] = new TFile("out_T5qqqqHg_2400_150.root");
    f[2] = new TFile("out_T5qqqqHg_2400_200.root");
    f[3] = new TFile("out_T5qqqqHg_2400_400.root");
    f[4] = new TFile("out_T5qqqqHg_2400_600.root");
    f[5] = new TFile("out_T5qqqqHg_2400_1000.root");

  }
   else if(Model.Contains("T5qqqqHg") && Gluino_m.Contains("2300")){
     Nlsp_m = {"127","150","200","400","600","1000"};
     n_files=6;
     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 127 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 150 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}, m_{#tilde{g}} = 2300 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 127 GeV","m_{#tilde{#chi}_{1}^{0}} = 150 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5qqqqHg_2300_127.root");
    f[1] = new TFile("out_T5qqqqHg_2300_150.root");
    f[2] = new TFile("out_T5qqqqHg_2300_200.root");
    f[3] = new TFile("out_T5qqqqHg_2300_400.root");
    f[4] = new TFile("out_T5qqqqHg_2300_600.root");
    f[5] = new TFile("out_T5qqqqHg_2300_1000.root");

  }

 else if(Model.Contains("T5qqqqHg") && Gluino_m.Contains("2200")){
     Nlsp_m = {"127","150","200","400","600","1000"};
     n_files=6;
     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 127 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 150 GeV)", "#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow q #bar{q} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G}, m_{#tilde{g}} = 2200 "};
 legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 127 GeV","m_{#tilde{#chi}_{1}^{0}} = 150 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5qqqqHg_2200_127.root");
    f[1] = new TFile("out_T5qqqqHg_2200_150.root");
    f[2] = new TFile("out_T5qqqqHg_2200_200.root");
    f[3] = new TFile("out_T5qqqqHg_2200_400.root");
    f[4] = new TFile("out_T5qqqqHg_2200_600.root");
    f[5] = new TFile("out_T5qqqqHg_2200_1000.root");

  }

   else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2100")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2100 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2100 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

     
    f[0] = new TFile("out_T5ttttZg_2100_10.root");
    f[1] = new TFile("out_T5ttttZg_2100_50.root");    
    f[2] = new TFile("out_T5ttttZg_2100_100.root");
    f[3] = new TFile("out_T5ttttZg_2100_200.root");
    f[4] = new TFile("out_T5ttttZg_2100_400.root");
    f[5] = new TFile("out_T5ttttZg_2100_600.root");
    f[6] = new TFile("out_T5ttttZg_2100_1000.root");    

  }

 else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2200")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2200 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2200 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5ttttZg_2200_10.root");
    f[1] = new TFile("out_T5ttttZg_2200_50.root");    
    f[2] = new TFile("out_T5ttttZg_2200_100.root");
    f[3] = new TFile("out_T5ttttZg_2200_200.root");
    f[4] = new TFile("out_T5ttttZg_2200_400.root");
    f[5] = new TFile("out_T5ttttZg_2200_600.root");
    f[6] = new TFile("out_T5ttttZg_2200_1000.root");    

  }
   else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2300")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2300 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2300 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5ttttZg_2300_10.root");
    f[1] = new TFile("out_T5ttttZg_2300_50.root");    
    f[2] = new TFile("out_T5ttttZg_2300_100.root");
    f[3] = new TFile("out_T5ttttZg_2300_200.root");
    f[4] = new TFile("out_T5ttttZg_2300_400.root");
    f[5] = new TFile("out_T5ttttZg_2300_600.root");
    f[6] = new TFile("out_T5ttttZg_2300_1000.root");    

  }
 else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2400")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2400 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2400 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5ttttZg_2400_10.root");
    f[1] = new TFile("out_T5ttttZg_2400_50.root");    
    f[2] = new TFile("out_T5ttttZg_2400_100.root");
    f[3] = new TFile("out_T5ttttZg_2400_200.root");
    f[4] = new TFile("out_T5ttttZg_2400_400.root");
    f[5] = new TFile("out_T5ttttZg_2400_600.root");
    f[6] = new TFile("out_T5ttttZg_2400_1000.root");    

  }

   else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2500")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2500 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2500 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5ttttZg_2500_10.root");
    f[1] = new TFile("out_T5ttttZg_2500_50.root");    
    f[2] = new TFile("out_T5ttttZg_2500_100.root");
    f[3] = new TFile("out_T5ttttZg_2500_200.root");
    f[4] = new TFile("out_T5ttttZg_2500_400.root");
    f[5] = new TFile("out_T5ttttZg_2500_600.root");
    f[6] = new TFile("out_T5ttttZg_2500_1000.root");    

  }
   else if(Model.Contains("T5ttttZg") && Gluino_m.Contains( "2600")){
     Nlsp_m = {"10","50","100","200","400","600","1000"};
     n_files=7;
     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 50 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 100 GeV)", "#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 200 GeV)" ,"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 400 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 600 GeV)","#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/H #tilde{G} (m_{#tilde{g}} = 2600 GeV, m_{#tilde {#chi}_{1}^{0}} = 1000 GeV)"};

     legend_title={"#tilde{g} #rightarrow t #bar{t} #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G}, m_{#tilde{g}} = 2600 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T5ttttZg_2600_10.root");
    f[1] = new TFile("out_T5ttttZg_2600_50.root");    
    f[2] = new TFile("out_T5ttttZg_2600_100.root");
    f[3] = new TFile("out_T5ttttZg_2600_200.root");
    f[4] = new TFile("out_T5ttttZg_2600_400.root");
    f[5] = new TFile("out_T5ttttZg_2600_600.root");
    f[6] = new TFile("out_T5ttttZg_2600_1000.root");    

  }

   else if(Model.Contains("T6ttZg") && Gluino_m.Contains( "1300")){
     Nlsp_m = {"10","50","100","200","400","600"};
     n_files=6;
       legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)","#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)"};
       legend_title={"#tilde{t} #rightarrow t#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1300 "};
       legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

     f[0] = new TFile("out_T6ttZg_1300_10.root");
     f[1] = new TFile("out_T6ttZg_1300_50.root");    
     f[2] = new TFile("out_T6ttZg_1300_100.root");
     f[3] = new TFile("out_T6ttZg_1300_200.root");
     f[4] = new TFile("out_T6ttZg_1300_400.root");
     f[5] = new TFile("out_T6ttZg_1300_600.root");

  }
  else if(Model.Contains("T6ttZg") && Gluino_m.Contains( "1400")){
     Nlsp_m = {"10","50","100","200","400","600"};
     n_files=6;
     legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)", " #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)"};

     legend_title={"#tilde{t} #rightarrow t#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1400 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};


    f[0] = new TFile("out_T6ttZg_1400_10.root");
    f[1] = new TFile("out_T6ttZg_1400_50.root");    
    f[2] = new TFile("out_T6ttZg_1400_100.root");
    f[3] = new TFile("out_T6ttZg_1400_200.root");
    f[4] = new TFile("out_T6ttZg_1400_400.root");
    f[5] = new TFile("out_T6ttZg_1400_600.root");

  }

  else if(Model.Contains("T6ttZg") && Gluino_m.Contains( "1600")){
     Nlsp_m = {"10","50","100","200","400","600"};
     n_files=6;
     legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)", " #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)"};

     legend_title={"#tilde{t} #rightarrow t#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1600 "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T6ttZg_1600_10.root");
    f[1] = new TFile("out_T6ttZg_1600_50.root");    
    f[2] = new TFile("out_T6ttZg_1600_100.root");
    f[3] = new TFile("out_T6ttZg_1600_200.root");
    f[4] = new TFile("out_T6ttZg_1600_400.root");
    f[5] = new TFile("out_T6ttZg_1600_600.root");

  }

  else if(Model.Contains("T6ttZg") && Gluino_m.Contains( "1800")){
     Nlsp_m = {"10","50","100","200","400","600"};
     n_files=6;
     legend_title={"#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 10 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 200 GeV)", " #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 400 GeV)", "#tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)"};
   legend_title={"#tilde{t} #rightarrow t#tilde{#chi}_{1}^{0}, #tilde{#chi}_{1}^{0} #rightarrow #gamma/Z #tilde{G} (m_{#tilde{g}} = 1800 "};
   legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 10 GeV","m_{#tilde{#chi}_{1}^{0}} = 50 GeV","m_{#tilde{#chi}_{1}^{0}} = 100 GeV","m_{#tilde{#chi}_{1}^{0}} = 200 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_T6ttZg_1800_10.root");
    f[1] = new TFile("out_T6ttZg_1800_50.root");    
    f[2] = new TFile("out_T6ttZg_1800_100.root");
    f[3] = new TFile("out_T6ttZg_1800_200.root");
    f[4] = new TFile("out_T6ttZg_1800_400.root");
    f[5] = new TFile("out_T6ttZg_1800_600.root");

  }

  //electroweakino production
else if(Model.Contains("TChiNG")){
     Nlsp_m = {"300","400","500","700"};
     n_files=4;
     legend_title={"TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 300 GeV","TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 400 GeV","TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 500 GeV", "TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 700 GeV"};

     legend_title={"TChiNg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 300 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 500 GeV","m_{#tilde{#chi}_{1}^{0}} = 700 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_TChiNG_0_300.root");
    f[1] = new TFile("out_TChiNG_0_400.root");    
    f[2] = new TFile("out_TChiNG_0_500.root");
    f[3] = new TFile("out_TChiNG_0_700.root");
    
  }

  else if(Model.Contains("TChiWG")){
     Nlsp_m = {"300","400","500","700"};
     n_files=4;
     legend_title={"TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 300 GeV","TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 400 GeV","TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 500 GeV", "TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} = 700 GeV"};

     legend_title={"TChiWg , M#tilde{#chi}_{1}^{0} = M#tilde{#chi}_{1}^{+/-} "};
     legend_texts_v1 = {"m_{#tilde{#chi}_{1}^{0}} = 300 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 500 GeV","m_{#tilde{#chi}_{1}^{0}} = 700 GeV","m_{#tilde{#chi}_{1}^{0}} = 400 GeV","m_{#tilde{#chi}_{1}^{0}} = 600 GeV","m_{#tilde{#chi}_{1}^{0}} = 1000 GeV"};

    f[0] = new TFile("out_TChiWG_0_300.root");
    f[1] = new TFile("out_TChiWG_0_400.root");    
    f[2] = new TFile("out_TChiWG_0_500.root");
    f[3] = new TFile("out_TChiWG_0_700.root");    
  }

  else if(Model.Contains("ttbarG")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};//{"300","400","500","700"};
    n_files=5;
    legend_title={"t #bar{t} + #gamma Jets : 2018","t #bar{t} + #gamma Jets : 2017","t #bar{t} + #gamma Jets : 2016postVFP","t #bar{t} + #gamma Jets : 2016preVFP","t #bar{t} + #gamma Jets : Run2"};

    f[0] = new TFile("Summer20UL18_TTGJets_inc_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_TTGJets_inc_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_TTGJets_inc_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_TTGJets_inc_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_TTGJets_inc_PhoIdloose_phopt40_BL.root");
  }

  else if(Model.Contains("ttbarJets")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};//{"300","400","500","700"};                                                                                           
    n_files=5;
    legend_title={"t #bar{t} + Jets : 2018","t #bar{t} + Jets : 2017","t #bar{t} + Jets : 2016postVFP","t #bar{t} + Jets : 2016preVFP","t #bar{t} + Jets : Run2"};

    f[0] = new TFile("Summer20UL18_TTJets_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_TTJets_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_TTJets_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_TTJets_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_TTJets_PhoIdloose_phopt40_BL.root");
  }

  else if(Model.Contains("WlnuJets")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};
    n_files=5;
    legend_title={"W(l#Nu) + Jets : 2018","W(l#Nu) + Jets : 2017","W(l#Nu) + Jets : 2016postVFP","W(l#Nu) + Jets : 2016preVFP","W(l#Nu) + Jets : Run2"};

    f[0] = new TFile("Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
  }

  else if(Model.Contains("WGJets")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};
    n_files=5;
    legend_title={"W(l#Nu) + #gamma Jets : 2018","W(l#Nu) + #gamma Jets : 2017","W(l#Nu) + #gamma Jets : 2016postVFP","W(l#Nu) + #gamma Jets : 2016preVFP","W(l#Nu) + #gamma Jets : Run2"};

    f[0] = new TFile("Summer20UL18_WGJets_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_WGJets_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_WGJets_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_WGJets_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_WGJets_PhoIdloose_phopt40_BL.root");
  }

  else if(Model.Contains("GJets")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};
    n_files=5;
    legend_title={"#gamma + Jets : 2018","#gamma + Jets : 2017","#gamma + Jets : 2016postVFP","#gamma + Jets : 2016preVFP","#gamma + Jets : Run2"};

    f[0] = new TFile("Summer20UL18_GJets_QCD_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_GJets_QCD_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_GJets_QCD_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_GJets_QCD_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_GJets_QCD_PhoIdloose_phopt40_BL.root");
  }

  else if(Model.Contains("ZnunuGJets")){
    Nlsp_m = {"59.74","41.53","16.5","19.5","137.19"};
    n_files=5;
    legend_title={"Z(#nu#nu) + #gamma Jets : 2018","Z(#nu#nu) + #gamma Jets : 2017","Z(#nu#nu) + #gamma Jets : 2016postVFP","Z(#nu#nu) + #gamma Jets : 2016preVFP","Z(#nu#nu) + #gamma Jets : Run2"};

    f[0] = new TFile("Summer20UL18_ZNuNu_PhoIdloose_phopt40_BL.root");
    f[1] = new TFile("Summer20UL17_ZNuNu_PhoIdloose_phopt40_BL.root");
    f[2] = new TFile("Summer20UL16_ZNuNu_PhoIdloose_phopt40_BL.root");
    f[3] = new TFile("Summer20UL16APV_ZNuNu_PhoIdloose_phopt40_BL.root");
    f[4] = new TFile("FullRun2_ZNuNu_PhoIdloose_phopt40_BL.root");
  }
  f1[0] = new TFile("FullRun2_TTGJets_inc_PhoIdloose_phopt40_BL.root");
  f1[4] = new TFile("FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_BL.root");
  f1[3] =  new TFile("FullRun2_WGJets_PhoIdloose_phopt40_BL.root");
  f1[1] =  new TFile("FullRun2_TTJets_PhoIdloose_phopt40_BL.root");
  f1[2] = new TFile("FullRun2_ZNuNu_PhoIdloose_phopt40_BL.root");
  f1[5] = new TFile("FullRun2_GJets_QCD_PhoIdloose_phopt40_BL.root");
int  i =6; 
int  n_var = varName.size();
 int   n_cut =baseline.size();
 int n_final = i+n_files;
 // for(int i_file=0; i_file<n_files;i_file++)
 //   {      
 //     for(int i_var=0; i_var<n_var;i_var++)
 // 	{
 // 	  vector<TH1D*> hist_list_Njets;
 // 	  vector<TH1D*> hist_list_Bjets;
 for(int i_cut=0; i_cut<n_cut;i_cut++){
   // vector<TH1D*> hist_list_Njets;
   // vector<TH1D*> hist_list_Bjets;
   for(int i_var=0; i_var<n_var;i_var++)
         {
	   vector<TH1D*> hist_list_Njets;
	     vector<TH1D*> hist_list_Bjets; 
	   for(int i_file=0; i_file<n_final;i_file++)
	     {
	       //	       cout<<i_file<<"\t"<<n_final<<endl;
	       sprintf(hist_name,"%s_%s",varName[i_var].c_str(),baseline[i_cut].c_str());
	       cout<<hist_name<<"\t"<<i_cut<<"\t"<<i_var<<"\t"<<i_file<<"\t"<<i_file-6<<endl;
	       TString temp = hist_name;
	       if(temp.Contains("Sbins") && temp.Contains("overlap")) sprintf(hist_name,"h_Sbins_LL_v1_%s",baseline[i_cut].c_str());
	       cout<< "updated "<<hist_name<<endl;
	       TH1D* h_resp = (TH1D*)f1[0]->Get(hist_name);
	       if(i_file<6){
		  h_resp = (TH1D*)f1[i_file]->Get(hist_name);
		  cout<<f1[i_file]->GetName()<<endl;
	       }
	       else
		 {
		   cout<<"onsdie  "<<hist_name<<"\t"<<i_cut<<"\t"<<i_var<<"\t"<<i_file-6<<"\t"<<endl;
		  h_resp = (TH1D*)f[i_file-6]->Get(hist_name);

		 }
	       cout<<"resp "<<h_resp->Integral()<<"\t"<<rebin[i_var]<<"\t"<<xmin[i_var]<<"\t"<<xmax[i_var]<<endl;
	       // h_resp->Rebin(rebin[i_var]);
	       // if(rebin[i_var]!=1){
	       // 	 h_resp->Rebin(2);
	       // }
	        h_resp= setMyRange(h_resp,xmin[i_var],xmax[i_var]+0.01*xmax[i_var]);
	       setLastBinAsOverFlow(h_resp);
	       h_resp->Rebin(rebin[i_var]);
	       if(rebin[i_var]!=1){
                 h_resp->Rebin(2);
               }
	       //h_resp= setMyRange(h_resp,xmin[i_var],xmax[i_var]+0.01*xmax[i_var]);

	       hist_list_Njets.push_back(h_resp); 
	       double factor=1.0;
	       // h_resp->Scale(factor/h_resp->Integral());
	       // h_resp2->Scale(factor/h_resp2->Integral());
	     }	
	  cout<<" hist_list_Njets.size() "<<hist_list_Njets.size()<<"\t "<<"baseline.size()  "<<baseline.size()<<endl;//hist_list_Bjets.size() "<<hist_list_Bjets.size()<<endl;
	  float energy=energyy[0];
	  int xrange=0.0;
	  sprintf(full_path,"%s/%s_%s_%s_%s_stackedWithbkg",pathname.c_str(),varName[i_var].c_str(),baseline[i_cut].c_str(),model.c_str(),gluino_m.c_str());
	  cout<<"varName "<< varName[i_var].c_str() <<"\t"<<i_var<<"\t"<<xlabel[i_var].c_str()<<"\t"<<rebin[i_var]<<"\t"<<ymin[i_var]<<"\t"<<ymax[i_var]<<"\t"<<xmin[i_var]<<"\t"<<xmax[i_var]<<"\t"<<legend_title[0].c_str()<<"\t"<<legend_texts_v1.size()<<endl;
	  if(i_var==2 || i_var==3)
	    generate_1Dplot(hist_list_Njets,full_path,xlabel[i_var].c_str(),"Entries",energy,rebin[i_var],ymin[i_var],ymax[i_var],xmin[i_var],xmax[i_var],legend_texts[i_cut].c_str(),false,true,false,true,legend_texts_v1,legend_title[0].c_str(), Model);
	  else
	    generate_1Dplot(hist_list_Njets,full_path,xlabel[i_var].c_str(),"Entries",energy,rebin[i_var],ymin[i_var],ymax[i_var],xmin[i_var],xmax[i_var],legend_texts[i_cut].c_str(),false,true,false,true,legend_texts_v1,legend_title[0].c_str(), Model);

	}

    }
}
      // 	  sprintf(full_path,"%s/%s_%s_%s_normalize",pathname.c_str(),string_png,varName[i_cut].c_str(),filetag[i_file].c_str());

	  
      // 	  // hist_list_Njets.at(0)->GetXaxis()->SetLabelSize(0.045);
      //     // hist_list_Njets.at(1)->GetXaxis()->SetLabelSize(0.045);
      // 	  // double factor=1.0;
      // 	  // hist_list_Njets.at(0)->Scale(factor/hist_list_Njets.at(0)->Integral());
      // 	  // hist_list_Njets.at(1)->Scale(factor/hist_list_Njets.at(1)->Integral());
      // 	  // hist_list_Njets.at(0)->GetXaxis()->SetLabelSize(0.045);
      // 	  // hist_list_Njets.at(1)->GetXaxis()->SetLabelSize(0.045);
      // 	  //	  TH1D* h_temp = (TH1D*)hist_list_Bjets.at(0)->Clone();	  
      // 	  hist_list_Bjets.at(0)->Scale(factor/hist_list_Bjets.at(0)->Integral());//(factor/hist_list_Njets.at(0)->Integral())*hist_list_Njets.at(0);//->Scale(factor/hist_list_Njets.at(0)->Integral());
      // 	  //	  hist_list_Bjets.push_back(h_temp);
      // 	  //	  TH1D*	h_temp1 = (TH1D*)hist_list_Bjets.at(1)->Clone();
      // 	  hist_list_Bjets.at(1)->Scale(factor/hist_list_Bjets.at(1)->Integral());//(factor/hist_list_Njets.at(1)->Integral())*hist_list_Njets.at(1);
      //     //h_temp1->Scale(factor/hist_list_Njets.at(1)->Integral());
      // 	  //          hist_list_Bjets.push_back(h_temp1);

      // 	  hNjets_ratio =(TH1D*)hist_list_Njets.at(1)->Clone();//(TH1D*)hist_list_Bjets.at(1)->Clone();
      // 	  //hNjets_ratio = hist_list_Bjets.at(1);
      // 	  hNjets_ratio->Divide(hist_list_Bjets.at(0));
      // 	  // hist_list_Bjets.at(1)->Draw("Hist");
      // 	  // hist_list_Bjets.at(1)->SetLineColor(kRed);
      // 	  // hist_list_Bjets.at(0)->SetLineColor(kGreen);
      // 	  // hist_list_Bjets.at(1)->SetLineStyle(3);
      // 	  // hist_list_Bjets.at(0)->Draw("Hist sames");
      // 	  //	  hNjets_ratio->Draw("Hist");
      // 	  // setLastBinAsOverFlow(hNjets_ratio);
      // 	  if(i_cut==2 || i_cut==3)
      //       generate_1Dplot(hist_list_Bjets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Normalize",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);
      //     else
      //       generate_1Dplot(hist_list_Bjets,hNjets_ratio,full_path,xLabel[i_cut].c_str(),"Normalize",energy,rebin[i_cut],ymin[i_cut],ymax[i_cut],xmin[i_cut],xmax[i_cut],leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, which_Lept);

      // 	}
      // //fout->Close();
      






#include <stdio.h>
#include<fstream>
#include <vector>
#include <string>
#include <map>
#include<iostream>
using namespace std;
void twofiles_overlay (string file1, string file2, string file3, string file4)                                                                                                                 
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
  char* hist_name8 = new char[200];
  char* hist_name9 = new char[200];
  char* title=new char[2000];
  char* full_path1 = new char[2000];
  char* full_path2 = new char[2000];
  char* full_path3 = new char[2000];
  char* full_path4 = new char[2000];
  char* path2 = new char[2000];

  char* hname1 = new char[2000];
  char* hname2 = new char[200];
  char* hname3 = new char[200];
 sprintf(hname,"%s",file1.c_str());
 sprintf(hname1,"%s",file2.c_str());
 sprintf(hname2,"%s",file3.c_str());
sprintf(hname3,"%s",file4.c_str());
 sprintf(path2,"./");//,file3.c_str());                                                                                                                     
  TFile * inputfile = new TFile(hname,"READ");
  TFile * inputfile1 = new TFile(hname1,"READ");
  TFile * inputfile2 = new TFile(hname2,"READ");
  TFile * inputfile3 = new TFile(hname3,"READ");
  sprintf(hist_name,"gr1d_xseclimit");
  sprintf(hist_name1,"gr1d_xseclimit");

  TGraph* h_Ht = (TGraph*)inputfile->Get(hist_name);
  TGraph* h_wHT = (TGraph*)inputfile1->Get(hist_name);
  TGraph* h_wHT1 = (TGraph*)inputfile2->Get(hist_name);
  TGraph* h_Ht1 = (TGraph*)inputfile3->Get(hist_name);

  sprintf(full_path1,"%s/overlay_twoBL_withHT_wHT_pt100_v2.png",path2);
  h_Ht->SetLineColor(kBlack);
  h_wHT->SetLineColor(kBlue);
  h_wHT1->SetLineColor(kMagenta);
  h_Ht1->SetLineColor(kRed);
  // h_wHT->Draw("");
  // h_Ht->Draw("same");
  TMultiGraph* mg = new TMultiGraph();
 TCanvas *canvas_n1 = new TCanvas(hist_name, hist_name,600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  mg->SetTitle(" ");
  h_Ht->GetXaxis()->SetTitle("mGluino(GeV)");
   h_wHT->GetXaxis()->SetTitle("mGluino(GeV)");
  h_wHT->GetYaxis()->SetTitleOffset(1.4);
  h_Ht->GetYaxis()->SetTitle("mNLSP (GeV)");
   h_wHT->GetYaxis()->SetTitle("mNLSP (GeV)");
   h_wHT1->GetXaxis()->SetTitle("mGluino(GeV)");
  h_wHT1->GetYaxis()->SetTitleOffset(1.4);

  // TAxis *axis5=  mg->GetYaxis();
  // axis5->SetRangeUser(0.02,0.3);
  canvas_n1->SetGrid();
  auto legend = new TLegend(0.15,0.5,0.5,0.7);
  legend->SetHeader("","C");
  legend->SetLineColor(kWhite);
  // legend->AddEntry(h_Ht,"MET>200 && pho_pt>100","l");
  // legend->AddEntry(h_wHT,"MET>200 && pho_pt>40 && mvacut","l");
  TLegendEntry* l3 = legend->AddEntry(h_Ht,"bin v3","l");//#chi^{2} method","ep");  
  l3->SetTextColor(kBlack);
  TLegendEntry* l2 = legend->AddEntry(h_Ht,"bin in BDT score","l");
  l2->SetTextColor(kBlue);
  TLegendEntry* l1 = legend->AddEntry(h_Ht,"bin v3 && mva cut","l");
  l1->SetTextColor(kMagenta);
  TLegendEntry* l = legend->AddEntry(h_Ht1,"bin in BDT score && mva cut","l");
  l->SetTextColor(kRed);

  gStyle->SetLegendTextSize(0.03);
  canvas_n1->cd();
  mg->Add(h_wHT1);

  mg->Add(h_Ht);
  mg->Add(h_wHT);
  //mg->Add(h_wHT1);
  mg->Add(h_Ht1);

  mg->GetYaxis()->SetTitle("mNLSP (GeV)");
  mg->GetXaxis()->SetTitle("mGluino(GeV)");
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetXaxis()->SetRangeUser(1400,2800);
  gPad->Modified();
  gPad->Update();
  // h_Ht->SetTitle("");
  // h_wHT->SetTitle("");
  //  h_wHT->Draw("");
  mg->Draw("ALP");

 legend->Draw("sames");
 gPad->Modified();
 gPad->Update();
 canvas_n1->Modified();
 canvas_n1->cd();
 canvas_n1->SetSelected(canvas_n1);
  canvas_n1->SaveAs(full_path1);
  

}

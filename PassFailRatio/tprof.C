#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <THStack.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TStyle.h>
#include "setNCUStyle.C"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

void tprofile() {
  TFile* file  = TFile::Open("JetHT-Run2015_passfailratio.4.root");
 
  TProfile* h_data1     = (TProfile*)(file->Get("tau21vsPR_700to1500_before"));
  TProfile* h_data2     = (TProfile*)(file->Get("tau21vsPR_700to1500_before"));
 
  setNCUStyle();
  TCanvas* a = new TCanvas("a","a",0,0,750,550);
  
  h_data1->SetXTitle("h1 jet pruned mass");
  //h_data1->SetYTitle(Form("%s",title.data()));
  h_data1->GetYaxis()->SetTitleOffset(1.4);
  h_data1->SetMarkerStyle(20);
  h_data1->SetMarkerSize(0.8); 
  h_data1->SetLineColor(kGreen+2);
  h_data1->SetMarkerColor(kGreen+2);
  h_data2->SetMarkerStyle(21);
  h_data2->SetMarkerSize(0.8);
  h_data2->SetLineColor(kBlue+2);
  h_data2->SetMarkerColor(kBlue+2);
  h_data1->Draw();
  h_data2->Draw("same");

  TLegend *legendx = new TLegend(0.5,0.5,0.8,0.65);
  //legendx->SetHeader(Form("%s",header.data()));
  legendx->AddEntry(h_data2,"700 GeV < dijet mass < 1500 GeV","lp");
  legendx->AddEntry(h_data1,"dijet mass > 1500 GeV", "lp");
  legendx->SetFillStyle(0);
  legendx->SetTextSize(0.035);
  legendx->Draw();
}


double quadfunc(double*v, double* p) {
  double x= v[0];
  return p[0]+p[1]*x+p[2]*x*x;
}			    

void pfratio() {
  //TFile* file  = TFile::Open("/Users/gdeleoniii/Documents/GREG/ROOT/ssss/JetHT-Run2015_bkgest.2.1.root");
  TFile* file  = TFile::Open("JetHT-Run2015_passfailratio.7.0.2.root");
  

  TH1F* h1 = (TH1F*)(file->Get("pass"));
  TH1F* h2 = (TH1F*)(file->Get("fail"));
  TH1F* h4 = (TH1F*)(file->Get("pfratio"));

  TProfile* h3 = (TProfile*)h1->Clone("h3");
  h3->Reset();
  h3->Divide(h1,h2);
  
  TF1* f_fitprmass = new TF1("f_fitprmass", quadfunc,-75,75,3);

  setNCUStyle();
  
  TCanvas* c4 = new TCanvas("c4","",0,0,600,600);                                                                                                  
  c4->cd();
  h3->Fit("f_fitprmass","","",-75,75);
  h3->SetLineColor(kBlack);                                                                                                                  
  h3->GetYaxis()->SetTitle("R_{p/f}");                                                                                                                           
  h3->GetXaxis()->SetTitle("m_{J} - m_{H} (GeV)");                                                                                                      
  h3->GetYaxis()->SetTitleOffset(1.6);                                                                                                                          
  h3->SetStats(0);                                                                                                                                              
  h3->Draw("e1");
  //h1->SetLineColor(kOrange-4);
  //h1->Draw("e1same");
  /*
  TLegend *legend = new TLegend(0.5,0.5,0.8,0.65);
  legend->AddEntry(h1,"NCU ntuples","lp");
  legend->AddEntry(h2,"DBT ntuples", "lp");
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);
  legend->Draw();
  */
  cout<<h3->GetBinContent(1)<<" "<<h4->GetBinContent(1)<<endl;
  cout<<h3->GetBinContent(2)<<" "<<h4->GetBinContent(2)<<endl;
  cout<<h3->GetBinContent(3)<<" "<<h4->GetBinContent(3)<<endl;
  cout<<h3->GetBinContent(4)<<" "<<h4->GetBinContent(4)<<endl;
  cout<<h3->GetBinContent(5)<<" "<<h4->GetBinContent(5)<<endl;
  TCanvas* c2 = new TCanvas("c2","",0,0,600,600);
  c2->cd();
  h4->Fit("f_fitprmass","","",-75,75);
  h4->GetYaxis()->SetTitle("R_{p/f}");
  h4->GetXaxis()->SetTitle("m_{J} - m_{H} (GeV)");
  h4->GetYaxis()->SetTitleOffset(1.6);
  h4->Draw("e1");

}

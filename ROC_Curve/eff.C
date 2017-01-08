#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"

void eff() {
  Float_t DBTsig[21] ={0};
  Float_t DBTbkg[21] ={0};

  Float_t FATsig[21] ={0};
  Float_t FATbkg[21] ={0};

  Float_t SUBsig[21] ={0};
  Float_t SUBbkg[21] ={0};

  ifstream fin_dbt;
  ifstream fin_fat;
  ifstream fin_sub;
  fin_dbt.open("ROC_Curve_DBT_JetCount.dat");
  fin_fat.open("ROC_Curve_FAT_JetCount.dat");
  fin_sub.open("ROC_Curve_SUB_JetCount.dat");
  for(int i=0;i<21;i++) {
    fin_dbt >> DBTsig[i];
    fin_fat >> FATsig[i];
    fin_sub >> SUBsig[i]; 
  }
  fin_dbt.close();
  fin_fat.close();
  fin_sub.close();

  std::string HT_name[]={"HT500to700","HT700to1000","HT1000to1500","HT1500to2000","HT2000toInf"};
  ifstream fin_ht_dbt[5];
  ifstream fin_ht_fat[5];
  ifstream fin_ht_sub[5];
  Float_t ht_dbt[5][21];
  Float_t ht_fat[5][21];
  Float_t ht_sub[5][21];
  for(int i=0;i<5;i++) {
    fin_ht_dbt[i].open(Form("ROC_Curve_DBT_JetCount_QCD_%s.dat",HT_name[i].data()));
    fin_ht_fat[i].open(Form("ROC_Curve_FAT_JetCount_QCD_%s.dat",HT_name[i].data()));
    fin_ht_sub[i].open(Form("ROC_Curve_SUB_JetCount_QCD_%s.dat",HT_name[i].data()));
    for(int j=0;j<21;j++) {
      fin_ht_dbt[i] >> ht_dbt[i][j];
      fin_ht_fat[i] >> ht_fat[i][j];
      fin_ht_sub[i] >> ht_sub[i][j];
    }
    fin_ht_dbt[i].close();
    fin_ht_fat[i].close();
    fin_ht_sub[i].close();
  }
    
  for(int i=0;i<21;i++) {
    DBTbkg[i] = ht_dbt[0][i] + ht_dbt[1][i] + ht_dbt[2][i] + ht_dbt[3][i] + ht_dbt[4][i];
    FATbkg[i] = ht_fat[0][i] + ht_fat[1][i] + ht_fat[2][i] + ht_fat[3][i] + ht_fat[4][i];
    SUBbkg[i] = ht_sub[0][i] + ht_sub[1][i] + ht_sub[2][i] + ht_sub[3][i] + ht_sub[4][i];
  }

  Float_t EffDBTsig[21] = {0};
  Float_t EffFATsig[21] = {0};
  Float_t EffSUBsig[21] = {0};
  Float_t EffDBTbkg[21] = {0};
  Float_t EffFATbkg[21] = {0};
  Float_t EffSUBbkg[21] = {0};

  Float_t DBTsigDenom = DBTsig[0];
  Float_t FATsigDenom = DBTsig[0];
  Float_t SUBsigDenom = DBTsig[0];
  Float_t DBTbkgDenom = DBTbkg[0];
  Float_t FATbkgDenom = DBTbkg[0];
  Float_t SUBbkgDenom = DBTbkg[0];


  for(int i=0;i<21;i++) {
    EffDBTsig[i] = DBTsig[i]/DBTsigDenom;
    EffDBTbkg[i] = DBTbkg[i]/DBTbkgDenom;
    EffFATsig[i] = FATsig[i]/FATsigDenom;
    EffFATbkg[i] = FATbkg[i]/FATbkgDenom;
    EffSUBsig[i] = SUBsig[i]/SUBsigDenom;
    EffSUBbkg[i] = SUBbkg[i]/SUBbkgDenom;

    std::cout<<EffDBTsig[i]<<" "<<EffDBTbkg[i]<<" "<<EffFATsig[i]<<" "<<EffFATbkg[i]<<" "<<EffSUBsig[i]<<" "<<EffSUBbkg[i]<<std::endl;
  }

  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);
  c3->cd();
  c3->SetLogy();
  TGraph *roc3 = new TGraph(21,EffDBTsig,EffDBTbkg);
  TGraph *roc2 = new TGraph(21,EffFATsig,EffFATbkg);
  TGraph *roc_curve = new TGraph(21,EffSUBsig,EffSUBbkg);
  roc_curve->SetLineWidth(2);
  roc_curve->SetLineStyle(1);
  roc_curve->SetMarkerStyle(0);
  roc_curve->SetLineColor(kRed-5);
  roc_curve->SetTitle("DoubleSV");
  roc_curve->GetXaxis()->SetLimits(0,1);
  roc_curve->SetMinimum(0);
  roc_curve->SetMaximum(1);
  roc_curve->GetYaxis()->SetTitle("#varepsilon_{Bkg}");
  roc_curve->GetXaxis()->SetTitle("#varepsilon_{Sig}");
  roc_curve->Draw("acp");
  roc2->SetLineWidth(2);
  roc2->SetLineStyle(1);
  roc2->SetMarkerStyle(0);
  roc2->SetLineColor(kGreen-5);
  roc2->Draw("cp");
  roc3->SetLineWidth(2);
  roc3->SetLineStyle(1);
  roc3->SetMarkerStyle(0);
  roc3->SetLineColor(kBlue-5);
  roc3->Draw("cp");

  TLegend *legend2 = new TLegend(0.16,0.75,0.41,0.89);
  legend2->SetHeader("70 GeV<M_{p}<200 GeV, p_{T}>300 GeV");
  legend2->AddEntry(roc3,"double-b-tag","l");
  legend2->AddEntry(roc_curve,"Subjet CSV", "l");
  legend2->AddEntry(roc2,"Fatjet CSV", "l");
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.04);
  legend2->Draw();
}



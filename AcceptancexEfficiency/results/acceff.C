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
#include <TGraphErrors.h>
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

using namespace std;
void acceff() {

  std::string bulkg_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  float bulkg_mass[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};

  std::string radion_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4500"};
  float radion_mass[10]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4500};
 
  float bulkg_accxeff[8][10],bulkg_err[8][10], radion_accxeff[8][10],radion_err[8][10];
  
  for(int i=0;i<10;i++) {
    ifstream fin;
    fin.open(Form("BulkGravTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root.dat",bulkg_name[i].data()));
    float denom[10], num[8][10];
    fin >> denom[i] >> num[0][i] >> num[1][i]>> num[2][i] >> num[3][i] >> num[4][i] >> num[5][i] >> num[6][i];
    fin.close();

    bulkg_accxeff[0][i]=num[0][i]/denom[i];
    bulkg_accxeff[1][i]=num[1][i]/denom[i];
    bulkg_accxeff[2][i]=num[2][i]/denom[i];
    bulkg_accxeff[3][i]=num[3][i]/denom[i];
    bulkg_accxeff[4][i]=num[4][i]/denom[i];
    bulkg_accxeff[5][i]=num[5][i]/denom[i];
    bulkg_accxeff[6][i]=num[6][i]/denom[i];

    bulkg_err[0][i] = sqrt((1-bulkg_accxeff[0][i])*bulkg_accxeff[0][i]/denom[i]);
    bulkg_err[1][i] = sqrt((1-bulkg_accxeff[1][i])*bulkg_accxeff[1][i]/denom[i]);
    bulkg_err[2][i] = sqrt((1-bulkg_accxeff[2][i])*bulkg_accxeff[2][i]/denom[i]);
    bulkg_err[3][i] = sqrt((1-bulkg_accxeff[3][i])*bulkg_accxeff[3][i]/denom[i]);
    bulkg_err[4][i] = sqrt((1-bulkg_accxeff[4][i])*bulkg_accxeff[4][i]/denom[i]);
    bulkg_err[5][i] = sqrt((1-bulkg_accxeff[5][i])*bulkg_accxeff[5][i]/denom[i]);
    bulkg_err[6][i] = sqrt((1-bulkg_accxeff[6][i])*bulkg_accxeff[6][i]/denom[i]);

  }

  for(int i=0;i<10;i++) {
    ifstream fin;
    fin.open(Form("RadionTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root.dat",radion_name[i].data()));
    float denom[10], num[8][10];
    fin >> denom[i] >> num[0][i] >> num[1][i]>> num[2][i] >> num[3][i] >> num[4][i] >> num[5][i] >> num[6][i];
    fin.close();

    radion_accxeff[0][i]=num[0][i]/denom[i];
    radion_accxeff[1][i]=num[1][i]/denom[i];
    radion_accxeff[2][i]=num[2][i]/denom[i];
    radion_accxeff[3][i]=num[3][i]/denom[i];
    radion_accxeff[4][i]=num[4][i]/denom[i];
    radion_accxeff[5][i]=num[5][i]/denom[i];
    radion_accxeff[6][i]=num[6][i]/denom[i];

    radion_err[0][i] = sqrt((1-radion_accxeff[0][i])*radion_accxeff[0][i]/denom[i]);
    radion_err[1][i] = sqrt((1-radion_accxeff[1][i])*radion_accxeff[1][i]/denom[i]);
    radion_err[2][i] = sqrt((1-radion_accxeff[2][i])*radion_accxeff[2][i]/denom[i]);
    radion_err[3][i] = sqrt((1-radion_accxeff[3][i])*radion_accxeff[3][i]/denom[i]);
    radion_err[4][i] = sqrt((1-radion_accxeff[4][i])*radion_accxeff[4][i]/denom[i]);
    radion_err[5][i] = sqrt((1-radion_accxeff[5][i])*radion_accxeff[5][i]/denom[i]);
    radion_err[6][i] = sqrt((1-radion_accxeff[6][i])*radion_accxeff[6][i]/denom[i]);
  }

  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,800,600);
  c3->SetLogy();
  c3->cd();
  TGraphErrors *BulkG1  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[0],0,bulkg_err[0]);
  TGraphErrors *BulkG2  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[1],0,bulkg_err[1]);
  TGraphErrors *BulkG3  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[2],0,bulkg_err[2]);
  TGraphErrors *BulkG4  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[3],0,bulkg_err[3]);
  TGraphErrors *BulkG5  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[4],0,bulkg_err[4]);
  TGraphErrors *BulkG6  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[5],0,bulkg_err[5]);
  TGraphErrors *BulkG7  = new TGraphErrors(10,bulkg_mass,bulkg_accxeff[6],0,bulkg_err[6]);


  BulkG1->GetYaxis()->SetRangeUser(0.01,1.05);
  BulkG1->GetYaxis()->SetTitle("Acceptance x Efficiency");
  BulkG1->GetXaxis()->SetTitle("M_{Bulk Graviton}");
  BulkG1->GetYaxis()->SetTitleOffset(0.9);

  BulkG1->SetLineWidth(2);
  BulkG1->SetMarkerSize(1.2);
  BulkG2->SetLineWidth(2);
  BulkG2->SetMarkerSize(1.2);
  BulkG3->SetLineWidth(2);
  BulkG3->SetMarkerSize(1.2);
  BulkG4->SetLineWidth(2);
  BulkG4->SetMarkerSize(1.2);
  BulkG5->SetLineWidth(2);
  BulkG5->SetMarkerSize(1.2);
  BulkG6->SetLineWidth(2);
  BulkG6->SetMarkerSize(1.2);
  BulkG7->SetLineWidth(2);
  BulkG7->SetMarkerSize(1.2);

  BulkG1->SetMarkerColor(kRed);
  BulkG1->SetLineColor(kRed);
  BulkG2->SetMarkerColor(kBlue);
  BulkG2->SetLineColor(kBlue);
  BulkG3->SetMarkerColor(kYellow);
  BulkG3->SetLineColor(kYellow);
  BulkG4->SetMarkerColor(kViolet);
  BulkG4->SetLineColor(kViolet);
  BulkG5->SetMarkerColor(kGreen);
  BulkG5->SetLineColor(kGreen);
  BulkG6->SetMarkerColor(kBlack);
  BulkG6->SetLineColor(kBlack);
  BulkG7->SetMarkerColor(kGray);
  BulkG7->SetLineColor(kGray);

  BulkG1->Draw();                                                                                                                                                  
  BulkG2->Draw("eplsame");
  BulkG3->Draw("eplsame");
  BulkG4->Draw("eplsame");
  BulkG5->Draw("eplsame");
  BulkG6->Draw("eplsame");
  BulkG7->Draw("eplsame");

  TLegend *legend = new TLegend(0.47,0.4,0.8,0.75);
  legend->AddEntry(BulkG1,"After Trigger","lp");                                                                                                   
  legend->AddEntry(BulkG2,"After Pre-selection", "lp");
  legend->AddEntry(BulkG3,"After #Delta#eta cut","lp");
  legend->AddEntry(BulkG4,"After Reduced Mass", "lp");
  legend->AddEntry(BulkG5,"After correctedPRmass", "lp");
  legend->AddEntry(BulkG6,"After #tau_{21} cut","lp");
  legend->AddEntry(BulkG7,"After double-b tagger", "lp");
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);
  legend->Draw();

  //c3->Print("BulkGrav_AccxEff.pdf");

  TCanvas* c =  new TCanvas("c","c",0,0,800,600);
  c->SetLogy();
  c->cd();
  TGraphErrors *Radion1 = new TGraphErrors(10,radion_mass,radion_accxeff[0],0,radion_err[0]);
  TGraphErrors *Radion2 = new TGraphErrors(10,radion_mass,radion_accxeff[1],0,radion_err[1]);
  TGraphErrors *Radion3 = new TGraphErrors(10,radion_mass,radion_accxeff[2],0,radion_err[2]);
  TGraphErrors *Radion4 = new TGraphErrors(10,radion_mass,radion_accxeff[3],0,radion_err[3]);
  TGraphErrors *Radion5 = new TGraphErrors(10,radion_mass,radion_accxeff[4],0,radion_err[4]);
  TGraphErrors *Radion6 = new TGraphErrors(10,radion_mass,radion_accxeff[5],0,radion_err[5]);
  TGraphErrors *Radion7 = new TGraphErrors(10,radion_mass,radion_accxeff[6],0,radion_err[6]);

  Radion1->GetYaxis()->SetRangeUser(0.01,1.05);
  Radion1->GetYaxis()->SetTitle("Acceptance x Efficiency");                                                                                                   
  Radion1->GetXaxis()->SetTitle("M_{Radion}");
  Radion1->GetYaxis()->SetTitleOffset(0.9);

  Radion1->SetLineWidth(2);
  Radion1->SetMarkerSize(1.2);
  Radion2->SetLineWidth(2);
  Radion2->SetMarkerSize(1.2);
  Radion3->SetLineWidth(2);
  Radion3->SetMarkerSize(1.2);
  Radion4->SetLineWidth(2);
  Radion4->SetMarkerSize(1.2);
  Radion5->SetLineWidth(2);
  Radion5->SetMarkerSize(1.2);
  Radion6->SetLineWidth(2);
  Radion6->SetMarkerSize(1.2);
  Radion7->SetLineWidth(2);
  Radion7->SetMarkerSize(1.2);

  Radion1->SetMarkerColor(kRed);
  Radion1->SetLineColor(kRed);
  Radion2->SetMarkerColor(kBlue);
  Radion2->SetLineColor(kBlue);
  Radion3->SetMarkerColor(kYellow);
  Radion3->SetLineColor(kYellow);
  Radion4->SetMarkerColor(kViolet);
  Radion4->SetLineColor(kViolet);
  Radion5->SetMarkerColor(kGreen);
  Radion5->SetLineColor(kGreen);
  Radion6->SetMarkerColor(kBlack);
  Radion6->SetLineColor(kBlack);
  Radion7->SetMarkerColor(kGray);
  Radion7->SetLineColor(kGray);
  
  Radion1->SetMarkerStyle(21);
  Radion2->SetMarkerStyle(21);
  Radion3->SetMarkerStyle(21);
  Radion4->SetMarkerStyle(21);
  Radion5->SetMarkerStyle(21);
  Radion6->SetMarkerStyle(21);
  Radion7->SetMarkerStyle(21);
  
  Radion1->Draw();
  Radion2->Draw("eplsame");
  Radion3->Draw("eplsame");
  Radion4->Draw("eplsame");
  Radion5->Draw("eplsame");
  Radion6->Draw("eplsame");
  Radion7->Draw("eplsame");
  
  TLegend *leg = new TLegend(0.47,0.4,0.8,0.75);
  leg->AddEntry(Radion1,"After Trigger","lp");
  leg->AddEntry(Radion2,"After Pre-selection", "lp");
  leg->AddEntry(Radion3,"After #Delta#eta cut","lp");
  leg->AddEntry(Radion4,"After Reduced Mass", "lp");
  leg->AddEntry(Radion5,"After correctedPRmass", "lp");
  leg->AddEntry(Radion6,"After #tau_{21} cut","lp");
  leg->AddEntry(Radion7,"After double-b tagger", "lp");
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->Draw();

  //  c->Print("Radion_AccxEff.pdf");
  
  /*std::string cut[]={"Trigger","Pre-Selection","correctedPRmass","deltaeta","reducedmass","Tau21","Double-b"};
  for(int i=0;i<7;i++) {
    ofstream fout1;
    ofstream fout2;
    fout1.open(Form("AcceptancexEfficiency_BulkGraviton_%s.dat",cut[i].data()),ios::out | ios::app);
    fout2.open(Form("AcceptancexEfficiency_Radion_%s.dat",cut[i].data()),ios::out | ios::app);
    for(int j=0;j<10;j++) {

      fout1<<bulkg_mass[j]<<"= "<<bulkg_accxeff[i][j]<<" +- "<<bulkg_err[i][j]<<endl;
      fout2<<radion_mass[j]<<"= "<<radion_accxeff[i][j]<<" +- "<<radion_err[i][j]<<endl;
      
    }
    fout1.close();
    fout2.close();
    }
  */
}

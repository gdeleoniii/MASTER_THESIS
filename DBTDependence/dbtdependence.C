#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "readSample.h"
#include "setNCUStyle.C" 

using namespace std;
void dbtdependence(std::string inputFile) {


  int total = 0;
  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);
    

  TProfile* p_tau21_700to1500_before = new TProfile("","",14,60,200,0,1);
  TProfile* p_tau21_1500_before = new TProfile("","",14,60,200,0,1);
  TProfile* p_tau21_700to1500_after = new TProfile("","",14,60,200,0,1);
  TProfile* p_tau21_1500_after = new TProfile("","",14,60,200,0,1);
  TProfile* p_dbt_700to1500_before = new TProfile("","",14,60,200,-1,1);
  TProfile* p_dbt_1500_before = new TProfile("","",14,60,200,-1,1);
  TProfile* p_dbt_700to1500_after = new TProfile("","",14,60,200,-1,1);
  TProfile* p_dbt_1500_after = new TProfile("","",14,60,200,-1,1);

  total += data.GetEntriesFast();
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    Float_t  mcWeight  = data.GetFloat("mcWeight");
    Float_t HT         = data.GetFloat("HT");
    //Float_t*  fatjet_doubleSV = data.GetPtrFloat("FATjet_DoubleSV");
    
    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
    

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	// std::cout << thisTrig << " : " << results << std::endl;
	
	if( (thisTrig.find("HLT_PFHT800")!= std::string::npos && results==1)
	    )
	  {
	    passTrigger=true;
	    break;
	  }


      }


    if(!passTrigger)continue;

    //3. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;

    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      //if(!passFatJetLooseID[ij])continue;
      if(!FATjetPassIDTight[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<60 || fatjetPRmassL2L3Corr[ij]>200)continue;

      //Double_t tau21 = (fatjetTau2[ij]/fatjetTau1[ij]);
      //if(tau21>0.6)continue;
      
      fatjet.push_back(ij);     
    }
    
    if(fatjet.size()<2)continue;
    
    int aa = fatjet[0]; //Mjj[0].second;
    int ee = fatjet[1]; //Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
 
    Double_t dEta = fabs(Jet1->Eta() - Jet2->Eta());
    if(dEta>1.3)continue;

    Float_t mff=(*Jet1+*Jet2).M();
    Float_t msubt = mff-(fatjetPRmassL2L3Corr[aa]-125)-(fatjetPRmassL2L3Corr[ee]-125);
    if(msubt<800)continue;

    Double_t leadtau21 = (fatjetTau2[aa]/fatjetTau1[aa]);
    Double_t subltau21 = (fatjetTau2[ee]/fatjetTau1[ee]);
    
    //if(subltau21>0.6)continue;

    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;

    //if(addjet_doubleSV[addJetIndex[1]]<0.6)continue;

    if(mff>700 && mff<1500) {
      p_tau21_700to1500_before->Fill(fatjetPRmassL2L3Corr[aa],leadtau21);
      p_dbt_700to1500_before->Fill(fatjetPRmassL2L3Corr[aa],addjet_doubleSV[addJetIndex[0]]);
    }

    if(mff>1500) { 
      p_tau21_1500_before->Fill(fatjetPRmassL2L3Corr[aa],leadtau21);
      p_dbt_1500_before->Fill(fatjetPRmassL2L3Corr[aa],addjet_doubleSV[addJetIndex[0]]);
    }

    if(addjet_doubleSV[addJetIndex[0]]>0.6) {
      if(mff>700 && mff<1500)p_tau21_700to1500_after->Fill(fatjetPRmassL2L3Corr[aa],leadtau21);
      if(mff>1500)p_tau21_1500_after->Fill(fatjetPRmassL2L3Corr[aa],leadtau21);
    }

    if(leadtau21<0.6) {
      if(mff>700 && mff<1500)p_dbt_700to1500_after->Fill(fatjetPRmassL2L3Corr[aa],addjet_doubleSV[addJetIndex[0]]);
      if(mff>1500)p_dbt_1500_after->Fill(fatjetPRmassL2L3Corr[aa],addjet_doubleSV[addJetIndex[0]]);
    }



  } //end of the event loop
  cout<<"entries="<<total<<endl;

  setNCUStyle();

  TCanvas* c4 = new TCanvas("c4","",0,0,600,600);                                                                                                  
  c4->cd();
  p_tau21_1500_before->SetXTitle("h1 jet pruned mass");
  p_tau21_1500_before->SetYTitle("#tau_{2}/#tau_{1}");
  p_tau21_1500_before->GetYaxis()->SetTitleOffset(1.4);
  //p_tau21_1500_before->GetYaxis()->SetRangeUser(0.61,0.35);
  p_tau21_1500_before->SetMarkerStyle(20);
  p_tau21_1500_before->SetMarkerSize(0.8); 
  p_tau21_1500_before->SetLineColor(kCyan-5);
  p_tau21_1500_before->SetMarkerColor(kCyan-5);
  p_tau21_700to1500_before->SetMarkerStyle(21);
  p_tau21_700to1500_before->SetMarkerSize(0.8);
  p_tau21_700to1500_before->SetLineColor(kBlue-5);
  p_tau21_700to1500_before->SetMarkerColor(kBlue-5);
  p_tau21_1500_before->Draw();
  p_tau21_700to1500_before->Draw("same");

  TLegend *legendx = new TLegend(0.36,0.76,0.63,0.9);
  legendx->SetHeader("before DBT cut");
  legendx->AddEntry(p_tau21_700to1500_before,"700 GeV < dijet mass < 1500 GeV","lp");
  legendx->AddEntry(p_tau21_1500_before,"dijet mass > 1500 GeV", "lp");
  legendx->SetFillStyle(0);
  legendx->SetTextSize(0.035);
  legendx->Draw();

  c4->Print("DBTdependence_tau21_beforedbtcut.pdf");

  TCanvas* c3 = new TCanvas("c3","",0,0,600,600);
  c3->cd();
  p_tau21_1500_after->SetXTitle("h1 jet pruned mass");
  p_tau21_1500_after->SetYTitle("#tau_{2}/#tau_{1}");
  //p_tau21_1500_after->GetYaxis()->SetRangeUser(0.6,0.37);
  p_tau21_1500_after->GetYaxis()->SetTitleOffset(1.4);
  p_tau21_1500_after->SetMarkerStyle(20);
  p_tau21_1500_after->SetMarkerSize(0.8);
  p_tau21_1500_after->SetLineColor(kCyan-5);
  p_tau21_1500_after->SetMarkerColor(kCyan-5);
  p_tau21_700to1500_after->SetMarkerStyle(21);
  p_tau21_700to1500_after->SetMarkerSize(0.8);
  p_tau21_700to1500_after->SetLineColor(kBlue-5);
  p_tau21_700to1500_after->SetMarkerColor(kBlue-5);
  p_tau21_1500_after->Draw();
  p_tau21_700to1500_after->Draw("same");

  TLegend *legendq = new TLegend(0.36,0.76,0.63,0.9);
  legendq->SetHeader("after DBT cut");
  legendq->AddEntry(p_tau21_700to1500_after,"700 GeV < dijet mass < 1500 GeV","lp");
  legendq->AddEntry(p_tau21_1500_after,"dijet mass > 1500 GeV", "lp");
  legendq->SetFillStyle(0);
  legendq->SetTextSize(0.035);
  legendq->Draw();

  c3->Print("DBTdependence_tau21_afterdbtcut.pdf");

  TCanvas* c2 = new TCanvas("c2","",0,0,600,600);
  c2->cd();
  p_dbt_1500_before->SetXTitle("h1 jet pruned mass");
  p_dbt_1500_before->SetYTitle("h1 bbtag");
  p_dbt_1500_before->GetYaxis()->SetTitleOffset(1.7);
  p_dbt_1500_before->SetMarkerStyle(20);
  p_dbt_1500_before->SetMarkerSize(0.8);
  p_dbt_1500_before->SetLineColor(kCyan-5);
  p_dbt_1500_before->SetMarkerColor(kCyan-5);
  p_dbt_700to1500_before->SetMarkerStyle(21);
  p_dbt_700to1500_before->SetMarkerSize(0.8);
  p_dbt_700to1500_before->SetLineColor(kBlue-5);
  p_dbt_700to1500_before->SetMarkerColor(kBlue-5);
  p_dbt_1500_before->Draw();
  p_dbt_700to1500_before->Draw("same");

  TLegend *legendt = new TLegend(0.36,0.76,0.63,0.9);
  legendt->SetHeader("before #tau_{2}/#tau_{1} cut");
  legendt->AddEntry(p_dbt_700to1500_before,"700 GeV < dijet mass < 1500 GeV","lp");
  legendt->AddEntry(p_dbt_1500_before,"dijet mass > 1500 GeV", "lp");
  legendt->SetFillStyle(0);
  legendt->SetTextSize(0.035);
  legendt->Draw();

  c2->Print("DBTdependence_dbt_beforetau21cut.pdf");

  TCanvas* c1 = new TCanvas("c1","",0,0,600,600);
  c1->cd();
  p_dbt_1500_after->SetXTitle("h1 jet pruned mass");
  p_dbt_1500_after->SetYTitle("h1 bbtag");
  p_dbt_1500_after->GetYaxis()->SetTitleOffset(1.7);
  p_dbt_1500_after->SetMarkerStyle(20);
  p_dbt_1500_after->SetMarkerSize(0.8);
  p_dbt_1500_after->SetLineColor(kCyan-5);
  p_dbt_1500_after->SetMarkerColor(kCyan-5);
  p_dbt_700to1500_after->SetMarkerStyle(21);
  p_dbt_700to1500_after->SetMarkerSize(0.8);
  p_dbt_700to1500_after->SetLineColor(kBlue-5);
  p_dbt_700to1500_after->SetMarkerColor(kBlue-5);
  p_dbt_1500_after->Draw();
  p_dbt_700to1500_after->Draw("same");

  TLegend *legenda = new TLegend(0.36,0.76,0.63,0.9);
  legenda->SetHeader("after #tau_{2}/#tau_{1} cut");
  legenda->AddEntry(p_dbt_700to1500_after,"700 GeV < dijet mass < 1500 GeV","lp");
  legenda->AddEntry(p_dbt_1500_after,"dijet mass > 1500 GeV", "lp");
  legenda->SetFillStyle(0);
  legenda->SetTextSize(0.035);
  legenda->Draw();

  c1->Print("DBTdependence_dbt_aftertau21cut.pdf");

  /*  TFile* outfile = new TFile(Form("%s_dbtdependence.4.root",name.data()),"recreate");
  p_tau21_700to1500_before->Write("tau21vsPR_700to1500_before");
  p_tau21_1500_before->Write("tau21vsPR_1500_before");
  p_tau21_700to1500_after->Write("tau21vsPR_700to1500_after");
  p_tau21_1500_after->Write("tau21vsPR_1500_after");
  p_dbt_700to1500_before->Write("DBTvsPR_700to1500_before");
  p_dbt_1500_before->Write("DBTvsPR_1500_before");
  p_dbt_700to1500_after->Write("DBTvsPR_700to1500_after");
  p_dbt_1500_after->Write("DBTvsPR_1500_after");
  outfile->Write();*/
}

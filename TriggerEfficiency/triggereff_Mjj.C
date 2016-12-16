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
#include <TGraphAsymmErrors.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"
 
using namespace std;
void triggereff_Mjj(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  TGraphAsymmErrors* hAE_Mjj = new TGraphAsymmErrors();
  TGraphAsymmErrors* hAE_Mjjred = new TGraphAsymmErrors();

  TH1F* h_Mjj1 =  new TH1F("","",80,400,5000);
  TH1F* h_Mjj2 =  new TH1F("","",80,400,5000);
  TH1F* h_Mjj3 =  new TH1F("","",80,400,5000);
  TH1F* h_Mjjred1 =  new TH1F("","",80,400,5000);
  TH1F* h_Mjjred2 =  new TH1F("","",80,400,5000);
  TH1F* h_Mjjred3 =  new TH1F("","",80,400,5000);
  h_Mjj1->Sumw2();
  h_Mjj2->Sumw2();
  h_Mjj3->Sumw2();
  h_Mjjred1->Sumw2();
  h_Mjjred2->Sumw2();
  h_Mjjred3->Sumw2();
  
  Long64_t nPass[20]={0};

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
     
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

       
    nPass[0]++;

    vector<int> fatjet;
    for(int ij=0; ij<nFJets; ij++) {
      fatjet.push_back(ij);
    }
    if(fatjet.size()<2)continue;

    int aa = fatjet[0]; 
    int ee = fatjet[1]; 
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    Float_t mff=(*Jet1+*Jet2).M();
    Float_t msubt = mff-(fatjetPRmassL2L3Corr[aa]-125)-(fatjetPRmassL2L3Corr[ee]-125);
 
    h_Mjj1->Fill(mff);
    h_Mjjred1->Fill(msubt);

    nPass[1]++;

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];

	//std::cout << thisTrig << " : " << results << std::endl;
        
        if( (thisTrig.find("HLT_PFHT800")!= std::string::npos && results==1)
            )
          {
            passTrigger=true;
            break;
          }


      }


    if(!passTrigger)continue;

    h_Mjj2->Fill(mff);
    h_Mjjred2->Fill(msubt);
    nPass[2]++;

  } //end of the event loop
  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,800,600);
  c3->cd();

  hAE_Mjj->BayesDivide(h_Mjj2,h_Mjj1,"v");
  hAE_Mjj->GetXaxis()->SetRangeUser(400,5000);
  hAE_Mjj->GetYaxis()->SetTitle("Trigger Efficiency");                                                                                                   
  hAE_Mjj->GetXaxis()->SetTitle("dijet mass");
  hAE_Mjj->GetXaxis()->SetTitleSize(0.05);
  hAE_Mjj->GetYaxis()->SetTitleSize(0.05);

  hAE_Mjj->SetLineWidth(2);
  hAE_Mjj->SetLineColor(kRed+3);
  //h_Mjj3->Divide(h_Mjj2,h_Mjj1);
  hAE_Mjj->Draw();

  TCanvas* c =  new TCanvas("c","c",0,0,800,600);
  c->cd();

  hAE_Mjjred->BayesDivide(h_Mjjred2,h_Mjjred1,"v");
  hAE_Mjjred->GetXaxis()->SetRangeUser(400,5000);
  hAE_Mjjred->GetYaxis()->SetTitle("Trigger Efficiency");
  hAE_Mjjred->GetXaxis()->SetTitle("reduced dijet mass");
  hAE_Mjjred->GetXaxis()->SetTitleSize(0.05);
  hAE_Mjjred->GetYaxis()->SetTitleSize(0.05);
  
  hAE_Mjjred->SetLineWidth(2);
  hAE_Mjjred->SetLineColor(kGreen+3);
  //hAE_Mjjred->Divide(h_Mjjred2,h_Mjjred1);
  hAE_Mjjred->Draw();

  std::string bulkg_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  std::string radion_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4500"};

  bool BulkGrav=(inputFile.find("BulkGrav")!= std::string::npos);
  if(BulkGrav) {
    for(int i=0;i<10;i++){
      bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos); 
      if(bulkmass) {
	c3->Print(Form("TriggerEff_Mjj_AsymmError_BulkGrav_%s.pdf",bulkg_name[i].data()));
	c->Print(Form("TriggerEff_Mjjred_AsymmError_BulkGrav_%s.pdf",bulkg_name[i].data()));
      }
    }
  }
  else {
    for(int i=0;i<10;i++){
      bool radionmass=(inputFile.find(Form("%s",radion_name[i].data()))!= std::string::npos);
      if(radionmass) {
	c3->Print(Form("TriggerEff_Mjj_AsymmError_Radion_%s.pdf",radion_name[i].data()));
	c->Print(Form("TriggerEff_Mjjred_AsymmError_Radion_%s.pdf",radion_name[i].data()));
      }
    }
  }
  

}

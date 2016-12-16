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
void triggereff_HT(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  TGraphAsymmErrors* h_triggeff = new TGraphAsymmErrors();

  TH1F* h_HT1 =  new TH1F("","",80,400,5000);
  TH1F* h_HT2 =  new TH1F("","",80,400,5000);
  TH1F* h_HT3 =  new TH1F("","",80,400,5000);
  h_HT1->Sumw2();
  h_HT2->Sumw2();
  h_HT3->Sumw2();

  Long64_t nPass[20]={0};

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nTHINJet         = data.GetInt("THINnJet");
    const int nTJets=nTHINJet;
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
     
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

       
    nPass[0]++;

    Float_t HT=0;
    for(int ij=0; ij<nTJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
      if(thisJet->Pt()<20)continue;
      if(fabs(thisJet->Eta())>3.0)continue;
       
      HT += fabs(thisJet->Pt());
    }
    h_HT1->Fill(HT);
 
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

    h_HT2->Fill(HT);
    nPass[2]++;

  } //end of the event loop
  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,800,600);
  c3->cd();

  h_triggeff->BayesDivide(h_HT2,h_HT1,"v");

  h_triggeff->GetXaxis()->SetRangeUser(400,5000);
  h_triggeff->GetYaxis()->SetTitle("Trigger Efficiency");                                                                                                   
  h_triggeff->GetXaxis()->SetTitle("HT");
  h_triggeff->GetXaxis()->SetTitleSize(0.05);
  h_triggeff->GetYaxis()->SetTitleSize(0.05);
  h_triggeff->SetLineWidth(2);
  h_triggeff->SetLineColor(kBlue+3);
  h_triggeff->Draw();
  
  std::string bulkg_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  std::string radion_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4500"};

  bool BulkGrav=(inputFile.find("BulkGrav")!= std::string::npos);
  if(BulkGrav) {
    for(int i=0;i<10;i++){
      bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos); 
      if(bulkmass)c3->Print(Form("TriggerEff_HT_AsymmError_BulkGrav_%s.pdf",bulkg_name[i].data()));
    }
  }
  else {
    for(int i=0;i<10;i++){
      bool radionmass=(inputFile.find(Form("%s",radion_name[i].data()))!= std::string::npos);
      if(radionmass)c3->Print(Form("TriggerEff_HT_AsymmError_Radion_%s.pdf",radion_name[i].data()));
    }
  }
  

}

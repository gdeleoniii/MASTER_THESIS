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
 
using namespace std;
void passfailratio(std::string inputFile) {

  int total = 0;
  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);
    


  float xbins1[] = {-75,-45,-20,10,40,75};
  const int nbins1=sizeof(xbins1)/sizeof(xbins1[0])-1;
  TH1F* hpass_lead = new TH1F("","",nbins1,xbins1);
  TH1F* hfail_lead = new TH1F("","",nbins1,xbins1);
  TH1F* hratio_lead = new TH1F("","",nbins1,xbins1);
  TH1F* hprmass = new TH1F("","",34,40,210);
  hpass_lead->Sumw2();
  hfail_lead->Sumw2();
  hratio_lead->Sumw2();
    
  Long64_t DENOM = 0;

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
      if(fatjetPRmassL2L3Corr[ij]<50 || fatjetPRmassL2L3Corr[ij]>200)continue;

      
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
    if(leadtau21<0.6)continue;
    if(subltau21>0.6)continue;    
    
    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;


    if(fatjetPRmassL2L3Corr[ee]<105 || fatjetPRmassL2L3Corr[ee]>135)continue; // sub-leading pass the purned mass cut
    if(addjet_doubleSV[addJetIndex[1]]>0.6)continue; // let the sub-leading pass the double-b tagger cut

    if(fatjetPRmassL2L3Corr[aa]>105 && fatjetPRmassL2L3Corr[aa]<135)continue;
    if(addjet_doubleSV[addJetIndex[0]]>0.6) {
      hpass_lead->Fill(fatjetPRmassL2L3Corr[aa]-125);
    }
    else if(addjet_doubleSV[addJetIndex[0]]<0.6) {
      hfail_lead->Fill(fatjetPRmassL2L3Corr[aa]-125);
    }



  } //end of the event loop
  cout<<"entries="<<total<<endl;
  cout<<"events="<<DENOM<<endl;

  hratio_lead->Divide(hpass_lead,hfail_lead);
 
  cout<<hratio_lead->GetBinContent(1)<<endl; 
  cout<<hratio_lead->GetBinContent(2)<<endl;
  cout<<hratio_lead->GetBinContent(3)<<endl;
  cout<<hratio_lead->GetBinContent(4)<<endl;
  cout<<hratio_lead->GetBinContent(5)<<endl;
  
  TFile* outfile = new TFile("JetHT-Run2015_passfailratio.7.0.2.root","recreate");
  hpass_lead->Write("pass");
  hfail_lead->Write("fail");
  hratio_lead->Write("pfratio");
  outfile->Write();

}

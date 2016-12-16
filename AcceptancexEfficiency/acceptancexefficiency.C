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
#include "setNCUStyle.C"
 
using namespace std;
void acceptancexefficiency(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

 Float_t nPass[20]={0};

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    nPass[0]++;

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    //Float_t*  fatjet_doubleSV = data.GetPtrFloat("FATjet_DoubleSV");
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));    

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

    nPass[1]++;   

    //3. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
       
    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!FATjetPassIDTight[ij])continue;
      //if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      
      fatjet.push_back(ij);	
    }
    
    if(fatjet.size()<2)continue;

    nPass[2]++;

    int aa = fatjet[0]; //Mjj[0].second;
    int ee = fatjet[1]; //Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    //if(fatjetPRmassL2L3Corr[aa]<105 || fatjetPRmassL2L3Corr[aa]>135)continue;
    //if(fatjetPRmassL2L3Corr[ee]<105 || fatjetPRmassL2L3Corr[ee]>135)continue;

    //nPass[3]++;

    Double_t dEta = fabs(Jet1->Eta() - Jet2->Eta());
    if(dEta>1.3)continue;

    nPass[3]++;

    Float_t mff=(*Jet1+*Jet2).M();
    Float_t msubt = mff-(fatjetPRmassL2L3Corr[aa]-125)-(fatjetPRmassL2L3Corr[ee]-125);
    if(msubt<800)continue;

    nPass[4]++;

    if(fatjetPRmassL2L3Corr[aa]<105 || fatjetPRmassL2L3Corr[aa]>135)continue;
    if(fatjetPRmassL2L3Corr[ee]<105 || fatjetPRmassL2L3Corr[ee]>135)continue;

    nPass[5]++;

    Double_t leadtau21 = (fatjetTau2[aa]/fatjetTau1[aa]);
    Double_t subltau21 = (fatjetTau2[ee]/fatjetTau1[ee]);
    if(leadtau21>0.6)continue;
    if(subltau21>0.6)continue;

    nPass[6]++;

    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;

    if(addjet_doubleSV[addJetIndex[0]]<0.6)continue;
    if(addjet_doubleSV[addJetIndex[1]]<0.6)continue;

    nPass[7]++;


  } //end of the event loop

  cout<<"Trigger = "<<nPass[1]/nPass[0]<<endl;
  cout<<"Pre-selection = "<<nPass[2]/nPass[1]<<endl;
  cout<<"Delta eta = "<<nPass[3]/nPass[2]<<endl;
  cout<<"Reduced mass = "<<nPass[4]/nPass[3]<<endl;
  cout<<"correctedPRmass = "<<nPass[5]/nPass[4]<<endl;
  cout<<"Tau21 = "<<nPass[6]/nPass[5]<<endl;
  cout<<"bbtag = "<<nPass[7]/nPass[6]<<endl;
  cout<<"Total Eff = "<<nPass[7]/nPass[0]<<endl;

  /*  ofstream fout;
  fout.open(Form("%s.dat",inputFile.data()),ios::out | ios::app);
  fout<<nPass[0]<<" "<<nPass[1]<<" "<<nPass[2]<<" "<<nPass[3]<<" "<<nPass[4]<<" "<<nPass[5]<<" "<<nPass[6]<<" "<<nPass[7]<<endl;
  fout.close();*/
}

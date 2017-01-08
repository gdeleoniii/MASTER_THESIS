#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"
#include "readSample.h"

void roc_curve(std::string inputFile) {

  int total = 0;
  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
   
  // Declare the histogram

  Long64_t DENOM = 0;
  Long64_t DBTnum[21] = {0};
  Long64_t FATnum[21] = {0};
  Long64_t SUBnum[21] = {0};

  total += data.GetEntriesFast();

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 10000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetCISVV2  = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);

    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      //if(thisJet->Pt()<200)continue;
      if(thisJet->Pt()<300)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      //if(!passFatJetLooseID[ij])continue;
      if(!FATjetPassIDTight[ij])continue;
      //if(fatjetPRmass[ij]<70||fatjetPRmass[ij]>200)continue;
      //if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      if(fatjetPRmassL2L3Corr[ij]<70 || fatjetPRmassL2L3Corr[ij]>200)continue;

      fatjet.push_back(ij);
    }

    if(fatjet.size()<2)continue;

    int aa = fatjet[0]; 
    int ee = fatjet[1]; 
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    Double_t dEta = fabs(Jet1->Eta() - Jet2->Eta());
    if(dEta>1.3)continue;

    Float_t mff=(*Jet1+*Jet2).M();
    if(mff<1000)continue;

    int addJetIndex[2]={-1,-1};    
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;}
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;}
    }
    
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;

    int fatjetIndex[2]={aa,ee};

    //---------------per jet counting------------------

    Long64_t nDBTPass[21] = {0};
    Long64_t nFATPass[21] = {0};
    DENOM += 2;

    Float_t add = -1;
    Float_t fat = 0;
    for(int i=0;i<21;i++) {
      for(int j=0;j<2;j++) {
	int ajet = addJetIndex[j];
	int fjet = fatjetIndex[j];
        if(addjet_doubleSV[ajet]>add)nDBTPass[i]++;
        if(fatjetCISVV2[fjet]>fat)nFATPass[i]++;
      }

      DBTnum[i] += nDBTPass[i];
      FATnum[i] += nFATPass[i];
      add += 0.1;
      fat += 0.05;
    }

    Long64_t nSUB1Pass[21] = {0};
    Long64_t nSUB2Pass[21] = {0};

    Float_t sub = 0;
    int fjet1 = fatjetIndex[0];
    int fjet2 = fatjetIndex[1];
    for(int i=0;i<21;i++) {
      for(int sub1=0;sub1<FATnSubSDJet[fjet1];sub1++) {
        if(FATsubjetSDCSV[fjet1][sub1]>sub)nSUB1Pass[i]++;
      }
      for(int sub2=0;sub2<FATnSubSDJet[fjet2];sub2++) {
        if(FATsubjetSDCSV[fjet2][sub2]>sub)nSUB2Pass[i]++;
      }

      if(nSUB1Pass[i] == 1) nSUB1Pass[i] = 0;
      if(nSUB2Pass[i] == 1) nSUB2Pass[i] = 0;
      SUBnum[i] += (nSUB1Pass[i]/2) + (nSUB2Pass[i]/2);
      sub +=0.05;
    }
    //-----------------------------------------------------------

  } // end of event loop
  fprintf(stderr, "Processed all events\n");

  std::string HT_name[]={"HT500to700","HT700to1000","HT1000to1500","HT1500to2000","HT2000toInf"};
  Double_t xsec[5]={32100.0,6831.0,1207.0,119.9,25.24};
  Float_t bbtag[21] = {0};
  Float_t fatcsv[21] = {0};
  Float_t subcsv[21] = {0};
  std::cout << "total= "<< total << std::endl;
   std::cout << "Denominator= "<< DENOM << std::endl;
  

  for(int g=0;g<5;g++) {
    bool HT=(inputFile.find(Form("%s",HT_name[g].data()))!= std::string::npos);
    if(HT) {
      std::cout<< "xsec = " << xsec[g]<< std::endl;
      std::cout << "DBT" << std::endl; 
     ofstream fout;
      fout.open(Form("ROC_Curve_DBT_JetCount_QCD_%s.dat",HT_name[g].data()),ios::out | ios::app);
      for(int i=0;i<21;i++) {
	bbtag[i] = (DBTnum[i]*xsec[g])/total;
	std::cout <<  DBTnum[i] << " " << bbtag[i] << std::endl;
	fout<<bbtag[i]<<endl;
      }
      fout.close();
      
      std::cout << "FAT" << std::endl;
      ofstream fout1;
      fout1.open(Form("ROC_Curve_FAT_JetCount_QCD_%s.dat",HT_name[g].data()),ios::out | ios::app);
      for(int i=0;i<21;i++) {
	fatcsv[i]=(FATnum[i]*xsec[g])/total;
	std::cout <<  FATnum[i] << " " << fatcsv[i] << std::endl;
	fout1<<fatcsv[i]<<endl;
      }
      fout1.close();
      
      std::cout << "SUB" << std::endl;
      ofstream fout2;
      fout2.open(Form("ROC_Curve_SUB_JetCount_QCD_%s.dat",HT_name[g].data()),ios::out | ios::app);
      for(int i=0;i<21;i++) {
	subcsv[i]=(SUBnum[i]*xsec[g])/total;
	std::cout <<  SUBnum[i] << " " << subcsv[i] << std::endl;
        fout2<<subcsv[i]<<endl;
      }
      fout2.close();
    }
  }
  

}

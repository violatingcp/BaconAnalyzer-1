//================================================================================================
//
// Perform preselection for W/Zprime(->qq)+jets events and produce bacon bits 
//
// Input arguments
//   argv[1] => lName = input bacon file name
//   argv[2] => lOption = dataset type: "mc", "data"
//   argv[3] => lJSON = JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[4] => lXS = cross section (pb), ignored for data 
//   argv[5] => weight = total weight, ignored for data
//________________________________________________________________________________________________

#include "../include/GenPartLoader.hh"
#include "../include/GenJetLoader.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TError.h>
#include <string>
#include <iostream>

// Object Processors
GenPartLoader   *fGen        = 0; 
GenJetLoader    *fGenJet     = 0; 
TH1F *fHist                  = 0;

const int NUM_PDF_WEIGHTS = 60;
// Load tree and return infile
TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  fHist        = (TH1F* ) lFile->FindObjectAny("TotalEvents");
  return lTree;
}
// === MAIN =======================================================================================================
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");
  const std::string lName        = argv[1];
  const std::string lOption      = argv[2];
  const std::string lJSON        = argv[3];
  const double      lXS          = atof(argv[4]);
  //const double      weight       = atof(argv[5]);
  const int      iSplit          = atoi(argv[6]);
  const int      maxSplit        = atoi(argv[7]);

  TTree *lTree = load(lName); 
  fGen      = new GenPartLoader(lTree);                     // 
  fGenJet   = new GenJetLoader (lTree,"AK8Puppi",2,200.);   // 2 Leading jets with gen pt > 200 
  
  TFile *lFile = TFile::Open("Output.root","RECREATE");
  TTree *lOut  = new TTree("Events","Events");

  //Setup histograms containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  // Setup Tree
  fGen ->setupTree (lOut,float(lXS),0.8);

  // Loop over events i0 = iEvent
  int neventsTotal = int(lTree->GetEntriesFast());
  int minEventsPerJob = neventsTotal / maxSplit;
  //int leftoverEvents = neventsTotal % maxSplit;
  int minEvent = iSplit * minEventsPerJob;
  int maxEvent = (iSplit+1) * minEventsPerJob;
  if (iSplit + 1 == maxSplit) maxEvent = neventsTotal;
  std::cout << neventsTotal << " total events" << std::endl;
  std::cout << iSplit << " iSplit " << std::endl;
  std::cout << maxSplit << " maxSplit " << std::endl;
  std::cout << minEvent << " min event" << std::endl;
  std::cout << maxEvent << " max event" << std::endl;  
  for(int i0 = minEvent; i0 < maxEvent; i0++) {
    //for(int i0 = 0; i0 < int(10000); i0++){ // for testing
    if (i0%1000 == 0) std::cout << i0 << " events processed " << std::endl;
    fGenJet->load(i0);
    std::vector<TJet*> pGenJets;
    fGenJet->selectVJets(pGenJets);
    // Check GenInfo
    float lWeight = 1;
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    if(pGenJets.size() == 0) continue;
    fGen->load(i0);
    if(fGen->leptonVeto()) continue;
    for(unsigned int i1 = 0; i1 < pGenJets.size(); i1++) {  
      fGen->chain(pGenJets[i1],0.1); //zcut 0.1
      lWeight = fGen->fWeight;
      lOut->Fill();
    }
   }
  std::cout << lTree->GetEntriesFast() << std::endl;
  lFile->cd();
  lOut->Write();  
  NEvents->Write();
  lFile->Close();
}

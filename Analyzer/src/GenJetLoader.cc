#include "../include/GenJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

GenJetLoader::GenJetLoader(TTree *iTree,std::string iJet,int iN,float iPt) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fN  = iN;
  fPt = iPt;
}
GenJetLoader::~GenJetLoader() { 
  delete fVJets;
  delete fVJetBr;
}
void GenJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
}
void GenJetLoader::selectVJets(std::vector<TJet*> &iGen) { 
  int lN = 0; 
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    if(lN > fN) continue;
    TJet *pVJet = (TJet*)((*fVJets)[i0]);    
    if(pVJet->genpt < fPt) continue;
    iGen.push_back(pVJet);
    lN++;
  }
}

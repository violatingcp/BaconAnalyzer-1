#ifndef GenJetLoader_H
#define GenJetLoader_H
#include "Utils.hh"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"

using namespace baconhep;

class GenJetLoader { 
public:
  GenJetLoader(TTree *iTree,std::string iJet,int iN,float iPt);
  ~GenJetLoader();
  void load (int iEvent);
  void selectVJets(std::vector<TJet*> &iGen);

  TClonesArray  *fVJets;
  TBranch       *fVJetBr;

protected: 
  TTree         *fTree;
  int            fN;
  float          fPt;
};
#endif

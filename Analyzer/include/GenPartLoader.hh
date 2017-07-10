#ifndef GenPartLoader_H
#define GenPartLoader_H
#include "Utils.hh"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace baconhep;

class GenPartLoader { 
public:
  GenPartLoader(TTree *iTree,bool iHadrons=true);
  ~GenPartLoader();
  void reset();
  void setupTree(TTree *iTree,float iXSIn,float iRadius);
  void load (int iEvent);
  void chain(TJet *iJet, float iZCut);
  float deltaR(float iPhi0,float iEta0,float iPhi1,float iEta1);
  float mother(TGenParticle* iPart, std::vector<TGenParticle*> &iPartons);
  int  find(std::vector<TGenParticle*> &iPartons,TGenParticle* iPart);
  bool isNeutrino(int iPdgId);
  float  phi(float iPhi0,float iPhi1);
  int  simplifiedPdg(int iPdgId);
  int  parentid(TGenParticle *iPart,bool iOutput=false);
  bool leptonVeto();
  //Debug
  void parentage(TGenParticle* iPart,TJet *iJet,std::vector<TGenParticle*> &iPartons);
  void printVtx(TGenParticle* iPart,TJet *iJet);  

  TClonesArray  *fGens;
  TBranch       *fGenBr;
  TGenEventInfo *fGenInfo;
  TBranch       *fGenInfoBr;

  float fWeight;
  std::vector<std::string> fLabels;
  std::vector<std::vector<float> > fVars;
  std::unordered_set<unsigned> fPartons;

protected: 
  TTree         *fTree;
  int fPartonBase;
  int fParticleBase;
  float fXS;
  float fRadius;
  float fDRHeavy;
};
#endif

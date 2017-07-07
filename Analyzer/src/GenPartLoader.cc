#include <iostream>
#include <assert.h> 
#include <string> 
#include "../include/GenPartLoader.hh"

using namespace baconhep;

GenPartLoader::GenPartLoader(TTree *iTree,bool iHadrons) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");

  fDRHeavy = 0.5;
  //  fGenVtx   = new TClonesArray("baconhep::TGenVtx");
  //  iTree->SetBranchAddress("GenVtx",       &fGenVtx);
  //  fGenVtxBr = iTree->GetBranch("GenVtx");
  if(iHadrons)  fPartons = {1, 2, 3, 4, 5, 21,15,11,12,111,211,221,113,213,130,310,311,321,411,421,433,441,443,511,521,513,523,531,541};
  if(!iHadrons) fPartons = {1, 2, 3, 4, 5, 21,15,11,12};
}
GenPartLoader::~GenPartLoader() { 
  delete fGenInfo;
  delete fGenInfoBr;

  delete fGens;
  delete fGenBr;

  //delete fGenVtx;
  //delete fGenVtxBr;
}
void GenPartLoader::reset() { 
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0].clear();
}
void GenPartLoader::setupTree(TTree *iTree,float iXSIn,float iRadius) { 
  reset();
  fRadius = iRadius;
  fTree   = iTree;
  fLabels.push_back("weight");
  fLabels.push_back("jet_pt");
  fLabels.push_back("jet_eta");
  fLabels.push_back("jet_phi");
  fLabels.push_back("jet_mass");
  fPartonBase = fLabels.size();
  fLabels.push_back("parton_pt");
  fLabels.push_back("parton_eta");
  fLabels.push_back("parton_phi");
  fLabels.push_back("parton_mass");
  fLabels.push_back("parton_pdgId");
  fParticleBase = fLabels.size();
  fLabels.push_back("part_pt");
  fLabels.push_back("part_eta");
  fLabels.push_back("part_phi");
  fLabels.push_back("part_pdgId");
  fLabels.push_back("part_d0");
  fLabels.push_back("part_mthrIndx");
  for(unsigned int i0 = 0; i0 < fLabels.size(); i0++) { 
    std::vector<float> lVars; 
    fVars.push_back(lVars); 
  }
 setupNtupleArr(iTree,fVars,fLabels);
 fXS = iXSIn;
}
void GenPartLoader::load(int iEvent) { 
  fGens     ->Clear();
  //fGenVtx   ->Clear();
  fGenBr    ->GetEntry(iEvent);
  fGenInfoBr->GetEntry(iEvent);
  //fGenVtxBr ->GetEntry(iEvent);
  fWeight = fXS*(fGenInfo->weight);
}
void GenPartLoader::chain(TJet *iJet, float iZCut) {
  reset();
  fVars[0].push_back(fWeight);
  fVars[1].push_back(iJet->genpt);
  fVars[2].push_back(iJet->geneta);
  fVars[3].push_back(iJet->genphi);
  fVars[4].push_back(iJet->genm);
  const std::unordered_set<unsigned> partonIDs = fPartons;
  auto validPID = [partonIDs](int id) -> bool {
    return (partonIDs.find(abs(id)) != partonIDs.end());
  };
  const std::unordered_set<unsigned> resonanceIDs = {6, 23, 24, 25, 54, 55}; 
  auto isResonance = [resonanceIDs](int id) -> bool {
    return ( (abs(id) > 10000) || // e.g. Z'
             (resonanceIDs.find(abs(id)) != resonanceIDs.end()) );
  };
  std::vector<TGenParticle*> lParton;
  std::vector<TGenParticle*> lChildren;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(!validPID(int(pGen->pdgId))) continue;        // Skip the intermediate stuff
    if(pGen->pt < iJet->genpt*iZCut) continue; 
    lChildren.clear();
    for(int i1 = 0; i1 < fGens->GetEntriesFast(); i1++) { 
      TGenParticle *pGen1 = (TGenParticle*)((*fGens)[i1]);
      //if (parentid(pGen1) != i0 || !validPID(int(pGen1->pdgId)) || isResonance(pGen1->pdgId)) continue;
      if (parentid(pGen1) != i0) continue;
      lChildren.push_back(pGen1);
    }
    int pPass = 0; 
    for(unsigned int i1 = 0; i1 < lChildren.size(); i1++) if(lChildren[i1]->pt > iJet->genpt*iZCut  && validPID(int(lChildren[i1]->pdgId)) && !isResonance(lChildren[i1]->pdgId) ) pPass++; 
    if(pPass > 1 && pPass == int(lChildren.size()) ) continue; //!!!!!!!! 2nd Left is to account for N daughters > 2 where one daughter fails condition and thus we lose assignment or in the case wherethe parton starts going to Hadrons in which case we would lose the leg
    bool pFoundParent      = false;
    bool pFoundResonance  = false;
    TGenParticle *iter = pGen;
    int pParent = parentid(iter);
    while(pParent > 1) { 
      iter = (TGenParticle*)((*fGens)[pParent]);
      //if (pFoundResonance) pFoundResonance = isResonance(iter->pdgId); //I really don't understand this so I commented it out
      if (find(lParton,iter) != -1) {pFoundParent = true; break;}
      pParent = parentid(iter);
    }
    if(pFoundParent && !pFoundResonance) continue;
    lParton.push_back(pGen);
    pParent = parentid(iter);
  }
  int lBase=fPartonBase;
  //std::cout <<" ===> Partons " << lParton.size() << std::endl;
  for(unsigned int i0 = 0; i0 < lParton.size(); i0++) { 
    fVars[lBase+0].push_back(lParton[i0]->pt);
    fVars[lBase+1].push_back(lParton[i0]->eta-iJet->geneta);
    fVars[lBase+2].push_back(phi(lParton[i0]->phi,iJet->genphi));
    fVars[lBase+3].push_back(lParton[i0]->mass);
    fVars[lBase+4].push_back(float(lParton[i0]->pdgId));
    //std::cout <<" ===> Parton " << i0 << " -id- " << lParton[i0]->pdgId << " - " << lParton[i0]->pt << " -dR- " << deltaR(lParton[i0]->phi,lParton[i0]->eta,iJet->genphi,iJet->geneta) << " === " << lParton[i0]->status << std::endl;
  }
  //std::cout << " ===> End Partons " << std::endl;
  lChildren.clear();
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(pGen->status != 1 || isNeutrino(pGen->pdgId)) continue;
    if(deltaR(pGen->phi,pGen->y,iJet->genphi,iJet->geneta) > fRadius) continue; 
    if(pGen->pt     < 0.2) continue;
    lChildren.push_back(pGen);
  }
  lBase=fParticleBase;
  for(unsigned int i0 = 0; i0 < lChildren.size(); i0++) { 
    fVars[lBase+0].push_back(lChildren[i0]->pt);
    fVars[lBase+1].push_back(lChildren[i0]->eta-iJet->geneta);
    fVars[lBase+2].push_back(phi(lChildren[i0]->phi,iJet->genphi));
    fVars[lBase+3].push_back(float(simplifiedPdg(lChildren[i0]->pdgId)));
    fVars[lBase+4].push_back(float(lChildren[i0]->d0));
    fVars[lBase+5].push_back(float(mother(lChildren[i0],lParton)));
    //if(mother(lChildren[i0],lParton) == -1 && parentid(lChildren[i0]) >  1) parentage(lChildren[i0],iJet,lParton);
    //if(mother(lChildren[i0],lParton) == -1) parentage(lChildren[i0],iJet,lParton);
  }
}
float GenPartLoader::deltaR(float iPhi0,float iEta0,float iPhi1,float iEta1) { 
  float pDEta = fabs(iEta0-iEta1);
  float pDPhi = fabs(iPhi0-iPhi1); if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}
float GenPartLoader::mother(TGenParticle* iPart, std::vector<TGenParticle*> &iPartons) { 
  int iId = -1;
  iId = find(iPartons,iPart);
  if(iId > -1) return float(iId); 
  int lParentId = parentid(iPart);
  TGenParticle *pParent =  (TGenParticle*)((*fGens)[lParentId]);
  while(find(iPartons,pParent) == -1 && lParentId > 1) {pParent =  (TGenParticle*)((*fGens)[lParentId]); lParentId = parentid(pParent);}
  if(lParentId > -1) iId = find(iPartons,pParent);
  pParent =  (TGenParticle*)((*fGens)[lParentId]);
  return float(iId); 
}
int GenPartLoader::find( std::vector<TGenParticle*> &iPartons,TGenParticle* iPart) { 
  for(unsigned int i0 = 0; i0 < iPartons.size(); i0++) if(iPartons[i0] == iPart) return i0;
  return -1;
}
bool GenPartLoader::isNeutrino(int iPdgId) { 
  return (abs(iPdgId) == 12 || abs(iPdgId) == 14 || abs(iPdgId) == 16);
}
int GenPartLoader::simplifiedPdg(int iPdgId) { 
  if(abs(iPdgId) == 130)  return 2112; //kLong  (neutral)
  if(abs(iPdgId) == 310)  return 2112; //kShort (neutral)
  if(abs(iPdgId) == 2112) return 2112; //Neutron
  if(abs(iPdgId) == 3322 || abs(iPdgId) == 3212) return 2112; //Neutral Strange Baryons
  if(abs(iPdgId) == 321)  return (iPdgId/abs(iPdgId))*211; //Charged Kaon
  if(abs(iPdgId) == 2212) return (iPdgId/abs(iPdgId))*211; //Proton
  if(abs(iPdgId) > 3000 && abs(iPdgId) < 4000 && abs(iPdgId) != 3322 && abs(iPdgId) != 3212) return (iPdgId/abs(iPdgId))*211; //Strange Baryons
  return iPdgId;
}
float  GenPartLoader::phi(float iPhi0,float iPhi1) { 
  double pDPhi = iPhi0-iPhi1;
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi) && pDPhi < 0) pDPhi =   2.*TMath::Pi()+pDPhi;
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi) && pDPhi > 0) pDPhi =  -2.*TMath::Pi()+pDPhi;
  return pDPhi;
}
//Debug Method
void GenPartLoader::parentage(TGenParticle* iPart,TJet *iJet,std::vector<TGenParticle*> &iPartons) { 
  const std::unordered_set<unsigned> partonIDs = fPartons;
  auto validPID = [partonIDs](int id) -> bool {
    return (partonIDs.find(abs(id)) != partonIDs.end());
  };
  const std::unordered_set<unsigned> resonanceIDs = {6, 23, 24, 25, 54, 55}; 
  auto isResonance = [resonanceIDs](int id) -> bool {
    //return ( (resonanceIDs.find(abs(id)) != resonanceIDs.end()) );
    return ( (abs(id) > 10000) || // e.g. Z'
	     (resonanceIDs.find(abs(id)) != resonanceIDs.end()) );
  };
  const std::unordered_set<unsigned> bIDs = {511,513,521,523,531,533,541,543,5122,5112,5114,5132,5212,5212,5222,5214,5224,5232,5332,5142};
  const std::unordered_set<unsigned> cIDs = {411,413,421,423,431,433,441,443,4122,4112,5114,4132,4212,4212,4222,4214,4224,4232,4332,4142,20433};
  auto isB = [bIDs](int id) -> bool {
    return (bIDs.find(abs(id)) != bIDs.end());
  };
  auto isC = [cIDs](int id) -> bool {
    return (cIDs.find(abs(id)) != cIDs.end());
  };
  int lParentId = parentid(iPart);
  TGenParticle *pParent =  (TGenParticle*)((*fGens)[lParentId]);
  int lNCount = 0; 
  bool lFromRes = false;
  //bool lFromUE    = true;
  bool lPassPtCut = false;
  //bool lNoHF      = true;
  while(lParentId > 1) {
    //if(isResonance(pParent->pdgId)) lFromUE = false;
    double pZ     = pParent->pt/iJet->genpt;
    if(pZ > 0.1 && validPID(pParent->pdgId) && !lFromRes) lPassPtCut = true;
    if(isResonance(pParent->pdgId) && !lPassPtCut) lFromRes = true;
    pParent =  (TGenParticle*)((*fGens)[lParentId]);
    //if(pZ > 0.1 && (isC(pParent->pdgId) || isB(pParent->pdgId)) ) lNoHF = false;
    lParentId = parentid(pParent); 
    lNCount++; 
  }
  //if(lFromRes || !lNoHF) return;
  //if(!lPassPtCut) return;
  if(iPart->pt > 100) { 
    for(int i0 = 0; i0 < TMath::Max(fGens->GetEntriesFast(),50); i0++) {   
      TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
      if(pGen == iPart) break;
      std::cout <<" Scanning ====> -id- " << i0  << " -pdgid- " << pGen->pdgId << " -pt- " << pGen->pt  << " -eta- " << pGen->eta << " -phi- " << pGen->phi <<  " -parent- " << pGen->parent << " -dr- " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) << std::endl;
    }
  }
  lParentId = parentid(iPart); 
  pParent =  (TGenParticle*)((*fGens)[lParentId]);
  if(iPart->pt > 100) std::cout << "===> Missing Parent -parent- " << iPart->parent << " -id- " << iPart->pdgId << " -pt- " << iPart->pt  << " -- " << lNCount << " -- " << lFromRes << " -- " << lPassPtCut << std::endl; 
  while(lParentId > 1) {
    pParent =  (TGenParticle*)((*fGens)[lParentId]);
    if(iPart->pt > 100) std::cout <<" ====> Searching Parent " << lParentId << " -id- " << pParent->pdgId << " -pt- " << pParent->pt << " -eta- " << pParent->eta << " -phi- " << pParent->phi << " - jet " << iJet->geneta << " - " << iJet->genphi  << " -dr- " <<  deltaR(pParent->phi,pParent->y,iJet->genphi,iJet->geneta) << " -chec- " << find(iPartons,pParent) << std::endl;
    lParentId = parentid(pParent,false); 
  }
  if(iPart->pt > 100) std::cout <<" ===> Done " << pParent->parent << std::endl;
}
int GenPartLoader::parentid(TGenParticle *iPart,bool iOutput) { 
  const std::unordered_set<unsigned> bIDs = {511,513,521,523,531,533,541,543,5122,5112,5114,5132,5212,5212,5222,5214,5224,5232,5332,5142};
  const std::unordered_set<unsigned> cIDs = {411,413,421,423,431,433,441,443,4122,4112,5114,4132,4212,4212,4222,4214,4224,4232,4332,4142,20433};
  auto isB = [bIDs](int id) -> bool {
    return (bIDs.find(abs(id)) != bIDs.end());
  };
  auto isC = [cIDs](int id) -> bool {
    return (cIDs.find(abs(id)) != cIDs.end());
  };
  /*
  if(iOutput) { 
    for(int i0 = 0; i0 < TMath::Max(fGens->GetEntriesFast(),50); i0++) {   
      TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
      //if(pGen == iPart) break;
      std::cout <<" Scanning ====> -id- " << i0  << " -pdgid- " << pGen->pdgId << " -pt- " << pGen->pt  << " -eta- " << pGen->eta << " -phi- " << pGen->phi <<  " -parent- " << pGen->parent << " -dr- " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) << std::endl;
    }
  }
  */
  int lId = iPart->parent;
  if(isC(iPart->pdgId) && lId > 1) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[lId]);
    if((fabs(pGen->pdgId) < 10 && fabs(pGen->pdgId) > 3 && pGen->pt > 5) || fabs(pGen->pdgId) > 10) return lId;
    lId = 0; 
  }
  /*
  if(isB(iPart->pdgId) && lId > 1) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[lId]);
    if((fabs(pGen->pdgId) < 10 && fabs(pGen->pdgId) > 4) || fabs(pGen->pdgId) > 10) return lId;
    lId = 0; 
  }
  */
  if(lId > 1 && lId != -2 && !isB(iPart->pdgId)) return lId;
  int lPdgId = 0;
  if(isC(iPart->pdgId)) lPdgId = 4;
  if(isB(iPart->pdgId)) lPdgId = 5;
  if(lPdgId == 0) return 0;
  lId = 0; 
  bool lUseCharge = true; if(abs(iPart->pdgId) > 1000) lUseCharge = false; // Ignore the charge for the baryons  b/c it is complicated
  int  lCharge    = 1;    if(lPdgId == 5) lCharge = -1; //Correct for the charge assignement of the mesons
  if(iOutput) std::cout << "====> A " << std::endl;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {   
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(iOutput) std::cout <<" Scanning ====> -id- " << i0  << " -pdgid- " << pGen->pdgId << " -pt- " << pGen->pt  << " -eta- " << pGen->eta << " -phi- " << pGen->phi <<  " -parent- " << pGen->parent << " -dr- " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) << std::endl;
    if(pGen == iPart)  break;
    if(lPdgId != fabs(pGen->pdgId) || (int(iPart->pdgId/abs(iPart->pdgId)) != int(lCharge*(pGen->pdgId/abs(pGen->pdgId))) && lUseCharge)) continue;
    if(deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) > fDRHeavy || pGen->pt < 0.3*iPart->pt) continue;
    if(iOutput) std::cout <<"   Link ====> -id- " << i0 << " -s- " << lPdgId << " -pdgid- " << pGen->pdgId << " -b- " << iPart->pdgId << " -- " << pGen->pt  << " -pt- " << iPart->pt << " -- " << pGen->eta << " -- " << pGen->phi <<  " -- " << pGen->parent << " -dr- " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) << std::endl;
    lId = i0;
  }
  if(iOutput) std::cout <<" ===> seraching done " << lId << std::endl;
  return lId;
}
void GenPartLoader::printVtx(TGenParticle* iPart,TJet *iJet) { 
  std::cout << " ===> Vtx -id- " << iPart->pdgId << " -vtxid- " << iPart->vtxId << " -- " << iPart->vtxFlav << std::endl;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {   
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(pGen->vtxId == iPart->vtxId) std::cout << " ----> " << i0 << " -- id -- " << pGen->pdgId << " -- " << pGen->pt << " -dr- " <<  deltaR(pGen->phi,pGen->y,iJet->genphi,iJet->geneta)  << " -vtxflav- " << pGen->vtxFlav << std::endl;
  }
  std::cout << " ===> Vtx end " << std::endl;
}
bool GenPartLoader::leptonVeto() { 
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {   
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(fabs(pGen->pdgId) != 11 && fabs(pGen->pdgId) != 13) continue;
    if(pGen->pt < 15) continue;
    TGenParticle *pParent = (TGenParticle*)((*fGens)[pGen->parent]);
    if(fabs(pParent->pdgId) == 24 || fabs(pParent->pdgId) == 23) return true;
  }
  return false;
}

#include <iostream>
#include <assert.h> 
#include <string> 
#include "../include/GenPartLoader.hh"

using namespace baconhep;

GenPartLoader::GenPartLoader(TTree *iTree) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");

  fDRHeavy = 0.5;
  //  fGenVtx   = new TClonesArray("baconhep::TGenVtx");
  //  iTree->SetBranchAddress("GenVtx",       &fGenVtx);
  // fGenVtxBr = iTree->GetBranch("GenVtx");
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
  //const std::unordered_set<unsigned> partonIDs = {1, 2, 3, 4, 5, 21,15,11,12,411,421,433,441,443,511,521,513,523,531,541};
  const std::unordered_set<unsigned> partonIDs = {1, 2, 3, 4, 5, 21,15,11,12};
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
  //std::cout <<" ====> Gen " << std::endl;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    //std::cout <<" ====> -id- " << i0 << " -pdgid- " << pGen->pdgId << " -- " << pGen->pt  << " -- " << pGen->eta << " -- " << pGen->phi <<  " -- " << pGen->parent << std::endl;
    if(!validPID(int(pGen->pdgId))) continue;        // Skip the intermediate stuff
    if(pGen->pt < iJet->genpt*iZCut) continue; 
    lChildren.clear();
    for(int i1 = 0; i1 < fGens->GetEntriesFast(); i1++) { 
      TGenParticle *pGen1 = (TGenParticle*)((*fGens)[i1]);
      if (parentid(pGen1) != i0 || !validPID(int(pGen1->pdgId)) || isResonance(pGen1->pdgId)) continue;
      lChildren.push_back(pGen1);
    }
    int pPass = 0; 
    for(unsigned int i1 = 0; i1 < lChildren.size(); i1++) if(lChildren[i1]->pt > iJet->genpt*iZCut) pPass++; 
    //if(pPass > 1) continue; //!!!!!!!!
    bool pFoundParent      = false;
    bool pFoundResonance  = false;
    TGenParticle *iter = pGen;
    int pParent = parentid(iter);
    while(pParent > 1) { 
      iter = (TGenParticle*)((*fGens)[pParent]);
      if (pFoundResonance) pFoundResonance = isResonance(iter->pdgId); //I really don't understand this but it won't hurt
      if (find(lParton,iter) != -1) {pFoundParent = true; break;}
      pParent = parentid(iter);
    }
    if(pFoundParent && !pFoundResonance) continue;
    //TGenParticle *pParent =  (TGenParticle*)((*fGens)[pGen->parent]);
    //std::cout << " Parton Adding ==> " << lParton.size() << " --- > " << pGen->pt << " -- " << pGen->pdgId << " ====> " << pParent->pdgId << " -- " << pParent->eta << " -- " << pParent->phi << " ---> " << deltaR(pParent->phi,pParent->eta,iJet->genphi,iJet->geneta)  << std::endl;
    lParton.push_back(pGen);
    pParent = parentid(iter);
  }
  std::cout <<" ====> Gen end" << std::endl;
  int lBase=fPartonBase;
  std::cout <<" ===> Partons " << lParton.size() << std::endl;
  for(unsigned int i0 = 0; i0 < lParton.size(); i0++) { 
    fVars[lBase+0].push_back(lParton[i0]->pt);
    fVars[lBase+1].push_back(lParton[i0]->eta-iJet->geneta);
    fVars[lBase+2].push_back(phi(lParton[i0]->phi,iJet->genphi));
    fVars[lBase+3].push_back(lParton[i0]->mass);
    fVars[lBase+4].push_back(float(lParton[i0]->pdgId));
    std::cout <<" ===> Parton " << i0 << " -id- " << lParton[i0]->pdgId << " - " << lParton[i0]->pt << " -dR- " << deltaR(lParton[i0]->phi,lParton[i0]->eta,iJet->genphi,iJet->geneta) << std::endl;
  }
  std::cout << " ===> End Partons " << std::endl;
  lChildren.clear();
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(pGen->status != 1 || isNeutrino(pGen->pdgId)) continue;
    if(deltaR(pGen->phi,pGen->y,iJet->genphi,iJet->geneta) > fRadius) continue; 
    if(pGen->pt     < 0.2) continue;
    lChildren.push_back(pGen);
    //    fPdgIds.push_back(simplifiedPdg(pGen->pdgId));
  }

  //std::sort(fPdgIds.begin(), fPdgIds.end());
  //auto last = std::unique(fPdgIds.begin(), fPdgIds.end());
  //fPdgIds.erase(last, fPdgIds.end()); 
  lBase=fParticleBase;
  for(unsigned int i0 = 0; i0 < lChildren.size(); i0++) { 
    fVars[lBase+0].push_back(lChildren[i0]->pt);
    fVars[lBase+1].push_back(lChildren[i0]->eta-iJet->geneta);
    fVars[lBase+2].push_back(phi(lChildren[i0]->phi,iJet->genphi));
    fVars[lBase+3].push_back(float(simplifiedPdg(lChildren[i0]->pdgId)));
    fVars[lBase+4].push_back(float(lChildren[i0]->d0));
    fVars[lBase+5].push_back(float(mother(lChildren[i0],lParton)));
    if(mother(lChildren[i0],lParton) == -1 && parentid(lChildren[i0]) >  1) parentage(lChildren[i0],iJet,lParton);
  }
}
float GenPartLoader::deltaR(float iPhi0,float iEta0,float iPhi1,float iEta1) { 
  float pDEta = fabs(iEta0-iEta1);
  float pDPhi = fabs(iPhi0-iPhi1); if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}
float GenPartLoader::mother(TGenParticle* iPart, std::vector<TGenParticle*> &iPartons) { 
  int iId = -1;
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
  if(fabs(iPhi0-iPhi1) > 2.*TMath::Pi()-fabs(iPhi0-iPhi1) && iPhi0-iPhi1 < 0) return 2.*TMath::Pi() -iPhi0+iPhi1; 
  if(fabs(iPhi0-iPhi1) > 2.*TMath::Pi()-fabs(iPhi0-iPhi1) && iPhi0-iPhi1 > 0) return -2.*TMath::Pi()-iPhi0+iPhi1; 
  return 0;
}
//Debug Method
void GenPartLoader::parentage(TGenParticle* iPart,TJet *iJet,std::vector<TGenParticle*> &iPartons) { 
  //const std::unordered_set<unsigned> partonIDs = {1, 2, 3, 4, 5, 21,15,11,12,411,421,433,441,443,511,521,513,523,531,541};
  const std::unordered_set<unsigned> partonIDs = {1, 2, 3, 4, 5, 21,15,11,12};
  auto validPID = [partonIDs](int id) -> bool {
    return (partonIDs.find(abs(id)) != partonIDs.end());
  };
  const std::unordered_set<unsigned> resonanceIDs = {6, 23, 24, 25, 54, 55}; 
  auto isResonance = [resonanceIDs](int id) -> bool {
    return ( (resonanceIDs.find(abs(id)) != resonanceIDs.end()) );
    //return ( (abs(id) > 10000) || // e.g. Z'
    //         (resonanceIDs.find(abs(id)) != resonanceIDs.end()) );
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
  bool lNoHF      = true;
  while(lParentId > 1) {
    //if(isResonance(pParent->pdgId)) lFromUE = false;
    double pZ     = pParent->pt/iJet->genpt;
    if(pZ > 0.1 && validPID(pParent->pdgId) && !lFromRes) lPassPtCut = true;
    if(isResonance(pParent->pdgId) && !lPassPtCut) lFromRes = true;
    pParent =  (TGenParticle*)((*fGens)[lParentId]);
    if(pZ > 0.1 && (isC(pParent->pdgId) || isB(pParent->pdgId)) ) lNoHF = false;
    lParentId = parentid(pParent); 
    lNCount++; 
  }
  if(lFromRes || lNoHF) return;
  lParentId = parentid(iPart); 
  pParent =  (TGenParticle*)((*fGens)[lParentId]);
  std::cout << "===> Missing Parent -parent- " << iPart->parent << " -id- " << iPart->pdgId << " -pt- " << iPart->pt  << " -- " << lNCount << " -- " << lFromRes << std::endl; 
  while(lParentId > 1) {
    pParent =  (TGenParticle*)((*fGens)[lParentId]);
    std::cout <<" ====> Searching Parent " << lParentId << " -id- " << pParent->pdgId << " -pt- " << pParent->pt << " -eta- " << pParent->eta << " -phi- " << pParent->phi << " - jet " << iJet->geneta << " - " << iJet->genphi  << " -dr- " <<  deltaR(pParent->phi,pParent->y,iJet->genphi,iJet->geneta) << " -chec- " << find(iPartons,pParent) << std::endl;
    lParentId = parentid(pParent,true); 
    //printVtx(pParent,iJet); 
  }
  std::cout <<" ===> Done " << pParent->parent << std::endl;
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
  int lId = iPart->parent;
  if(iOutput) std::cout <<" ===> seraching " << iPart->pdgId  << " -parent- " << lId << std::endl;
  if(lId > 1 && lId != -2 && !isB(iPart->pdgId) && !isC(iPart->pdgId)) return lId;
  /*
  if(isB(iPart->pdgId) && lId > 1) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[lId]);
    if(fabs(pGen->pdgId) == 5) return lId;
    }
  */
  int lPdgId = 0;
  if(isC(iPart->pdgId)) lPdgId = 4;
  if(isB(iPart->pdgId)) lPdgId = 5;
  if(lPdgId == 0) return 0;
  lId = 0; 
  bool lUseCharge = true; if(abs(iPart->pdgId) > 1000) lUseCharge = false; // Ignore the charge for the baryons  b/c it is complicated
  int  lCharge    = 1;    if(lPdgId == 5) lCharge = -1; //Correct for the charge assignement of the mesons
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {   
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(pGen == iPart)  break;
    if(abs(pGen->pdgId) != lPdgId || (int(iPart->pdgId/abs(iPart->pdgId)) != int(lCharge*(pGen->pdgId/abs(pGen->pdgId))) && lUseCharge) ) continue;
    if(iOutput) std::cout <<" ====> -id- " << i0 << " -s- " << lPdgId << " -pdgid- " << pGen->pdgId << " -b- " << iPart->pdgId << " -- " << pGen->pt  << " -pt- " << iPart->pt << " -- " << pGen->eta << " -- " << pGen->phi <<  " -- " << pGen->parent << " -dr- " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) << std::endl;
    if(deltaR(pGen->phi,pGen->y,iPart->phi,iPart->y) > fDRHeavy || pGen->pt < 0.7*iPart->pt) continue;
    //std::cout <<" ======> pass " << deltaR(pGen->phi,pGen->y,iPart->phi,iPart->eta) << std::endl;
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

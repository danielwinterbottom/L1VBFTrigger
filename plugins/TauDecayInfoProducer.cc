#include "TriggerStudies/L1VBFTrigger/plugins/TauDecayInfoProducer.hh"
#include <string>
#include <vector>
#include <bitset>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

using namespace edm;
using namespace std;

TauDecayInfoProducer::TauDecayInfoProducer(const edm::ParameterSet& config)
 { 
  ps = config;
  
  edm::InputTag inputTag_GenParticleCollection = ps.getUntrackedParameter<edm::InputTag>("inputTag_GenParticleCollection",edm::InputTag("genParticles")); 
  m_inputTag_GenParticleCollection = consumes<reco::GenParticleCollection>(inputTag_GenParticleCollection); 

}

TauDecayInfoProducer::~TauDecayInfoProducer() {}

void TauDecayInfoProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {

  Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(m_inputTag_GenParticleCollection, genParticles);
  
  higgsTauTauDecayMode = 0;
  
  unsigned tauCount  = 0;
  unsigned elecCount = 0;
  unsigned muCount   = 0;
  
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    
    if(fabs(p.pdgId())==15){
    
        
      if(p.mother()->pdgId() == 25){
          
        bool isMu   = false;
        bool isElec = false;
          
        unsigned n = p.numberOfDaughters();
        for(size_t j = 0; j < n; ++ j) {
          const reco::GenParticle *d = (reco::GenParticle*) p.daughter( j );

          
          if(fabs(d->pdgId())==11)      isElec = true;
          else if(fabs(d->pdgId())==13) isMu = true;
        }
        
        if(isElec) elecCount++;
        if(isMu) muCount++;
        tauCount++;
      }
    }
  }
  
  if     (elecCount == 2 && muCount == 0 && tauCount == 2) higgsTauTauDecayMode = 1;
  else if(elecCount == 1 && muCount == 1 && tauCount == 2) higgsTauTauDecayMode = 2;
  else if(elecCount == 1 && muCount == 0 && tauCount == 2) higgsTauTauDecayMode = 3;
  else if(elecCount == 0 && muCount == 2 && tauCount == 2) higgsTauTauDecayMode = 4;
  else if(elecCount == 0 && muCount == 1 && tauCount == 2) higgsTauTauDecayMode = 5;
  else if(elecCount == 0 && muCount == 0 && tauCount == 2) higgsTauTauDecayMode = 6;
  else                                                     higgsTauTauDecayMode = 7;
  
}

void TauDecayInfoProducer::beginJob() {
  ic::StaticTree::tree_->Branch("higgsTauTauDecayMode", &higgsTauTauDecayMode);

}

void TauDecayInfoProducer::endJob() {}

// define this as a plug-in
DEFINE_FWK_MODULE(TauDecayInfoProducer);

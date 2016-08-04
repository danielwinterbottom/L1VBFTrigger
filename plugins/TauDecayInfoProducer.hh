#ifndef L1VBFTrigger_TauDecayInfoProducer_hh
#define L1VBFTrigger_TauDecayInfoProducer_hh
#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TTree.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "../interface/StaticTree.hh"

/**
 * @brief See documentation [here](\ref objs-genparticle)
 */
class TauDecayInfoProducer : public edm::EDProducer {
 public:
  explicit TauDecayInfoProducer(const edm::ParameterSet &);
  ~TauDecayInfoProducer();

 private:
  virtual void beginJob();
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  
  edm::ParameterSet ps;
  
  TTree * tree_;
  
  // Input tags

  edm::EDGetTokenT< reco::GenParticleCollection > m_inputTag_GenParticleCollection;
  
  unsigned higgsTauTauDecayMode;
  
};

#endif

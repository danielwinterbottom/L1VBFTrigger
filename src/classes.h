#include <vector>
#include <map>
#include <utility>
#include <string>
#include "TriggerStudies/L1VBFTrigger/interface/Candidate.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TObject.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TEGamma.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TMuon.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TTau.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TJet.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TSum.hh"

namespace { struct dictionary {
  ic::Candidate dummy1;
  std::vector<ic::Candidate> dummy2;
  ic::L1TObject              dictL1TObject;
  std::vector<ic::L1TObject> dictL1TObjectCollection;
  ic::L1TEGamma              dictL1TEGamma;
  std::vector<ic::L1TEGamma> dictL1TEGammaCollection;
  ic::L1TMuon                dictL1TMuon;
  std::vector<ic::L1TMuon>   dictL1TMuonCollection;
  ic::L1TTau                 dictL1TTau;
  std::vector<ic::L1TTau>    dictL1TTauCollection;
  ic::L1TJet                 dictL1TJet;
  std::vector<ic::L1TJet>    dictL1TJetCollection;
  ic::L1TSum                 dictL1TSum;
  std::vector<ic::L1TSum>    dictL1TSumCollection;
  
};
}


#ifndef TriggerStudies_L1VBFTrigger_L1TJet_hh
#define TriggerStudies_L1VBFTrigger_L1TJet_hh

#include "TriggerStudies/L1VBFTrigger/interface/L1TObject.hh"

#include <vector>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"

namespace ic {
  
  class L1TJet : public ic::L1TObject {
  public:
    L1TJet();
    virtual ~L1TJet();
    virtual void Print() const;
    
    #ifndef SKIP_CINT_DICT
  public:
    ClassDef(L1TJet, 1);
    #endif
  };
  
  typedef std::vector<ic::L1TJet> L1TJetCollection;
}

#endif

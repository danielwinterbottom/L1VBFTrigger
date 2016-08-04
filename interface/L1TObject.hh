#ifndef TriggerStudies_L1VBFTrigger_L1TObject_hh
#define TriggerStudies_L1VBFTrigger_L1TObject_hh

#include "TriggerStudies/L1VBFTrigger/interface/Candidate.hh"

#include <vector>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"

namespace ic {

  class L1TObject : public ic::Candidate {
  public:
    L1TObject();
    virtual ~L1TObject();
    
  #ifndef SKIP_CINT_DICT
  public:
    ClassDef(L1TObject, 1);
  #endif
};

typedef std::vector<ic::L1TObject> L1TObjectCollection;
}

#endif

#ifdef __CINT__
#include <utility>
#include <string>

#pragma link off all globals;
#pragma link off all classes;
// #pragma link C++ all typedef;
#pragma link off all functions;
// #pragma link off nestedtypedef;
#pragma link C++ nestedclasses;

#pragma link C++ class std::pair<std::string, bool>+;
#pragma link C++ class std::pair<unsigned long, float>+;

#pragma link C++ class std::pair<std::string,ic::Candidate>+;
#pragma link C++ class std::vector<std::pair<std::string,ic::Candidate> >+;

#pragma link C++ class ic::Candidate+;
#pragma link C++ class std::vector<ic::Candidate>+;

#pragma link C++ class ic::L1TObject+;
#pragma link C++ class std::vector<ic::L1TObject>+;

#pragma link C++ class ic::L1TMuon+;
#pragma link C++ class std::vector<ic::L1TMuon>+;

#pragma link C++ class ic::L1TTau+;
#pragma link C++ class std::vector<ic::L1TTau>+;

#pragma link C++ class ic::L1TEGamma+;
#pragma link C++ class std::vector<ic::L1TEGamma>+;

#pragma link C++ class ic::L1TSum+;
#pragma link C++ class std::vector<ic::L1TSum>+;

#pragma link C++ class ic::L1TJet+;
#pragma link C++ class std::vector<ic::L1TJet>+;


#endif

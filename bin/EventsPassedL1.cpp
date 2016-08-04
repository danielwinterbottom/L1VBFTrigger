// ROOT includes
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TBranch.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <string>
#include <cstdlib>
#include <math.h> 

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"

// ICHiggsTauTau Objects
#include "TriggerStudies/L1VBFTrigger/interface/L1TObject.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TEGamma.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TMuon.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TTau.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TJet.hh"
#include "TriggerStudies/L1VBFTrigger/interface/L1TSum.hh"
#include "TriggerStudies/L1VBFTrigger/interface/Candidate.hh"

struct greater_Candidate{
  bool operator() (const ic::Candidate *a,const ic::Candidate *b) {
    return (a->pt() > b->pt());
  }
};

struct ProgramOptions{
  ProgramOptions(){
    isData         = false;
    useEmu         = true;
    input          = "input.dat";
    minEvent       = 0;
    maxEvent       = -1;
    maxEvents      = -1;
    cutsInput      = "CutsInput.txt";
    cutsInputLine  = 1;
    outputFilename = "output.root";
    decayType      = 0;
  }

  bool                  isData;
  bool                  useEmu;
  std::string           input;
  int                   minEvent;
  int                   maxEvent;
  std::vector<unsigned> runs;
  int                   maxEvents;
  std::string           cutsInput;
  unsigned              cutsInputLine;
  std::string           outputFilename;
  unsigned              decayType;
  
};

bool processArgs(int argc, char* argv[], ProgramOptions &opt){
  
  if(argc < 7){
    std::cout << "Wrong number of arguments, expected minimum 7 arguments." << std::endl;
    return 0;
  }
  
  for (int i=1; i < argc; ++i) {
    std::string arg = argv[i];
    if     (arg == "--isData") {opt.isData  = true;}
    else if(arg == "--useEmu") {opt.useEmu  = true;}
    else if(arg == "--input"){
      if(i+1<argc){i++; opt.input = argv[i];}
    } else if(arg == "--minEvent"){
      if(i+1<argc){i++; opt.minEvent = std::atoi(argv[i]);}
    }else if(arg == "--maxEvent"){
      if(i+1<argc){i++; opt.maxEvent = std::atoi(argv[i]);}
    }else if(arg == "--cutsInput"){
      if(i+1<argc){i++; opt.cutsInput = argv[i];}
    }else if(arg == "--cutsInputLine"){
      if(i+1<argc){i++; opt.cutsInputLine = std::atoi(argv[i]);}
    }else if(arg == "--maxEvents"){
      if(i+1<argc){i++; opt.maxEvents = atoi(argv[i]);}
    }else if(arg == "--runs"){
      if(i+1<argc){
        i++; 
        std::stringstream inVar(argv[i]);
        unsigned i;
        while(inVar >> i){
          opt.runs.push_back(i);
          if(inVar.peek() == ','){inVar.ignore();}
        }
      }
    }else if(arg == "--outputFilename"){
      if(i+1<argc){i++; opt.outputFilename = argv[i];}
    }else if(arg == "--decay"){
        if(i+1<argc){i++; 
          if     (std::string(argv[i]) == "EleEle"){options.decayType = 1}
          else if(std::string(argv[i]) == "EleMuo"){options.decayType = 2}
          else if(std::string(argv[i]) == "EleHad"){options.decayType = 3}
          else if(std::string(argv[i]) == "MuoMuo"){options.decayType = 4}
          else if(std::string(argv[i]) == "MuoHad"){options.decayType = 5}
          else if(std::string(argv[i]) == "HadHad"){options.decayType = 6}        
        }
  }
  
  return 1;
}

int main(int argc, char* argv[]){
    
  ProgramOptions options;
  if(!processArgs(argc, argv, options)){return -1;}
  
  std::cout << "===== program options ====="         << std::endl;
  std::cout << "isData         = " << options.isData         << std::endl;
  std::cout << "useEmu         = " << options.useEmu         << std::endl;
  std::cout << "input          = " << options.input          << std::endl;
  std::cout << "minEvent       = " << options.minEvent       << std::endl;
  std::cout << "maxEvent       = " << options.maxEvent       << std::endl;
  std::cout << "maxEvents      = " << options.maxEvents      << std::endl;
  std::cout << "cutsInput      = " << options.cutsInput      << std::endl;
  std::cout << "cutsInputLine  = " << options.cutsInputLine  << std::endl;
  std::cout << "outputFilename = " << options.outputFilename << std::endl;
  
  bool SpecifiedRun = false;
  
  if(argc > 8 && options.isData){
    std::cout << "===== runs =====" << std::endl;
    for(unsigned i=0; i<options.runs.size(); ++i){
      std::cout << options.runs[i] << std::endl;
      SpecifiedRun = true;
    }
  }
  
  TFile *fOut = new TFile(options.outputFilename.c_str(),"RECREATE");
  
  std::vector<TH1D*> myHists;
  
  if(options.isData){
    for(unsigned i=0; i<options.runs.size(); ++i){
      std::string histName = Form("h_eventspassed_%u",options.runs[i]);
      TH1D *hist = new TH1D(histName.c_str(),histName.c_str(),3,0,3);
      myHists.push_back(hist);
    }
  }
  else{
    TH1D *hist;
    hist->SetName("h_eventspassed");
    myHists.push_back(hist);
  }
  
  TString filename1;
  TChain *chIn1 = new TChain("l1EventTree/L1EventTree");
  TChain *chIn2;
  if(options.useEmu) chIn2 = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  else chIn2 = new TChain("l1UpgradeTree/L1UpgradeTree");
  TChain *chIn3 new TChain("icEventProducer/EventTree");
  
  std::ifstream infile1;
  infile1.open(options.input);
  while (infile1 >> filename1){
    chIn1->Add(filename1);
    chIn2->Add(filename1);
  }
  infile1.close();
  
  L1Analysis::L1AnalysisEventDataFormat     *productL1Event   = 0;
  L1Analysis::L1AnalysisL1UpgradeDataFormat *productL1Upgrade = 0;
  
  productL1Event = new L1Analysis::L1AnalysisEventDataFormat();
  chIn1->SetBranchAddress("Event",&productL1Event);
  
  productL1Upgrade= new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  chIn2->SetBranchAddress("L1Upgrade",&productL1Upgrade);
  
  unsigned decay;
  chIn3->SetBranchAddress("higgsTauTauDecayMode",&decay);
  
  double EGPtCut; double MuPtCut; double Tau1PtCut; double Tau2PtCut; double Jet1PtCut; double Jet2PtCut; double Jet3PtCut; bool EGIso; bool MuIso; bool TauIso; double EGEta; double MuEta; double TauEta; double MjjPtMin1; double MjjPtMin2; double MjjCut;
  std::ifstream infile2;
  infile2.open (options.cutsInput);
  unsigned linecount = 1;
  while(infile2 >> EGPtCut >> MuPtCut >> Tau1PtCut >> Tau2PtCut >> Jet1PtCut >> Jet2PtCut >> Jet3PtCut >> EGIso >> MuIso >> TauIso >> EGEta >> MuEta >> TauEta >> MjjPtMin1 >> MjjPtMin2 >> MjjCut){
    if(linecount == options.cutsInputLine){ 
      break;
    } else{
      ++linecount;
    }
  }
  infile2.close();
  std::cout << "===== L1 cuts =====" << std::endl;
  
  std::cout << "Leading jet Pt cut = " << Jet1PtCut << std::endl;
  std::cout << "Sub-leading jet Pt = " << Jet2PtCut << std::endl;
  std::cout << "Mjj jet1 min Pt =    " << MjjPtMin1 << std::endl;
  std::cout << "Mjj jet2 min Pt =    " << MjjPtMin2 << std::endl;
  std::cout << "Mjj cut =            " << MjjCut    << std::endl;
  
  if(EGPtCut > 0){
    std::cout << "EG Pt cut =          "               << EGPtCut << std::endl;
    if(EGEta > 0) std::cout << "EG max eta =         " << EGEta   << std::endl;
    std::cout << "Isolated EG's =                    ";
    if(EGIso > 0) std::cout << "true" << std::endl;
    else          std::cout << "false" << std::endl;
  }
  if(MuPtCut > 0){
    std::cout << "Muon Pt cut =        "               << MuPtCut << std::endl;
    if(MuEta > 0) std::cout << "Muons max eta =      " << MuEta   << std::endl;
    std::cout << "Isolated muons's =   ";
    if(MuIso > 0) std::cout << "true" << std::endl;
    else          std::cout << "false" << std::endl;
  }
  if(Tau1PtCut > 0 || Tau2PtCut > 0){
    std::cout << "Leading tau Pt cut =     "                << Tau1PtCut << std::endl;
    std::cout << "Sub-leading tau Pt cut = "                << Tau2PtCut << std::endl;
    if(TauEta > 0) std::cout << "Tau max eta =            " << TauEta    << std::endl;
    std::cout << "Isolated taus's =                       ";
    if(TauIso > 0) std::cout << "true" << std::endl;
    else std::cout << "false" << std::endl;
  }
  
  unsigned nentries1;
  if(options.maxEvent != -1){
    nentries1 = options.maxEvent;  
  }else if(options.maxEvents == -1){
    nentries1 = chIn1->GetEntries();
  }else{
    nentries1 = (unsigned)options.maxEvents;
  }
  if(nentries1 > chIn1->GetEntries()) nentries1 = chIn1->GetEntries();
  
  std::cout << "Running over " << nentries1 << " events" << std::endl;
  
  unsigned n_report = 10000;
  unsigned totaleventcount = 0;
  
  for(unsigned event=(unsigned)options.minEvent; event<nentries1; event++){
    
    chIn1->GetEntry(event);
    chIn2->GetEntry(event);
    totaleventcount++;
    
    unsigned histIndex = 0;
    
    if(SpecifiedRun && options.isData){
      unsigned run = productL1Event->run;
      bool isSpecifiedRun = false;
      for(unsigned i=0; i<options.runs.size(); ++i){
        if(run == options.runs[i]){
          histIndex = i;
          isSpecifiedRun = true;
        }
      }
      if(!isSpecifiedRun){
        continue;
      }
    }
    
    if(!isData && options.decayType !=0){
      if(options.decayType != decay){
        continue;
      }
    }
    
    myHists[histIndex]->Fill(0);
      
    std::vector<ic::L1TEGamma*> l1electrons;
    unsigned short int nL1TEGamma = productL1Upgrade->nEGs;
    for(unsigned i=0; i<nL1TEGamma; i++){
      
      if(productL1Upgrade->egBx[i] != 0){continue;}
      
      ic::L1TEGamma *thisEG = new ic::L1TEGamma;
      thisEG->set_pt    (productL1Upgrade->egEt    [i]);
      thisEG->set_eta   (productL1Upgrade->egEta   [i]);
      thisEG->set_phi   (productL1Upgrade->egPhi   [i]);
      double theta = 2*atan(exp(-1*productL1Upgrade->egEta[i]));
      double Energy = productL1Upgrade->egEt[i]/sin(theta);
      thisEG->set_energy(Energy);
      thisEG->isolation = productL1Upgrade->egIso[i];
      
      l1electrons.push_back(thisEG);
    }
    std::vector<ic::L1TEGamma*> l1electrons2 = std::vector<ic::L1TEGamma*>(l1electrons);
    
    std::vector<ic::L1TMuon*> l1muons;
    unsigned short int nL1TMuon = productL1Upgrade->nMuons;
    for(unsigned i=0; i<nL1TMuon; i++){
      
      if(productL1Upgrade->muonBx[i] != 0){continue;}
      
      ic::L1TMuon *thisMuon = new ic::L1TMuon();
      thisMuon->set_pt    (productL1Upgrade->muonEt    [i]);
      thisMuon->set_eta   (productL1Upgrade->muonEta   [i]);
      thisMuon->set_phi   (productL1Upgrade->muonPhi   [i]);
      thisMuon->quality = productL1Upgrade->muonQual [i];
      double theta = 2*atan(exp(-1*productL1Upgrade->muonEta[i]));
      double Energy = productL1Upgrade->muonEt[i]/sin(theta);
      thisMuon->set_energy(Energy);
      thisMuon->isolation = productL1Upgrade->muonIso  [i];
      
      l1muons.push_back(thisMuon);
    
    }
    std::vector<ic::L1TMuon*> l1muons2 = std::vector<ic::L1TMuon*>(l1muons);
    
    std::vector<ic::L1TTau*> l1taus;
    unsigned short int nL1TTaus = productL1Upgrade->nTaus;
    for(unsigned i=0; i<nL1TTaus; i++){
      
      if(productL1Upgrade->tauBx[i] != 0){continue;}
    
      ic::L1TTau *thisTau = new ic::L1TTau;
      thisTau->set_pt    (productL1Upgrade->tauEt    [i]);
      thisTau->set_eta   (productL1Upgrade->tauEta   [i]);
      thisTau->set_phi   (productL1Upgrade->tauPhi   [i]);
      double theta = 2*atan(exp(-1*productL1Upgrade->tauEta[i]));
      double Energy = productL1Upgrade->tauEt[i]/sin(theta);
      thisTau->set_energy(Energy);
      thisTau->isolation = productL1Upgrade->tauIso  [i];
      
      l1taus.push_back(thisTau);
    
    }
    
    std::vector<ic::L1TJet*> l1jets;
    unsigned short int nL1TJets = productL1Upgrade->nJets;
    for(unsigned i=0; i<nL1TJets; i++){
      
      if(productL1Upgrade->jetBx[i] != 0){continue;}
      
      ic::L1TJet *thisJet = new ic::L1TJet();
      thisJet->set_pt    (productL1Upgrade->jetEt    [i]);
      thisJet->set_eta   (productL1Upgrade->jetEta   [i]);
      thisJet->set_phi   (productL1Upgrade->jetPhi   [i]);
      double theta = 2*atan(exp(-1*productL1Upgrade->jetEta[i]));
      double Energy = productL1Upgrade->jetEt[i]/sin(theta);
      thisJet->set_energy(Energy);
      
      l1jets.push_back(thisJet);
    }    
    
    std::vector<ic::L1TSum*> l1sums;
    unsigned short int nL1TSums = productL1Upgrade->nSums;
    for(unsigned i=0; i<nL1TSums; i++){
      
      if(productL1Upgrade->sumBx[i] != 0){continue;}
      
      ic::L1TSum *thisSum = new ic::L1TSum;
      thisSum->set_pt    (productL1Upgrade->sumEt [i]);
      thisSum->set_eta   (0);
      thisSum->set_phi   (productL1Upgrade->sumPhi[i]);
      thisSum->set_energy(productL1Upgrade->sumEt [i]);
    
      short int type = productL1Upgrade->sumType[i];
      if     (type==L1Analysis::EtSumType::kTotalEt)  {thisSum->sumType = ic::L1TSum::SumType::kTotalEt;}
      else if(type==L1Analysis::EtSumType::kTotalHt)  {thisSum->sumType = ic::L1TSum::SumType::kTotalHt;}
      else if(type==L1Analysis::EtSumType::kMissingEt){thisSum->sumType = ic::L1TSum::SumType::kMissingEt;}
      else if(type==L1Analysis::EtSumType::kMissingHt){thisSum->sumType = ic::L1TSum::SumType::kMissingHt;}
      else if(type==L1Analysis::EtSumType::kTotalEtx) {thisSum->sumType = ic::L1TSum::SumType::kTotalEtx;}
      else if(type==L1Analysis::EtSumType::kTotalEty) {thisSum->sumType = ic::L1TSum::SumType::kTotalEty;}
      else if(type==L1Analysis::EtSumType::kTotalHtx) {thisSum->sumType = ic::L1TSum::SumType::kTotalHtx;}
      else if(type==L1Analysis::EtSumType::kTotalHty) {thisSum->sumType = ic::L1TSum::SumType::kTotalHty;}
    
      l1sums.push_back(thisSum);
    }
    
    unsigned n_l1jets_ = l1jets.size();
    unsigned n_l1electrons_ = l1electrons.size();
    unsigned n_l1muons_ = l1muons.size();
    unsigned n_l1taus_ = l1taus.size();
      
    std::sort(l1electrons.begin(),l1electrons.end(),greater_Candidate());
    std::sort(l1muons.begin(),l1muons.end(),greater_Candidate());
    std::sort(l1taus.begin(),l1taus.end(),greater_Candidate());
    std::sort(l1jets.begin(),l1jets.end(),greater_Candidate());
    
    for(int i= n_l1muons_-1; i>=0; i--){
        if(l1muons[i]->quality < 12) l1muons.erase(l1muons.begin()+i);
        n_l1muons_ = l1muons.size();
    }
    
    if(TauIso){
        for(int i= n_l1taus_-1; i>=0; i--) if(l1taus[i]->isolation == 0) l1taus.erase(l1taus.begin()+i);
        n_l1taus_ = l1taus.size();
    }
    
    if(EGIso){
        for(int i= n_l1electrons_-1; i>=0; i--) if(l1electrons[i]->isolation == 0) l1electrons.erase(l1electrons.begin()+i);
        n_l1electrons_ = l1electrons.size();
    }
    if(MuIso){
        for(int i= n_l1muons_-1; i>=0; i--) if(l1muons[i]->isolation == 0) l1muons.erase(l1muons.begin()+i);
        n_l1muons_ = l1muons.size();
    }
    
    if(TauEta > 0){
        for(int i= n_l1taus_-1; i>=0; i--) if(std::fabs(l1taus[i]->vector().Rapidity()) > TauEta) l1taus.erase(l1taus.begin()+i);
        n_l1taus_ = l1taus.size();
    }
    
    if(EGEta > 0){
        for(int i= n_l1electrons_-1; i>=0; i--) if(std::fabs(l1electrons[i]->vector().Rapidity()) > EGEta) l1electrons.erase(l1electrons.begin()+i);
        n_l1electrons_ = l1electrons.size();
    }
    if(MuEta > 0){
        for(int i= n_l1muons_-1; i>=0; i--) if(std::fabs(l1muons[i]->vector().Rapidity()) > MuEta) l1muons.erase(l1muons.begin()+i);
        n_l1muons_ = l1muons.size();
    }
    
    bool Filter = false;
    
    if(EGPtCut > 0){
      if(n_l1electrons_ < 1) Filter = true;
      else if (l1electrons[0]->vector().Pt() < EGPtCut) Filter = true;
    }
    if(MuPtCut > 0){
      if(n_l1muons_ < 1) Filter = true;
      else if (l1muons[0]->vector().Pt() < MuPtCut) Filter = true;    
    }
    if(Tau1PtCut > 0){
      if(n_l1taus_ < 1) Filter = true;
      else if (l1taus[0]->vector().Pt() < Tau1PtCut) Filter = true;    
    }
    if(Tau2PtCut > 0){
      if(n_l1taus_ < 2) Filter = true;
      else if (l1taus[1]->vector().Pt() < Tau2PtCut) Filter = true;    
    }
    if(Jet1PtCut > 0){
      if(n_l1jets_ < 1) Filter = true;
      else if (l1jets[0]->vector().Pt() < Jet1PtCut) Filter = true;    
    }
    if(Jet2PtCut > 0){
      if(n_l1jets_ < 2) Filter = true;
      else if (l1jets[1]->vector().Pt() < Jet2PtCut) Filter = true;    
    }
    if(Jet3PtCut > 0){
      if(n_l1jets_ < 3) Filter = true;
      else if (l1jets[2]->vector().Pt() < Jet3PtCut) Filter = true;    
    }
    
    double LargestMjj = -1;
    if(MjjCut > 0){
      for(unsigned i=0; i<n_l1jets_; ++i){
        for(unsigned j=i+1; j < n_l1jets_; ++j){
          if(l1jets[i]->vector().Pt() >= MjjPtMin1 && l1jets[j]->vector().Pt() >= MjjPtMin2){
            double Mjj = (l1jets[i]->vector()+l1jets[j]->vector()).M();
            if(Mjj > LargestMjj) LargestMjj = Mjj;    
          }
        }
      }
      if(LargestMjj < MjjCut) Filter = true;  
    }
    
    if(!Filter){
    
      myHists[histIndex]->Fill(1);
     
      //Check if event fired any other trigger
     
      //SingleJet170
      bool PassedSingleJet170 = false;
      if(n_l1jets_ >= 1){ 
        if(l1jets[0]->vector().Pt() >= 170) {
          PassedSingleJet170 = true;  
        }
      }
      //DoubleJetC100
      bool PassedDoubleJet100 = false;
      unsigned Cjetcount1=0;
      for(unsigned i = 0; i<n_l1jets_; i++){
        if(std::fabs(l1jets[i]->vector().Rapidity())<=3 && l1jets[i]->vector().Pt() >= 100) Cjetcount1++;    
      }
      if(Cjetcount1 >= 2) {
        PassedDoubleJet100 = true;
      }
     
      //QuadJetC40
      bool PassedQuadJetC50 = false;
      unsigned Cjetcount=0;
      for(unsigned i = 0; i<n_l1jets_; i++){
        if(std::fabs(l1jets[i]->vector().Rapidity())<=3 && l1jets[i]->vector().Pt() >= 50) Cjetcount++;    
      }
      if(Cjetcount >= 4) {
        PassedQuadJetC50 = true;
      }
      
      //DoubleIsoTau28er
      bool PassedDoubleIsoTau28er = false;
      unsigned taucount=0;
      for(unsigned i = 0; i< n_l1taus_; i++){
        if(std::fabs(l1taus[i]->vector().Rapidity())<=2.1 && l1taus[i]->vector().Pt() >= 28 && l1taus[i]->isolation !=0) taucount++;    
      }
      if(taucount >=2) {
        PassedDoubleIsoTau28er = true;  
      }
     
      //HTT300
      
      double HTT = -1;
      
      for(unsigned i=0; i< l1sums.size(); i++){  
        if(l1sums[i]->sumType==1) HTT = l1sums[i]->vector().Pt();
      }
      
      bool PassedHTT300 = false;
      if(HTT >= 300) {
        PassedHTT300 = true; 
      }
     
      //TripleJet_84_68_48
      bool PassedTripleJet_88_72_56 = false;
      std::vector<ic::L1TJet> centraljets;
      std::vector<ic::L1TJet> forwardsjets;
      for(unsigned i = 0; i<n_l1jets_; i++){
        if(std::fabs(l1jets[i]->vector().Rapidity()) <= 3.) centraljets.push_back(*l1jets[i]);
        else forwardsjets.push_back(*l1jets[i]);
      }
      if(centraljets.size() >=3){
        if(centraljets[0].vector().Pt() >= 88 && centraljets[1].vector().Pt() > 72 && centraljets[2].vector().Pt() >= 56) PassedTripleJet_88_72_56 = true;
      }
      if(centraljets.size() >=2 && forwardsjets.size() >=1){
        if     (centraljets[0].vector().Pt() >= 88 && centraljets[1].vector().Pt() > 72 && forwardsjets[0].vector().Pt() >= 56) PassedTripleJet_88_72_56 = true;
        else if(centraljets[0].vector().Pt() >= 88 && centraljets[1].vector().Pt() > 56 && forwardsjets[0].vector().Pt() >= 72) PassedTripleJet_88_72_56 = true;
        else if(centraljets[0].vector().Pt() >= 72 && centraljets[1].vector().Pt() > 56 && forwardsjets[0].vector().Pt() >= 88) PassedTripleJet_88_72_56 = true;
      }
      
      //SingleMu22, SingleMu20er, Mu6_HTT200
      bool PassedSingleMu22 = false;
      bool PassedSingleMu20er = false;
      bool Mu6;
      for(unsigned i = 0; i<l1muons2.size(); i++){
        if(l1muons2[i]->vector().Pt() >= 22 && l1muons2[i]->quality>=12) PassedSingleMu22 = true;
        if(l1muons2[i]->vector().Pt() >= 20 && l1muons2[i]->quality>=12 && std::fabs(l1muons2[i]->vector().Rapidity()) <= 2.1) PassedSingleMu20er = true;
        if(l1muons2[i]->vector().Pt() >= 6 && l1muons2[i]->quality>=12) Mu6 = true;
      }
      
      //Mu6HTT200
      bool PassedMu6HTT200 = false;
      if(HTT >= 200 && Mu6) {
        PassedMu6HTT200 = true; 
      }
     
      //SingleEG34, SingleIsoEG28, SingleIsoEG26er, 
      
      bool PassedEG34 = false;
      bool PassedIsoEG28 = false;
      bool PassedIsoEG26er = false;
      
      for(unsigned i = 0; i<l1electrons2.size(); i++){
        if(l1electrons2[i]->vector().Pt() >= 34 ) PassedEG34 = true;
        if(l1electrons2[i]->vector().Pt() >= 26 && l1electrons2[i]->isolation != 0 ) PassedIsoEG28 = true;
        if(l1electrons2[i]->vector().Pt() >= 26 && l1electrons2[i]->isolation != 0 && std::fabs(l1electrons2[i]->vector().Rapidity()) <= 2.1) PassedIsoEG26er = true;
      }
      
      //MET80
      bool PassedMET80 = false;
      for(unsigned i=0; i<l1sums.size(); ++i){
        double MET = 0;  
        if(l1sums[i]->sumType==2) MET = l1sums[i]->vector().Pt();
        if(MET >= 80) PassedMET80 = true;
      }
      
      if(!(PassedSingleJet170 || PassedDoubleJet100 || PassedQuadJetC50 || PassedDoubleIsoTau28er || PassedHTT300 || PassedTripleJet_88_72_56 ||  PassedSingleMu22 || PassedSingleMu20er || PassedMu6HTT200 || PassedSingleMu22 || PassedSingleMu20er || PassedMu6HTT200 || PassedEG34 || PassedIsoEG26er || PassedIsoEG28 || PassedMET80)) myHists[histIndex]->Fill(2);
      
    }
    
    if(event % n_report == 0) std::cout << "Finished processing " << event << "th event." << std::endl;
    
    for(unsigned obj=0; obj<l1electrons.size();    obj++){delete l1electrons[obj];}
    for(unsigned obj=0; obj<l1muons.size();        obj++){delete l1muons[obj];}
    for(unsigned obj=0; obj<l1taus.size();         obj++){delete l1taus[obj];}
    for(unsigned obj=0; obj<l1jets.size();         obj++){delete l1jets[obj];}
    for(unsigned obj=0; obj<l1sums.size();         obj++){delete l1sums[obj];}
  }
  
  std::cout << "events processed " << totaleventcount << std::endl;

  fOut->Write();  
  fOut->Close();

  return 0;  

}

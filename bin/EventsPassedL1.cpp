// ROOT includes
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TBranch.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <string>
#include <cstdlib>
#include <math.h> 
#include <set>

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

struct info {
  unsigned run;
  std::vector<double> instlumis;
  std::vector<unsigned> lumisections;
};

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
    minFile        = 1;
    maxFile        = -1;
    maxEvents      = -1;
    cutsInput      = "CutsInput.txt";
    cutsInputLine  = 1;
    outputFilename = "output.root";
    decayType      = 0;
    doRunSummary   = false;
    makeHistograms = false;
    doBinnedRates  = true;
    
  }

  bool                  isData;
  bool                  useEmu;
  std::string           input;
  int                   minFile;
  int                   maxFile;
  std::vector<unsigned> runs;
  int                   maxEvents;
  std::string           cutsInput;
  unsigned              cutsInputLine;
  std::string           outputFilename;
  unsigned              decayType;
  bool makeHistograms;
  bool doRunSummary;
  bool doBinnedRates;
  
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
    else if(arg == "--makeHistograms") {opt.makeHistograms  = true;}
    else if(arg == "--doRunSummary") {opt.doRunSummary  = true;}
    else if(arg == "--doBinnedRates") {opt.doBinnedRates = true;}
    else if(arg == "--input"){
      if(i+1<argc){i++; opt.input = argv[i];}
    } else if(arg == "--minFile"){
      if(i+1<argc){i++; opt.minFile = std::atoi(argv[i]);}
    }else if(arg == "--maxFile"){
      if(i+1<argc){i++; opt.maxFile = std::atoi(argv[i]);}
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
    }else if(arg == "--decayType"){
        if(i+1<argc){i++; 
          if     (std::string(argv[i]) == "EleEle"){opt.decayType = 1;}
          else if(std::string(argv[i]) == "EleMuo"){opt.decayType = 2;}
          else if(std::string(argv[i]) == "EleHad"){opt.decayType = 3;}
          else if(std::string(argv[i]) == "MuoMuo"){opt.decayType = 4;}
          else if(std::string(argv[i]) == "MuoHad"){opt.decayType = 5;}
          else if(std::string(argv[i]) == "HadHad"){opt.decayType = 6;}        
        }
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
  std::cout << "minFile        = " << options.minFile       << std::endl;
  std::cout << "maxFile        = " << options.maxFile       << std::endl;
  std::cout << "maxEvents      = " << options.maxEvents      << std::endl;
  std::cout << "cutsInput      = " << options.cutsInput      << std::endl;
  std::cout << "cutsInputLine  = " << options.cutsInputLine  << std::endl;
  std::cout << "outputFilename = " << options.outputFilename << std::endl;
  std::cout << "decayType      = " << options.decayType      << std::endl;
  std::cout << "makeHistograms = " << options.makeHistograms << std::endl;
  std::cout << "doRunSummary   = " << options.doRunSummary   << std::endl;
  std::cout << "doBinnedRates  = " << options.doBinnedRates  << std::endl;
  
  TFile *fOut = new TFile(options.outputFilename.c_str(),"RECREATE");

  bool SpecifiedRun = false;
  
  if(argc > 8 && options.isData){
    std::cout << "===== runs =====" << std::endl;
    for(unsigned i=0; i<options.runs.size(); ++i){
      std::cout << options.runs[i] << std::endl;
      SpecifiedRun = true;
    }
  }
  
  std::vector<TH1D*> myHists;
  std::vector<TH1D*> leadJetVary;
  std::vector<TH1D*> subleadJetVary;
  std::vector<TH1D*> MjjVary;
  std::vector<TH1D*> MjjMinVary;
  std::vector<TH1D*> MjjMaxVary;
  std::vector<TH1D*> muVary;
  TH1D *lumibinned_num = new TH1D("lumibinned_num","lumibinned_num",20,0,2);
  TH1D *lumibinned_denum = new TH1D("lumibinned_denum","lumibinned_denum",20,0,2);
  TH1D *gains = new TH1D("gains","gains",12,0,12);
  lumibinned_num->Sumw2();
  //lumibinned_denum->Sumw2();
  
  TH1D *lumibinned_num_pure = new TH1D("lumibinned_num_pure","lumibinned_num_pure",20,0,2);
  lumibinned_num_pure->Sumw2();
  
  TH2D *hist1_2D = new TH2D("leadJetVary",    "leadJetVary"   ,20,0,2,201,-1,200);
  TH2D *hist2_2D = new TH2D("hsubleadJetVary","subleadJetVary",20,0,2,201,-1,200);
  TH2D *hist3_2D = new TH2D("hMjjVary",       "MjjVary"       ,20,0,2,901,-1,900);
  TH2D *hist4_2D = new TH2D("hMjjMinVary",    "MjjMinVary"    ,20,0,2,201,-1,200);
  TH2D *hist5_2D = new TH2D("hMjjMaxVary",    "MjjMaxVary"    ,20,0,2,201,-1,200);
  TH2D *hist6_2D = new TH2D("hMuVary",        "MuVar"         ,20,0,2,31,-1,30);
  
  hist1_2D->Sumw2();
  hist2_2D->Sumw2();
  hist3_2D->Sumw2();
  hist4_2D->Sumw2();
  hist5_2D->Sumw2();
  hist6_2D->Sumw2();
  
  if(options.isData && SpecifiedRun){
    for(unsigned i=0; i<options.runs.size(); ++i){
      std::string histName = Form("h_eventspassed_%u",options.runs[i]);
      TH1D *hist = new TH1D(histName.c_str(),histName.c_str(),3,0,3);
      myHists.push_back(hist);
      
      if(options.makeHistograms){
        std::string histName1 = Form("h_leadJetVary_%u",options.runs[i]);
        std::string histName2 = Form("h_subleadJetVary_%u",options.runs[i]);
        std::string histName3 = Form("h_MjjVary_%u",options.runs[i]);
        std::string histName4 = Form("h_MjjMinVary_%u",options.runs[i]);
        std::string histName5 = Form("h_MjjMaxVary_%u",options.runs[i]);
        std::string histName6 = Form("h_MuVary_%u",options.runs[i]);
        
        TH1D *hist1 = new TH1D(histName1.c_str(),histName1.c_str(),200,0,200);
        TH1D *hist2 = new TH1D(histName2.c_str(),histName2.c_str(),200,0,200);
        TH1D *hist3 = new TH1D(histName3.c_str(),histName3.c_str(),900,0,900);
        TH1D *hist4 = new TH1D(histName4.c_str(),histName4.c_str(),200,0,200);
        TH1D *hist5 = new TH1D(histName5.c_str(),histName5.c_str(),200,0,200);
        TH1D *hist6 = new TH1D(histName6.c_str(),histName5.c_str(),30,0,30);
        
        leadJetVary.push_back(hist1);
        subleadJetVary.push_back(hist2);
        MjjVary.push_back(hist3);
        MjjMinVary.push_back(hist4);
        MjjMaxVary.push_back(hist5);
        muVary.push_back(hist6);
        
      }
  
    }
  }
  else {
    TH1D *hist = new TH1D("h_eventspassed","h_eventspassed",3,0,3);
    myHists.push_back(hist);
    
    if(options.makeHistograms){
      TH1D *hist1 = new TH1D("h_leadJetVary","h_leadJetVary",200,0,200);
      TH1D *hist2 = new TH1D("h_subleadJetVary","h_subleadJetVary",200,0,200);
      TH1D *hist3 = new TH1D("h_MjjVary","h_MjjVary",900,0,900);
      TH1D *hist4 = new TH1D("h_MjjMinVary","h_MjjMinVary",200,0,200);
      TH1D *hist5 = new TH1D("h_MjjMaxVary","h_MjjMaxVary",200,0,200);
      TH1D *hist6 = new TH1D("h_MuVary","h_MuVar",30,0,30);
      
      leadJetVary.push_back(hist1);
      subleadJetVary.push_back(hist2);
      MjjVary.push_back(hist3);
      MjjMinVary.push_back(hist4);
      MjjMaxVary.push_back(hist5);
      muVary.push_back(hist6);
    }
  }
  
  std::vector<info> runInfos;
  if(options.doBinnedRates){
    for(std::vector<unsigned>::iterator it = options.runs.begin(); it != options.runs.end(); ++it){
      info runInfo;
      runInfo.run = *it;
      std::string lumifilename = Form("XMLLumis/%u_Lumi.dat",*it);
      std::cout<< lumifilename << std::endl;
      std::ifstream lumifile;
      lumifile.open(lumifilename);
      unsigned lumiSection;
      double instlumi;
      while (lumifile >> lumiSection >> instlumi){
        runInfo.lumisections.push_back(lumiSection);
        runInfo.instlumis.push_back(instlumi);
      }
      runInfos.push_back(runInfo);
      lumifile.close();
    }
  }
  
  TString filename1;
  TChain *chIn1 = new TChain("l1EventTree/L1EventTree");
  TChain *chIn2;
  if(options.useEmu) chIn2 = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  else chIn2 = new TChain("l1UpgradeTree/L1UpgradeTree");
  TChain *chIn3 = new TChain();
  if(!options.isData) chIn3 = new TChain("icEventProducer/EventTree");
  
  std::ifstream infile1;
  infile1.open(options.input);
  int fileCount=1;
  while (infile1 >> filename1){
    if(fileCount >= options.minFile && fileCount < options.maxFile){
      std::cout << "Adding file: " << filename1 << std::endl;  
      chIn1->Add(filename1);
      chIn2->Add(filename1);
      if(!options.isData) chIn3->Add(filename1);
    }
    fileCount++;
  }
  infile1.close();
  
  L1Analysis::L1AnalysisEventDataFormat     *productL1Event   = 0;
  L1Analysis::L1AnalysisL1UpgradeDataFormat *productL1Upgrade = 0;
  
  productL1Event = new L1Analysis::L1AnalysisEventDataFormat();
  chIn1->SetBranchAddress("Event",&productL1Event);
  
  productL1Upgrade= new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  chIn2->SetBranchAddress("L1Upgrade",&productL1Upgrade);
  
  unsigned decayValue = 0;
  if(!options.isData) chIn3->SetBranchAddress("higgsTauTauDecayMode",&decayValue);
  
  double EGPtCut; double MuPtCut; double Tau1PtCut; double Tau2PtCut; double Jet1PtCut; double Jet2PtCut; double Jet3PtCut; bool EGIso; bool MuIso; bool TauIso; double EGEta; double MuEta; double TauEta; double Jet1Eta; double MjjPtMin1; double MjjPtMin2; double MjjCut;
  std::ifstream infile2;
  infile2.open (options.cutsInput);
  unsigned linecount = 1;
  while(infile2 >> EGPtCut >> MuPtCut >> Tau1PtCut >> Tau2PtCut >> Jet1PtCut >> Jet2PtCut >> Jet3PtCut >> EGIso >> MuIso >> TauIso >> EGEta >> MuEta >> TauEta >> Jet1Eta >> MjjPtMin1 >> MjjPtMin2 >> MjjCut){
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
  if(options.maxEvents == -1){
    nentries1 = chIn1->GetEntries();
  }else{
    nentries1 = (unsigned)options.maxEvents;
  }
  if(nentries1 > chIn1->GetEntries()) nentries1 = chIn1->GetEntries();
  
  std::cout << "Running over " << nentries1 << " events" << std::endl;
  
  unsigned n_report = 10000;
  unsigned totaleventcount = 0;
  
  std::set<unsigned> summary_runs;
  
  for(unsigned event=0; event<nentries1; event++){
      
    if(options.doRunSummary){
      unsigned run = productL1Event->run;
      summary_runs.insert(run);
    }
      
    if(event % n_report == 0) std::cout << "Processing " << event << "th event." << std::endl;
    
    chIn1->GetEntry(event);
    chIn2->GetEntry(event);
    if(!options.isData) chIn3->GetEntry(event);
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
    
    if(!options.isData && options.decayType !=0){
      if(options.decayType != decayValue){
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
    
    std::vector<ic::L1TTau*> l1taus_filtered(l1taus);
    std::vector<ic::L1TJet*> l1jets_filtered(l1jets);
    std::vector<ic::L1TEGamma*> l1electrons_filtered(l1electrons);
    std::vector<ic::L1TMuon*> l1muons_filtered(l1muons);
    std::vector<ic::L1TJet*> l1jets_etafiltered(l1jets);
    
    
    for(int i= l1muons_filtered.size()-1; i>=0; i--){
        if(l1muons_filtered[i]->quality < 12) l1muons_filtered.erase(l1muons_filtered.begin()+i);
    }
    
    if(TauIso){
        for(int i= l1taus_filtered.size()-1; i>=0; i--) if(l1taus_filtered[i]->isolation == 0) l1taus_filtered.erase(l1taus_filtered.begin()+i);
    }
    
    if(EGIso){
        for(int i= l1electrons_filtered.size()-1; i>=0; i--) if(l1electrons_filtered[i]->isolation == 0) l1electrons_filtered.erase(l1electrons_filtered.begin()+i);
    }
    if(MuIso){
        for(int i= l1muons_filtered.size()-1; i>=0; i--) if(l1muons_filtered[i]->isolation == 0) l1muons_filtered.erase(l1muons_filtered.begin()+i);
    }
    
    if(TauEta > 0){
        for(int i= l1taus_filtered.size()-1; i>=0; i--) if(std::fabs(l1taus_filtered[i]->vector().Rapidity()) > TauEta) l1taus_filtered.erase(l1taus_filtered.begin()+i);
    }
    
    if(EGEta > 0){
        for(int i= l1electrons_filtered.size()-1; i>=0; i--) if(std::fabs(l1electrons_filtered[i]->vector().Rapidity()) > EGEta) l1electrons_filtered.erase(l1electrons_filtered.begin()+i);
    }
    if(MuEta > 0){
        for(int i= l1muons_filtered.size()-1; i>=0; i--) if(std::fabs(l1muons_filtered[i]->vector().Rapidity()) > MuEta) l1muons_filtered.erase(l1muons_filtered.begin()+i);
    }
    if(Jet1Eta > 0){
      for(int i= l1jets_etafiltered.size()-1; i>=0; i--) if(std::fabs(l1jets_etafiltered[i]->vector().Rapidity()) > Jet1Eta) l1jets_etafiltered.erase(l1jets_etafiltered.begin()+i);
    }
    
    n_l1electrons_ = l1electrons_filtered.size();
    n_l1muons_     = l1muons_filtered.size();
    n_l1taus_      = l1taus_filtered.size();
    n_l1jets_      = l1jets_filtered.size();
    
    bool Filter = false;
    
    if(EGPtCut > 0){
      if(n_l1electrons_ < 1){ Filter = true;}
      else if (l1electrons_filtered[0]->vector().Pt() < EGPtCut){ Filter = true;}
    }
    if(MuPtCut > 0){
      if(n_l1muons_ < 1){ Filter = true;}
      else if (l1muons_filtered[0]->vector().Pt() < MuPtCut){
          //std::cout << l1muons_filtered[0]->vector().Pt() << std::endl;
          //std::cout << l1muons_filtered[0]->quality << std::endl;
          //std::cout << l1muons_filtered[0]->vector().Rapidity() << std::endl;
          Filter = true;}    
    }
    if(Tau1PtCut > 0){
      if(n_l1taus_ < 1){ Filter = true;}
      else if (l1taus_filtered[0]->vector().Pt() < Tau1PtCut){ Filter = true;}   
    }
    if(Tau2PtCut > 0){
      if(n_l1taus_ < 2){ Filter = true;}
      else if (l1taus_filtered[1]->vector().Pt() < Tau2PtCut){ Filter = true;}    
    }
    if(Jet1PtCut > 0){
      if(l1jets_filtered.size() < 1){ Filter = true;}
      else if (l1jets_etafiltered[0]->vector().Pt() < Jet1PtCut){ Filter = true;}
    }
    if(Jet2PtCut > 0){
      if(n_l1jets_ < 2){ Filter = true;}
      else if (l1jets_filtered[1]->vector().Pt() < Jet2PtCut){ Filter = true;}
      if(Jet1Eta>0){
        Filter=true;
        for(unsigned i=0; i<l1jets_filtered.size(); ++i){
          for(unsigned j=i+1; j<l1jets_filtered.size(); ++j){
             if((l1jets_filtered[i]->vector().Pt() > Jet1PtCut && std::fabs(l1jets_filtered[i]->vector().Rapidity()) >= Jet1Eta && l1jets_filtered[j]->vector().Pt() > Jet2PtCut && std::fabs(l1jets_filtered[j]->vector().Rapidity()) < Jet1Eta) || (l1jets_filtered[i]->vector().Pt() > Jet1PtCut && std::fabs(l1jets_filtered[i]->vector().Rapidity()) < Jet1Eta && l1jets_filtered[j]->vector().Pt() > Jet2PtCut && std::fabs(l1jets_filtered[j]->vector().Rapidity()) >= Jet1Eta) || (l1jets_filtered[i]->vector().Pt() > Jet1PtCut && std::fabs(l1jets_filtered[i]->vector().Rapidity()) < Jet1Eta && l1jets_filtered[j]->vector().Pt() > Jet2PtCut && std::fabs(l1jets_filtered[j]->vector().Rapidity()) < Jet1Eta)) Filter = false;   
          }
        }    
      }
    }
    if(Jet3PtCut > 0){
      if(n_l1jets_ < 3){ Filter = true;}
      else if (l1jets_filtered[2]->vector().Pt() < Jet3PtCut){ Filter = true;}    
    }
    
    double LargestMjj = -1;
    double LargestMjjJet1Pt = -1;
    double LargestMjjJet2Pt = -1;
    if(MjjCut > 0){
      for(unsigned i=0; i<n_l1jets_; ++i){
        for(unsigned j=i+1; j < n_l1jets_; ++j){
          if(l1jets_filtered[i]->vector().Pt() >= MjjPtMin1 && l1jets_filtered[j]->vector().Pt() >= MjjPtMin2){
            double Mjj = (l1jets_filtered[i]->vector()+l1jets_filtered[j]->vector()).M();
            if(Mjj > LargestMjj){
                LargestMjj = Mjj;
                LargestMjjJet1Pt = l1jets_filtered[i]->vector().Pt();
                LargestMjjJet2Pt = l1jets_filtered[j]->vector().Pt();
            }
          }
        }
      }
      if(LargestMjj < MjjCut){ Filter = true;} 
      
      bool PassedDoubleIsoTau36er = false;
      bool PassedSingleEG = false;
      bool PassedSingleMu = false;
      bool PassedMET90 = false;
      
      unsigned taucount = 0;
      for(unsigned i = 0; i< l1taus.size(); i++){
        if(std::fabs(l1taus[i]->vector().Rapidity())<=2.1 && l1taus[i]->vector().Pt() >= 36 && l1taus[i]->isolation !=0){
          taucount++;
        }
      }
      if(taucount >=2) {
        PassedDoubleIsoTau36er = true;  
      }
      
      bool PassedEG36 = false;
      bool PassedIsoEG34 = false;
      bool PassedIsoEG32er = false;
      bool PassedTripleJet_92_76_64 = false;
      
      for(unsigned i = 0; i<l1electrons.size(); i++){
        if(l1electrons[i]->vector().Pt() >= 36 ) PassedEG36 = true;
        if(l1electrons[i]->vector().Pt() >= 32 && l1electrons[i]->isolation != 0 ) PassedIsoEG34 = true;
        if(l1electrons[i]->vector().Pt() >= 32 && l1electrons[i]->isolation != 0 && std::fabs(l1electrons[i]->vector().Rapidity()) <= 2.1) PassedIsoEG32er = true;
      }
      if(PassedIsoEG32er || PassedIsoEG34 || PassedEG36) PassedSingleEG = true;
      
      bool PassedSingleMu22 = false;
      bool PassedSingleMu20er = false;
      for(unsigned i = 0; i<l1muons.size(); i++){
        if(l1muons[i]->vector().Pt() >= 22 && l1muons[i]->quality>=12) PassedSingleMu22 = true;
        if(l1muons[i]->vector().Pt() >= 20 && l1muons[i]->quality>=12 && std::fabs(l1muons[i]->vector().Rapidity()) <= 2.1) PassedSingleMu20er = true;
      }
      
      if(PassedSingleMu20er || PassedSingleMu22) PassedSingleMu = true;
      
      for(unsigned i=0; i<l1sums.size(); ++i){
        double MET = 0;  
        if(l1sums[i]->sumType==2) MET = l1sums[i]->vector().Pt();
        if(MET >= 90) PassedMET90 = true;
      }
      
      std::vector<ic::L1TJet> centraljets;
      std::vector<ic::L1TJet> forwardsjets;
      for(unsigned i = 0; i<l1jets.size(); i++){
        if(std::fabs(l1jets[i]->vector().Rapidity()) <= 3.) centraljets.push_back(*l1jets[i]);
        else forwardsjets.push_back(*l1jets[i]);
      }
      if(centraljets.size() >=3){
        if(centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() > 76 && centraljets[2].vector().Pt() >= 64) PassedTripleJet_92_76_64 = true;
      }
      if(centraljets.size() >=2 && forwardsjets.size() >=1){
        if     (centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() >= 76 && forwardsjets[0].vector().Pt() >= 64) PassedTripleJet_92_76_64 = true;
        else if(centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() >= 64 && forwardsjets[0].vector().Pt() >= 76) PassedTripleJet_92_76_64 = true;
        else if(centraljets[0].vector().Pt() >= 76 && centraljets[1].vector().Pt() >= 64 && forwardsjets[0].vector().Pt() >= 88) PassedTripleJet_92_76_64 = true;
      }
      
      gains->Fill(0);
      if(!Filter) gains->Fill(1);
      if(PassedDoubleIsoTau36er){
        gains->Fill(2);
      } else if(!Filter) gains->Fill(3);
      if(PassedSingleEG){
        gains->Fill(4); 
      } else if(!Filter) gains->Fill(5);
      if(PassedSingleMu){
        gains->Fill(6);
      } else if(!Filter) gains->Fill(7);
      if(PassedMET90){
        gains->Fill(8);
      } else if(!Filter) gains->Fill(9);
      if(PassedTripleJet_92_76_64){
        gains->Fill(10);
      } else if(!Filter) gains->Fill(11);
    }

    /*
    bool FailedOR = true;
    for(unsigned i=0; i<l1muons_filtered.size(); ++i){
      if(l1muons_filtered[0]->vector().Pt()>=5) FailedOR = false;
    }
    for(unsigned i=0; i<l1electrons_filtered.size(); ++i){
      if(l1electrons_filtered[0]->vector().Pt()>=10) FailedOR = false;
    }
    for(unsigned i=0; i<l1taus_filtered.size(); ++i){
      if(l1taus_filtered[0]->vector().Pt()>=30) FailedOR = false;
    }
    if(n_l1jets_>3){
      if(l1jets[0]->vector().Pt()>=Jet1PtCut && l1jets[1]->vector().Pt() >=60 && l1jets[2]->vector().Pt()>=Jet2PtCut) FailedOR = false;    
    }
    for(unsigned i=0; i<l1sums.size(); ++i){
        double MET = 0;  
        if(l1sums[i]->sumType==2) MET = l1sums[i]->vector().Pt();
        if(MET >= 50) FailedOR = false;
      }
    if(FailedOR) Filter = true;
    */
    if(options.doBinnedRates){
      unsigned run = productL1Event->run;
      unsigned lumiSection = productL1Event->lumi;
      double lumi = -1;
      for(unsigned r=0; r<runInfos.size(); ++r){
        if(runInfos[r].run == run){
          for(unsigned l=0; l<runInfos[r].lumisections.size(); ++l){
              //std::cout << runInfos[r].lumisections[l] << std::endl;
            if(runInfos[r].lumisections[l] == lumiSection){
              lumi = runInfos[r].instlumis[l]/10000;
              //std::cout << lumi << std::endl;
              break;
            }
          }
          break;
        }
      }
      int NBunch = -1;
      if(run >= 273158 && run <= 273411) NBunch = 589;
      else if (run >= 277069 && run <= 277148) NBunch = 2064;
      double w = NBunch*11.247;
      if(NBunch != -1 && lumi != -1){
        lumibinned_denum->Fill(lumi);
        hist1_2D->Fill(lumi, -1);
        hist2_2D->Fill(lumi, -1);
        hist3_2D->Fill(lumi, -1);
        hist4_2D->Fill(lumi, -1);
        hist5_2D->Fill(lumi, -1);
        hist6_2D->Fill(lumi, -1);
        if(!Filter){
          lumibinned_num->Fill(lumi,w);
          hist1_2D->Fill(lumi, l1jets_filtered[0]->vector().Pt(),w);
          hist2_2D->Fill(lumi, l1jets_filtered[1]->vector().Pt(),w);
          hist3_2D->Fill(lumi, LargestMjj,w);
          hist4_2D->Fill(lumi, LargestMjjJet2Pt,w);
          hist5_2D->Fill(lumi, LargestMjjJet1Pt,w);
        if(MuPtCut > 0) hist6_2D->Fill(lumi, l1muons_filtered[0]->vector().Pt(),w);
        }
      }
    }
    
    if(!Filter){
    
      myHists[histIndex]->Fill(1);
      
      if(options.makeHistograms){
        leadJetVary   [histIndex]->Fill(l1jets_filtered[0]->vector().Pt());
        subleadJetVary[histIndex]->Fill(l1jets_filtered[1]->vector().Pt());
        MjjVary       [histIndex]->Fill(LargestMjj);
        MjjMinVary    [histIndex]->Fill(LargestMjjJet2Pt);
        MjjMaxVary    [histIndex]->Fill(LargestMjjJet1Pt);
        if(MuPtCut > 0) muVary        [histIndex]->Fill(l1muons_filtered[0]->vector().Pt());
      }
      
      //Check if event fired any other trigger to estimate pure rate
     
      //SingleJet180
      //std::cout << "SingleJet" << std::endl;
      bool PassedSingleJet180 = false;
      if(l1jets.size()>= 1){ 
        if(l1jets[0]->vector().Pt() >= 180) {
          PassedSingleJet180 = true;
        //  std::cout<< l1jets[0]->vector().Pt() << std::endl;;
        }
      }
      //DoubleJetC112
      //std::cout << "DoubleJet" << std::endl;
      bool PassedDoubleJet112 = false;
      unsigned Cjetcount1=0;
      for(unsigned i = 0; i<l1jets.size(); i++){
        if(std::fabs(l1jets[i]->vector().Rapidity())<=3 && l1jets[i]->vector().Pt() >= 112){
          Cjetcount1++;
         //std::cout << l1jets[i]->vector().Rapidity() << std::endl;
         //std::cout << l1jets[i]->vector().Pt() << std::endl;
        }
      }
      if(Cjetcount1 >= 2) {
        PassedDoubleJet112 = true;
      }
     
      //QuadJetC50
      bool PassedQuadJetC50 = false;
      unsigned Cjetcount=0;
      for(unsigned i = 0; i<l1jets.size(); i++){
        if(std::fabs(l1jets[i]->vector().Rapidity())<=3 && l1jets[i]->vector().Pt() >= 50) Cjetcount++;    
      }
      if(Cjetcount >= 4) {
        PassedQuadJetC50 = true;
      }
      
      //DoubleIsoTau36er
      bool PassedDoubleIsoTau36er = false;
      bool PassedSingleTau120er = false;
      bool PassedDoubleTau50er = false;
      unsigned taucount=0;
      unsigned taucount2=0;
      //std::cout << "DoubleIsoTau" << std::endl;
      for(unsigned i = 0; i< l1taus.size(); i++){
        if(std::fabs(l1taus[i]->vector().Rapidity())<=2.1 && l1taus[i]->vector().Pt() >= 36 && l1taus[i]->isolation !=0){
          taucount++;
         // std::cout << l1taus[i]->vector().Rapidity() << std::endl;
         // std::cout << l1taus[i]->isolation << std::endl;
         // std::cout << l1taus[i]->vector().Pt() << std::endl;
        }
        if(std::fabs(l1taus[i]->vector().Rapidity())<=2.1 && l1taus[i]->vector().Pt() >= 50) taucount2++;
        if(std::fabs(l1taus[i]->vector().Rapidity())<=2.1 && l1taus[i]->vector().Pt() >= 120) PassedSingleTau120er = true;
      }
      if(taucount >=2) {
        PassedDoubleIsoTau36er = true;  
      }
      if(taucount2 >=2) PassedDoubleTau50er = true;
     
      //HTT300
      
      double HTT = -1;
      
      for(unsigned i=0; i< l1sums.size(); i++){  
        if(l1sums[i]->sumType==1) HTT = l1sums[i]->vector().Pt();
      }
      
      bool PassedHTT300 = false;
      if(HTT >= 300) {
        PassedHTT300 = true; 
      }
     
      //TripleJet_92_76_64
      bool PassedTripleJet_92_76_64 = false;
      std::vector<ic::L1TJet> centraljets;
      std::vector<ic::L1TJet> forwardsjets;
      for(unsigned i = 0; i<l1jets.size(); i++){
        if(std::fabs(l1jets[i]->vector().Rapidity()) <= 3.) centraljets.push_back(*l1jets[i]);
        else forwardsjets.push_back(*l1jets[i]);
      }
      if(centraljets.size() >=3){
        if(centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() > 76 && centraljets[2].vector().Pt() >= 64) PassedTripleJet_92_76_64 = true;
      }
      if(centraljets.size() >=2 && forwardsjets.size() >=1){
        if     (centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() >= 76 && forwardsjets[0].vector().Pt() >= 64) PassedTripleJet_92_76_64 = true;
        else if(centraljets[0].vector().Pt() >= 92 && centraljets[1].vector().Pt() >= 64 && forwardsjets[0].vector().Pt() >= 76) PassedTripleJet_92_76_64 = true;
        else if(centraljets[0].vector().Pt() >= 76 && centraljets[1].vector().Pt() >= 64 && forwardsjets[0].vector().Pt() >= 88) PassedTripleJet_92_76_64 = true;
      }
      
      //std::cout << "Passed Triple Jet " << std::endl;
      //if(PassedTripleJet_92_76_64){
      //  std::cout << "Central Jets:" << std::endl;
      //  std::cout << centraljets[0].vector().Pt() << "  " << centraljets[1].vector().Pt() << centraljets[2].vector().Pt() << std::endl;
      //  std::cout << forwardsjets[0].vector().Pt() << "  " << forwardsjets[1].vector().Pt() << forwardsjets[2].vector().Pt() << std::endl;
      //}
      
      //SingleMu22, SingleMu20er, Mu6_HTT200
      bool PassedSingleMu22 = false;
      bool PassedSingleMu20er = false;
      bool Mu6 = false;
      bool PassedDoubleMu_12_5 = false;
      bool PassedDoubleMu0er_MaxdEta = false;
      for(unsigned i = 0; i<l1muons.size(); i++){
        if(l1muons[i]->vector().Pt() >= 22 && l1muons[i]->quality>=12) PassedSingleMu22 = true;
        if(l1muons[i]->vector().Pt() >= 20 && l1muons[i]->quality>=12 && std::fabs(l1muons[i]->vector().Rapidity()) <= 2.1) PassedSingleMu20er = true;
        if(l1muons[i]->vector().Pt() >= 6 && l1muons[i]->quality>=12) Mu6 = true;
        if(l1muons[i]->vector().Pt() >= 12 && l1muons[i]->quality>=8){
          for(unsigned j=i+1; j<l1muons.size(); j++){
            if(l1muons[j]->vector().Pt() >= 5 && l1muons[i]->quality>=8) PassedDoubleMu_12_5 = true; 
          }
        }
        
        if(l1muons[i]->vector().Pt() >= 0 && l1muons[i]->quality>=8 && std::fabs(l1muons[i]->vector().Rapidity()) <= 1.4){
          for(unsigned j=i+1; j<l1muons.size(); j++){
            if(l1muons[j]->vector().Pt() >= 0 && l1muons[i]->quality>=8 && std::fabs(l1muons[i]->vector().Rapidity()) <= 1.4){
              int charge = l1muons[i]->charge*l1muons[j]->charge;
              double dEta = std::fabs(l1muons[i]->vector().Rapidity() - l1muons[j]->vector().Rapidity());
              if (charge <0 && dEta<=1.8) PassedDoubleMu0er_MaxdEta = true;
            }
          }
        }
      }
      
      //Mu6HTT200
      bool PassedMu6HTT200 = false;
      if(HTT >= 200 && Mu6) {
        PassedMu6HTT200 = true; 
      }
     
      //SingleEG36, SingleIsoEG36, SingleIsoEG32er, 
      
      bool PassedEG36 = false;
      bool PassedIsoEG34 = false;
      bool PassedIsoEG32er = false;
      bool PassedDoubleEG_24_17 = false;
      bool PassedDoubleEG_25_12 = false;
      
      for(unsigned i = 0; i<l1electrons.size(); i++){
        if(l1electrons[i]->vector().Pt() >= 36 ) PassedEG36 = true;
        if(l1electrons[i]->vector().Pt() >= 32 && l1electrons[i]->isolation != 0 ) PassedIsoEG34 = true;
        if(l1electrons[i]->vector().Pt() >= 32 && l1electrons[i]->isolation != 0 && std::fabs(l1electrons[i]->vector().Rapidity()) <= 2.1) PassedIsoEG32er = true;
        if(l1electrons[i]->vector().Pt() >= 24 ){
          for(unsigned j = i+1; j<l1electrons.size(); j++){
            if(l1electrons[j]->vector().Pt() >= 17) PassedDoubleEG_24_17 = true; 
          }
        }
        if(l1electrons[i]->vector().Pt() >= 25 ){
          for(unsigned j = i+1; j<l1electrons.size(); j++){
            if(l1electrons[j]->vector().Pt() >= 12) PassedDoubleEG_25_12 = true;    
          }
        }
      }
      
      //MET90
      bool PassedMET90 = false;
      bool MET55 = false;
      for(unsigned i=0; i<l1sums.size(); ++i){
        double MET = 0;  
        if(l1sums[i]->sumType==2) MET = l1sums[i]->vector().Pt();
        if(MET >= 90) PassedMET90 = true;
        if(MET >= 55) MET55 = true;
      }
      
      bool Mu20 = false;
      bool IsoEG6 = false;
      bool PassedEG22Tau26 = false;
      bool PassedDoubleMuMET = false;
      bool PassedMuTau = false;
      bool PassedJetMu = false;
      bool PassedEGMu = false;
      unsigned mucount1 = 0;
      for(unsigned i=0; i<l1muons.size(); ++i){
        if(l1muons[i]->quality >=12 && l1muons[i]->vector().Pt() >= 20) Mu20 = true; 
        if(l1muons[i]->quality >=8) mucount1++;
        if(l1muons[i]->quality >=12 && l1muons[i]->vector().Pt() >= 18 && std::fabs(l1muons[i]->vector().Rapidity())<=2.1){
          for(unsigned j=0; j< l1taus.size(); ++j){
            if(std::fabs(l1taus[j]->vector().Rapidity())<=2.1 && l1taus[j]->vector().Pt() >=20) PassedMuTau = true;
          }
        }
        if(l1muons[i]->quality >=12 && l1muons[i]->vector().Pt() >= 5 && std::fabs(l1muons[i]->vector().Rapidity())<=2.1){
          for(unsigned j=0; j< l1electrons.size(); ++j){
            if((std::fabs(l1electrons[j]->isolation)!=0 && l1electrons[j]->vector().Pt() >=20) || l1electrons[j]->vector().Pt() >=23) PassedEGMu = true;
          }
        }
        if(l1muons[i]->quality >=12 && l1muons[i]->vector().Pt() >= 3){
          for(unsigned j=0; j< l1jets.size(); ++j){
            if(l1jets[j]->vector().Pt() >=120){
                double dEta = std::fabs(l1jets[j]->vector().Rapidity() - l1muons[i]->vector().Rapidity());
                double dPhi = std::fabs(l1jets[j]->vector().Phi() - l1muons[i]->vector().Phi());
                if(dEta<=0.4 && dPhi<=0.4) PassedJetMu = true;
            }
          }
        }
      }
      if(mucount1 >=2 && MET55) PassedDoubleMuMET = true;
      for(unsigned i=0; i<l1electrons.size(); ++i){
        if(l1electrons[i]->isolation !=0 && l1electrons[i]->vector().Pt() >=6) IsoEG6 = true;
        if(l1electrons[i]->isolation !=0 && l1electrons[i]->vector().Pt() >=22 && std::fabs(l1electrons[i]->vector().Rapidity())<=2.1 ){
          for(unsigned j=0; j< l1taus.size(); ++j){
            if(l1taus[j]->isolation != 0 && std::fabs(l1taus[j]->vector().Rapidity())<=2.1 && l1taus[j]->vector().Pt() >=26){
              double dEta = std::fabs(l1taus[j]->vector().Rapidity() - l1electrons[i]->vector().Rapidity());
              if(dEta >= 0.2) PassedEG22Tau26 = true;   
            }
          }
            
        }
      }
      
      bool PassedMu20IsoEG6 = Mu20 && IsoEG6;
      bool PassedOthers = false;

      if(PassedSingleJet180 || PassedDoubleJet112 || PassedQuadJetC50 || PassedDoubleIsoTau36er || PassedHTT300 || PassedTripleJet_92_76_64 ||  PassedSingleMu22 || PassedSingleMu20er || PassedMu6HTT200 || PassedSingleMu22 || PassedSingleMu20er || PassedMu6HTT200 || PassedEG36 || PassedIsoEG32er || PassedIsoEG34 || PassedMET90 || PassedDoubleEG_24_17 || PassedDoubleEG_25_12 || PassedSingleTau120er || PassedDoubleTau50er || PassedDoubleMu_12_5 || PassedDoubleMu0er_MaxdEta || PassedMu20IsoEG6 || PassedEG22Tau26 || PassedDoubleMuMET || PassedMuTau || PassedJetMu || PassedEGMu){
        PassedOthers = true;
      }
      if (!PassedOthers) myHists[histIndex]->Fill(2);
      
      if(options.doBinnedRates){
        unsigned run = productL1Event->run;
        unsigned lumiSection = productL1Event->lumi;
        double lumi = -1;
        for(unsigned r=0; r<runInfos.size(); ++r){
          if(runInfos[r].run == run){
            for(unsigned l=0; l<runInfos[r].lumisections.size(); ++l){
              if(runInfos[r].lumisections[l] == lumiSection){
                lumi = runInfos[r].instlumis[l]/10000;
                break;
              }
            }
            break;
          }
        }
        int NBunch = -1;
        if(run >= 273158 && run <= 273411) NBunch = 589;
        else if (run >= 277069 && run <= 277148) NBunch = 2064;
        double w = NBunch*11.247;
        if(NBunch != -1 && lumi != -1){
          if(!PassedOthers){
            lumibinned_num_pure->Fill(lumi,w);
          }
        }
      }
      
      /*if(!Filter){
        std::cout << "Iso Tau ";
        for(unsigned i =0; i<l1taus.size(); ++i){
          if(l1taus[i]->isolation !=0 && std::fabs(l1taus[i]->vector().Rapidity()) >=2.1){
            std::cout << l1taus[i]->vector().Pt() << std::endl;
            break;
          }
        }
        std::cout << "Tau ";
        for(unsigned i =0; i<l1taus.size(); ++i){
          if(std::fabs(l1taus[i]->vector().Rapidity()) >=2.1){
            std::cout << l1taus[i]->vector().Pt() << std::endl;
            break;
          }
        }
        std::cout << "Muon "; 
        for(unsigned i =0; i<l1muons.size(); ++i){
          if(l1muons[i]->quality>=12){
            std::cout << l1muons[i]->vector().Pt() << std::endl;
            break;
          }
        }
        std::cout << "EG "; 
        for(unsigned i =0; i<l1electrons.size(); ++i){
            std::cout << l1electrons[i]->vector().Pt() << std::endl;
            break;
        }
        std::cout << "Iso EG "; 
        for(unsigned i =0; i<l1electrons.size(); ++i){
          if(l1electrons[i]->isolation !=0){
            std::cout << l1electrons[i]->vector().Pt() << std::endl;
            break;
          }
        }
      }*/
      
    }
    

  
    for(unsigned obj=0; obj<l1electrons.size();    obj++){delete l1electrons[obj];}
    for(unsigned obj=0; obj<l1muons.size();        obj++){delete l1muons[obj];}
    for(unsigned obj=0; obj<l1taus.size();         obj++){delete l1taus[obj];}
    for(unsigned obj=0; obj<l1jets.size();         obj++){delete l1jets[obj];}
    for(unsigned obj=0; obj<l1sums.size();         obj++){delete l1sums[obj];}
    l1electrons.clear();
    l1muons.clear();
    l1jets.clear();
    l1taus.clear();
    l1sums.clear();
    l1electrons_filtered.clear();
    l1muons_filtered.clear();
    l1jets_filtered.clear();
    l1taus_filtered.clear();
  
  }
  
  std::cout << "events processed " << totaleventcount << std::endl;
  
  /*if(options.makeHistograms){
    for(unsigned i=0; i<myHists.size(); ++i){
      for(unsigned j=1; j<=(unsigned)leadJetVary[i]->GetNbinsX(); ++j){
        leadJetVary[i]->SetBinContent(j,leadJetVary[i]->Integral(j, leadJetVary[i]->GetNbinsX()+1));
      }
      for(unsigned j=1; j<=(unsigned)subleadJetVary[i]->GetNbinsX(); ++j){
        subleadJetVary[i]->SetBinContent(j,subleadJetVary[i]->Integral(j, subleadJetVary[i]->GetNbinsX()+1));
      }
      for(unsigned j=1; j<=(unsigned)MjjVary[i]->GetNbinsX(); ++j){
        MjjVary[i]->SetBinContent(j,MjjVary[i]->Integral(j, MjjVary[i]->GetNbinsX()+1));
      }
      for(unsigned j=1; j<=(unsigned)MjjMinVary[i]->GetNbinsX(); ++j){
        MjjMinVary[i]->SetBinContent(j,MjjMinVary[i]->Integral(j, MjjMinVary[i]->GetNbinsX()+1));
      }
      for(unsigned j=1; j<=(unsigned)MjjMaxVary[i]->GetNbinsX(); ++j){
        MjjMaxVary[i]->SetBinContent(j,MjjMaxVary[i]->Integral(j, MjjMaxVary[i]->GetNbinsX()+1));
      }
      for(unsigned j=1; j<=(unsigned)muVary[i]->GetNbinsX(); ++j){
        muVary[i]->SetBinContent(j,muVary[i]->Integral(j, muVary[i]->GetNbinsX()+1));
      }
      leadJetVary   [i]->Scale(1/myHists[i]->GetBinContent(2));
      subleadJetVary[i]->Scale(1/myHists[i]->GetBinContent(2));
      MjjVary       [i]->Scale(1/myHists[i]->GetBinContent(2));
      MjjMinVary    [i]->Scale(1/myHists[i]->GetBinContent(2));
      MjjMaxVary    [i]->Scale(1/myHists[i]->GetBinContent(2));
      muVary        [i]->Scale(1/myHists[i]->GetBinContent(2));
      
    }
  }*/
  
  if(options.doRunSummary){
    std::set<unsigned>::iterator it;
    std::ofstream outfile;
    outfile.open("runSummary.txt");
    for (it = summary_runs.begin(); it != summary_runs.end(); ++it){ 
      outfile << *it << "\n";  
    }
    outfile.close();
  }

  fOut->Write();  
  fOut->Close();

  return 0;  

}

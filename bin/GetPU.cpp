#include <iostream>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
    
int main(int argc, char* argv[]){

  if(argc != 2){
    std::cout << "Wrong number of arguments, expected 1 argument." << std::endl;
    return 0;
  }

  std::string run = argv[1];
  
  std::string line;
  std::string infilename = run+"_Lumi.dat";
  std::ifstream in (infilename);
  
  std::string outfilename = run+"_LumiAndPU.dat";
  std::ifstream outfile (outfilename);

  unsigned lumisection;
  double lumi;
  
  while (in >> lumisection >> lumi)
  { 
    std::ofstream out ("temp.txt");  
    std::string line = Form("{\"273158\": [[%u, %u]]}", lumisection);
    out << line << std::endl;
    out.close();
    
    std::string command = "pileupCalc.py -i temp.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80  MyDataPileupHistogram.root";
    system(command.c_str());
    
    TFile *f0 = new TFile("MyDataPileupHistogram.root");
    TH1D *h0 = (TH1D*)f0->Get(hist1);
    double PU = h0->GetMean();
    
    outfile << lumisection << "  " << lumi << "  " << PU << std::endl;
    
  }
  std::string rmcommand = "rm temp.txt";
  system(rmcommand.c_str());
  in.close();
  outfile.close();
}

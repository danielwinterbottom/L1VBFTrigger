#include <iostream>
#include <string>
#include <fstream>
    
int main(int argc, char* argv[]){

  if(argc != 2){
    std::cout << "Wrong number of arguments, expected 1 argument." << std::endl;
    return 0;
  }

  std::string run = argv[1];
  
  std::string line;
  std::string infilename = run+".xml";
  std::string outfilename = run+"_Lumi.dat";
  std::ifstream in(infilename);
  std::ofstream out(outfilename);
  
  std::cout << "Starting scan of XML file.." << std::endl;

  while (std::getline(in,line))
  { 
    std::string lumistart = "<TD TITLE=\"instLumi\" ALIGN=CENTER><FONT SIZE=-3>";
    std::string lumiend = "</FONT></TD><TD TITLE=\"delivLumi\" ALIGN=CENTER><FONT SIZE=-3>";
    int l = line.find(lumiend);
    int f = line.find(lumistart);
    
    if(f!=-1 && l!=-1){
      std::string lumi = line.substr(f + lumistart.length(), l - f - lumistart.length());
      
      std::string lumisectionstart = "LUMISECTION=";
      std::string lumisectionend = "</A>";
      int i = line.find(lumisectionend);
      int j = line.find(lumisectionstart);

      if(i!=-1 && j!=-1){
        std::string lumisection = line.substr(j + lumisectionstart.length(), i - j - lumisectionstart.length());
        lumisection = lumisection.substr(0,(lumisection.length()-1)/2);
        out << lumisection << "  " << lumi << "\n";
        std::cout << lumisection << "    " << lumi << std::endl;
      }
    }

  }
  in.close();
  out.close();
}

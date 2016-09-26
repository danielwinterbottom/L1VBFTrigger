{
//std::vector<std::string> samples = {"HToTauTau_HadHad", "HToInvisible","Run2016E"};//,
//std::vector<std::string> bases = {"DoubleJet_90_30_Mj30j30_580","DoubleJet_90_30_Mj30j30_620","DoubleJet_94_30_Mj30j30_620","DoubleJet_97_30_Mj30j30_620"};
std::vector<std::string> samples = {"HToTauTau_HadHad","Run2016E"};
std::vector<std::string> bases = {"DoubleJet30_Mj30j30_360_Mu6","DoubleJet30_Mj30j30_400_Mu6"};
    
for(unsigned j=0; j<samples.size(); ++j){
  for(unsigned k=0; k<bases.size(); ++k){
      
    std::vector<std::string> histnames;
    std::vector<std::string> xlabels;
    std::string temp;
    std::string templabels;
    
    bool isData = false;
    if(samples[j] == "Run2016E") isData = true;
    
    std::string filename = "output/"+bases[k]+"/"+samples[j];
    if(samples[j] == "Run2016E"){
      filename+="_output";
      //if(bases[k] == "DoubleJet_90_30_Mj30j30_580") filename += "2";
      filename += "2";
    }
    filename += ".root";
    unsigned run = 277072;
    
    TFile *f0 = new TFile(filename.c_str());
    
    temp = "h_leadJetVary"   ;
    templabels = "lead jet Pt [GeV]";
    if(isData) temp += Form("_%u",run);
    histnames.push_back(temp);
    xlabels.push_back(templabels);
    temp = "h_subleadJetVary";
    templabels = "sub-lead jet Pt [GeV]";
    if(isData) temp += Form("_%u",run);
    histnames.push_back(temp);
    xlabels.push_back(templabels);
    temp = "h_MjjVary";
    templabels = "Mjj [GeV]";
    if(isData) temp += Form("_%u",run);
    histnames.push_back(temp);
    xlabels.push_back(templabels);
    temp = "h_MjjMinVary"    ;
    templabels = "Mjj jets min Pt [GeV]";
    if(isData) temp += Form("_%u",run);
    histnames.push_back(temp);
    xlabels.push_back(templabels);
    temp = "h_MjjMaxVary"    ;
    templabels = "Mjj jets max Pt [GeV]";
    if(isData) temp += Form("_%u",run);
    histnames.push_back(temp);
    xlabels.push_back(templabels);
    
    
    for(unsigned i=0; i<histnames.size(); ++i){
      TH1D *h1 = (TH1D*)f0->Get(histnames[i].c_str());
      h1->SetFillStyle(0);
      if(!isData) h1->GetYaxis()->SetTitle("Efficiency");
      else        h1->GetYaxis()->SetTitle("Rate");
      
      h1->GetXaxis()->SetTitle(xlabels[i].c_str());
      
      h1->Draw();
      
      TLatex *   tex = new TLatex(0.15,0.92,(samples[j]+", Base Algo: "+bases[k]).c_str()); 
      tex->SetNDC();
      tex->SetTextFont(44);
      tex->SetTextSize(23);
      tex->SetLineWidth(2);
      tex->Draw();

    
      c1->Print(("output/"+bases[k]+"/"+samples[j]+"_"+bases[k]+"_"+histnames[i]+".png").c_str());
    }
  }
}

}
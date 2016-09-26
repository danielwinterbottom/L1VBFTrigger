int SignalEff(unsigned k=1, unsigned l=0){
  std::string dirname;
  if(l == 1) dirname = "MuTau_Mu9";
  else if (l == 0) dirname = "TauTau_Jet97";
  TCanvas c1;
  std::string filename = Form("output/DoubleJet_TauTau%u.root",k);
  TFile *f0 = new TFile(filename.c_str());
  TH1D *h1 = (TH1D*)f0->Get("h_eventspassed");
  TH1D *leadJetVary = (TH1D*)f0->Get("h_leadJetVary");
  TH1D *subleadJetVary = (TH1D*)f0->Get("h_subleadJetVary");
  TH1D *MjjVary = (TH1D*)f0->Get("h_MjjVary");
  TH1D *muVary = (TH1D*)f0->Get("h_MuVary");
  
  for(unsigned j=1; j<=(unsigned)leadJetVary->GetNbinsX(); ++j){
    leadJetVary->SetBinContent(j,leadJetVary->Integral(j, leadJetVary->GetNbinsX()+1));
  }
  for(unsigned j=1; j<=(unsigned)subleadJetVary->GetNbinsX(); ++j){
    subleadJetVary->SetBinContent(j,subleadJetVary->Integral(j, subleadJetVary->GetNbinsX()+1));
  }
  for(unsigned j=1; j<=(unsigned)MjjVary->GetNbinsX(); ++j){
    MjjVary->SetBinContent(j,MjjVary->Integral(j, MjjVary->GetNbinsX()+1));
  }
  for(unsigned j=1; j<=(unsigned)muVary->GetNbinsX(); ++j){
    muVary->SetBinContent(j,muVary->Integral(j, muVary->GetNbinsX()+1));
  }
  leadJetVary   ->Scale(1/h1->GetBinContent(1));
  subleadJetVary->Scale(1/h1->GetBinContent(1));
  MjjVary       ->Scale(1/h1->GetBinContent(1));
  muVary        ->Scale(1/h1->GetBinContent(1));
  //subleadJetVary->Draw();
  
  TLatex *   tex1;
  
  double eff;
  std::string eff_string;
  
  muVary->SetFillStyle(0);
  muVary->GetXaxis()->SetTitle("L1 muon p_{T} cut [GeV]");
  muVary->GetYaxis()->SetTitle("Efficiency");
  eff = muVary->GetBinContent(muVary->GetXaxis()->FindBin(11));
  eff_string = Form ("Efficiency @ 11 GeV = %.2f", eff );
  muVary->Draw();
  tex1 = new TLatex(0.50,0.92, eff_string.c_str());
  tex1->SetNDC();   
  tex1->SetTextFont(44);   
  tex1->SetTextSize(23);   
  tex1->SetLineWidth(2);   
  tex1->Draw();
  c1.Print((dirname+"/muPt.pdf").c_str());
  c1.Print((dirname+"/muPt.png").c_str());
  
  MjjVary->SetFillStyle(0);
  MjjVary->GetXaxis()->SetTitle("L1 M_{jj} cut [GeV]");
  MjjVary->GetYaxis()->SetTitle("Efficiency");
  eff = MjjVary->GetBinContent(MjjVary->GetXaxis()->FindBin(700));
  eff_string = Form ("Efficiency @ 700 GeV = %.2f", eff );
  MjjVary->Draw();
  tex1 = new TLatex(0.50,0.92, eff_string.c_str());
  tex1->SetNDC();   
  tex1->SetTextFont(44);   
  tex1->SetTextSize(23);   
  tex1->SetLineWidth(2);   
  tex1->Draw();
  c1.Print((dirname+"/Mjj.pdf").c_str());
  c1.Print((dirname+"/Mjj.png").c_str());
  
  leadJetVary->SetFillStyle(0);
  leadJetVary->GetXaxis()->SetTitle("L1 lead jet p_{T} cut [GeV]");
  leadJetVary->GetYaxis()->SetTitle("Efficiency");
  eff = leadJetVary->GetBinContent(leadJetVary->GetXaxis()->FindBin(105));
  eff_string = Form ("Efficiency @ 105 GeV = %.2f", eff );
  leadJetVary->Draw();
  tex1 = new TLatex(0.50,0.92, eff_string.c_str());
  tex1->SetNDC();   
  tex1->SetTextFont(44);   
  tex1->SetTextSize(23);   
  tex1->SetLineWidth(2);   
  tex1->Draw();
  c1.Print((dirname+"/Jet1Pt.pdf").c_str());
  c1.Print((dirname+"/Jet1Pt.png").c_str());
  
  subleadJetVary->SetFillStyle(0);
  subleadJetVary->GetXaxis()->SetTitle("L1 sub-lead jet p_{T} cut [GeV]");
  subleadJetVary->GetYaxis()->SetTitle("Efficiency");
  eff = subleadJetVary->GetBinContent(subleadJetVary->GetXaxis()->FindBin(58));
  eff_string = Form ("Efficiency @ 58 GeV = %.2f", eff );
  subleadJetVary->Draw();
  tex1 = new TLatex(0.50,0.92, eff_string.c_str());
  tex1->SetNDC();   
  tex1->SetTextFont(44);   
  tex1->SetTextSize(23);   
  tex1->SetLineWidth(2);   
  tex1->Draw();
  c1.Print((dirname+"/Jet2Pt.pdf").c_str());
  c1.Print((dirname+"/Jet2Pt.png").c_str());
  
  bool Proceed = true; 
  while(Proceed){
    std::string temp;
    cin >> temp;
    if(temp == "y") Proceed = false;    
  }

  return k;
}

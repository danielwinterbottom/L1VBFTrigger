int RateVsVar(unsigned k=1, unsigned l=0){
  TCanvas c1;
  std::string filename = Form("output/DoubleJet%u.root",k);
  std::string histname;
  std::string outname;
  std::string Xname;
  double xmin;
  double xmax;
  double step;
  if(true){
    histname = "hMuVary";
    outname = "RateOut/Mu.pdf";
    Xname = "L1 muon p_{T} cut [GeV]";
    xmin = 8;
    xmax = 18;
    step = 1;
  }
  if(true){
    histname = "leadJetVary";
    outname = "RateOut/JetPt.pdf";
    Xname = "L1 lead jet p_{T} cut [GeV]";
    xmin = 97;
    xmax = 120;
    step = 1;
  }
  if(false){
    histname = "hsubleadJetVary";
    outname = "RateOut/Jet2Pt.pdf";
    Xname = "L1 sub-lead jet p_{T} cut [GeV]";
    xmin = 30;
    xmax = 60;
    step = 1;
  }
  if(false){
    histname = "hMjjVary";
    outname = "RateOut/Mjj.pdf";
    Xname = "L1 M_{jj} cut [GeV]";
    xmin = 620;
    xmax = 750;
    step = 10;
  }
  
  TFile *f0 = new TFile(filename.c_str());    
  TH2D *h1 = (TH2D*)f0->Get(histname.c_str());
  TH1D *hout = new TH1D("hout", "hout",20,0,2);
  TH1D *hdenum =  new TH1D("hdenum", "hdenum",20,0,2);
  hout->Sumw2();
  TH1D *h_rate = new TH1D("h_rate","h_rate",(int)(xmax-xmin)/step,xmin,xmax);
  h_rate->SetFillStyle(0);
  h_rate->GetXaxis()->SetTitle(Xname.c_str());
  h_rate->GetYaxis()->SetTitle("L1 Rate @ 1.5e34cm^{-2}s^{-1} [kHz]");
  for(unsigned i=(unsigned)xmin; i<=xmax; i+=(unsigned)step){
    cout << i << endl;
    for(unsigned j=1; j<=(unsigned)h1->GetNbinsX(); ++j){
        hout->SetBinContent(j,h1->Integral(j,j,h1->GetYaxis()->FindBin(i),h1->GetNbinsY()+1));
        double error=0;
        h1->IntegralAndError(j,j,h1->GetYaxis()->FindBin(i),h1->GetNbinsY()+1,error);
        hout->SetBinError(j,error);
        hdenum->SetBinContent(j,h1->GetBinContent(j,1));
    }
    hout->Divide(hdenum);
    

    //hout->Draw();
    double rate1pt2=0;
    double rate1pt5=0;
    //hout->GetXaxis()->SetRangeUser(0.4,1.6);
    
    if(l==0){
      TF1  *f1 = new TF1("f1","pol1",0.4,1.6);
      hout->Fit("f1","R");
      rate1pt2 = f1->GetParameter(0) + f1->GetParameter(1)*1.2;
      rate1pt5 = f1->GetParameter(0) + f1->GetParameter(1)*1.5;
    }
    if(l==1){
      TF1  *f1 = new TF1("f1","pol2",0.4,1.6);
      hout->Fit("f1","R");
      rate1pt2 = f1->GetParameter(0) + f1->GetParameter(1)*1.2 + f1->GetParameter(2)*1.2*1.2;
      rate1pt5 = f1->GetParameter(0) + f1->GetParameter(1)*1.5 + f1->GetParameter(2)*1.5*1.5;
    }
    h_rate->SetBinContent(h_rate->GetXaxis()->FindBin(i),rate1pt5);
    //std::cout << "Rate @ 1.2e34: " << rate1pt2 << std::endl;
    //std::cout << "Rate @ 1.5e34: " << rate1pt5 << std::endl;
  }
  h_rate->Draw();
  c1.Print(outname.c_str());
  c1.Print(outname.c_str());
  bool Proceed = true;
  while(Proceed){
    std::string temp;
    cin >> temp;
    if(temp == "y") Proceed = false;    
  }
  
  return k;
}

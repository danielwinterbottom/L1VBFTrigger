{
TCanvas c1;
std::string filename = "output/RateBinnedByLumiMu8.root";
std::string num_name = "lumibinned_num";
std::string denum_name = "lumibinned_denum";
TFile *f0 = new TFile(filename.c_str());    
TH1D *h1 = (TH1D*)f0->Get(num_name.c_str());    
TH1D *h2 = (TH1D*)f0->Get(denum_name.c_str());    
h1->Divide(h2);
h1->Draw();
h1->GetXaxis()->SetTitle("Instantaneous Lumi [cm^{-2}s^{-1}]");
h1->GetYaxis()->SetTitle("L1 Rate [kHz]");
h1->GetXaxis()->SetRangeUser(0.4,1.6);
TF1  *f1 = new TF1("f1","pol1",0.4,1.6);
//f1->FixParameter(0,1);
h1->Fit("f1","R");
//f1->Draw();
std::string outputname = "DoubleJet30_Mj30j30_400_Mu8";
c1.Print(("RateVsLumiPlots/"+outputname+".png").c_str());
c1.Print(("RateVsLumiPlots/"+outputname+".pdf").c_str());
}

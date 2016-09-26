int PlotRateVsLumi(unsigned k=1, unsigned l = 0, unsigned opt = 0){
    
std::string outputname;
if(k == 1) outputname = "DoubleJet_90_30_Mj30j30_580";
else if(k == 2) outputname = "DoubleJet_90_50_Mj30j30_580";
else if(k == 3) outputname = "DoubleJet_90_30_Mj30j30_620";
else if(k == 4) outputname = "DoubleJet_94_30_Mj30j30_620";
else if(k == 5) outputname = "DoubleJet_97_30_Mj30j30_620";
else if(k == 6) outputname = "DoubleJet_97_30_Mj30j30_630";
else if(k == 7) outputname = "DoubleJet30_Mj30j30_360_Mu6";
else if(k == 8) outputname = "DoubleJet30_Mj30j30_360_Mu7";
else if(k == 9) outputname = "DoubleJet30_Mj30j30_360_Mu8";
else if(k == 10) outputname = "DoubleJet30_Mj30j30_400_Mu6";
else if(k == 11) outputname = "DoubleJet30_Mj30j30_400_Mu7";
else if(k == 12) outputname = "DoubleJet30_Mj30j30_400_Mu8";
else if(k == 13) outputname = "DoubleJet30_Mj30j30_360_Mu6er";
else if(k == 14) outputname = "DoubleJet30_Mj30j30_360_Mu7er";
else if(k == 15) outputname = "DoubleJet30_Mj30j30_360_Mu8er";
else if(k == 16) outputname = "DoubleJet30_Mj30j30_400_Mu6er";
else if(k == 17) outputname = "DoubleJet30_Mj30j30_400_Mu7er";
else if(k == 18) outputname = "DoubleJet30_Mj30j30_400_Mu8er";
else if(k == 19) outputname = "DoubleJet30_Mj30j30_360_Mu10";
else if(k == 20) outputname = "DoubleJet30_Mj30j30_360_Mu10er";
else if(k == 21) outputname = "DoubleJet30_Mj30j30_360_Mu9";
else if(k == 22) outputname = "DoubleJet30_Mj30j30_400_Mu9";
else if(k == 23) outputname = "DoubleJet30_Mj30j30_400_Mu10";
else if(k == 29) outputname = "DoubleJet_90_30_Mj30j30_620";
else if(k == 30) outputname = "DoubleJet_97_30_Mj30j30_620";
else if(k == 31) outputname = "DoubleJet_105_30_Mj30j30_620";
else if(k == 32) outputname = "DoubleJet30_Mj30j30_400_Mu8";
else if(k == 33) outputname = "DoubleJet30_Mj30j30_400_Mu9";
else if(k == 34) outputname = "DoubleJet30_Mj30j30_440_Mu9";
else outputname = "output";

TCanvas c1;
std::string filename = Form("output/DoubleJet%u.root",k);
std::string num_name = "lumibinned_num";
std::string denum_name = "lumibinned_denum";
if(opt == 1){
  num_name = "lumibinned_num_pure"; 
  //denum_name = "lumibinned_denum_pure";
  filename = Form("output/PureRate%u.root",k);
}
TFile *f0 = new TFile(filename.c_str());    
TH1D *h1 = (TH1D*)f0->Get(num_name.c_str());    
TH1D *h2 = (TH1D*)f0->Get(denum_name.c_str());    
h1->Divide(h2);
h1->Draw();
h1->GetXaxis()->SetTitle("Instantaneousi Lumi [e34 cm^{-2}s^{-1}]");
h1->GetYaxis()->SetTitle("L1 Rate [kHz]");
h1->GetXaxis()->SetRangeUser(0.4,1.6);

double rate1pt2=0;
double rate1pt5=0;
if(l==0){
  TF1  *f1 = new TF1("f1","pol1",0.4,1.6);
  h1->Fit("f1","R");
  rate1pt2 = f1->GetParameter(0) + f1->GetParameter(1)*1.2;
  rate1pt5 = f1->GetParameter(0) + f1->GetParameter(1)*1.5;
}
if(l==1){
  TF1  *f1 = new TF1("f1","pol2",0.4,1.6);
  h1->Fit("f1","R");
  rate1pt2 = f1->GetParameter(0) + f1->GetParameter(1)*1.2 + f1->GetParameter(2)*1.2*1.2;
  rate1pt5 = f1->GetParameter(0) + f1->GetParameter(1)*1.5 + f1->GetParameter(2)*1.5*1.5;
}

std::string rate1pt2_string = Form ("Rate @ 1.2 e34 cm^{-2}s^{-1} = %.1f kHz", rate1pt2 );
std::string rate1pt5_string = Form ("Rate @ 1.5 e34 cm^{-2}s^{-1} = %.1f kHz", rate1pt5 );
TLatex *   tex1 = new TLatex(0.45,0.2, rate1pt2_string.c_str());
tex1->SetNDC();   
tex1->SetTextFont(44);   
tex1->SetTextSize(23);   
tex1->SetLineWidth(2);   
tex1->Draw();

TLatex *   tex2 = new TLatex(0.45,0.15, rate1pt5_string.c_str());
tex2->SetNDC();   
tex2->SetTextFont(44);   
tex2->SetTextSize(23);   
tex2->SetLineWidth(2);   
tex2->Draw();

TLatex *   tex3 = new TLatex(0.15,0.92, outputname.c_str());
tex3->SetNDC();   
tex3->SetTextFont(44);   
tex3->SetTextSize(23);   
tex3->SetLineWidth(2);   
tex3->Draw();

c1.Print(("RateVsLumiPlots/"+outputname+".png").c_str());
c1.Print(("RateVsLumiPlots/"+outputname+".pdf").c_str());

return k;
}

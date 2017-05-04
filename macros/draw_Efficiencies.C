#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

Double_t fitfunc(Double_t *v, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = v[0]/par[2];
  Double_t fitval;
  if (arg>1)
    fitval = par[0]*(1.- (1./(par[1]*log(arg)+1))) + par[3]*(1-par[0]*(1.- (1./(par[1]*log(arg)+1))));
  else
    fitval = 0.;
  
  return fitval;
}

void draw_Efficiencies() {

  // inputs
  TFile *inputBinps = new TFile("myOutputFileBinps.root");  
  TFile *inputCamer = new TFile("myOutputFileCameroniBtf.root");  
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *binp1 = (TGraphAsymmErrors*)inputBinps->Get("binp1");
  TGraphAsymmErrors *binp3 = (TGraphAsymmErrors*)inputBinps->Get("binp3");
  TGraphAsymmErrors *binp4 = (TGraphAsymmErrors*)inputBinps->Get("binp4");
  
  TGraphAsymmErrors *effMib25Corr_1200 = (TGraphAsymmErrors*)inputCamer->Get("effMib25Corr_1200");
  TGraphAsymmErrors *effRm8Corr_1200   = (TGraphAsymmErrors*)inputCamer->Get("effRm8Corr_1200");
  TGraphAsymmErrors *effRm5Corr_1200   = (TGraphAsymmErrors*)inputCamer->Get("effRm5Corr_1200");

  TGraphAsymmErrors *effZS2Corr  = (TGraphAsymmErrors*)inputH4->Get("effZS2Corr");
  TGraphAsymmErrors *effMib3Corr = (TGraphAsymmErrors*)inputH4->Get("effMib3Corr");

  TF1* f1 = new TF1("f1",fitfunc,0.,2700.,4);
  f1->SetParameters(1.,1.,300.,0.);  
  f1->SetParLimits(0,0.,1.);
  f1->SetParLimits(1,0.,1.);
  f1->SetParLimits(2,0.,1000.);
  f1->SetParLimits(3,0.,1.);

  //binp3->Fit("f1","B");

  TF1* f2 = new TF1("f2",fitfunc,0.,2700.,4);
  f2->SetParameters(0.5,0.5,1400.,0.5,0.5);
  f2->SetParLimits(0,0.,1.);
  f2->SetParLimits(1,0.,1.);
  f2->SetParLimits(2,1300.,1500.);
  f2->SetParLimits(3,0.,1.);
  f2->SetParLimits(4,0.,1.);
  effMib25Corr_1200->Fit("f2","B");

  TF1* f3 = new TF1("f3",fitfunc,0.,2700.,4);

  // cosmetics
  binp1->SetMarkerStyle(20);
  binp1->SetMarkerColor(1);
  binp1->SetLineColor(1);
  binp1->SetLineStyle(1);
  binp1->SetLineWidth(2);
  binp1->SetMarkerSize(1.1);
  binp3->SetMarkerStyle(21);
  binp3->SetMarkerColor(kRed+2);
  binp3->SetLineColor(kRed+2);
  binp3->SetLineStyle(2);
  binp3->SetLineWidth(2);
  binp3->SetMarkerStyle(21);
  binp4->SetMarkerStyle(22);
  binp4->SetMarkerColor(kGreen+3);
  binp4->SetLineColor(kGreen+3);
  binp4->SetLineStyle(3);
  binp4->SetLineWidth(2);
  binp4->SetMarkerSize(1.4);

  effMib25Corr_1200->SetMarkerStyle(23);
  effMib25Corr_1200->SetMarkerColor(kRed+3);
  effMib25Corr_1200->SetLineColor(kRed+3);
  effMib25Corr_1200->SetLineStyle(4);
  effMib25Corr_1200->SetLineWidth(2);
  effMib25Corr_1200->SetMarkerSize(1.4);
  effRm8Corr_1200->SetMarkerStyle(20);
  effRm8Corr_1200->SetMarkerColor(5);
  effRm8Corr_1200->SetLineColor(5);
  effRm8Corr_1200->SetLineStyle(1);
  effRm8Corr_1200->SetLineWidth(2);
  effRm8Corr_1200->SetMarkerSize(1);

  effRm5Corr_1200->SetMarkerStyle(33);
  effRm5Corr_1200->SetMarkerColor(kCyan+2);
  effRm5Corr_1200->SetLineColor(kCyan+2);
  effRm5Corr_1200->SetLineStyle(5);
  effRm5Corr_1200->SetLineWidth(2);
  effRm5Corr_1200->SetMarkerSize(1.7);

  effZS2Corr->SetMarkerStyle(24);
  effZS2Corr->SetMarkerColor(kMagenta+3);
  effZS2Corr->SetLineColor(kMagenta+3);
  effZS2Corr->SetLineStyle(6);
  effZS2Corr->SetLineWidth(2);
  effZS2Corr->SetMarkerSize(1.1);
  effMib3Corr->SetMarkerStyle(25);
  effMib3Corr->SetMarkerColor(kBlue-1);
  effMib3Corr->SetLineColor(kBlue-1);
  effMib3Corr->SetLineStyle(7);
  effMib3Corr->SetLineWidth(2);

  // plotting
  setStyle();

  //TH2F* H2 = new TH2F("H2","",3200,0.,3200.,100,0.,1.4);
  TH2F* H2 = new TH2F("H2","",2700,0.,2700.,100,0.,1.4);
  H2->GetXaxis()->SetTitle("MCP-stack bias (V)");
  H2->GetYaxis()->SetTitle("Efficiency");
  TLine *line = new TLine(0,1,3200,1);

  TLegend *leg;
  leg = new TLegend(0.18,0.70,0.38,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(binp1, "40x2+40x2, Amp bias=2100V", "pl");
  leg->AddEntry(binp3, "90x1+40x2, Amp bias=2100V", "pl");
  leg->AddEntry(binp4, "90x2+40x2, Amp bias=2100V", "pl");
  //leg->AddEntry(effZS2Corr,  "40x3, electrons at 50 GeV", "pl");

  TLegend *legA;
  legA = new TLegend(0.58,0.70,0.78,0.90);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.03);
  legA->SetFillColor(0);
  legA->AddEntry(effMib25Corr_1200, "80x2+40x1, Amp bias=1200V", "pl");
  //legA->AddEntry(effRm8Corr_1200,   "80x1+40x1, Amp bias=1200V", "pl");
  legA->AddEntry(effRm5Corr_1200,   "80x1+40x1, Amp bias=1200V", "pl");
  //legA->AddEntry(effMib3Corr, "40x2, electrons at 50 GeV", "pl");

  TLegend *leg2;
  leg2 = new TLegend(0.18,0.57,0.38,0.65);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetFillColor(0);
  leg2->AddEntry(binp1,      "electrons at 491 MeV", "p");
  leg2->AddEntry(effZS2Corr, "electrons at 50 GeV",  "p");

  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  binp1->Draw("PCsame");
  binp3->Draw("PCsame");
  binp4->Draw("PCsame");
  effMib25Corr_1200->Draw("PCsame");
  //effRm8Corr_1200->Draw("PCsame");
  effRm5Corr_1200->Draw("PCsame");
  //effZS2Corr->Draw("PCsame");
  //effMib3Corr->Draw("PCsame");
  leg->Draw();
  legA->Draw();
  //leg2->Draw();
  TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw(); 
  c1->SaveAs("summaryEff_v3.png");
  c1->SaveAs("summaryEff_v3.pdf");


}

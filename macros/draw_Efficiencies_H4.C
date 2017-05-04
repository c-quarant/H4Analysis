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
  if (par[3] != 0) arg = v[0]/par[3];
  Double_t fitval;
  if (arg>1)
    fitval = par[0] + (par[1]*(1.- (1./(par[2]*log(arg)+1))));
  else
    fitval = 0.;
  
  return fitval;
}

Double_t fitfunc2(Double_t *v, Double_t *par) {

  Double_t fitval;
  fitval = 1./(1. + par[0]*exp(par[1]*(v[0]+par[2])));
 
  return fitval;
}

void draw_Efficiencies_H4() {

  gStyle->SetOptTitle(0); 
  //gStyle->SetOptStat(1110); 
  gStyle->SetOptStat(0000); 
  //gStyle->SetOptFit(1); 
  gStyle->SetOptFit(0); 
  //gStyle->SetErrorX(0); 

  // inputs
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *effZS2Corr  = (TGraphAsymmErrors*)inputH4->Get("effZS2Corr");
  TGraphAsymmErrors *effMib3Corr = (TGraphAsymmErrors*)inputH4->Get("effMib3Corr");
  TGraphAsymmErrors *effZS2PMTCorr = (TGraphAsymmErrors*)inputH4->Get("effZS2PMTCorr");

  TF1* f1 = new TF1("f1",fitfunc,1300.,3200.,4);
  f1->SetLineColor(kGreen+2);
  f1->SetLineStyle(7);
  f1->SetLineWidth(2);

  TF1* f2 = new TF1("f2",fitfunc,1300.,3200.,4);
  f2->SetLineColor(kRed+1);
  f2->SetLineStyle(8);
  f2->SetLineWidth(2);

  TF1* f4 = new TF1("f4",fitfunc2,1300.,3200.,3);
  f4->SetLineColor(kGray+1);
  f4->SetLineStyle(6);
  f4->SetLineWidth(2);

  f1->SetParameters(0.2,1.,1.,2250.);  
  f1->SetParLimits(1,0.,1.);
  f1->SetParLimits(3,1800.,2250.);
  effZS2Corr->Fit("f1","B");

  f2->SetParameters(0.2,1.,1.,2000.);  
  f2->SetParLimits(1,0.,1.);
  f2->SetParLimits(3,1800.,2050.);
  effMib3Corr->Fit("f2","B");

  f4->SetParameters(1.,1.,1.,1.); 
  f4->SetParLimits(0,0.,1.); 
  //f4->SetParLimits(1,0.,1.); 
  f4->SetParLimits(2,-3500,-1000.); 
  effZS2PMTCorr->Fit("f4");

  // cosmetics
  effZS2Corr->SetMarkerStyle(20);
  effZS2Corr->SetMarkerColor(kMagenta+3);
  effZS2Corr->SetLineColor(kMagenta+3);
  effZS2Corr->SetLineColor(kGreen+2);
  effZS2Corr->SetLineStyle(7);
  effZS2Corr->SetLineWidth(2);
  effZS2Corr->SetMarkerSize(1.1);

  effMib3Corr->SetMarkerStyle(21);
  effMib3Corr->SetMarkerColor(kBlue-1);
  effMib3Corr->SetLineColor(kBlue-1);
  effMib3Corr->SetLineColor(kRed+1);
  effMib3Corr->SetLineStyle(8);
  effMib3Corr->SetLineWidth(2);
  effMib3Corr->SetMarkerSize(1.1);

  effZS2PMTCorr->SetMarkerStyle(23);
  effZS2PMTCorr->SetMarkerColor(kOrange+1);
  effZS2PMTCorr->SetLineColor(kOrange+1);  
  effZS2PMTCorr->SetLineColor(kGray+1);
  effZS2PMTCorr->SetLineStyle(6);
  effZS2PMTCorr->SetLineWidth(2);
  effZS2PMTCorr->SetMarkerSize(1.3);

  // plotting
  setStyle();

  TH2F* H2 = new TH2F("H2","",2000,1300.,3300.,100,0.,1.05);
  H2->GetXaxis()->SetTitle("MCP-stack bias (V)");
  H2->GetYaxis()->SetTitle("Efficiency");
  H2->GetXaxis()->SetNdivisions(505);
  TLine *line = new TLine(1300,1,3300,1);

  TLegend *leg;
  leg = new TLegend(0.18,0.65,0.55,0.87);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->AddEntry(effZS2PMTCorr,  "40x3 PMT-MCP", "pl");
  leg->AddEntry(effZS2Corr,  "40x3 iMCP", "pl");
  leg->AddEntry(effMib3Corr,  "40x2 iMCP", "pl");


  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  effZS2Corr->Draw("Psame");
  effMib3Corr->Draw("Psame");
  effZS2PMTCorr->Draw("Psame");
  leg->Draw("same");
  TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 50 GeV}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw(); 
  c1->SaveAs("summaryEff_H4.png");
  c1->SaveAs("summaryEff_H4.pdf");
}


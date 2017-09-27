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

void draw_Efficiencies_H4Abs() {

  gStyle->SetOptTitle(0); 
  //gStyle->SetOptStat(1110); 
  gStyle->SetOptStat(0000); 
  //gStyle->SetOptFit(1); 
  gStyle->SetOptFit(0); 
  //gStyle->SetErrorX(0); 

  // inputs
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *effAbsZS2G  = (TGraphAsymmErrors*)inputH4->Get("effAbsZS2G");
  TGraphAsymmErrors *effAbsZS1G = (TGraphAsymmErrors*)inputH4->Get("effAbsZS1G");
  TGraphAsymmErrors *effAbsMiB3G = (TGraphAsymmErrors*)inputH4->Get("effAbsMiB3G");

  for(int ii = 0; ii < effAbsZS2G->GetN();ii++)
  {
        double x,y;
        double y_errorUp,y_errorDown;
        effAbsZS2G->GetPoint(ii,x,y);
        if(x != 0) continue;
        effAbsZS2G->SetPoint(ii,x,0.661674);
  }

  for(int ii = 0; ii < effAbsMiB3G->GetN();ii++)
  {
        double x,y;
        double y_errorUp,y_errorDown;
        effAbsMiB3G->GetPoint(ii,x,y);
        if(x != 0) continue;
        effAbsMiB3G->SetPoint(ii,x,0.5);
  }
  
  // cosmetics
  effAbsZS2G->SetMarkerStyle(24);
  effAbsZS2G->SetMarkerColor(kMagenta+3);
  effAbsZS2G->SetLineColor(kMagenta+3);
  effAbsZS2G->SetLineStyle(6);
  effAbsZS2G->SetLineWidth(2);
  effAbsZS2G->SetMarkerSize(1.1);
  effAbsZS1G->SetMarkerStyle(26);
  effAbsZS1G->SetMarkerColor(kViolet+2);
  effAbsZS1G->SetLineColor(kViolet+2);
  effAbsZS1G->SetLineStyle(7);
  effAbsZS1G->SetLineWidth(2);
  effAbsZS1G->SetMarkerSize(1.1);
  effAbsMiB3G->SetMarkerStyle(25);
  effAbsMiB3G->SetMarkerColor(kBlue-1);
  effAbsMiB3G->SetLineColor(kBlue-1);
  effAbsMiB3G->SetLineStyle(8);
  effAbsMiB3G->SetLineWidth(2);
  effAbsMiB3G->SetMarkerSize(1.1);

  // plotting
  setStyle();

  TH2F* H2 = new TH2F("H2","",55,0.,5.5,100,0.,1.05);
  H2->GetXaxis()->SetTitle("Radiation lengths (X_{0})");
  H2->GetYaxis()->SetTitle("Efficiency");
  //H2->GetXaxis()->SetNdivisions(505);
  TLine *line = new TLine(0.,1,5.5,1);

  TLegend *leg;
  leg = new TLegend(0.68,0.55,0.85,0.67);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->AddEntry(effAbsZS2G,  "40x3", "pl");
  //leg->AddEntry(effAbsZS1G,  "40x3 PMT-MCP", "pl");
  leg->AddEntry(effAbsMiB3G,  "40x2", "pl");

  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  effAbsZS2G->Draw("PLsame");
  //effAbsZS1G->Draw("PCsame");
  effAbsMiB3G->Draw("PLsame");
  leg->Draw("same");
  TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 20 GeV}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw(); 
  c1->SaveAs("effVsXo_v2.png");
  c1->SaveAs("effVsXo_v2.pdf");
}


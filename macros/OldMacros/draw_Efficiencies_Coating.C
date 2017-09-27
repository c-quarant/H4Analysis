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

Double_t fitfunc(Double_t *v, Double_t *par);
Double_t fitfunc2(Double_t *v, Double_t *par);
std::pair<float,std::pair<float,float> > computeStat( std::vector<double>, std::vector<double>, std::vector<double>); 

void draw_Efficiencies_Coating() {

  gStyle->SetOptTitle(0); 
  //gStyle->SetOptStat(1110); 
  gStyle->SetOptStat(0000); 
  //gStyle->SetOptFit(1); 
  gStyle->SetOptFit(0); 
  //gStyle->SetErrorX(0); 

  // inputs
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *effSee = (TGraphAsymmErrors*)inputH4->Get("effSee");
  TGraphAsymmErrors *effMgO_PadA = (TGraphAsymmErrors*)inputH4->Get("effMgO_PadA");
  TGraphAsymmErrors *effMgO_PadB = (TGraphAsymmErrors*)inputH4->Get("effMgO_PadB");
  TGraphAsymmErrors *effMgO_PadC = (TGraphAsymmErrors*)inputH4->Get("effMgO_PadC");
  TGraphAsymmErrors *effMgO_PadD = (TGraphAsymmErrors*)inputH4->Get("effMgO_PadD");

  TGraphAsymmErrors *effMgO = new TGraphAsymmErrors();
  
  std::vector<std::vector<double> > mean;
  mean.resize(effMgO_PadA->GetN());
  std::vector<std::vector<double> > error_up;
  error_up.resize(effMgO_PadA->GetN());
  std::vector<std::vector<double> > error_down;
  error_down.resize(effMgO_PadA->GetN());

  for(int ii = 0; ii<effMgO_PadA->GetN(); ii++)
  {
      double x,y;
      double errorUp, errorDown;

      effMgO_PadA->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadA->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadA->GetErrorYlow(ii));

      effMgO_PadB->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadB->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadB->GetErrorYlow(ii));

      effMgO_PadC->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadC->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadC->GetErrorYlow(ii));

      effMgO_PadD->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadD->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadD->GetErrorYlow(ii));
  }

  for(int ii = 0; ii<effMgO_PadA->GetN(); ii++)
  {
      float Mean = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).first;
      float Error_up = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).second.first;
      float Error_down = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).second.second;
        
      double x,y;
      effMgO_PadA->GetPoint(ii,x,y);
      effMgO->SetPoint(ii,x,Mean);
  }
 
  TF1* f1 = new TF1("f1",fitfunc,1000.,3200.,4);
  f1->SetLineColor(kRed+1);
  f1->SetLineStyle(7);
  f1->SetLineWidth(2);

  TF1* f2 = new TF1("f2",fitfunc,1000.,3200.,4);
  f2->SetLineColor(kGreen+1);
  f2->SetLineStyle(8);
  f2->SetLineWidth(2);

  f1->SetParameters(0.2,0.1,1.,1450.);  
  f1->SetParLimits(1,0.,1.);
  f1->SetParLimits(3,1400.,1450.);
  effSee->Fit("f1","B,0");

  TF1* f1_post = new TF1("f2_post",fitfunc,1100.,2500.,4);
  f1_post->FixParameter(0,f1->GetParameter(0));
  f1_post->FixParameter(1,f1->GetParameter(1));
  f1_post->FixParameter(2,f1->GetParameter(2));
  f1_post->FixParameter(3,f1->GetParameter(3));
  f1_post->SetLineColor(kRed+1);
  f1_post->SetLineStyle(7);
  f1_post->SetLineWidth(2);

  f2->SetParameters(0.2,0.1,1.,2100.);  
  f2->SetParLimits(1,0.,1.);
  //f2->SetParLimits(2,0.,1.);
  f2->SetParLimits(3,2000.,2100.);
  effMgO->Fit("f2","B,0");

  TF1* f2_post = new TF1("f2_post",fitfunc,1900.,3000.,4);
  f2_post->FixParameter(0,f2->GetParameter(0));
  f2_post->FixParameter(1,f2->GetParameter(1));
  f2_post->FixParameter(2,f2->GetParameter(2));
  f2_post->FixParameter(3,f2->GetParameter(3));
  f2_post->SetLineColor(kGreen+1);
  f2_post->SetLineStyle(8);
  f2_post->SetLineWidth(2);

  // cosmetics
  effSee->SetMarkerStyle(20);
  effSee->SetMarkerColor(kCyan+2);
  effSee->SetLineColor(kCyan+2);  
  effSee->SetLineColor(kRed+1);
  effSee->SetLineStyle(7);
  effSee->SetLineWidth(2);
  effSee->SetMarkerSize(1.3);

  effMgO->SetMarkerStyle(21);
  effMgO->SetMarkerColor(kOrange+1);
  effMgO->SetLineColor(kOrange+1);  
  effMgO->SetLineColor(kGreen+1);
  effMgO->SetLineStyle(8);
  effMgO->SetLineWidth(2);
  effMgO->SetMarkerSize(1.1);

  // plotting
  setStyle();

  TH2F* H2 = new TH2F("H2","",2500,800.,3300.,100,0.,1.05);
  H2->GetXaxis()->SetTitle("MCP-stack bias (V)");
  H2->GetYaxis()->SetTitle("Efficiency");
  TLine *line = new TLine(800,1,3300,1);

  TLegend *leg;
  leg = new TLegend(0.22,0.75,0.55,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(effSee,  "40x2 SEE-iMCP, electrons at 50 GeV", "pl");
  leg->AddEntry(effMgO,  "40x2 MgO-iMCP, muons at 150 GeV", "pl");


  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  effSee->Draw("Psame");
  effMgO->Draw("Psame");
  f1_post->Draw("Lsame");
  f2_post->Draw("Lsame");
  leg->Draw("same");
  c1->SaveAs("summaryEff_Coating.png");
  c1->SaveAs("summaryEff_Coating.pdf");
}

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
  fitval = 1./(par[0] + par[1]*exp(par[2]*v[0]));
 
  return fitval;
}

std::pair<float,std::pair<float,float> > computeStat( std::vector<double> mean, std::vector<double> error_up, std::vector<double> error_down)
{
  float Mean = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Mean = Mean +mean.at(ii);

  Mean = Mean/mean.size();
  
  float Error_up = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Error_up = Error_up + (mean.at(ii)-Mean)*(mean.at(ii)-Mean);

  Error_up = Error_up/(mean.size()-1);
  Error_up = sqrt(Error_up);

  float Error_down = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Error_down = Error_down + (mean.at(ii)-Mean)*(mean.at(ii)-Mean);

  Error_down = Error_down/(mean.size()-1);
  Error_down = sqrt(Error_down);

  std::pair<float,std::pair<float,float> > outpair = std::make_pair(Mean,std::make_pair(Error_up,Error_down));
  return outpair ;
}


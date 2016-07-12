#include "TFile.h" 
#include "TTree.h" 
#include "TH1F.h" 
#include "TGraphAsymmErrors.h"
#include "TCanvas.h" 
#include "TLegend.h" 
#include "TROOT.h"
#include "TStyle.h"

#include <iostream>
#include <fstream> 

void SetEx(TGraphAsymmErrors* gae, Double_t Ex);

void DrawEfficiency(std::string inputs)
{

    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());

    TGraphAsymmErrors* eff_MiB_10m_HVMID_1200 = (TGraphAsymmErrors*)inputFile->Get("eff_MiB_10m_4567_4575_HVMID_1200"); 
    TGraphAsymmErrors* eff_MiB_10m_HVMID_1000 = (TGraphAsymmErrors*)inputFile->Get("eff_MiB_10m_4560_4564_and_4576_4582_HVMID_1000");
    TGraphAsymmErrors* eff_MiB_25m_HVMID_1000 = (TGraphAsymmErrors*)inputFile->Get("eff_MiB_25m_4521_4531_HVMID_1000");
    TGraphAsymmErrors* eff_MiB_25m_HVMID_1300 = (TGraphAsymmErrors*)inputFile->Get("eff_MiB_25m_4537_4553_HVMID_1300");

    SetEx(eff_MiB_10m_HVMID_1200,0.);
    SetEx(eff_MiB_10m_HVMID_1000,0.);
    SetEx(eff_MiB_25m_HVMID_1000,0.);
    SetEx(eff_MiB_25m_HVMID_1300,0.);
    
    eff_MiB_25m_HVMID_1300->GetXaxis()->SetTitle("HV (V)");
    eff_MiB_25m_HVMID_1300->GetYaxis()->SetTitle("#epsilon");

    eff_MiB_25m_HVMID_1300->GetYaxis()->SetRangeUser(0.,1.);

    eff_MiB_25m_HVMID_1000->SetMarkerStyle(20);
    eff_MiB_25m_HVMID_1000->SetMarkerSize(0.9);
    eff_MiB_25m_HVMID_1000->SetMarkerColor(kRed+1);
    eff_MiB_25m_HVMID_1300->SetMarkerStyle(24);
    eff_MiB_25m_HVMID_1300->SetMarkerSize(0.9);
    eff_MiB_25m_HVMID_1300->SetMarkerColor(kRed+1);

    eff_MiB_10m_HVMID_1000->SetMarkerStyle(20);
    eff_MiB_10m_HVMID_1000->SetMarkerSize(0.9);
    eff_MiB_10m_HVMID_1000->SetMarkerColor(kGreen+1);
    eff_MiB_10m_HVMID_1200->SetMarkerStyle(24);
    eff_MiB_10m_HVMID_1200->SetMarkerSize(0.9);
    eff_MiB_10m_HVMID_1200->SetMarkerColor(kGreen+1);

    TLegend* legend = new TLegend(0.47, 0.22, 0.94, 0.44);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(eff_MiB_25m_HVMID_1000,"MiB 25 #mum, HVMID = 1 kV","P");
    legend -> AddEntry(eff_MiB_25m_HVMID_1300,"MiB 25 #mum, HVMID = 1.3 kV","P");
    legend -> AddEntry(eff_MiB_10m_HVMID_1000,"MiB 10 #mum, HVMID = 1 kV","P");
    legend -> AddEntry(eff_MiB_10m_HVMID_1200,"MiB 10 #mum, HVMID = 1.2 kV","P");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetGrid();
    eff_MiB_25m_HVMID_1300->Draw("AP");
    eff_MiB_25m_HVMID_1000->Draw("P,same");
    eff_MiB_10m_HVMID_1200->Draw("P,same");
    eff_MiB_10m_HVMID_1000->Draw("P,same");
    legend -> Draw("same");

    c1 -> Print("T9_testbeam_Summer2016_MiB_Efficiency.png","png");
    c1 -> Print("T9_testbeam_Summer2016_MiB_Efficiency.pdf","pdf");
    
}

void SetEx(TGraphAsymmErrors* gae, Double_t Ex)
{
   Int_t np = gae->GetN();
   for (Int_t i=0; i<np; i++) {
      gae->SetPointEXhigh(i,Ex);
      gae->SetPointEXlow(i,Ex);
   }
}
	

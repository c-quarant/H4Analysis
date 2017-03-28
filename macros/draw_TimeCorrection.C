#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

#include<iostream>
#include<string>
#include<fstream>

void draw_TimeCorrection()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
    gStyle->SetErrorX(0);

    
    TFile* inputFile = TFile::Open("timingCorrection_CFD50_thres20.root");
     
    TH1F* timingCorrection_wrtMiB2 = (TH1F*)inputFile->Get("timingCorrection_wrtMiB2");
    TH1F* timingCorrection_wrtRm2 = (TH1F*)inputFile->Get("timingCorrection_wrtRm2");
    TF1* fit_corr_wrtMiB2 = (TF1*)inputFile->Get("timingCorrection_wrtMiB2_f");
    TF1* fit_corr_wrtRm2 = (TF1*)inputFile->Get("timingCorrection_wrtRm2_f");

    timingCorrection_wrtMiB2->SetMarkerColor(kBlue+1); 
    timingCorrection_wrtMiB2->SetLineColor(kBlue+1);
    timingCorrection_wrtMiB2->SetMarkerStyle(20);
    timingCorrection_wrtMiB2->SetMarkerSize(1.5);
    fit_corr_wrtMiB2->SetLineColor(kRed+1);

    timingCorrection_wrtRm2->SetMarkerColor(kGreen+1); 
    timingCorrection_wrtRm2->SetLineColor(kGreen+1);
    timingCorrection_wrtRm2->SetMarkerStyle(21);
    timingCorrection_wrtRm2->SetMarkerSize(1.5);
    fit_corr_wrtRm2->SetLineColor(kViolet+1);

    //std::cout << timingCorrection_wrtMiB2->GetNbinsX() << std::endl;
    float ampBinning[9] = {0.};
    for(int ii = 1; ii<= timingCorrection_wrtMiB2->GetNbinsX(); ii++)
        ampBinning[ii] = timingCorrection_wrtMiB2->GetBinCenter(ii)+timingCorrection_wrtMiB2->GetBinWidth(ii)/2.;

    float corrBinning[300] = {0.};
    for(int ii = 0; ii<= 300; ii++)
        corrBinning[ii] = -0.2+ii*0.00075;
 
    setStyle();

    TH2F* H2 = new TH2F("H2","",8,ampBinning,300,corrBinning);
    H2->GetXaxis()->SetTitle("amplitude (ADC counts)");
    H2->GetYaxis()->SetTitle("#Deltat (ns)");

    TLegend* legend = new TLegend(0.52, 0.52, 0.69, 0.64);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.05);
    
    legend -> AddEntry(timingCorrection_wrtMiB2,"i-MCP","P");
    legend -> AddEntry(timingCorrection_wrtRm2,"PMT-MCP","P");

    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    H2->Draw();
    //c1->cd();
    //c1->SetLogx();
    timingCorrection_wrtMiB2->Draw("P,same");
    timingCorrection_wrtRm2->Draw("P.same");
    fit_corr_wrtMiB2->Draw("L,same");
    fit_corr_wrtRm2->Draw("L,same");
    legend->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw();
    c1 -> Print(std::string("timing_correction_Final.png").c_str(),"png");
    c1 -> Print(std::string("timing_correction_Final.pdf").c_str(),"pdf");
}

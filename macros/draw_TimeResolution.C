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

void draw_TimeResolution()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
    gStyle->SetErrorX(0);

    
    TFile* inputFile = TFile::Open("Data_TimeResolution_vs_amp_BINP3_CFD50_thres20_onlyWrtMiB2_Inclusive.root");
     
    TH1F* time_wrtMiB2 = (TH1F*)inputFile->Get("time_wrtMiB2");  
    //TH1F* time_wrtMiB2 = (TH1F*)inputs->Get("time_wrtRm2");  
    time_wrtMiB2->SetLineColor(kCyan+1);
    time_wrtMiB2->SetFillColor(kCyan+1);

    TF1* g_res = (TF1*)inputFile->Get("g_res");
    g_res->SetLineColor(kRed);
    g_res->SetNpx(10000);

    setStyle();

    TH2F* H2 = new TH2F("H2","",200,-0.5,0.5,700,0.,700.);
    H2->GetXaxis()->SetTitle("#Deltat (ns)");
    H2->GetYaxis()->SetTitle("events/10 ps");

    TLegend* legend = new TLegend(0.63, 0.75, 0.70, 0.82);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.03);
    
    legend -> AddEntry(time_wrtMiB2,"90x1+40x2: #sigma_{t} = 25 #pm 3 ps","F");
    
    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.04);
    latexLabel->SetTextColor(kBlack);
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica

    TLatex *latexLabel2 = new TLatex();
    latexLabel2->SetTextSize(0.04);
    latexLabel2->SetTextColor(kBlack);
    latexLabel2->SetNDC();
    latexLabel2->SetTextFont(42); // helvetica
    
    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    H2->Draw();
    //c1->cd();
    time_wrtMiB2->Draw("H,same");
    //g_res->Draw("same");
    //legend ->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw(); 
    c1 -> Print("TimeResolution_inclusive.png","png");
    c1 -> Print("TimeResolution_inclusive.pdf","pdf");   

}

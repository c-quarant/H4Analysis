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

void draw_TimeResolution_H4_Abs()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    ifstream infile;    

    TGraphAsymmErrors* g1 = new TGraphAsymmErrors();

    bool isF1_open = false;

    infile.open("res_vs_HV_ZS2_abs_final_v2.dat");    
    if(!infile.fail())
    {
       int raw = 0;
       isF1_open = true;
       while(!infile.eof())
       {
          if(infile.eof())
             break;
          float X0, res, resError;
          infile >> X0 >> res >> resError;
          std::cout << X0 << " " << res << " " << resError << " " << raw << std::endl;
          g1->SetPoint(raw,X0,res);
          g1->SetPointEYlow(raw,resError); 
          g1->SetPointEYhigh(raw,resError); 
          raw++;
       }    
    }
    infile.close(); 

    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.8);
    g1->SetMarkerColor(kBlue+1);
    g1->SetLineColor(kBlue+1);
    g1->SetLineStyle(4);
    g1->SetLineWidth(2);

    setStyle(); 

    TH2F* H2 = new TH2F("H2","",6,0.,6.,100,10.,70.);
    H2->GetXaxis()->SetTitle("Radiation length (X_{0})");
    H2->GetYaxis()->SetTitle("#sigma_{t} (ps)");

    TLegend* legend = new TLegend(0.62, 0.72, 0.79, 0.80);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    if(isF1_open) legend -> AddEntry(g1,"40x3","PEL");

    TLine *line = new TLine(0,20,6,20);
    line->SetLineColor(kBlack);
    line->SetLineStyle(4);
    line->SetLineWidth(2);
    
    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    H2->Draw();
    g1->Draw("PC,same");
    legend->Draw("same");
    line->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 50 GeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw(); 
    c1 -> Print("Resolution_AbsScan_H4.png","png");
    c1 -> Print("Resolution_AbsScan_H4.pdf","pdf");
}


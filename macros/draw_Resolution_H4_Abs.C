#include "TFile.h" 
#include "TTree.h" 
#include "TH1F.h" 
#include "TH2F.h" 
#include "TF1.h" 
#include "TGraphAsymmErrors.h"
#include "TCanvas.h" 
#include "TLegend.h" 
#include "TROOT.h"
#include "TStyle.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TLine.h"

#include <iostream>
#include <fstream> 

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
    g1->GetXaxis()->SetRangeUser(0.,6.);
    g1->GetYaxis()->SetRangeUser(0,100);
    g1->GetXaxis()->SetTitle("Radiation length (X_{0})");
    g1->GetYaxis()->SetTitle("#sigma_{t} (ps)");
    g1->GetXaxis()->SetTitleSize(0.045);
    g1->GetYaxis()->SetTitleSize(0.045);;

    TLegend* legend = new TLegend(0.62, 0.72, 0.79, 0.84);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    if(isF1_open) legend -> AddEntry(g1,"i-MCP, 3 layers","PEL");

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

    TLine *line = new TLine(0,20,5.55,20);
    line->SetLineColor(kBlack);
    line->SetLineStyle(4);
    line->SetLineWidth(2);
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    g1->Draw("APC");
    legend->Draw("same");
    latexLabel->DrawLatex(0.10, 0.93,"Preliminary");
    latexLabel2->DrawLatex(0.65, 0.93,"Electrons at 20 GeV");
    line->Draw("same");
    c1 -> Print(std::string("Resolution_AbsScan_H4.png").c_str(),"png");
    c1 -> Print(std::string("Resolution_AbsScan_H4.pdf").c_str(),"pdf");
}


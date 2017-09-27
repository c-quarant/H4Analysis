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

void ReverseXAxis (TGraph *g);
void ReverseXGraph (TGraph *g);

void draw_TimeResolution_vs_fraction()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
    gStyle->SetErrorX(0);

    
    TFile* inputFile = TFile::Open("Final_TimeResolution_vs_eff_BINP3_CFD50_thres20_onlyWrtMiB2.root");
     
    TGraphAsymmErrors* TimeResolution_vs_frac_wrtMiB2 = (TGraphAsymmErrors*)inputFile->Get("TimeResolution_vs_eff_wrtMiB2_BINP3_CFD50_thres20_onlyWrtMiB2");   

    TimeResolution_vs_frac_wrtMiB2->SetMarkerStyle(20);
    TimeResolution_vs_frac_wrtMiB2->SetMarkerSize(0.9);
    TimeResolution_vs_frac_wrtMiB2->SetMarkerColor(kBlack);
    TimeResolution_vs_frac_wrtMiB2->SetLineColor(kBlack);
    TimeResolution_vs_frac_wrtMiB2->SetLineWidth(1);

    for(int ii = 0; ii < TimeResolution_vs_frac_wrtMiB2->GetN();ii++)
    {
        double x,y;
        double x_errorUp,y_errorUp, x_errorDown,y_errorDown;
        TimeResolution_vs_frac_wrtMiB2->GetPoint(ii,x,y);
        x_errorDown = TimeResolution_vs_frac_wrtMiB2->GetErrorXlow(ii);
        x_errorUp = TimeResolution_vs_frac_wrtMiB2->GetErrorXhigh(ii);
        y_errorDown = TimeResolution_vs_frac_wrtMiB2->GetErrorYlow(ii);
        y_errorUp = TimeResolution_vs_frac_wrtMiB2->GetErrorYhigh(ii); 

        TimeResolution_vs_frac_wrtMiB2->SetPoint(ii,x,y);
        TimeResolution_vs_frac_wrtMiB2->SetPointEXlow(ii,x_errorDown); 
        TimeResolution_vs_frac_wrtMiB2->SetPointEXhigh(ii,x_errorUp); 
        TimeResolution_vs_frac_wrtMiB2->SetPointEYlow(ii,y_errorDown); 
        TimeResolution_vs_frac_wrtMiB2->SetPointEYhigh(ii,y_errorUp); 
    }

    bool useBestResolution = true;
    bool isF1_open = false;
    if(useBestResolution){
       ifstream infile;    
       infile.open("resolution_vs_amp_BINP3_Fraction_v2.txt");    
       if(!infile.fail())
       {
        int raw = 0;
        isF1_open = true;
        while(!infile.eof())
        {
          if(infile.eof())
             break;
          float res, resError;
          infile >> res >> resError;
          std::cout << res << " " << resError << " " << raw << std::endl;

          double x,y;
          TimeResolution_vs_frac_wrtMiB2->GetPoint(raw,x,y);
          TimeResolution_vs_frac_wrtMiB2->SetPoint(raw,x,res);
          TimeResolution_vs_frac_wrtMiB2->SetPointEYlow(raw,resError); 
          TimeResolution_vs_frac_wrtMiB2->SetPointEYhigh(raw,resError); 
          raw++;
        }    
     }
     infile.close(); 
    }

    setStyle(); 
    
    TH2F* H2 = new TH2F("H2","",11,0.,110.,90,10.,40.);
    H2->GetXaxis()->SetTitle("events fraction (%)");
    //H2->GetXaxis()->SetTitleSize(0.045);
    H2->GetYaxis()->SetTitle("#sigma_{t} (ps)");

    TLegend* legend = new TLegend(0.58, 0.65, 0.72, 0.85);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.05);

    //legend -> AddEntry(TimeResolution_vs_frac_wrtMiB2_SIM_shifted,"Toy Simulation","PL");
    legend -> AddEntry(TimeResolution_vs_frac_wrtMiB2,"90x1+40x2","PL");
    
    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    H2->Draw();
    TimeResolution_vs_frac_wrtMiB2->Draw("P,same");
    //ReverseXAxis(TimeResolution_vs_frac_wrtMiB2);
    //ReverseXGraph(TimeResolution_vs_frac_wrtMiB2);
    legend->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw(); 
    c1 -> Print("BINP3_Resolution_vs_fraction_CFD50_wrtMiB2.png","png");
    c1 -> Print("BINP3_Resolution_vs_fraction_CFD50_wrtMiB2.pdf","pdf");  
}

void ReverseXAxis (TGraph *g)
{
   // Remove the current axis
   g->GetXaxis()->SetLabelOffset(999);
   g->GetXaxis()->SetTickLength(0);
   // Redraw the new axis
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
                                gPad->GetUymin(),
                                gPad->GetUxmin(),      
                                gPad->GetUymin(),
                                g->GetXaxis()->GetXmin(),
                                g->GetXaxis()->GetXmax(),
                                510,"-SDH");  
   newaxis->SetLabelOffset(-0.03);    
   newaxis->Draw();
}

void ReverseXGraph (TGraph *g)
{
   // Create a new graph
   Int_t n = g->GetN();
   Double_t *x = g->GetX();
   Double_t *y = g->GetY();
   Double_t xr[100];
   Double_t dx = g->GetXaxis()->GetXmin()+g->GetXaxis()->GetXmax();
   for (Int_t i=0; i<n; i++) {
      xr[i] = -x[i]+dx;
   }
   gr = new TGraph(n,xr,y);  
   gr->SetMarkerStyle(20);
   gr->SetLineColor(kRed);
   gr->SetMarkerColor(kRed);
   gr->Draw("PL");
}

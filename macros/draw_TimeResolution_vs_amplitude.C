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

void draw_TimeResolution_vs_amplitude()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
    gStyle->SetErrorX(0);

    
    TFile* inputFile = TFile::Open("Final_TimeResolution_vs_amp_BINPs_CFD50_thres20_onlyWrtMiB2.root");
     
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_SIM = (TGraphAsymmErrors*)inputFile->Get("TimeResolution_vs_amp_SIM_thres20._CFD50_normalNoise");
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_SIM_shifted = new TGraphAsymmErrors(); 
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputFile->Get("TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20_onlyWrtMiB2");   

    TimeResolution_vs_amp_wrtMiB2->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerSize(0.9);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerColor(kRed+1);
    TimeResolution_vs_amp_wrtMiB2->SetLineColor(kRed+1);
    TimeResolution_vs_amp_wrtMiB2->SetLineWidth(1);

    TF1* g_res;
    g_res = new TF1("g_res","sqrt([0]*[0]/(x*x)+[1]*[1])",0.,3500.); 
    g_res->SetParName(0,"a");   
    g_res->SetParName(1,"b"); 
    g_res->SetParameters(0,900.);
    g_res->SetParameters(1,5);
    g_res->SetParLimits(1,0.,40.);      
    TimeResolution_vs_amp_wrtMiB2->Fit("g_res","B");

    float p_const_data = 0.;
    float p_const_SIM = 0.;
    float delta = fabs(g_res->Eval(3000.)-TimeResolution_vs_amp_wrtMiB2_SIM->Eval(3000.));
    p_const_data = g_res->GetParameter(1);
    
    for(int ii = 1; ii < TimeResolution_vs_amp_wrtMiB2_SIM->GetN();ii++)
    {
        double x,y;
        double x_errorUp,y_errorUp, x_errorDown,y_errorDown;
        TimeResolution_vs_amp_wrtMiB2_SIM->GetPoint(ii,x,y);
        //std::cout << ii << " " << x << " " << y << std::endl;
        x_errorDown = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorXlow(ii);
        x_errorUp = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorXhigh(ii);
        y_errorDown = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorYlow(ii);
        y_errorUp = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorYhigh(ii);  

        float res = y+delta;
        if(ii == 0) TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPoint(ii,x,99999.);
        else{
           TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPoint(ii-1,x,res);
           TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEXlow(ii-1,x_errorDown); 
           TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEXhigh(ii-1,x_errorUp); 
           TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEYlow(ii-1,y_errorDown); 
           TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEYhigh(ii-1,y_errorUp); 
        }
    }

    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerStyle(24);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerSize(0.9);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerColor(kBlue+1);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetLineColor(kBlue+1);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetLineWidth(1);

    g_res->SetLineWidth(2);
    g_res->SetLineColor(kViolet+1);
    
    setStyle(); 
    
    float ampBinning[9] = {20.};
    for(int ii = 0; ii < TimeResolution_vs_amp_wrtMiB2->GetN();ii++)
    {
        double x,y;
        TimeResolution_vs_amp_wrtMiB2->GetPoint(ii,x,y);
        //std::cout << ii << " " << x-TimeResolution_vs_amp_wrtMiB2->GetErrorXlow(ii) <<  " " << x << " " << x+TimeResolution_vs_amp_wrtMiB2->GetErrorXhigh(ii) << " " << y << std::endl;
        ampBinning[ii+1] = x+TimeResolution_vs_amp_wrtMiB2->GetErrorXhigh(ii);
    }
    
    float resBinning[701] = {0.};
    for(int ii = 0; ii<= 700; ii++)
        resBinning[ii] = 10+ii*0.1;

    TH2F* H2 = new TH2F("H2","",8,ampBinning,700,resBinning);
    H2->GetXaxis()->SetTitle("amplitude (ADC counts)");
    H2->GetYaxis()->SetTitle("#sigma_{t} (ps)");

    TLegend* legend = new TLegend(0.58, 0.45, 0.72, 0.62);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);

    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2_SIM_shifted,"Toy Simulation","PL");
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2,"90x1+40x2","PL");
    
    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.04);
    latexLabel->SetTextColor(kViolet+1);
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica

    char LatexText[1000];
    sprintf(LatexText,"#sigma_{t} = #frac{a}{amplitude} #oplus b",(int)g_res->GetParameter(1),(int)g_res->GetParError(1));

    TLatex *latexLabel2 = new TLatex();
    latexLabel2->SetTextSize(0.04);
    latexLabel2->SetTextColor(kViolet+1);
    latexLabel2->SetNDC();
    latexLabel2->SetTextFont(42); // helvetica

    char LatexText2[1000];
    sprintf(LatexText2,"b = %.0f #pm %.0f ps",g_res->GetParameter(1),g_res->GetParError(1)); 

    TLatex *latexLabel3 = new TLatex();
    latexLabel3->SetTextSize(0.04);
    latexLabel3->SetTextColor(kViolet+1);
    latexLabel3->SetNDC();
    latexLabel3->SetTextFont(42); // helvetica

    char LatexText3[1000];
    //sprintf(LatexText3,"a = %.0f #pm %.0f ADC x ps",g_res->GetParameter(0),g_res->GetParError(0)); 
    sprintf(LatexText3,"a = 2.8 #pm 0.9 ADC x ns",g_res->GetParameter(0),g_res->GetParError(0)); 

    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    c1->SetLogx();
    H2->Draw();
    TimeResolution_vs_amp_wrtMiB2->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->Draw("P,same");  
    latexLabel->DrawLatex(0.58, 0.83,LatexText);
    latexLabel3->DrawLatex(0.58, 0.74,LatexText3);
    latexLabel2->DrawLatex(0.58, 0.65,LatexText2);
    //g_res->SetStats(0);
    g_res->Draw("same");
    c1->Update();
    TPaveStats *st = (TPaveStats*)TimeResolution_vs_amp_wrtMiB2->FindObject("stats");
    st->SetY1NDC(0.); //new x start position
    st->SetY2NDC(0.); //new x end position
    st->Draw();
    legend->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw(); 
    c1 -> Print("BINP3_ResolutionFit_CFD50_wrtMiB2.png","png");
    c1 -> Print("BINP3_ResolutionFit_CFD50_wrtMiB2.pdf","pdf");  
}

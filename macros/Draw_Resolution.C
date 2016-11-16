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

#include <iostream>
#include <fstream> 

void drawSIM(TFile* inputs);
void drawBINPS_wrtMiB2(TFile* inputs, std::string iTiming);
void drawBINPS_MiB2_compare(TFile* inputs, std::string iMCP, std::string iTiming);
void fitBINPs_wrtMiB2(TFile* inputs, std::string iMCP, std::string iTiming, std::string Fit);

void Draw_Resolution(std::string inputs,std::string iMCP,std::string iTiming,std::string Fit,int nPlot)
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    
    if(nPlot == 0) drawSIM(inputFile);
    else if(nPlot == 1) drawBINPS_MiB2_compare(inputFile,iMCP, iTiming);
    else if(nPlot == 2) drawBINPS_wrtMiB2(inputFile, iTiming);
    else if(nPlot == 3) fitBINPs_wrtMiB2(inputFile,iMCP, iTiming, Fit);
    
}

void drawSIM(TFile* inputs)
{
    TGraphAsymmErrors* TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise");
    TGraphAsymmErrors* TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise");
    TGraphAsymmErrors* TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres50_LED50_normalNoise");
    TGraphAsymmErrors* TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres50_LED50_doubleNoise");

    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetXaxis()->SetTitle("amp_max");
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetYaxis()->SetRangeUser(5E-4,0.3);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerStyle(20);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerSize(0.6);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerColor(kBlack);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetLineColor(kBlack); 
    
    TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise->SetMarkerStyle(20);
    TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise->SetMarkerSize(0.6);
    TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise->SetMarkerColor(kBlue+1);
    TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise->SetLineColor(kBlue+1); 

    TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise->SetMarkerStyle(20);
    TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise->SetMarkerSize(0.6);
    TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise->SetMarkerColor(kRed+1);
    TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise->SetLineColor(kRed+1); 

    TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise->SetMarkerStyle(24);
    TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise->SetMarkerSize(0.6);
    TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise->SetMarkerColor(kGreen+1);
    TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise->SetLineColor(kGreen+1); 

    TLegend* legend = new TLegend(0.52, 0.52, 0.69, 0.64);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise,"Noise: 3.5 ADC","P");
    legend -> AddEntry(TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise,"Noise: 7.0 ADC","P");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogy();
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->Draw("AP");
    //TimeResolution_vs_amp_SIM_thres20_LED50_normalNoise->Draw("P,same");
    TimeResolution_vs_amp_SIM_thres20_CFD50_doubleNoise->Draw("P,same");
    //TimeResolution_vs_amp_SIM_thres20_LED50_doubleNoise->Draw("P,same");
    legend ->Draw("same");
    c1 -> Print("SIM_Resolution.png","png");
    c1 -> Print("SIM_Resolution.pdf","pdf");
    
}

void drawBINPS_wrtMiB2(TFile* inputs, std::string iTiming)
{
    TGraphAsymmErrors* TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise");
    
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20;
    if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_"+iTiming+"_thres20").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_"+iTiming+"_thres50").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_"+iTiming+"_thres100").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_"+iTiming+"_thres150").c_str());

    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20;
    if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_"+iTiming+"_thres20_onlyWrtMiB2").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_"+iTiming+"_thres50_onlyWrtMiB2").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_"+iTiming+"_thres100_onlyWrtMiB2").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_"+iTiming+"_thres150_onlyWrtMiB2").c_str());

    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20;
    if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_"+iTiming+"_thres20_onlyWrtMiB2").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_"+iTiming+"_thres50_onlyWrtMiB2").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_"+iTiming+"_thres100_onlyWrtMiB2").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_"+iTiming+"_thres150_onlyWrtMiB2").c_str());

    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20;
    if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_"+iTiming+"_thres20").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_"+iTiming+"_thres50").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_"+iTiming+"_thres100").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_"+iTiming+"_thres150").c_str());

    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetXaxis()->SetTitle("amp_max");
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetYaxis()->SetTitle("#sigma_{t}(ps)");
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetXaxis()->SetRangeUser(0.,3100.);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->GetYaxis()->SetRangeUser(1.,200.);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerStyle(20);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerSize(0.6);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetMarkerColor(kBlack);
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->SetLineColor(kBlack);

    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->SetMarkerColor(kRed+1);
    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->SetLineColor(kRed+1);  
    
    TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20->SetMarkerColor(kViolet+1);
    TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20->SetLineColor(kViolet+1);  

    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->SetMarkerColor(kBlue+1);
    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->SetLineColor(kBlue+1);  

    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->SetMarkerColor(kGreen+1);
    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->SetLineColor(kGreen+1); 
    
    TLegend* legend = new TLegend(0.32, 0.52, 0.49, 0.74);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise,"Simulation","PL");
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20,std::string("BINP1 "+iTiming).c_str(),"PL");
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20,std::string("BINP2 "+iTiming).c_str(),"PL");
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20,std::string("BINP3 "+iTiming).c_str(),"PL");
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20,std::string("BINP4 "+iTiming).c_str(),"PL");

    /*TLegend* legendLog = new TLegend(0.52, 0.37, 0.69, 0.59);
    legendLog -> SetFillColor(kWhite);
    legendLog -> SetFillStyle(1000);
    legendLog -> SetLineWidth(0);
    legendLog -> SetLineColor(kWhite);
    legendLog -> SetTextFont(42);  
    legendLog -> SetTextSize(0.04);
    
    legendLog -> AddEntry(TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise,"Simulation","PL");
    legendLog -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20,std::string("BINP1 "+iTiming).c_str(),"PL");
    legendLog -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20,std::string("BINP3 "+iTiming).c_str(),"PL");
    legendLog -> AddEntry(TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20,std::string("BINP4 "+iTiming).c_str(),"PL");*/
    
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    //c1->SetLogy();
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->Draw("AP");
    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_BINP2_CFD50_thres20->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->Draw("P,same");
    legend ->Draw("same");
    c1 -> Print(std::string("BINPs_Resolution_"+iTiming+"_wrtMiB2.png").c_str(),"png");
    c1 -> Print(std::string("BINPs_Resolution_"+iTiming+"_wrtMiB2.pdf").c_str(),"pdf");

    /*TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise->Draw("AP");
    TimeResolution_vs_amp_wrtMiB2_BINP1_CFD50_thres20->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_BINP3_CFD50_thres20->Draw("P,same");
    TimeResolution_vs_amp_wrtMiB2_BINP4_CFD50_thres20->Draw("P,same");
    legendLog ->Draw("same");
    c2 -> Print("BINPs_Resolution_wrtMiB2_logy.png","png");
    c2 -> Print("BINPs_Resolution_wrtMiB2_logy.pdf","pdf");*/
}

void drawBINPS_MiB2_compare(TFile* inputs, std::string iMCP, std::string iTiming)
{
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2;
    if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres20").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres50").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres100").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres150").c_str());

    TGraphAsymmErrors* TimeResolution_vs_amp_onlyWrtMiB2;
    if(iTiming == "CFD50") TimeResolution_vs_amp_onlyWrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres20_onlyWrtMiB2").c_str());
    else if(iTiming == "LED50") TimeResolution_vs_amp_onlyWrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres50_onlyWrtMiB2").c_str());
    else if(iTiming == "LED100") TimeResolution_vs_amp_onlyWrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres100_onlyWrtMiB2").c_str());
    else if(iTiming == "LED150") TimeResolution_vs_amp_onlyWrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_"+iTiming+"_thres150_onlyWrtMiB2").c_str());

    
    TimeResolution_vs_amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
    TimeResolution_vs_amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
    TimeResolution_vs_amp_wrtMiB2->GetXaxis()->SetRangeUser(0.,3100.);
    TimeResolution_vs_amp_wrtMiB2->GetYaxis()->SetRangeUser(1.,200.);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerColor(kBlue+1);
    TimeResolution_vs_amp_wrtMiB2->SetLineColor(kBlue+1);

    TimeResolution_vs_amp_onlyWrtMiB2->SetMarkerStyle(24);
    TimeResolution_vs_amp_onlyWrtMiB2->SetMarkerSize(0.6);
    TimeResolution_vs_amp_onlyWrtMiB2->SetMarkerColor(kRed+1);
    TimeResolution_vs_amp_onlyWrtMiB2->SetLineColor(kRed+1);  
    
    TLegend* legend = new TLegend(0.32, 0.62, 0.49, 0.74);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(TimeResolution_vs_amp_wrtMiB2,std::string(iMCP+"_"+iTiming+"_wrtMIB2").c_str(),"PL");
    legend -> AddEntry(TimeResolution_vs_amp_onlyWrtMiB2,std::string(iMCP+"_"+iTiming+"_onlyWrtMIB2").c_str(),"PL");

    TCanvas* c1 = new TCanvas();
    c1->cd();
    //c1->SetLogy();
    TimeResolution_vs_amp_wrtMiB2->Draw("AP");
    TimeResolution_vs_amp_onlyWrtMiB2->Draw("P,same");
    legend ->Draw("same");
    c1 -> Print(std::string(iMCP+"_Resolution_"+iTiming+"_wrtMiB2.png").c_str(),"png");
    c1 -> Print(std::string(iMCP+"_Resolution_"+iTiming+"_wrtMiB2.pdf").c_str(),"pdf");
}

void fitBINPs_wrtMiB2(TFile* inputs, std::string iMCP, std::string iTiming, std::string Fit)
{
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_SIM = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise");
    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2_SIM_shifted = new TGraphAsymmErrors();  

    TGraphAsymmErrors* TimeResolution_vs_amp_wrtMiB2;
    if(iMCP == "SIM")TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get("TimeResolution_vs_amp_SIM_thres20_CFD50_normalNoise");
    else if(iMCP != "BINP3"){
       if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_CFD50_thres20").c_str());
       else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED50_thres50").c_str());
       else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED100_thres100").c_str());
       else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED150_thres150").c_str());
    }else{
       if(iTiming == "CFD50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_CFD50_thres20_onlyWrtMiB2").c_str());
       else if(iTiming == "LED50") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED50_thres50_onlyWrtMiB2").c_str());
       else if(iTiming == "LED100") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED100_thres100_onlyWrtMiB2").c_str());
       else if(iTiming == "LED150") TimeResolution_vs_amp_wrtMiB2 = (TGraphAsymmErrors*)inputs->Get(std::string("TimeResolution_vs_amp_wrtMiB2_"+iMCP+"_LED150_thres150_onlyWrtMiB2").c_str());
    }  
    TimeResolution_vs_amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
    TimeResolution_vs_amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
    TimeResolution_vs_amp_wrtMiB2->GetXaxis()->SetRangeUser(0.,3100.);
    TimeResolution_vs_amp_wrtMiB2->GetYaxis()->SetRangeUser(1.,200.);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerStyle(20);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2->SetMarkerColor(kBlack);
    TimeResolution_vs_amp_wrtMiB2->SetLineColor(kBlack);


    TF1* g_res;
    TF1* g_res_SIM;
    if(Fit == "f1"){
       g_res = new TF1("g_res","sqrt([0]*[0]/(x*x)+[1]*[1])",100.,2400.); 
       g_res->SetParName(0,"a");   
       g_res->SetParName(1,"b");    

       if(iMCP != "SIM"){
          g_res->SetParameters(0,900.);
          g_res->SetParameters(1,20);
          g_res->SetParLimits(1,10.,40.);
       }else if(iMCP == "BINP2"){
          g_res->SetParameters(0,100.);
          g_res->SetParameters(1,5);
          g_res->SetParLimits(0,0.,10000.);
          g_res->SetParLimits(1,0.,40.);
       }else {
          g_res->SetParameters(0,900.);
          g_res->SetParameters(1,5);
          g_res->SetParLimits(1,0.,40.);
       }
    }else if(Fit == "f2"){
       g_res = new TF1("g_res","sqrt([0]*[0]/(x*x)+[1]*[1]/(log(x)*log(x))+[2]*[2])",100.,2400.); 
       g_res->SetParName(0,"a");   
       g_res->SetParName(1,"b");   
       g_res->SetParName(2,"c");    

       if(iMCP != "SIM"){
          //g_res->SetParameters(0,900.);
          g_res->SetParameters(2,20);
          g_res->SetParLimits(2,10.,40.);
       }else if(iMCP == "BINP2"){
          g_res->SetParameters(2,5);
          g_res->SetParLimits(0,0.,999999999.);
          g_res->SetParLimits(1,0.,999999999.);
          g_res->SetParLimits(2,0.,40.);
       }else{
          //g_res->SetParameters(0,900.);
          g_res->SetParameters(2,5);
          g_res->SetParLimits(2,0.,40.);
       }
    }
 
    if(Fit == "f1"){
       g_res_SIM = new TF1("g_res_SIM","sqrt([0]*[0]/(x*x)+[1]*[1])",100.,2400.); 
       g_res_SIM->SetParameters(0,900.);
       g_res_SIM->SetParameters(1,20);
       g_res_SIM->SetParLimits(1,10.,40.);
    }else if(Fit == "f2"){
       g_res_SIM = new TF1("g_res_SIM","sqrt([0]*[0]/(x*x)+[1]*[1]/(log(x)*log(x))+[2]*[2])",100.,2400.); 
       //g_res_SIM->SetParameters(0,900.);
       g_res_SIM->SetParameters(2,20);
       g_res_SIM->SetParLimits(2,10.,40.);
    }

    TimeResolution_vs_amp_wrtMiB2->Fit("g_res","B");
    TimeResolution_vs_amp_wrtMiB2_SIM->Fit("g_res_SIM","B");

    float p_const_data = 0.;
    float p_const_SIM = 0.;
    if(Fit == "f1"){
       p_const_data = g_res->GetParameter(1);
       p_const_SIM = g_res_SIM->GetParameter(1);
    }else{
       p_const_data = g_res->GetParameter(2);
       p_const_SIM = g_res_SIM->GetParameter(2);
    }

    for(int ii = 0; ii < TimeResolution_vs_amp_wrtMiB2_SIM->GetN();ii++)
    {
        double x,y;
        double x_errorUp,y_errorUp, x_errorDown,y_errorDown;
        TimeResolution_vs_amp_wrtMiB2_SIM->GetPoint(ii,x,y);
        x_errorDown = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorXlow(ii);
        x_errorUp = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorXhigh(ii);
        y_errorDown = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorYlow(ii);
        y_errorUp = TimeResolution_vs_amp_wrtMiB2_SIM->GetErrorYhigh(ii);  

        float res = sqrt(y*y-p_const_SIM*p_const_SIM+p_const_data*p_const_data);
        TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPoint(ii,x,res);
        TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEXlow(ii,x_errorDown); 
        TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEXhigh(ii,x_errorUp); 
        TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEYlow(ii,y_errorDown); 
        TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetPointEYhigh(ii,y_errorUp); 
    }

    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerStyle(24);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerSize(0.6);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetMarkerColor(kBlue+1);
    TimeResolution_vs_amp_wrtMiB2_SIM_shifted->SetLineColor(kBlue+1);

    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.04);
    latexLabel->SetTextColor(kRed);
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica

    TCanvas* c1 = new TCanvas();
    c1->cd();
    if(iMCP != "SIM"){
       TimeResolution_vs_amp_wrtMiB2_SIM_shifted->Draw("AP");  
       TimeResolution_vs_amp_wrtMiB2->Draw("P,same");
    }else{
       TimeResolution_vs_amp_wrtMiB2->Draw("AP");
    }
    if(Fit == "f1") latexLabel->DrawLatex(0.65, 0.65,"#sigma_{t} = #frac{a}{amp_max} #oplus b");
    else if(Fit == "f2") latexLabel->DrawLatex(0.55, 0.65,"#sigma_{t} = #frac{a}{amp_max} #oplus #frac{b}{log(amp_max)} #oplus c"); 
    g_res->Draw("same");
    c1 -> Print(std::string(iMCP+"_ResolutionFit_"+iTiming+"_wrtMiB2_"+Fit+".png").c_str(),"png");
    c1 -> Print(std::string(iMCP+"_ResolutionFit_"+iTiming+"_wrtMiB2_"+Fit+".pdf").c_str(),"pdf");
}
     

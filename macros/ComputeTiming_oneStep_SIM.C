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

int nBins = 75;
float stepSize = 50.;

std::string ShiftVar(TTree* h4, std::string Var);
void FinalTiming(TTree* h4, TFile* inputFile, std::string Timing, std::string thresMCP, bool doScan);

void ComputeTiming_oneStep_SIM(std::string inputs, std::string Timing, std::string thresMCP, bool doScan = false)
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    TTree* h4 = (TTree*)inputFile->Get("h4");
    
    FinalTiming(h4, inputFile, Timing, thresMCP, doScan);
}

void FinalTiming(TTree* h4, TFile* inputFile, std::string Timing, std::string thresMCP, bool doScan)
{
    TH1F* h_time = new TH1F("h_time","",400,-0.2,0.2);
    TH2F* h_time_vs_amp = new TH2F("h_time_vs_amp","",75,0.,3750.,40000,-20.,20.);
  
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    std::string shiftedVar = ShiftVar(h4, std::string("time[BINP4"+iTiming+"]"));

    h4->Draw(std::string(shiftedVar+" >> h_time").c_str(),std::string("amp_max[BINP4]>"+thresMCP+" && amp_max[BINP4]<5000.").c_str()); 
    h4->Draw(std::string(shiftedVar+":amp_max[BINP4] >> h_time_vs_amp").c_str(),std::string("amp_max[BINP4]>"+thresMCP+" && amp_max[BINP4]<5000.").c_str()); 

    TF1* g_res = new TF1("g_res","gaus",-1.,1.);
    g_res->SetParameters(0,h_time->GetEntries()/2.);
    g_res->SetParameters(1,0.);
    g_res->SetParameters(2,0.05);
    g_res->SetParLimits(0,0.,h_time->GetEntries());
    g_res->SetParLimits(1,h_time->GetBinCenter(h_time->GetMaximumBin())-0.1,h_time->GetBinCenter(h_time->GetMaximumBin())+0.1);
    g_res->SetParLimits(2,0.,0.1);
    h_time->Fit("g_res","B");
    
    h_time->GetXaxis()->SetTitle("t-t_{ref} (ns)");
  
    TCanvas* c1 = new TCanvas();
    c1->cd();
    h_time->Draw("hist");
    g_res->Draw("same");
    c1 -> Print(std::string("TimeResolution_SIM_thres"+thresMCP+"_"+Timing+".png").c_str(),"png");
    c1 -> Print(std::string("TimeResolution_SIM_thres"+thresMCP+"_"+Timing+".pdf").c_str(),"pdf");
    
    h_time_vs_amp->FitSlicesY();
    TH2F* time_vs_amp_2 = (TH2F*)inputFile->Get("h_time_vs_amp_2");
    time_vs_amp_2->GetXaxis()->SetTitle("amp_max");
    time_vs_amp_2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    time_vs_amp_2->SetAxisRange(time_vs_amp_2->GetMinimum()-0.01,time_vs_amp_2->GetMaximum()+0.01, "Y");
    time_vs_amp_2->SetMarkerStyle(20);
    time_vs_amp_2->SetMarkerSize(0.9);
    time_vs_amp_2->SetMarkerColor(kBlack);
    time_vs_amp_2->SetLineColor(kBlack); 

    TCanvas* c2 = new TCanvas();
    c2->cd();
    time_vs_amp_2->Draw();
    c2 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+"_auto.png").c_str(),"png");
    c2 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+"_auto.pdf").c_str(),"pdf");

    TGraphAsymmErrors* g_Res_vs_Amp = new TGraphAsymmErrors(); 

    if(doScan == true)
    {
       std::vector<TH1F*> resHist;
       resHist.resize(nBins);
       std::vector<TF1*> resFit;
       resFit.resize(nBins);

       int iPoint = 0;
       
       std::string Selection = "";
       for(int ii = 0; ii <nBins; ii++)
       {
           char Name [50];
           sprintf (Name,"h_Res_%d",ii);
           resHist[ii] = new TH1F(Name,Name,40000,-20.,20.);     

           char cutMin [10];
           char cutMax [10];
           sprintf (cutMin,"%f",ii*stepSize);
           sprintf (cutMax,"%f",(ii+1)*stepSize);
           if(ii == 0) Selection = "amp_max[BINP4]>=0. && amp_max<"+ std::string(cutMax);
           else Selection = "amp_max[BINP4]>="+std::string(cutMin) +" && amp_max<"+ std::string(cutMax);

           h4->Draw(std::string("time[BINP4] >> "+std::string(Name)).c_str(),Selection.c_str());

           char NameFit [50];
           sprintf (NameFit,"f_Res_%d",ii);
           resFit[ii] = new TF1(NameFit,"gaus",resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())-1.,resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())+1.);
           resFit[ii]->SetParameters(0,resHist[ii]->GetEntries()/2);
           resFit[ii]->SetParameters(1,0.);
           resFit[ii]->SetParameters(2,0.5);
           resFit[ii]->SetParLimits(0,0.,resHist[ii]->GetEntries());
           resFit[ii]->SetParLimits(1,resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())-0.1,resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())+0.1);
           resFit[ii]->SetParLimits(2,0.,1.);

           if(resHist[ii]->GetEntries() != 0){
              resHist[ii]->Fit(NameFit,"B");
              g_Res_vs_Amp->SetPoint(iPoint,(ii+0.5)*stepSize,resFit[ii]->GetParameter(2));
              g_Res_vs_Amp->SetPointError(iPoint,(ii+0.5)*stepSize,(ii+0.5)*stepSize,resFit[ii]->GetParError(2)/2.,resFit[ii]->GetParError(2)/2.);
              iPoint++;
           }
       }

       g_Res_vs_Amp->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp->GetYaxis()->SetTitle("#sigma_{t}(ns)");
       g_Res_vs_Amp->SetMarkerStyle(20);
       g_Res_vs_Amp->SetMarkerSize(0.7);
       g_Res_vs_Amp->SetMarkerColor(kBlack);
       g_Res_vs_Amp->SetLineColor(kBlack);
       //g_Res_vs_Amp->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-0.01,0.15);

       TCanvas* c3 = new TCanvas();
       c3->cd();
       g_Res_vs_Amp->Draw("AP");
       c3 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+".png").c_str(),"png");
       c3 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+".pdf").c_str(),"pdf");   
    }

    TFile* output = new TFile(std::string("SIM_TimeResolution_vs_ampMax_"+Timing+"_thres"+thresMCP+".root").c_str(),"RECREATE");
    output->cd();
    time_vs_amp_2->Write(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+"_auto").c_str());
    g_Res_vs_Amp->Write(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing).c_str());
    output->Close();   
}

std::string ShiftVar(TTree* h4, std::string Var)
{
    TH1F* h = new TH1F("h","h",40000,-20.,20.);
    h4->Draw(std::string(Var+std::string(" >> h")).c_str(),"amp_max[BINP4]>20. && amp_max[BINP4]<5000.");
    
    
    char Mean [100];
    sprintf(Mean,"%f",h->GetMean());
    std::string sMean = std::string(Mean);

    std::string shiftedVar = "";
    if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       shiftedVar = Var+std::string("+")+std::string(sMean);
    }else{
       shiftedVar = Var+std::string("-")+std::string(sMean);
    }

    return shiftedVar;
    delete h;
}


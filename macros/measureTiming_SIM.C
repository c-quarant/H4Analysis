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
#include "TRandom3.h"

#include <iostream>
#include <fstream> 

int nBins = 75;
float stepSize = 50.;

std::string ShiftVar(TTree* h4, std::string Var);
void FinalTiming(TTree* h4, TTree* digi, TFile* inputFile, std::string Timing, std::string thresMCP, bool doScan);

void measureTiming_SIM(std::string inputs, std::string Timing, std::string thresMCP, bool doScan = true)
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    TTree* h4 = (TTree*)inputFile->Get("h4");
    TTree* digi = (TTree*)inputFile->Get("digi");
    
    FinalTiming(h4, digi, inputFile, Timing, thresMCP, doScan);
}

void FinalTiming(TTree* h4, TTree* digi, TFile* inputFile, std::string Timing, std::string thresMCP, bool doScan)
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

       int           BINP4;
       int           CFD;
       int           LED50;
       int           LED100;
       int           LED150;
       //ULong64_t       index;
       //Uint          n_channels;
       //Uint          n_timetypes;
       float         b_charge[1];   //[n_channels]
       float         b_slope[1];   //[n_channels]
       float         b_rms[1];   //[n_channels]
       float         time[4];   //[n_timetypes]
       float         time_chi2[4];   //[n_timetypes]
       float         maximum[1];   //[n_channels]
       float         time_maximum[1];   //[n_channels]
       float         amp_max[1];   //[n_channels]
       float         time_max[1];   //[n_channels]
       float         chi2_max[1];   //[n_channels]
       float         charge_tot[1];   //[n_channels]
       float         charge_sig[1];   //[n_channels]
       float         fit_ampl[1];   //[n_channels]
       float         fit_time[1];   //[n_channels]
       float         fit_chi2[1];   //[n_channels]
       float         calibration[1];   //[n_channels]

       TBranch        *b_BINP4;   //!
       TBranch        *b_CFD;   //!
       TBranch        *b_LED50;   //!
       TBranch        *b_LED100;   //!
       TBranch        *b_LED150;   //!
       //TBranch        *b_index;   //!
       //TBranch        *b_n_channels;   //!
       //TBranch        *b_n_timetypes;   //!
       TBranch        *b_b_charge;   //!
       TBranch        *b_b_slope;   //!
       TBranch        *b_b_rms;   //!
       TBranch        *b_time;   //!
       TBranch        *b_time_chi2;   //!
       TBranch        *b_maximum;   //!
       TBranch        *b_time_maximum;   //!
       TBranch        *b_amp_max;   //!
       TBranch        *b_time_max;   //!
       TBranch        *b_chi2_max;   //!
       TBranch        *b_charge_tot;   //!
       TBranch        *b_charge_sig;   //!
       TBranch        *b_fit_ampl;   //!
       TBranch        *b_fit_time;   //!
       TBranch        *b_fit_chi2;   //!
       TBranch        *b_calibration;   //!

       digi->SetBranchAddress("BINP4", &BINP4, &b_BINP4);
       digi->SetBranchAddress("CFD", &CFD, &b_CFD);
       digi->SetBranchAddress("LED50", &LED50, &b_LED50);
       digi->SetBranchAddress("LED100", &LED100, &b_LED100);
       digi->SetBranchAddress("LED150", &LED150, &b_LED150);
       //digi->SetBranchAddress("index", &index, &b_index);
       //digi->SetBranchAddress("n_channels", &n_channels, &b_n_channels);
       //digi->SetBranchAddress("n_timetypes", &n_timetypes, &b_n_timetypes);
       digi->SetBranchAddress("b_charge", b_charge, &b_b_charge);
       digi->SetBranchAddress("b_slope", b_slope, &b_b_slope);
       digi->SetBranchAddress("b_rms", b_rms, &b_b_rms);
       digi->SetBranchAddress("time", time, &b_time);
       digi->SetBranchAddress("time_chi2", time_chi2, &b_time_chi2);
       digi->SetBranchAddress("maximum", maximum, &b_maximum);
       digi->SetBranchAddress("time_maximum", time_maximum, &b_time_maximum);
       digi->SetBranchAddress("amp_max", amp_max, &b_amp_max);
       digi->SetBranchAddress("time_max", time_max, &b_time_max);
       digi->SetBranchAddress("chi2_max", chi2_max, &b_chi2_max);
       digi->SetBranchAddress("charge_tot", charge_tot, &b_charge_tot);
       digi->SetBranchAddress("charge_sig", charge_sig, &b_charge_sig);
       digi->SetBranchAddress("fit_ampl", fit_ampl, &b_fit_ampl);
       digi->SetBranchAddress("fit_time", fit_time, &b_fit_time);
       digi->SetBranchAddress("fit_chi2", fit_chi2, &b_fit_chi2);
       digi->SetBranchAddress("calibration", calibration, &b_calibration);

       std::vector<TH1F*> resHist;
       resHist.resize(nBins);
       std::vector<TF1*> resFit;
       resFit.resize(nBins);

       int iPoint = 0;
       
       TRandom3 rGen(0);
       
       std::string Selection = "";
       for(int ii = 0; ii <nBins; ii++)
       {
           char Name [50];
           sprintf (Name,"h_Res_%d",ii);
           resHist[ii] = new TH1F(Name,Name,40000,-20.,20.);     
       }

       for(int entry = 0; entry < digi->GetEntries(); entry++)
       {

        if(entry%1000==0) std::cout<<"--- Reading entry = "<<entry<<std::endl;
        digi->GetEntry(entry);    

        float smear = rGen.Gaus(0.,5.E-3);
        //std::cout << amp_max[BINP4] << " " << time[BINP4] << " " << smear << std::endl; 
        for(int ii = 0; ii <nBins; ii++)
        {
            if(ii == 0 && (amp_max[BINP4]>=20. && amp_max[BINP4]<(ii+1)*50.)) resHist[ii]->Fill(time[BINP4]+smear);   
            else if(amp_max[BINP4]>=ii*50. && amp_max[BINP4]<(ii+1)*50.) resHist[ii]->Fill(time[BINP4]+smear);       
        }

       }

       for(int ii = 0; ii <nBins; ii++)
       {

           char NameFit [50];
           sprintf (NameFit,"f_Res_%d",ii);
           resFit[ii] = new TF1(NameFit,"gaus",resHist[ii]->GetMean()-3.*resHist[ii]->GetRMS(),resHist[ii]->GetMean()+3.*resHist[ii]->GetRMS());
           resFit[ii]->SetParameters(0,resHist[ii]->GetEntries()/2.);
           resFit[ii]->SetParameters(1,0.);
           resFit[ii]->SetParameters(2,0.5);
           resFit[ii]->SetParLimits(0,0.,resHist[ii]->GetEntries());
           resFit[ii]->SetParLimits(1,resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())-1.,resHist[ii]->GetBinCenter(resHist[ii]->GetMaximumBin())+1.);
           resFit[ii]->SetParLimits(2,0.,1.);

           if(resHist[ii]->GetEntries() != 0){
              resHist[ii]->Fit(resFit[ii]->GetName(),"B");
              float sigma = resFit[ii]->GetParameter(2);
              std::cout << "Sigma = " << sigma << std::endl;
              g_Res_vs_Amp->SetPoint(iPoint,(ii+0.5)*stepSize,sigma*1000.);
              g_Res_vs_Amp->SetPointError(iPoint,0.5*stepSize,0.5*stepSize,resFit[ii]->GetParError(2)*1000.,resFit[ii]->GetParError(2)*1000.);
              iPoint++;
           }
       }

       g_Res_vs_Amp->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp->GetYaxis()->SetTitle("#sigma_{t}(ps)");
       g_Res_vs_Amp->SetMarkerStyle(20);
       g_Res_vs_Amp->SetMarkerSize(0.7);
       g_Res_vs_Amp->SetMarkerColor(kBlack);
       g_Res_vs_Amp->SetLineColor(kBlack);
       g_Res_vs_Amp->GetYaxis()->SetRangeUser(1.,200.);

       TCanvas* c3 = new TCanvas();
       c3->cd();
       c3->SetLogy();
       g_Res_vs_Amp->Draw("AP");
       c3 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+".png").c_str(),"png");
       c3 -> Print(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+".pdf").c_str(),"pdf");   
    }

    TFile* output = new TFile(std::string("SIM_TimeResolution_vs_ampMax_"+Timing+"_thres"+thresMCP+"_normalNoise.root").c_str(),"RECREATE");
    output->cd();
    time_vs_amp_2->Write(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+"_auto_normalNoise").c_str());
    g_Res_vs_Amp->Write(std::string("TimeResolution_vs_amp_SIM_thres"+thresMCP+"_"+Timing+"_normalNoise").c_str());
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


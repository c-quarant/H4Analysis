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

#include <iostream>
#include <fstream> 

void Timing(TTree* h4, std::string iMCP, TFile* inputFile);
void TimeCorrection(TTree* h4, std::string iMCP, TFile* inputFile);
void Amp_vs_TimeRel(TTree* h4, std::string iMCP);
void TimeMax_vs_Amp(TTree* h4, std::string iMCP);

void ComputeTiming(std::string inputs, std::string iMCP)
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());

    TTree* h4 = (TTree*)inputFile->Get("h4");

    Timing(h4, iMCP, inputFile);
    //TimeCorrection(h4, iMCP, inputFile);
    //Amp_vs_TimeRel(h4, iMCP);
    //TimeMax_vs_Amp(h4, iMCP);

}

void Timing(TTree* h4, std::string iMCP, TFile* inputFile)
{
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",400,-1,1);
    TH1F* time_wrtRm2 = new TH1F("time_wrtRm2","",400,-1,1);

    //TF1 *g_res = new TF1("g_res","gaus",-0.6,0.6);
    TF1 *g_res = new TF1("g_res","gaus",-1.,1.);
    //TF1 *g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-1.,1.);
   
    //--- CFD 0.5 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> time_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[MiB2]-(6.626+0.240*log(amp_max[M10]-0.354)) >> time_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> time_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[MiB2]-(8.128-0.121*log(amp_max[M10]+12.36)) >> time_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> time_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[MiB2]-(7.232+322.2*1/(amp_max[M10]+354.)) >> time_wrtMiB2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");*/
    //--- LED 150 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> time_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[MiB2]-(7.308+390.*1/(amp_max[M10]+310.2)) >> time_wrtMiB2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");*/
    g_res->SetParameters(1,-0.2,0.2);
    g_res->SetParameters(2,0.,0.3);
    time_wrtMiB2->Fit("g_res","B"); 

    TCanvas* c1 = new TCanvas();
    c1->cd();
    time_wrtMiB2->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c1 -> Print("TimeResolution_MiB_25mu_wrtMiB2_CFD50.png","png");
       c1 -> Print("TimeResolution_MiB_25mu_wrtMiB2_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c1 -> Print("TimeResolution_MiB_10mu_wrtMiB2_CFD50.png","png");
       c1 -> Print("TimeResolution_MiB_10mu_wrtMiB2_CFD50.pdf","pdf");
    }
      
    //--- CFD 0.5 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> time_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[Rm2]-(1.016+0.240*log(amp_max[M10]-1.184)) >> time_wrtRm2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2]-time[Rm2]-(-5.62-0.00001726*amp_max[MiB2]) >> time_wrtRm2","amp_max[MiB2]>20 && amp_max[Rm2]>200 && fabs(time[MiB2]-time[Rm2]+5.604) < 0.2 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2]-(0.9249-0.000005989*amp_max[MiB2]))<0.25");
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> time_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[Rm2]-(2.513-0.121*log(amp_max[M10]+11.91)) >> time_wrtRm2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED50]-time[Rm2]-(-5.081-0.158*log(amp_max[MiB2]-24.41)) >> time_wrtRm2","amp_max[MiB2]>50 && amp_max[Rm2]>200 && fabs(time[MiB2+LED50]-time[Rm2]+6.001) < 0.2 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED50]-(0.123+0.199*log(amp_max[MiB2]+70.)))<0.25");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> time_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[Rm2]-(1.62+319.6*1/(amp_max[M10]+351.7)) >> time_wrtRm2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");
     else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED100]-time[Rm2]-(-4.863-0.177*log(amp_max[MiB2]-55.78)) >> time_wrtRm2","amp_max[MiB2]>100 && amp_max[Rm2]>200 && fabs(time[MiB2+LED100]-time[Rm2]+5.874) < 0.2 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED100]-(-0.348+0.251*log(amp_max[MiB2]+80.11)))<0.25");*/
    //--- LED 150 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> time_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[Rm2]-(1.696+387.7*1/(amp_max[M10]+307.3)) >> time_wrtRm2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED150]-time[Rm2]-(-4.529-0.214*log(amp_max[MiB2]-36.48)) >> time_wrtRm2","amp_max[MiB2]>150 && amp_max[Rm2]>200 && fabs(time[MiB2+LED150]-time[Rm2]+5.785) < 0.2 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED150]-(-0.910+0.317*log(amp_max[MiB2]+170.)))<0.25");*/
    //g_res->SetParameters(0,0.,time_wrtRm2->GetEntries()/2.);
    //g_res->SetParameters(1,-0.2,0.2);
    //g_res->SetParameters(2,0.,0.1);
    /*g_res->SetParameters(3,0.,time_wrtRm2->GetEntries()/2.);
    g_res->SetParameters(4,0.,0.5);
    g_res->SetParLimits(0,0.,time_wrtRm2->GetEntries()/2.);
    g_res->SetParLimits(1,-0.2,0.2);
    g_res->SetParLimits(2,0.,0.1);
    g_res->SetParLimits(3,0.,time_wrtRm2->GetEntries()/2.);
    g_res->SetParLimits(4,0.,0.5);*/
    g_res->SetParameters(1,-0.2,0.2);
    g_res->SetParameters(2,0.,0.3);
    time_wrtRm2->Fit("g_res","B");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    time_wrtRm2->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c2 -> Print("TimeResolution_MiB_25mu_wrtRm2_CFD50.png","png");
       c2 -> Print("TimeResolution_MiB_25mu_wrtRm2_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c2 -> Print("TimeResolution_MiB_10mu_wrtRm2_CFD50.png","png");
       c2 -> Print("TimeResolution_MiB_10mu_wrtRm2_CFD50.pdf","pdf");
    }else if(iMCP == "MiB2"){
       c2 -> Print("TimeResolution_MiB2_wrtRm2_CFD50.png","png");
       c2 -> Print("TimeResolution_MiB2_wrtRm2_CFD50.pdf","pdf");
    }

    int nBins = 10;
    float stepSize = 340.; 

    std::vector<TH1F*> resHist;
    resHist.resize(nBins);
    std::vector<TF1*> resFit;
    resFit.resize(nBins);
    std::vector<TF1*> resFitAlt;
    resFitAlt.resize(nBins);

    TGraphAsymmErrors* g_Res_vs_Amp_wrtMiB2 = new TGraphAsymmErrors(); 
    TGraphAsymmErrors* g_Res_vs_Amp_wrtRm2 = new TGraphAsymmErrors(); 

    for(int ii = 0; ii < nBins; ii++)
    {
        char Name [50];
        sprintf (Name,"h_Res_1_%d",ii);
        resHist[ii] = new TH1F(Name,Name,400,-2.,2.);

        char NameFit [50];
        sprintf (NameFit,"f_Res_1_%d",ii);
        resFit[ii] = new TF1(NameFit,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-0.5,0.5);

        char NameFitAlt [50];
        sprintf (NameFitAlt,"f_ResAlt_2_%d",ii);
        resFitAlt[ii] = new TF1(NameFitAlt,"gaus",-0.5,0.5);
        
        //--- CFD 0.5 ---
        char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10]-time[MiB2]-(6.626+0.240*log(amp_max[M10]-0.354)) >> %s",Name);
        //--- LED 50 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED50]-time[MiB2]-(8.128-0.121*log(amp_max[M10]+12.36)) >> %s",Name);*/
        //--- LED 100 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED100]-time[MiB2]-(7.232+322.2*1/(amp_max[M10]+354.)) >> %s",Name);*/
        //--- LED 150 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED150]-time[MiB2]-(7.308+390.*1/(amp_max[M10]+310.2)) >> %s",Name);*/

        if(iMCP == "M25")
           h4->Draw(drawM25,cutM25); 
        else if(iMCP == "M10")
           h4->Draw(drawM10,cutM10);

        char NameOutputM25_png [500];
        sprintf (NameOutputM25_png,"TimeResolution_MiB_25mu_wrtMiB2_%d_CFD50.png",ii+1);
        char NameOutputM25_pdf [500];
        sprintf (NameOutputM25_pdf,"TimeResolution_MiB_25mu_wrtMiB2_%d_CFD50.pdf",ii+1);
        char NameOutputM10_png [500];
        sprintf (NameOutputM10_png,"TimeResolution_MiB_10mu_wrtMiB2_%d_CFD50.png",ii+1);
        char NameOutputM10_pdf [500];
        sprintf (NameOutputM10_pdf,"TimeResolution_MiB_10mu_wrtMiB2_%d_CFD50.pdf",ii+1);

        resFit[ii]->SetParameters(0,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParameters(1,-0.2,0.2);
        resFit[ii]->SetParameters(2,0.,0.2);
        resFit[ii]->SetParameters(3,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParameters(4,0.,0.5);
        resFit[ii]->SetParLimits(0,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParLimits(1,-0.2,0.2);
        resFit[ii]->SetParLimits(2,0.,0.2);
        resFit[ii]->SetParLimits(3,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParLimits(4,0.,0.5);
        resHist[ii]->Fit(NameFit,"B");

        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){ 
           resFitAlt[ii]->SetParameters(1,-0.2,0.2);
           resFitAlt[ii]->SetParameters(2,0.,0.2);
           resFitAlt[ii]->SetParLimits(1,-0.2,0.2);
           resFitAlt[ii]->SetParLimits(2,0.,0.2);
           resHist[ii]->Fit(NameFitAlt,"B");
        }
   
        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){
           g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt[ii]->GetParameter(2));
           g_Res_vs_Amp_wrtMiB2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt[ii]->GetParError(2),resFitAlt[ii]->GetParError(2));
        }else{     
           g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ii*stepSize+stepSize/2.,resFit[ii]->GetParameter(2));
           g_Res_vs_Amp_wrtMiB2->SetPointError(ii,stepSize/2.,stepSize/2.,resFit[ii]->GetParError(2),resFit[ii]->GetParError(2));
        }
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        resHist[ii]->Draw("hist");
        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){
           resFitAlt[ii]->Draw("same");
        }else{     
           resFit[ii]->Draw("same");
        }
        if(iMCP == "M25"){
           c2 -> Print(NameOutputM25_png,"png");
           c2 -> Print(NameOutputM25_pdf,"pdf");
        }else if(iMCP == "M10"){
           c2 -> Print(NameOutputM10_png,"png");
           c2 -> Print(NameOutputM10_pdf,"pdf");
        }
        
        delete c2;
    }

    g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(0.05,0.16);

    TCanvas* c3 = new TCanvas();
    c3->cd();
    g_Res_vs_Amp_wrtMiB2->Draw("AP");
    if(iMCP == "M25"){
       c3 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_CFD50.png","png");
       c3 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c3 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_CFD50.png","png");
       c3 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_CFD50.pdf","pdf");
    }

    resHist.clear();
    resHist.resize(30);
    resFit.clear();
    resFit.resize(30);
    resFitAlt.clear();
    resFitAlt.resize(30);
    
    for(int ii = 0; ii < nBins; ii++)
    {
        char Name [50];
        sprintf (Name,"h_Res_2_%d",ii);
        resHist[ii] = new TH1F(Name,Name,400,-2.,2.);

        char NameFit [50];
        sprintf (NameFit,"f_Res_2_%d",ii);
        resFit[ii] = new TF1(NameFit,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-0.5,0.5);

        char NameFitAlt [50];
        sprintf (NameFitAlt,"f_ResAlt_2_%d",ii);
        resFitAlt[ii] = new TF1(NameFitAlt,"gaus",-0.5,0.5);
        
        //--- CFD 0.5 ---
        char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10]-time[Rm2]-(1.016+0.240*log(amp_max[M10]-1.184)) >> %s",Name);
        //--- LED 50 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED50]-time[Rm2]-(2.513-0.121*log(amp_max[M10]+11.91)) >> %s",Name);*/
        //--- LED 100 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED100]-time[Rm2]-(1.62+319.6*1/(amp_max[M10]+351.7)) >> %s",Name);*/
        //--- LED 150 ---
        /*char cutM25 [1000];
        sprintf (cutM25,"amp_max[M25]>%f && amp_max[M25]<%f && amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM25 [1000];
        sprintf (drawM25,"time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.71)) >> %s",Name);

        char cutM10 [1000];
        sprintf (cutM10,"amp_max[M10]>%f && amp_max[M10]<%f && amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25",ii*stepSize,(ii+1)*stepSize);
        char drawM10 [1000];
        sprintf (drawM10,"time[M10+LED150]-time[Rm2]-(1.696+387.7*1/(amp_max[M10]+307.3)) >> %s",Name);*/

        if(iMCP == "M25")
           h4->Draw(drawM25,cutM25); 
        else if(iMCP == "M10")
           h4->Draw(drawM10,cutM10);

        char NameOutputM25_png [500];
        sprintf (NameOutputM25_png,"TimeResolution_MiB_25mu_wrtRm2_%d_CFD50.png",ii+1);
        char NameOutputM25_pdf [500];
        sprintf (NameOutputM25_pdf,"TimeResolution_MiB_25mu_wrtRm2_%d_CFD50.pdf",ii+1);
        char NameOutputM10_png [500];
        sprintf (NameOutputM10_png,"TimeResolution_MiB_10mu_wrtRm2_%d_CFD50.png",ii+1);
        char NameOutputM10_pdf [500];
        sprintf (NameOutputM10_pdf,"TimeResolution_MiB_10mu_wrtRm2_%d_CFD50.pdf",ii+1);

        resFit[ii]->SetParameters(0,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParameters(1,-0.2,0.2);
        resFit[ii]->SetParameters(2,0.,0.2);
        resFit[ii]->SetParameters(3,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParameters(4,0.,0.5);
        resFit[ii]->SetParLimits(0,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParLimits(1,-0.2,0.2);
        resFit[ii]->SetParLimits(2,0.,0.2);
        resFit[ii]->SetParLimits(3,0.,resHist[ii]->GetEntries()/2.);
        resFit[ii]->SetParLimits(4,0.,0.5);
        resHist[ii]->Fit(NameFit,"B");
        
        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){ 
           resFitAlt[ii]->SetParameters(1,-0.2,0.2);
           resFitAlt[ii]->SetParameters(2,0.,0.2);
           resFitAlt[ii]->SetParLimits(1,-0.2,0.2);
           resFitAlt[ii]->SetParLimits(2,0.,0.2);
           resHist[ii]->Fit(NameFitAlt,"B");
        }
   
        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){
           g_Res_vs_Amp_wrtRm2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt[ii]->GetParameter(2));
           g_Res_vs_Amp_wrtRm2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt[ii]->GetParError(2),resFitAlt[ii]->GetParError(2));
        }else{     
           g_Res_vs_Amp_wrtRm2->SetPoint(ii,ii*stepSize+stepSize/2.,resFit[ii]->GetParameter(2));
           g_Res_vs_Amp_wrtRm2->SetPointError(ii,stepSize/2.,stepSize/2.,resFit[ii]->GetParError(2),resFit[ii]->GetParError(2));
        }
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        resHist[ii]->Draw("hist");
        if(resFit[ii]->GetParameter(0) == 0 || resFit[ii]->GetParameter(1) == 0 || resFit[ii]->GetParameter(2) == 0 || resFit[ii]->GetParameter(3) == 0 || resFit[ii]->GetParameter(4) == 0){
           resFitAlt[ii]->Draw("same");
        }else{     
           resFit[ii]->Draw("same");
        }
        if(iMCP == "M25"){
           c2 -> Print(NameOutputM25_png,"png");
           c2 -> Print(NameOutputM25_pdf,"pdf");
        }else if(iMCP == "M10"){
           c2 -> Print(NameOutputM10_png,"png");
           c2 -> Print(NameOutputM10_pdf,"pdf");
        }
        
        delete c2;
    }

    g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(0.05,0.16);

    TCanvas* c4 = new TCanvas();
    c4->cd();
    g_Res_vs_Amp_wrtRm2->Draw("AP");
    if(iMCP == "M25"){
       c4 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_25mu_CFD50.png","png");
       c4 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c4 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_10mu_CFD50.png","png");
       c4 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_10mu_CFD50.pdf","pdf");
    }
    
    TH2F* timeResol_vs_amp_wrtMiB2 = new TH2F("timeResol_vs_amp_wrtMiB2","",12,0.,3600.,500,-20.,20.);
    TH2F* timeResol_vs_amp_wrtRm2 = new TH2F("timeResol_vs_amp_wrtRm2","",12,0.,3600.,500,-20.,20.);
     
    //--- CFD 0.5 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[MiB2]-(6.626+0.240*log(amp_max[M10]-0.354)):amp_max[M10] >> timeResol_vs_amp_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[MiB2]-(8.128-0.121*log(amp_max[M10]+12.36)):amp_max[M10] >> timeResol_vs_amp_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[MiB2]-(7.232+322.2*1/(amp_max[M10]+354.)):amp_max[M10] >> timeResol_vs_amp_wrtMiB2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");*/
    //--- LED 150 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]-(8.611+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[MiB2]-(7.308+390.*1/(amp_max[M10]+310.2)):amp_max[M10] >> timeResol_vs_amp_wrtMiB2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");*/

    timeResol_vs_amp_wrtMiB2->FitSlicesY();
    TH2F* timeResol_vs_amp_wrtMiB2_2 = (TH2F*)inputFile->Get("timeResol_vs_amp_wrtMiB2_2");    timeResol_vs_amp_wrtMiB2_2->Draw(); 

    timeResol_vs_amp_wrtMiB2_2->GetXaxis()->SetTitle("amp_max");
    timeResol_vs_amp_wrtMiB2_2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    timeResol_vs_amp_wrtMiB2_2->SetAxisRange(0.05,0.16, "Y");
    //timeResol_vs_amp_wrtMiB2_2->SetAxisRange(0.,1200., "X");
    timeResol_vs_amp_wrtMiB2_2->SetMarkerStyle(20);
    timeResol_vs_amp_wrtMiB2_2->SetMarkerSize(0.9);
    timeResol_vs_amp_wrtMiB2_2->SetMarkerColor(kBlack);
    timeResol_vs_amp_wrtMiB2_2->SetLineColor(kBlack);

    TCanvas* c5 = new TCanvas();
    c5->cd();
    timeResol_vs_amp_wrtMiB2_2->Draw();
    if(iMCP == "M25"){
       c5 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_FitSlicesY_CFD50.png","png");
       c5 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_FitSlicesY_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c5 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_FitSlicesY_CFD50.png","png");
       c5 -> Print("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_FitSlicesY_CFD50.pdf","pdf");
    }

    //--- CFD 0.5 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(8.611+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[Rm2]-(1.016+0.240*log(amp_max[M10]-1.184)):amp_max[M10] >> timeResol_vs_amp_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 &&b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[Rm2]-(2.513-0.121*log(amp_max[M10]+11.91)):amp_max[M10] >> timeResol_vs_amp_wrtRm2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[Rm2]-(1.62+319.6*1/(amp_max[M10]+351.7)):amp_max[M10] >> timeResol_vs_amp_wrtRm2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");*/
    //--- LED 150 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]-(2.746+0.126*log(amp_max[M25]-28.78)):amp_max[M25] >> timeResol_vs_amp_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[M25]-time[MiB2]-9.434) < 1 && fabs(time[Rm2]-time[M25]+3.577) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25");
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[Rm2]-(1.696+387.7*1/(amp_max[M10]+307.3)):amp_max[M10] >> timeResol_vs_amp_wrtRm2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");*/

    timeResol_vs_amp_wrtRm2->FitSlicesY();
    TH2F* timeResol_vs_amp_wrtRm2_2 = (TH2F*)inputFile->Get("timeResol_vs_amp_wrtRm2_2");    timeResol_vs_amp_wrtRm2_2->Draw(); 

    timeResol_vs_amp_wrtRm2_2->GetXaxis()->SetTitle("amp_max");
    timeResol_vs_amp_wrtRm2_2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    timeResol_vs_amp_wrtRm2_2->SetAxisRange(0.05,0.16, "Y");
    //timeResol_vs_amp_wrtRm2_2->SetAxisRange(0.,1200., "X");
    timeResol_vs_amp_wrtRm2_2->SetMarkerStyle(20);
    timeResol_vs_amp_wrtRm2_2->SetMarkerSize(0.9);
    timeResol_vs_amp_wrtRm2_2->SetMarkerColor(kBlack);
    timeResol_vs_amp_wrtRm2_2->SetLineColor(kBlack);

    TCanvas* c6 = new TCanvas();
    c6->cd();
    timeResol_vs_amp_wrtRm2_2->Draw();
    if(iMCP == "M25"){
       c6 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_25mu_FitSlicesY_CFD50.png","png");
       c6 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_25mu_FitSlicesY_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c6 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_10mu_FitSlicesY_CFD50.png","png");
       c6 -> Print("TimeResolution_vs_amp_wrtRm2_MiB_10mu_FitSlicesY_CFD50.pdf","pdf");
    }
}

void TimeCorrection(TTree* h4, std::string iMCP, TFile* inputFile)
{
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",30,0.,3000.,500,-20.,20.);
    TH2F* timingCorrection_wrtRm2 = new TH2F("timingCorrection_wrtRm2","",30,0.,3000.,550,-20.,20.);

    //--- CFD 0.5 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[MiB2]:amp_max[M10] >> timingCorrection_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");*/
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[MiB2]:amp_max[M10] >> timingCorrection_wrtMiB2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[MiB2]:amp_max[M10] >> timingCorrection_wrtMiB2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");*/
    //--- LED 150 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[MiB2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[MiB2]:amp_max[M10] >> timingCorrection_wrtMiB2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");

    timingCorrection_wrtMiB2->FitSlicesY();
    TH2F* timingCorrection_wrtMiB2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtMiB2_1");    //timingCorrection_wrtMiB2_1->Draw(); 

    if(iMCP == "M25") timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle("amp_max[25]");
    else if(iMCP == "M10") timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle("amp_max[10]");
    timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time-time[MiB2]");
    timingCorrection_wrtMiB2_1->SetAxisRange(7.,9., "Y");
    //timingCorrection_wrtMiB2_1->SetAxisRange(0.,1200., "X");
    timingCorrection_wrtMiB2_1->SetMarkerStyle(20);
    timingCorrection_wrtMiB2_1->SetMarkerSize(0.9);
    timingCorrection_wrtMiB2_1->SetMarkerColor(kBlack);
    timingCorrection_wrtMiB2_1->SetLineColor(kBlack);

    //TF1* fit_corr1 = new TF1("fit_corr1","[0]+[1]*log(x+[2])",150,3000.);
    //timingCorrection_wrtMiB2_1->Fit("fit_corr1");
    TF1* fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",150.,3000.);
    timingCorrection_wrtMiB2_1->Fit("fit_corr1");

    
    TCanvas* c3 = new TCanvas();
    c3->cd();
    timingCorrection_wrtMiB2_1->Draw();
    fit_corr1->Draw("same");
    if(iMCP == "M25"){
       c3 -> Print("timingCorrection_wrtMiB2_MiB_25mu_CFD50.png","png");
       c3 -> Print("timingCorrection_wrtMiB2_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c3 -> Print("timingCorrection_wrtMiB2_MiB_10mu_CFD50.png","png");
       c3 -> Print("timingCorrection_wrtMiB2_MiB_10mu_CFD50.pdf","pdf");
    }

    //--- CFD 0.5 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]:amp_max[M25] >> timingCorrection_wrtRm2","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10]-time[Rm2]:amp_max[M10] >> timingCorrection_wrtRm2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10]-(-0.238+0.287*log(amp_max[M10]+70.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2]-time[Rm2]:amp_max[MiB2] >> timingCorrection_wrtRm2","amp_max[MiB2]>20 && amp_max[Rm2]>200 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2]-(0.925-0.000005989*amp_max[MiB2]))<0.25");*/
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED50]-time[Rm2]:amp_max[M10] >> timingCorrection_wrtRm2","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED50]-(-1.5821+0.621*log(amp_max[M10]+30.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED50]-time[Rm2]:amp_max[MiB2] >> timingCorrection_wrtRm2","amp_max[MiB2]>50 && amp_max[Rm2]>200 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED50]-(-0.124+0.199*log(amp_max[MiB2]+70.)))<0.25");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED100]-time[Rm2]:amp_max[M10] >> timingCorrection_wrtRm2","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED100]-(-3.024+0.796*log(amp_max[M10]+80.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED100]-time[Rm2]:amp_max[MiB2] >> timingCorrection_wrtRm2","amp_max[MiB2]>100 && amp_max[Rm2]>200 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED100]-(-0.348+0.251*log(amp_max[MiB2]+80.11)))<0.25");*/
    //--- LED 150 ---
    if(iMCP == "M25")
       h4->Draw("time[M25]-time[Rm2]:amp_max[M25] >> timingCorrection_wrtMiB2","amp_max[M25]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4 && fabs(time_max[M25]-time[M25]-(1.116+0.000275*amp_max[M25]))<0.25"); 
    else if(iMCP == "M10")
       h4->Draw("time[M10+LED150]-time[Rm2]:amp_max[M10] >> timingCorrection_wrtRm2","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4 && fabs(time_max[M10]-time[M10+LED150]-(-3.946+0.901*log(amp_max[M10]+130.)))<0.25");
    else if(iMCP == "MiB2")
       h4->Draw("time[MiB2+LED150]-time[Rm2]:amp_max[MiB2] >> timingCorrection_wrtRm2","amp_max[MiB2]>150 && amp_max[Rm2]>200 && b_rms[MiB2] <= 4 && fabs(time_max[MiB2]-time[MiB2+LED150]-(-0.910+0.317*log(amp_max[MiB2]+170.)))<0.25");
    
    timingCorrection_wrtRm2->FitSlicesY();
    TH2F* timingCorrection_wrtRm2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtRm2_1");    timingCorrection_wrtRm2_1->Draw(); 

    if(iMCP == "M25") timingCorrection_wrtRm2_1->GetXaxis()->SetTitle("amp_max[25]");
    else if(iMCP == "M10") timingCorrection_wrtRm2_1->GetXaxis()->SetTitle("amp_max[10]");
    else if(iMCP == "MiB2") timingCorrection_wrtRm2_1->GetXaxis()->SetTitle("amp_max[MiB2]");
    timingCorrection_wrtRm2_1->GetYaxis()->SetTitle("time-time[Rm2]");
    timingCorrection_wrtRm2_1->SetAxisRange(-6.5,-5., "Y");
    //timingCorrection_wrtRm2_1->SetAxisRange(0.,1200., "X");
    timingCorrection_wrtRm2_1->SetMarkerStyle(20);
    timingCorrection_wrtRm2_1->SetMarkerSize(0.9);
    timingCorrection_wrtRm2_1->SetMarkerColor(kBlack);
    timingCorrection_wrtRm2_1->SetLineColor(kBlack);

    TF1* fit_corr2 = new TF1("fit_corr2","[0]+[1]*log(x+[2])",150.,3000.);
    //TF1* fit_corr2 = new TF1("fit_corr2","pol1",20.,3000.);
    timingCorrection_wrtRm2_1->Fit("fit_corr2");
    //TF1* fit_corr2 = new TF1("fit_corr2","[0]+[1]*1/(x+[2])",20.,3000.);
    //timingCorrection_wrtRm2_1->Fit("fit_corr2");

    TCanvas* c4 = new TCanvas();
    c4->cd();
    timingCorrection_wrtRm2_1->Draw();
    fit_corr2->Draw("same");
    if(iMCP == "M25"){
       c4 -> Print("timingCorrection_wrtRm2_MiB_25mu_CFD50.png","png");
       c4 -> Print("timingCorrection_wrtRm2_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c4 -> Print("timingCorrection_wrtRm2_MiB_10mu_CFD50.png","png");
       c4 -> Print("timingCorrection_wrtRm2_MiB_10mu_CFD50.pdf","pdf");
    }else if(iMCP == "MiB2"){
       c4 -> Print("timingCorrection_wrtRm2_MiB2_CFD50.png","png");
       c4 -> Print("timingCorrection_wrtRm2_MiB2_CFD50.pdf","pdf");
    }

}

void Amp_vs_TimeRel(TTree* h4, std::string iMCP)
{
  
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",2250,-50,175,300,-1.,1.5,1000.,1500.);
    TProfile2D* p2D_amp_vs_time_zoom1 = new TProfile2D("p2D_amp_vs_time_zoom1","",300,-10,20,300,-1.,1.5,1000.,1500.);
    TProfile2D* p2D_amp_vs_time_zoom2 = new TProfile2D("p2D_amp_vs_time_zoom2","",300,-10,20,160,0.5,1.1,1000.,1500.);

    if(iMCP == "M25")
       h4->Draw("maximum[M25]:WF_val/maximum[M25]:WF_time-time[M25] >> p2D_amp_vs_time","amp_max[M25]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M25 && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("maximum[M10]:WF_val/maximum[M10]:WF_time-time[M10] >> p2D_amp_vs_time","amp_max[M10]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M10  && fabs(time[MiB2]-time[M10]+7.966) < 1 && fabs(time[Rm2]-time[M10]+2.354) < 1 && fabs(time[MiB2]-time[Rm2]+5.609) < 0.2 && b_rms[M10] <= 4","goff");
    
    if(iMCP == "M25") p2D_amp_vs_time->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    else if(iMCP == "M10") p2D_amp_vs_time->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    if(iMCP == "M25") p2D_amp_vs_time->GetYaxis()->SetTitle("WF_val/amp_max[M25]");
    else if(iMCP == "M10") p2D_amp_vs_time->GetYaxis()->SetTitle("WF_val/amp_max[M10]");
    p2D_amp_vs_time->SetAxisRange(450.,1750., "Z");
    
    TCanvas* c3 = new TCanvas();
    c3->cd();
    p2D_amp_vs_time->Draw("COLZ");
    if(iMCP == "M25"){
       c3 -> Print("amp_max_vs_time_MiB_25mu_CFD50.png","png");
       c3 -> Print("amp_max_vs_time_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c3 -> Print("amp_max_vs_time_MiB_10mu_CFD50.png","png");
       c3 -> Print("amp_max_vs_time_MiB_10mu_CFD50.pdf","pdf");
    }

    if(iMCP == "M25")
       h4->Draw("amp_max[M25]:WF_val/amp_max[M25]:WF_time-time[M25] >> p2D_amp_vs_time_zoom1","amp_max[M25]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M25 && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goof");
    else if(iMCP == "M10")
       h4->Draw("amp_max[M10]:WF_val/amp_max[M10]:WF_time-time[M10] >> p2D_amp_vs_time_zoom1","amp_max[M10]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M10  && fabs(time[MiB2]-time[M10]+7.966) < 1 && fabs(time[Rm2]-time[M10]+2.354) < 1 && fabs(time[MiB2]-time[Rm2]+5.609) < 0.2 && b_rms[M10] <= 4","goff");
    
    if(iMCP == "M25") p2D_amp_vs_time_zoom1->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    else if(iMCP == "M10") p2D_amp_vs_time_zoom1->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    if(iMCP == "M25") p2D_amp_vs_time_zoom1->GetYaxis()->SetTitle("WF_val/amp_max[M25]");
    else if(iMCP == "M10") p2D_amp_vs_time_zoom1->GetYaxis()->SetTitle("WF_val/amp_max[M10]");
    p2D_amp_vs_time_zoom1->SetAxisRange(450.,1750., "Z");
    
    TCanvas* c4 = new TCanvas();
    c4->cd();
    p2D_amp_vs_time_zoom1->Draw("COLZ");
    if(iMCP == "M25"){
       c4 -> Print("amp_max_vs_time_MiB_25mu_zoom1_CFD50.png","png");
       c4 -> Print("amp_max_vs_time_MiB_25mu_zoom1_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c4 -> Print("amp_max_vs_time_MiB_10mu_zoom1_CFD50.png","png");
       c4 -> Print("amp_max_vs_time_MiB_10mu_zoom1_CFD50.pdf","pdf");
    }

    if(iMCP == "M25")
       h4->Draw("amp_max[M25]:WF_val/amp_max[M25]:WF_time-time[M25] >> p2D_amp_vs_time_zoom2","amp_max[M25]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M25 && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("amp_max[M10]:WF_val/amp_max[M10]:WF_time-time[M10] >> p2D_amp_vs_time_zoom2","amp_max[M10]>500 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M10  && fabs(time[MiB2]-time[M10]+7.966) < 1 && fabs(time[Rm2]-time[M10]+2.354) < 1 && fabs(time[MiB2]-time[Rm2]+5.609) < 0.2 && b_rms[M10] <= 4","goff");
    
    if(iMCP == "M25") p2D_amp_vs_time_zoom2->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    else if(iMCP == "M10") p2D_amp_vs_time_zoom2->GetXaxis()->SetTitle("WF_time-time[M25] (ns)");
    if(iMCP == "M25") p2D_amp_vs_time_zoom2->GetYaxis()->SetTitle("WF_val/amp_max[M25]");
    else if(iMCP == "M10") p2D_amp_vs_time_zoom2->GetYaxis()->SetTitle("WF_val/amp_max[M10]");
    p2D_amp_vs_time_zoom2->SetAxisRange(450.,1750., "Z");
    
    TCanvas* c5 = new TCanvas();
    c5->cd();
    p2D_amp_vs_time_zoom2->Draw("COLZ");
    if(iMCP == "M25"){
       c5 -> Print("amp_max_vs_time_MiB_25mu_zoom2_CFD50.png","png");
       c5 -> Print("amp_max_vs_time_MiB_25mu_zoom2_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c5 -> Print("amp_max_vs_time_MiB_10mu_zoom2_CFD50.png","png");
       c5 -> Print("amp_max_vs_time_MiB_10mu_zoom2_CFD50.pdf","pdf");
    }
}

void TimeMax_vs_Amp(TTree* h4, std::string iMCP)
{
    TH2D* h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    TH2D* h2_time_maximum_vs_amp = new TH2D("h2_time_maximum_vs_amp","",300,0.,3000.,500,0.,5.);

    TF1* pol1_max = new TF1("pol1_max","pol1",0.,3000.);
    TF1* pol1_maximum = new TF1("pol1_maximum","pol1",0.,3000.);
    TF1* fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",150.,3000.);
    TF1* fit_corr_maximum = new TF1("fit_corr_maximum","[0]+[1]*log(x+[2])",150.,3000.);

    //--- CFD 0.5 ---
    /*if(iMCP == "M25")
       h4->Draw("time_max[M25]-time[M25]:amp_max[M25] >> h2_time_max_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M25 && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_max[M10]-time[M10]:amp_max[M10] >> h2_time_max_vs_amp","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_max[MiB2]-time[MiB2]:amp_max[MiB2] >> h2_time_max_vs_amp","(amp_max[MiB2]>20 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.604) < 0.2 && b_rms[MiB2] <= 4","goff");*/
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time_max[M25]-time[M25]:amp_max[M25] >> h2_time_max_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_max[M10]-time[M10+LED50]:amp_max[M10] >> h2_time_max_vs_amp","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_max[MiB2]-time[MiB2+LED50]:amp_max[MiB2] >> h2_time_max_vs_amp","(amp_max[MiB2]>50 && amp_max[Rm2]>200) && fabs(time[MiB2+LED50]-time[Rm2]+6.001) < 0.2 && b_rms[MiB2] <= 4","goff");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time_max[M25]-time[M25]:amp_max[M25] >> h2_time_max_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_max[M10]-time[M10+LED100]:amp_max[M10] >> h2_time_max_vs_amp","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_max[MiB2]-time[MiB2+LED100]:amp_max[MiB2] >> h2_time_max_vs_amp","(amp_max[MiB2]>100 && amp_max[Rm2]>200) && fabs(time[MiB2+LED100]-time[Rm2]+5.874) < 0.2 && b_rms[MiB2] <= 4","goff");*/  
    //--- LED 150 ---
    if(iMCP == "M25")
       h4->Draw("time_max[M25]-time[M25]:amp_max[M25] >> h2_time_max_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_max[M10]-time[M10+LED150]:amp_max[M10] >> h2_time_max_vs_amp","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_max[MiB2]-time[MiB2+LED150]:amp_max[MiB2] >> h2_time_max_vs_amp","(amp_max[MiB2]>150 && amp_max[Rm2]>200) && fabs(time[MiB2+LED150]-time[Rm2]+5.785) < 0.2 && b_rms[MiB2] <= 4","goff");
    

    fit_corr_max->SetParameters(2,130.,170.);
    fit_corr_max->SetParLimits(2,130.,170.);
    
    if(iMCP == "M25"){
       h2_time_max_vs_amp->GetYaxis()->SetTitle("time_max-time[M25] (ns)");
       h2_time_max_vs_amp->GetXaxis()->SetTitle("amp_max[M25]");
    }else if(iMCP == "M10"){
       h2_time_max_vs_amp->GetYaxis()->SetTitle("time_max-time[M10] (ns)");
       h2_time_max_vs_amp->GetXaxis()->SetTitle("amp_max[M10]");
    }else if(iMCP == "MiB2"){
       h2_time_max_vs_amp->GetYaxis()->SetTitle("time_max-time[MiB2] (ns)");
       h2_time_max_vs_amp->GetXaxis()->SetTitle("amp_max[MiB2]");
    }
    h2_time_max_vs_amp->SetAxisRange(0.,1750., "Z");
    if(iMCP == "M25") h2_time_max_vs_amp->Fit("");
    else if(iMCP == "M10") h2_time_max_vs_amp->Fit("fit_corr_max","B","",150.,3000.);
    else if(iMCP == "MiB2") h2_time_max_vs_amp->Fit("fit_corr_max","B","",150.,3000.);
    
    TCanvas* c6 = new TCanvas();
    c6->cd();
    h2_time_max_vs_amp->Draw("COLZ");
    if(iMCP == "M25") fit_corr_max->Draw("same");
    else if(iMCP == "M10") fit_corr_max->Draw("same");
    else if(iMCP == "MiB2") fit_corr_max->Draw("same");
    if(iMCP == "M25"){
       c6 -> Print("deltaT_max_vs_amp_max_MiB_25mu_CFD50.png","png");
       c6 -> Print("deltaT_max_vs_amp_max_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c6 -> Print("deltaT_max_vs_amp_max_MiB_10mu_CFD50.png","png");
       c6 -> Print("deltaT_max_vs_amp_max_MiB_10mu_CFD50.pdf","pdf");
    }else if(iMCP == "MiB2"){
       c6 -> Print("deltaT_max_vs_amp_max_MiB2_CFD50.png","png");
       c6 -> Print("deltaT_max_vs_amp_max_MiB2_CFD50.pdf","pdf");
    }

    //--- CFD 0.5 ---
    /*if(iMCP == "M25")
       h4->Draw("time_maximum[M25]-time[M25]:amp_max[M25] >> h2_time_maximum_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && WF_ch == M25 && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_maximum[M10]-time[M10]:amp_max[M10] >> h2_time_maximum_vs_amp","time_maximum[M10]-time[M10]<3.5 && amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10]+8.186) < 1 && fabs(time[Rm2]-time[M10]+2.576) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_maximum[MiB2]-time[MiB2]:amp_max[MiB2] >> h2_time_maximum_vs_amp","(amp_max[MiB2]>20 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[Rm2]+5.604) < 0.2 && b_rms[MiB2] <= 4","goff");*/
    //--- LED 50 ---
    /*if(iMCP == "M25")
       h4->Draw("time_maximum[M25]-time[M25]:amp_max[M25] >> h2_time_maximum_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_maximum[M10]-time[M10+LED50]:amp_max[M10] >> h2_time_maximum_vs_amp","amp_max[M10]>50 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED50]+7.326) < 1 && fabs(time[Rm2]-time[M10+LED50]+1.721) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_maximum[MiB2]-time[MiB2+LED50]:amp_max[MiB2] >> h2_time_maximum_vs_amp","(amp_max[MiB2]>50 && amp_max[Rm2]>200) && fabs(time[MiB2+LED50]-time[Rm2]+6.001) < 0.2 && b_rms[MiB2] <= 4","goff");*/
    //--- LED 100 ---
    /*if(iMCP == "M25")
       h4->Draw("time_maximum[M25]-time[M25]:amp_max[M25] >> h2_time_maximum_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_maximum[M10]-time[M10+LED100]:amp_max[M10] >> h2_time_maximum_vs_amp","amp_max[M10]>100 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED100]+7.617) < 1 && fabs(time[Rm2]-time[M10+LED100]+2.042) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_maximum[MiB2]-time[MiB2+LED100]:amp_max[MiB2] >> h2_time_maximum_vs_amp","(amp_max[MiB2]>100 && amp_max[Rm2]>200) && fabs(time[MiB2+LED100]-time[Rm2]+5.874) < 0.2 && b_rms[MiB2] <= 4","goff");*/ 
    //--- LED 150 ---
    if(iMCP == "M25")
       h4->Draw("time_maximum[M25]-time[M25]:amp_max[M25] >> h2_time_maximum_vs_amp","amp_max[M25]>20 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M25]+9.347) < 1 && fabs(time[Rm2]-time[M25]+3.800) < 1 && fabs(time[MiB2]-time[Rm2]+5.863) < 0.2 && b_rms[M25] <= 4","goff");
    else if(iMCP == "M10")
       h4->Draw("time_maximum[M10]-time[M10+LED150]:amp_max[M10] >> h2_time_maximum_vs_amp","amp_max[M10]>150 && (amp_max[MiB2]>200 && amp_max[Rm2]>200) && fabs(time[MiB2]-time[M10+LED150]+7.631) < 1 && fabs(time[Rm2]-time[M10+LED150]+2.026) < 1 && fabs(time[MiB2]-time[Rm2]+5.611) < 0.2 && b_rms[M10] <= 4","goff");
    else if(iMCP == "MiB2")
       h4->Draw("time_maximum[MiB2]-time[MiB2+LED150]:amp_max[MiB2] >> h2_time_maximum_vs_amp","(amp_max[MiB2]>150 && amp_max[Rm2]>200) && fabs(time[MiB2+LED150]-time[Rm2]+5.785) < 0.2 && b_rms[MiB2] <= 4","goff");

    fit_corr_maximum->SetParameters(2,130.,170.);
    fit_corr_maximum->SetParLimits(2,130.,170.);
    
    if(iMCP == "M25"){
       h2_time_maximum_vs_amp->GetYaxis()->SetTitle("time_maximum-time[M25] (ns)");
       h2_time_maximum_vs_amp->GetXaxis()->SetTitle("amp_max[M25]");
    }else if(iMCP == "M10"){
       h2_time_maximum_vs_amp->GetYaxis()->SetTitle("time_maximum-time[M10] (ns)");
       h2_time_maximum_vs_amp->GetXaxis()->SetTitle("amp_max[M10]");
    }else if(iMCP == "MiB2"){
       h2_time_maximum_vs_amp->GetYaxis()->SetTitle("time_maximum-time[MiB2] (ns)");
       h2_time_maximum_vs_amp->GetXaxis()->SetTitle("amp_max[MiB2]");
    }
    h2_time_maximum_vs_amp->SetAxisRange(0.,1750., "Z");
    if(iMCP == "M25") h2_time_maximum_vs_amp->Fit("pol1_maximum");
    else if(iMCP == "M10") h2_time_maximum_vs_amp->Fit("fit_corr_maximum","B","",150.,3000.);
    else if(iMCP == "MiB2") h2_time_maximum_vs_amp->Fit("fit_corr_maximum","B","",150.,3000.);
    TCanvas* c7 = new TCanvas();
    c7->cd();
    h2_time_maximum_vs_amp->Draw("COLZ");
    if(iMCP == "M25") fit_corr_maximum->Draw("same");
    else if(iMCP == "M10") fit_corr_maximum->Draw("same");
    else if(iMCP == "MiB2") fit_corr_maximum->Draw("same");
    if(iMCP == "M25"){
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB_25mu_CFD50.png","png");
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB_25mu_CFD50.pdf","pdf");
    }else if(iMCP == "M10"){
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB_10mu_CFD50.png","png");
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB_10mu_CFD50.pdf","pdf");
    }else if(iMCP == "MiB2"){
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB2_CFD50.png","png");
       c7 -> Print("deltaT_maximum_vs_amp_max_MiB2_CFD50.pdf","pdf");
    }
}
	

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


//int nBins = 26;
//float ampMin[27] = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650, 700., 750., 800., 850., 900., 950., 1000., 1250, 1500., 1750., 2000.,2500.,3000.};
int nBins = 7;
float ampMin[8] = {0., 400., 800., 1200., 1700., 2000.,2500.,3000.};
float timeMin[4001];

//BINP
std::string amp_max_MiB2 = "200.";
std::string time_max_MiB2 = "150.";
std::string time_max_Rm2 = "150.";
std::string amp_max_Rm2 = "200.";
std::string scintMin = "200.";
std::string scintMax = "700.";
std::string timeChi2 = "99999999.";
bool doScan_Corr = false;


void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan, bool doScanEff, bool doDoubleGauss, std::string HodoSelection, bool doOnlyWrtMiB2);
void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan_Corr, std::string HodoSelection, bool doOnlyWrtMiB2);
void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2);
void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2);
void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void CheckEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void HodoValidation(TTree* h4, TTree* hodo, TTree* adc, int nParticle);
std::string AddSelection(TTree*, std::string, std::string, std::string, bool);
std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel);

void ComputeTiming_oneStep(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2 = false, bool doFirstStep = true, bool doPulseShapes = false, bool doScan = false, bool doScanEff = false, bool doDoubleGauss = false, std::string hodoXmin="0", std::string hodoXmax="30", std::string hodoYmin="0", std::string hodoYmax="30")
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    TTree* h4 = (TTree*)inputFile->Get("h4");
    TTree* hodo = (TTree*)inputFile->Get("hodo");
    TTree* adc = (TTree*)inputFile->Get("adc");
    
    std::vector<float>* Params = new std::vector<float>;
    std::vector<float>* Params_wrtMiB2 = new std::vector<float>;
    std::vector<float>* Params_wrtRm2 = new std::vector<float>;

    for(int ii = 0; ii < 4001; ii++)
        timeMin[ii] = 0.01*ii-20.;

    std::string nameiMCP = iMCP;
    if(iMCP == "M25") nameiMCP = "MiB_25mu";
    else if(iMCP == "M10") nameiMCP = "MiB_10mu";
    else if(iMCP == "M8") nameiMCP = "Rm_8mu";
    else if(iMCP == "M5") nameiMCP = "Rm_5mu";

    std::string HodoSelection = " && X>"+hodoXmin+" && X<"+hodoXmax+" && Y>"+hodoYmin+" && Y<"+hodoYmax;
    
    //HodoValidation(h4, hodo, adc, 2);
    if(doFirstStep == true){
       TimeChi2(h4, iMCP, nameiMCP, thresMCP, maxMCP, doOnlyWrtMiB2);
       AmpMax(h4, iMCP, nameiMCP, thresMCP, maxMCP, doOnlyWrtMiB2);
       Hodoscope(h4, iMCP, nameiMCP, thresMCP, maxMCP, doOnlyWrtMiB2);
       CheckEfficiency(h4, inputs, iMCP, nameiMCP, thresMCP, maxMCP);
    }else if(doPulseShapes == true) PulseShapes(h4, iMCP, nameiMCP, thresMCP, maxMCP, HodoSelection, doOnlyWrtMiB2);
    else{
       AmpVsTime_Selection(h4, iMCP, nameiMCP,Timing, Params, thresMCP, maxMCP, HodoSelection, doOnlyWrtMiB2);
       TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, maxMCP, doScan_Corr, HodoSelection, doOnlyWrtMiB2);
       FinalTiming(h4, inputs, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, maxMCP, doScan, doScanEff, doDoubleGauss, HodoSelection, doOnlyWrtMiB2);
    }
}

void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan, bool doScanEff, bool doDoubleGauss, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",400,-1.,1.);
    TH1F* time_wrtRm2 = new TH1F("time_wrtRm2","",400,-1.,1.);
    TH1F* time_wrtMiB2_noCorrection = new TH1F("time_wrtMiB2_noCorrection","",1600,-10.,10.);
    TH1F* time_wrtRm2_noCorrection = new TH1F("time_wrtRm2_noCorrection","",1600,-10.,10.);
    TH2F* time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins,ampMin,4000,timeMin);
    TH2F* time_vs_amp_wrtRm2 = new TH2F("time_vs_amp_wrtRm2","",nBins,ampMin,4000,timeMin);

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    //TF1 *g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-5.,5.);
    TF1 *g_res;

    char Selection1 [1000];  
    char Selection2 [1000];
    char Selection3 [1000];
    char Selection4 [1000];
    char Selection5 [1000];
    
    //(time-time_max) vs amp selection
    //if(iMCP == "M25") sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*amp_max[")+iMCP+std::string("]))<0.25")).c_str(),Params->at(0),Params->at(1));
    sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    //time correction
    if(iMCP != "MiB2"){
       sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
    }
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
    }

    //no time correction
    std::string Selection6 = "time["+iMCP+iTiming+"]-time[MiB2] >>";
    std::string Selection7 = "time["+iMCP+iTiming+"]-time[Rm2] >>";
    std::string Selection8;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection8,"0.",false);
      Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection8,"0.",false);
      //Selection8 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection8,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
         //Selection8 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection8,"0.",false);
      }else{
         Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection8,"0.",false);
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
         //Selection8 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection8,"0.",false);
      }   
    }
    std::string Selection9 = Selection8+std::string(" && ")+std::string(Selection1)+HodoSelection; 

    std::string Selection10 = Selection2+std::string(" time_wrtMiB2");
    std::string Selection11 = Selection3+std::string(" time_vs_amp_wrtMiB2");
    std::string Selection12 = Selection4+std::string(" time_wrtRm2");
    std::string Selection13 = Selection5+std::string(" time_vs_amp_wrtRm2"); 
    std::string Selection14 = Selection6+std::string(" time_wrtMiB2_noCorrection");
    std::string Selection15 = Selection7+std::string(" time_wrtRm2_noCorrection");

    //std::cout << "wrt MiB2 = " << Selection14.c_str() << ", " << Selection9.c_str() << std::endl;
    //std::cout << "wrt Rm2 = " << Selection15.c_str() << ", " << Selection9.c_str() << std::endl;

    if(iMCP != "MiB2"){
       h4->Draw(Selection10.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str()); 
       h4->Draw(Selection11.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str()); 
       h4->Draw(Selection14.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str()); 
    }
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       h4->Draw(Selection12.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str());
       h4->Draw(Selection13.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str()); 
       h4->Draw(Selection15.c_str(),std::string(Selection9+" && amp_max["+iMCP+"]>"+thresMCP).c_str()); 
    }

    if(time_wrtMiB2->GetEntries() < 2000) time_wrtMiB2->Rebin(8);
    if(time_wrtRm2->GetEntries() < 2000) time_wrtRm2->Rebin(8);

    if(iMCP != "MiB2"){
       if(doDoubleGauss == false) g_res = new TF1("g_res","gaus",-1.,1.);
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",-1.,1.);  
       g_res->SetParameters(0,time_wrtMiB2->GetEntries()/2.);
       g_res->SetParameters(1,0.);
       g_res->SetParameters(2,0.1);
       if(doDoubleGauss == true) g_res->SetParameters(3,time_wrtMiB2->GetEntries()/2.);
       if(doDoubleGauss == true) g_res->SetParameters(4,0.);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.1);
       g_res->SetParLimits(0,0.,time_wrtMiB2->GetEntries());
       g_res->SetParLimits(1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       g_res->SetParLimits(2,0.,0.5);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtMiB2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,0.5);
       time_wrtMiB2->Fit("g_res","B"); 
       
       char Sigma[100];
       TLatex *latexLabel = new TLatex();
       float sigma_eff;
       float s_sigma;
       if(doDoubleGauss == true){
          float f1 = g_res->GetParameter(0)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float f2 = g_res->GetParameter(3)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          sigma_eff = sqrt(f1*g_res->GetParameter(2)*g_res->GetParameter(2)+f2*g_res->GetParameter(5)*g_res->GetParameter(5) + f1*g_res->GetParameter(1)*g_res->GetParameter(1)+f2*g_res->GetParameter(4)*g_res->GetParameter(4) - (f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4))*(f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4)));
          s_sigma = sqrt((f1*f1*g_res->GetParameter(2)*g_res->GetParameter(2)*g_res->GetParError(2)*g_res->GetParError(2)+f2*f2*g_res->GetParameter(4)*g_res->GetParameter(4)*g_res->GetParError(4)*g_res->GetParError(4)/(sigma_eff*sigma_eff))); 
       }else{
           sigma_eff = g_res->GetParameter(2); 
           s_sigma = g_res->GetParError(2);
       }

       if(iMCP != "MiB2" && iMCP != "Rm2") sigma_eff = sqrt(sigma_eff*sigma_eff-23.E-3*23.E-3);
       
       sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

       latexLabel->SetTextSize(0.05);
       latexLabel->SetNDC();
       latexLabel->SetTextFont(42); // helvetica

       time_wrtMiB2->GetXaxis()->SetTitle("t-t_{ref} (ns)");

       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";
  
       TCanvas* c1 = new TCanvas();
       c1->cd();
       time_wrtMiB2->Draw("hist");
       latexLabel->DrawLatex(0.72, 0.55,Sigma);
       g_res->Draw("same");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
       delete g_res; 
    }

    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){ 
       if(doDoubleGauss == false) g_res = new TF1("g_res","gaus",-1.,1.);
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",-1.,1.);
       g_res->SetParameters(0,time_wrtRm2->GetEntries()/2.);
       g_res->SetParameters(1,0.);
       g_res->SetParameters(2,0.1);
       if(doDoubleGauss == true) g_res->SetParameters(3,time_wrtRm2->GetEntries()/2.);
       if(doDoubleGauss == true) g_res->SetParameters(4,0.);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.1);
       g_res->SetParLimits(0,0.,time_wrtRm2->GetEntries());
       g_res->SetParLimits(1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       g_res->SetParLimits(2,0.,0.5);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtRm2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,0.5);
       time_wrtRm2->Fit("g_res","B"); 
       
       char Sigma[100];
       TLatex *latexLabel = new TLatex();
       float sigma_eff;
       float s_sigma; 
       if(doDoubleGauss == true){
          float f1 = g_res->GetParameter(0)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float f2 = g_res->GetParameter(3)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          sigma_eff = sqrt(f1*g_res->GetParameter(2)*g_res->GetParameter(2)+f2*g_res->GetParameter(5)*g_res->GetParameter(5) + f1*g_res->GetParameter(1)*g_res->GetParameter(1)+f2*g_res->GetParameter(4)*g_res->GetParameter(4) - (f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4))*(f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4)));
          s_sigma = sqrt((f1*f1*g_res->GetParameter(2)*g_res->GetParameter(2)*g_res->GetParError(2)*g_res->GetParError(2)+f2*f2*g_res->GetParameter(4)*g_res->GetParameter(4)*g_res->GetParError(4)*g_res->GetParError(4)/(sigma_eff*sigma_eff)));
       }else{
           sigma_eff = g_res->GetParameter(2); 
           s_sigma = g_res->GetParError(2);
       }

       if(iMCP != "MiB2" && iMCP != "Rm2") sigma_eff = sqrt(sigma_eff*sigma_eff-17.E-3*17.E-3);
       
       sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

       latexLabel->SetTextSize(0.05);
       latexLabel->SetNDC();
       latexLabel->SetTextFont(42); // helvetica

       time_wrtRm2->GetXaxis()->SetTitle("t-t_{ref} (ns)");
        
       TCanvas* c2 = new TCanvas();
       c2->cd();
       time_wrtRm2->Draw("hist");
       latexLabel->DrawLatex(0.72, 0.55,Sigma);
       g_res->Draw("same");
       c2 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
       delete g_res; 

    }   

    if(iMCP != "MiB2"){
       if(doDoubleGauss == false) g_res = new TF1("g_res","gaus",time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.1);
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.1);
       g_res->SetParameters(0,0.,time_wrtMiB2_noCorrection->GetEntries());
       g_res->SetParameters(1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
       g_res->SetParameters(2,0.,1.);
       if(doDoubleGauss == true) g_res->SetParameters(3,0.,time_wrtMiB2_noCorrection->GetEntries());
       if(doDoubleGauss == true) g_res->SetParameters(4,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.,1.);
       g_res->SetParLimits(0,0.,time_wrtMiB2_noCorrection->GetEntries());
       g_res->SetParLimits(1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
       g_res->SetParLimits(2,0.,1.);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtMiB2_noCorrection->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,1.);
       time_wrtMiB2_noCorrection->Fit("g_res","B"); 

       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       TCanvas* c3 = new TCanvas();
       c3->cd();
       time_wrtMiB2_noCorrection->Draw("hist");
       g_res->Draw("same");
       c3 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+"_noCorrection.png").c_str(),"png");
       c3 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+"_noCorrection.pdf").c_str(),"pdf");
       delete g_res; 
    }
    
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.1);
       g_res->SetParameters(0,0.,time_wrtRm2_noCorrection->GetEntries());
       g_res->SetParameters(1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
       g_res->SetParameters(2,0.,1.);
       if(doDoubleGauss == true) g_res->SetParameters(3,0.,time_wrtRm2_noCorrection->GetEntries());
       if(doDoubleGauss == true) g_res->SetParameters(4,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.,1.);
       g_res->SetParLimits(0,0.,time_wrtRm2_noCorrection->GetEntries());
       g_res->SetParLimits(1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
       g_res->SetParLimits(2,0.,1.);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtRm2_noCorrection->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,1.);
       time_wrtRm2_noCorrection->Fit("g_res","B");

       TCanvas* c4 = new TCanvas();
       c4->cd();
       time_wrtRm2_noCorrection->Draw("hist");
       g_res->Draw("same");
       c4 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtRm2_"+Timing+"_thres"+thresMCP+"_noCorrection.png").c_str(),"png");
       c4 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtRm2_"+Timing+"_thres"+thresMCP+"_noCorrection.pdf").c_str(),"pdf");
       delete g_res;
    }

    TH2F* time_vs_amp_wrtMiB2_2;
    TH2F* time_vs_amp_wrtRm2_2;

    if(iMCP != "MiB2"){
       time_vs_amp_wrtMiB2->FitSlicesY();
       time_vs_amp_wrtMiB2_2 = (TH2F*)inputFile->Get("time_vs_amp_wrtMiB2_2");
       time_vs_amp_wrtMiB2_2->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       time_vs_amp_wrtMiB2_2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
       time_vs_amp_wrtMiB2_2->SetAxisRange(time_vs_amp_wrtMiB2_2->GetMinimum()-0.01,time_vs_amp_wrtMiB2_2->GetMaximum()+0.01, "Y");
       time_vs_amp_wrtMiB2_2->SetMarkerStyle(20);
       time_vs_amp_wrtMiB2_2->SetMarkerSize(0.9);
       time_vs_amp_wrtMiB2_2->SetMarkerColor(kBlack);
       time_vs_amp_wrtMiB2_2->SetLineColor(kBlack); 

       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       TCanvas* c5 = new TCanvas();
       c5->cd();
       time_vs_amp_wrtMiB2_2->Draw();
       c5 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+"_auto.png").c_str(),"png");
       c5 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+"_auto.pdf").c_str(),"pdf");
    }

    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       time_vs_amp_wrtRm2->FitSlicesY();
       time_vs_amp_wrtRm2_2 = (TH2F*)inputFile->Get("time_vs_amp_wrtRm2_2");
       time_vs_amp_wrtRm2_2->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       time_vs_amp_wrtRm2_2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
       time_vs_amp_wrtRm2_2->SetAxisRange(time_vs_amp_wrtRm2_2->GetMinimum()-0.01,time_vs_amp_wrtRm2_2->GetMaximum()+0.01, "Y");
       time_vs_amp_wrtRm2_2->SetMarkerStyle(20);
       time_vs_amp_wrtRm2_2->SetMarkerSize(0.9);
       time_vs_amp_wrtRm2_2->SetMarkerColor(kBlack);
       time_vs_amp_wrtRm2_2->SetLineColor(kBlack); 

       TCanvas* c6 = new TCanvas();
       c6->cd();
       time_vs_amp_wrtRm2_2->Draw();
       c6 -> Print(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_auto.png").c_str(),"png");
       c6 -> Print(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_auto.pdf").c_str(),"pdf");
   }

   TGraphAsymmErrors* g_Res_vs_Amp_wrtMiB2 = new TGraphAsymmErrors(); 
   TGraphAsymmErrors* g_Res_vs_Amp_wrtRm2 = new TGraphAsymmErrors(); 

   TGraphAsymmErrors* g_Res_vs_eff_wrtMiB2 = new TGraphAsymmErrors(); 
   TGraphAsymmErrors* g_Res_vs_eff_wrtRm2 = new TGraphAsymmErrors(); 

   if(doScan == true){
    std::vector<float> points_wrtMiB2;
    std::vector<float> points_wrtRm2;
    
    std::vector<TH1F*> resHist_wrtMiB2;
    resHist_wrtMiB2.resize(nBins);
    std::vector<TF1*> resFit_wrtMiB2;
    resFit_wrtMiB2.resize(nBins);
    std::vector<TF1*> resFitAlt_wrtMiB2;
    resFitAlt_wrtMiB2.resize(nBins);

    std::vector<TH1F*> resHist_wrtRm2;
    resHist_wrtRm2.resize(nBins);
    std::vector<TF1*> resFit_wrtRm2;
    resFit_wrtRm2.resize(nBins);
    std::vector<TF1*> resFitAlt_wrtRm2;
    resFitAlt_wrtRm2.resize(nBins);

    int iPoint_MiB2 = 0;
    int iPoint_Rm2 = 0;

    for(int ii = 0; ii < nBins; ii++)
    {
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_Res_1_%d_wrtMiB2",ii);
        resHist_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,400,-2.,2.);

        char Name_wrtRm2 [50];
        sprintf (Name_wrtRm2,"h_Res_1_%d_wrtRm2",ii);
        resHist_wrtRm2[ii] = new TH1F(Name_wrtRm2,Name_wrtRm2,400,-2.,2.);

        char cutMin [10];
        char cutMax [10];
        sprintf (cutMin,"%f",ampMin[ii]);
        sprintf (cutMax,"%f",ampMin[ii+1]);

        std::string Selection16 = Selection9+std::string(" && amp_max[")+iMCP+std::string("]>=")+std::string(cutMin)+std::string(" && amp_max[")+iMCP+std::string("]<")+std::string(cutMax);
        std::string Selection17 = Selection2+std::string(Name_wrtMiB2);
        std::string Selection18 = Selection4+std::string(Name_wrtRm2);

        //std::cout <<  Selection2 << " " << std::string(Name_wrtMiB2) <<  std::endl;
       
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];

        std::string wrtMCP = "";
        if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

        sprintf (NameOutput_wrtMiB2_png,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(), wrtMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(), wrtMCP.c_str());
        sprintf (NameOutput_wrtRm2_png,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
         
        if(iMCP != "MiB2"){
           h4->Draw(Selection17.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",resHist_wrtMiB2[ii]->GetMean()-3*resHist_wrtMiB2[ii]->GetRMS(),resHist_wrtMiB2[ii]->GetMean()+3*resHist_wrtMiB2[ii]->GetRMS());
           else resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"gaus",resHist_wrtMiB2[ii]->GetMean()-3*resHist_wrtMiB2[ii]->GetRMS(),resHist_wrtMiB2[ii]->GetMean()+3*resHist_wrtMiB2[ii]->GetRMS());
           resFitAlt_wrtMiB2[ii]->SetParameters(0,resHist_wrtMiB2[ii]->GetEntries()/2);
           resFitAlt_wrtMiB2[ii]->SetParameters(1,0.);
           resFitAlt_wrtMiB2[ii]->SetParameters(2,0.5);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(3,resHist_wrtMiB2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(1,0.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(5,0.5);
           resFitAlt_wrtMiB2[ii]->SetParLimits(0,0.,resHist_wrtMiB2[ii]->GetEntries());
           resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,1.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(3,0.,resHist_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(4,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(5,0.,1.);
           resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");

           float sigma_eff = 0.;
           float s_sigma = 0.;
           if(doDoubleGauss == true){
              float f1 = resFitAlt_wrtMiB2[ii]->GetParameter(0)/(resFitAlt_wrtMiB2[ii]->GetParameter(0)+resFitAlt_wrtMiB2[ii]->GetParameter(3));
              float f2 = resFitAlt_wrtMiB2[ii]->GetParameter(3)/(resFitAlt_wrtMiB2[ii]->GetParameter(0)+resFitAlt_wrtMiB2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_wrtMiB2[ii]->GetParameter(2)*resFitAlt_wrtMiB2[ii]->GetParameter(2)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(5)*resFitAlt_wrtMiB2[ii]->GetParameter(5) + f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4)*resFitAlt_wrtMiB2[ii]->GetParameter(4) - (f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4))*(f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4)));
              s_sigma = sqrt((f1*f1*resFitAlt_wrtMiB2[ii]->GetParameter(2)*resFitAlt_wrtMiB2[ii]->GetParameter(2)*resFitAlt_wrtMiB2[ii]->GetParError(2)*resFitAlt_wrtMiB2[ii]->GetParError(2)+f2*f2*resFitAlt_wrtMiB2[ii]->GetParameter(4)*resFitAlt_wrtMiB2[ii]->GetParameter(4)*resFitAlt_wrtMiB2[ii]->GetParError(4)*resFitAlt_wrtMiB2[ii]->GetParError(4)/(sigma_eff*sigma_eff)));
           }else{
              sigma_eff = resFitAlt_wrtMiB2[ii]->GetParameter(2);
              s_sigma = resFitAlt_wrtMiB2[ii]->GetParError(2);
           }

           if(sigma_eff*1000. < 10.) continue;
           if(iMCP != "MiB2" && iMCP != "Rm2") sigma_eff = sqrt(sigma_eff*sigma_eff-23.E-3*23.E-3);  
           g_Res_vs_Amp_wrtMiB2->SetPoint(iPoint_MiB2,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2,sigma_eff*1000.);
           g_Res_vs_Amp_wrtMiB2->SetPointError(iPoint_MiB2,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,s_sigma*1000.,s_sigma*1000.);

           iPoint_MiB2++;

           points_wrtMiB2.push_back(sigma_eff*1000.);
        
           TCanvas* c2 = new TCanvas();
           c2->cd();
           resHist_wrtMiB2[ii]->Draw("hist");
           resFitAlt_wrtMiB2[ii]->Draw("same");
           c2 -> Print(NameOutput_wrtMiB2_png,"png");
           c2 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c2;
        }

        if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
           h4->Draw(Selection18.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",resHist_wrtRm2[ii]->GetMean()-3.*resHist_wrtRm2[ii]->GetRMS(),resHist_wrtRm2[ii]->GetMean()+3.*resHist_wrtRm2[ii]->GetRMS());
           else resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_wrtRm2[ii]->GetMean()-3.*resHist_wrtRm2[ii]->GetRMS(),resHist_wrtRm2[ii]->GetMean()+3.*resHist_wrtRm2[ii]->GetRMS());
           resFitAlt_wrtRm2[ii]->SetParameters(0,resHist_wrtRm2[ii]->GetEntries()/2.);
           resFitAlt_wrtRm2[ii]->SetParameters(1,0.);
           resFitAlt_wrtRm2[ii]->SetParameters(2,0.05);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(3,resHist_wrtRm2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(1,0.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(5,0.05);
           resFitAlt_wrtRm2[ii]->SetParLimits(0,0.,resHist_wrtRm2[ii]->GetEntries());
           resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,0.5);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(3,0.,resHist_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(4,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(5,0.,0.5);
           resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");

           float sigma_eff = 0.;
           float s_sigma = 0.;
           if(doDoubleGauss == true){
              float f1 = resFitAlt_wrtRm2[ii]->GetParameter(0)/(resFitAlt_wrtRm2[ii]->GetParameter(0)+resFitAlt_wrtRm2[ii]->GetParameter(3));
              float f2 = resFitAlt_wrtRm2[ii]->GetParameter(3)/(resFitAlt_wrtRm2[ii]->GetParameter(0)+resFitAlt_wrtRm2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_wrtRm2[ii]->GetParameter(2)*resFitAlt_wrtRm2[ii]->GetParameter(2)+f2*resFitAlt_wrtRm2[ii]->GetParameter(5)*resFitAlt_wrtRm2[ii]->GetParameter(5) + f1*resFitAlt_wrtRm2[ii]->GetParameter(1)*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4)*resFitAlt_wrtRm2[ii]->GetParameter(4) - (f1*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4))*(f1*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4)));
              s_sigma = sqrt((f1*f1*resFitAlt_wrtRm2[ii]->GetParameter(2)*resFitAlt_wrtRm2[ii]->GetParameter(2)*resFitAlt_wrtRm2[ii]->GetParError(2)*resFitAlt_wrtRm2[ii]->GetParError(2)+f2*f2*resFitAlt_wrtRm2[ii]->GetParameter(4)*resFitAlt_wrtRm2[ii]->GetParameter(4)*resFitAlt_wrtRm2[ii]->GetParError(4)*resFitAlt_wrtRm2[ii]->GetParError(4)/(sigma_eff*sigma_eff)));
           }else{
              sigma_eff = resFitAlt_wrtRm2[ii]->GetParameter(2);    
              s_sigma = resFitAlt_wrtRm2[ii]->GetParError(2);    
           }

           if(sigma_eff*1000. < 10.) continue;
           if(iMCP != "MiB2" && iMCP != "Rm2") sigma_eff = sqrt(sigma_eff*sigma_eff-17.E-3*17.E-3);  
           g_Res_vs_Amp_wrtRm2->SetPoint(iPoint_Rm2,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2.,sigma_eff*1000.);
           g_Res_vs_Amp_wrtRm2->SetPointError(iPoint_Rm2,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2.,s_sigma*1000.,s_sigma*1000.);

           iPoint_Rm2++;

           points_wrtRm2.push_back(sigma_eff*1000.);
        
           TCanvas* c3 = new TCanvas();
           c3->cd();
           resHist_wrtRm2[ii]->Draw("hist");
           resFitAlt_wrtRm2[ii]->Draw("same");
           c3 -> Print(NameOutput_wrtRm2_png,"png");
           c3 -> Print(NameOutput_wrtRm2_pdf,"pdf");
           delete c3;
       }
    }

    if(points_wrtMiB2.size() != 0) std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    if(points_wrtRm2.size() != 0) std::sort(points_wrtRm2.begin(),points_wrtRm2.end());
    
    if(iMCP != "MiB2"){
       g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
       g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
       g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
       g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-5,points_wrtMiB2.at(points_wrtMiB2.size()-1)+5);

       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       TCanvas* c7 = new TCanvas();
       c7->cd();
       g_Res_vs_Amp_wrtMiB2->Draw("AP");
       c7 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
       c7 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
    }
    
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
       g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
       g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
       g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
       g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
       g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-5,points_wrtRm2.at(points_wrtRm2.size()-1)+5);

       TCanvas* c8 = new TCanvas();
       c8->cd();
       g_Res_vs_Amp_wrtRm2->Draw("AP");
       c8 -> Print(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c8 -> Print(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }          
    
   }

   std::string wrtMCP = "";
   if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";
   TFile* output_ampMax = new TFile(std::string("Data_TimeResolution_vs_amp_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".root").c_str(),"RECREATE");
   output_ampMax->cd();
   g_Res_vs_Amp_wrtMiB2->Write(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP).c_str());
   if(doOnlyWrtMiB2 == false) g_Res_vs_Amp_wrtRm2->Write(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP).c_str());
   output_ampMax->Close();

   int nBinsROC = 75;
   float stepROC = 10.;

   /*if(Timing == "LED50") nBinsROC = 60;
   if(Timing == "LED100") nBinsROC = 55;
   if(Timing == "LED150") nBinsROC = 55; */
   
   if(doScanEff == true){
    
    std::vector<float> points_ROC_wrtMiB2;
    std::vector<float> points_ROC_wrtRm2;
    
    std::vector<TH1F*> resHist_ROC_wrtMiB2;
    resHist_ROC_wrtMiB2.resize(nBinsROC);
    std::vector<TF1*> resFit_ROC_wrtMiB2;
    resFit_ROC_wrtMiB2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_wrtMiB2;
    resFitAlt_ROC_wrtMiB2.resize(nBinsROC);

    std::vector<TH1F*> resHist_ROC_wrtRm2;
    resHist_ROC_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFit_ROC_wrtRm2;
    resFit_ROC_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_wrtRm2;
    resFitAlt_ROC_wrtRm2.resize(nBinsROC);

    float res_min_wrtMiB2 = 99999.;
    float amp_max_min_wrtMiB2 = 0.;
    float efficiency_min_wrtMiB2 = 1.1;
    float efficiency_wrtMin_wrtMiB2 = 0.;
    float error_wrtMiB2 = 0.;
    float res_min_wrtRm2 = 99999.;
    float amp_max_min_wrtRm2 = 0.;
    float efficiency_min_wrtRm2 = 1.1;
    float efficiency_wrtMin_wrtRm2 = 0.;
    float error_wrtRm2 = 0.;

    int iPoint_MiB2 = 0;
    int iPoint_Rm2 = 0;

    wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

    for(int ii = 0; ii < nBinsROC; ii++)
    {
        if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";
        else wrtMCP = "";
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_Res_%d_ROC_wrtMiB2%s",ii,wrtMCP.c_str());
        resHist_ROC_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,400,-2.,2.);

        char Name_wrtRm2 [50];
        sprintf (Name_wrtRm2,"h_Res_%d_ROC_wrtRm2",ii);
        resHist_ROC_wrtRm2[ii] = new TH1F(Name_wrtRm2,Name_wrtRm2,400,-2.,2.);

        std::string::size_type sz;     // alias of size_t

       float thres = std::stof (thresMCP,&sz);

        char cutMin [10];
        char cutMax [10];
        sprintf (cutMin,"%f",ii*stepROC+thres);

        std::cout << "CutMin = " << cutMin << std::endl;
        
        std::string Selection16 = Selection9+std::string(" && amp_max[")+iMCP+std::string("]>=")+std::string(cutMin);
        std::string Selection17 = Selection2+std::string(Name_wrtMiB2);
        std::string Selection18 = Selection4+std::string(Name_wrtRm2);

        if(iMCP != "MiB2"){
           h4->Draw(Selection17.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",resHist_ROC_wrtMiB2[ii]->GetMean()-3*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetMean()+3*resHist_ROC_wrtMiB2[ii]->GetRMS());
           else resFitAlt_ROC_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"gaus",resHist_ROC_wrtMiB2[ii]->GetMean()-3*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetMean()+3*resHist_ROC_wrtMiB2[ii]->GetRMS());
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(0,resHist_ROC_wrtMiB2[ii]->GetEntries()/2);
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(1,0.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(2,0.5);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(3,resHist_ROC_wrtMiB2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(1,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(5,0.5);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(0,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(1,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(2,0.,1.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(3,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(4,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(5,0.,1.);
           resHist_ROC_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");

           float sigma_eff = 0.;
           float s_sigma = 0.;
           if(doDoubleGauss == true){
              float f1 = resFitAlt_ROC_wrtMiB2[ii]->GetParameter(0)/(resFitAlt_ROC_wrtMiB2[ii]->GetParameter(0)+resFitAlt_ROC_wrtMiB2[ii]->GetParameter(3));
              float f2 = resFitAlt_ROC_wrtMiB2[ii]->GetParameter(3)/(resFitAlt_ROC_wrtMiB2[ii]->GetParameter(0)+resFitAlt_ROC_wrtMiB2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2)+f2*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(5)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(5) + f1*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4) - (f1*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4))*(f1*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4)));
              s_sigma = sqrt((f1*f1*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2)*resFitAlt_ROC_wrtMiB2[ii]->GetParError(2)*resFitAlt_ROC_wrtMiB2[ii]->GetParError(2)+f2*f2*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4)*resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4)*resFitAlt_ROC_wrtMiB2[ii]->GetParError(4)*resFitAlt_ROC_wrtMiB2[ii]->GetParError(4)/(sigma_eff*sigma_eff)));
           }else{
              sigma_eff = resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2);
              s_sigma = resFitAlt_ROC_wrtMiB2[ii]->GetParError(2);
           }
           
           //std::cout << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0) << " " << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(1) << " " << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(2) << std::endl;
           
           if(1000.*sigma_eff < 10) continue;
           g_Res_vs_eff_wrtMiB2->SetPoint(iPoint_MiB2,ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0),1000.*sigma_eff);
           g_Res_vs_eff_wrtMiB2->SetPointError(iPoint_MiB2,ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(1),ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(2),s_sigma*1000.,s_sigma*1000.);
           iPoint_MiB2++;

           points_ROC_wrtMiB2.push_back(1000.*sigma_eff);

           if(1000.*sigma_eff<res_min_wrtMiB2){
              res_min_wrtMiB2 = 1000.*sigma_eff;
              amp_max_min_wrtMiB2 = ii*stepROC+thres;
              efficiency_wrtMin_wrtMiB2 = ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0);
              error_wrtMiB2 = s_sigma*1000.;
           } 
           
           if(ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0) < efficiency_min_wrtMiB2) efficiency_min_wrtMiB2 = ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0); 
        }

        if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
           h4->Draw(Selection18.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",resHist_ROC_wrtRm2[ii]->GetMean()-3*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetMean()+3*resHist_ROC_wrtRm2[ii]->GetRMS());
           else resFitAlt_ROC_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_ROC_wrtRm2[ii]->GetMean()-3*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetMean()+3*resHist_ROC_wrtRm2[ii]->GetRMS());
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(0,resHist_ROC_wrtRm2[ii]->GetEntries()/2);
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(1,0.);
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(2,0.5);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(3,resHist_ROC_wrtRm2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(1,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(5,0.5);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(0,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(1,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(2,0.,1.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(3,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(4,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(5,0.,1.);
           resHist_ROC_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");

           float sigma_eff = 0.;
           float s_sigma = 0.;
           if(doDoubleGauss == true){
              float f1 = resFitAlt_ROC_wrtRm2[ii]->GetParameter(0)/(resFitAlt_ROC_wrtRm2[ii]->GetParameter(0)+resFitAlt_ROC_wrtRm2[ii]->GetParameter(3));
              float f2 = resFitAlt_ROC_wrtRm2[ii]->GetParameter(3)/(resFitAlt_ROC_wrtRm2[ii]->GetParameter(0)+resFitAlt_ROC_wrtRm2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_ROC_wrtRm2[ii]->GetParameter(2)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(2)+f2*resFitAlt_ROC_wrtRm2[ii]->GetParameter(5)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(5) + f1*resFitAlt_ROC_wrtRm2[ii]->GetParameter(1)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4) - (f1*resFitAlt_ROC_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4))*(f1*resFitAlt_ROC_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4)));
              s_sigma = sqrt((f1*f1*resFitAlt_ROC_wrtRm2[ii]->GetParameter(2)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(2)*resFitAlt_ROC_wrtRm2[ii]->GetParError(2)*resFitAlt_ROC_wrtRm2[ii]->GetParError(2)+f2*f2*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4)*resFitAlt_ROC_wrtRm2[ii]->GetParameter(4)*resFitAlt_ROC_wrtRm2[ii]->GetParError(4)*resFitAlt_ROC_wrtRm2[ii]->GetParError(4)/(sigma_eff*sigma_eff)));
           }else{
              sigma_eff = resFitAlt_ROC_wrtRm2[ii]->GetParameter(2);
              s_sigma = resFitAlt_ROC_wrtRm2[ii]->GetParError(2);
           }
           
           //std::cout << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0) << " " << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(1) << " " << ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(2) << std::endl;
           if(1000.*sigma_eff < 10) continue; 
           g_Res_vs_eff_wrtRm2->SetPoint(iPoint_Rm2,ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0),1000.*sigma_eff);
           g_Res_vs_eff_wrtRm2->SetPointError(iPoint_Rm2,ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(1),ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(2),s_sigma*1000.,s_sigma*1000.);
           iPoint_Rm2++;

           points_ROC_wrtRm2.push_back(1000.*sigma_eff);

           if(1000.*sigma_eff<res_min_wrtRm2){
              res_min_wrtRm2 = 1000.*sigma_eff;
              amp_max_min_wrtRm2 = ii*stepROC+thres;
              efficiency_wrtMin_wrtRm2 = ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0);
              error_wrtRm2 = s_sigma*1000.;
           } 
           
           if(ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0) < efficiency_min_wrtRm2) efficiency_min_wrtRm2 = ComputeEfficiency(h4, inputs, iMCP, Selection16, Selection9).at(0);
        }
    }

    if(iMCP != "MiB2") std::sort(points_ROC_wrtMiB2.begin(),points_ROC_wrtMiB2.end());
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false) std::sort(points_ROC_wrtRm2.begin(),points_ROC_wrtRm2.end());

   char latex1_Min[500];
   sprintf (latex1_Min,"Minimum: ");
   char latex1_Amp[500];
   sprintf (latex1_Amp,"amp_max > %.0f ",amp_max_min_wrtMiB2);
   char latex1_Eff[500];
   sprintf (latex1_Eff,"#epsilon = %.2f ",efficiency_wrtMin_wrtMiB2);
   char latex1_Sigma[500];
   sprintf (latex1_Sigma,"#sigma_{t} = %.0f +/- %.0f ps",res_min_wrtMiB2, error_wrtMiB2);
   
   char latex2_Min[500];
   sprintf (latex2_Min,"Minimum: ");
   char latex2_Amp[500];
   sprintf (latex2_Amp,"amp_max > %.0f ",amp_max_min_wrtRm2);
   char latex2_Eff[500];
   sprintf (latex2_Eff,"#epsilon = %.2f ",efficiency_wrtMin_wrtRm2);
   char latex2_Sigma[500];
   sprintf (latex2_Sigma,"#sigma_{t} = %.0f +/- %.0f ps",res_min_wrtRm2, error_wrtRm2);
    
   TLatex *latexLabel1_Min = new TLatex();      
   latexLabel1_Min->SetTextSize(0.04);
   latexLabel1_Min->SetNDC();
   latexLabel1_Min->SetTextFont(42); // helvetica
   latexLabel1_Min->SetTextColor(kRed);

   TLatex *latexLabel1_Amp = new TLatex();      
   latexLabel1_Amp->SetTextSize(0.04);
   latexLabel1_Amp->SetNDC();
   latexLabel1_Amp->SetTextFont(42); // helvetica
   latexLabel1_Amp->SetTextColor(kRed);

   TLatex *latexLabel1_Eff = new TLatex();      
   latexLabel1_Eff->SetTextSize(0.04);
   latexLabel1_Eff->SetNDC();
   latexLabel1_Eff->SetTextFont(42); // helvetica
   latexLabel1_Eff->SetTextColor(kRed);

   TLatex *latexLabel1_Sigma = new TLatex();      
   latexLabel1_Sigma->SetTextSize(0.04);
   latexLabel1_Sigma->SetNDC();
   latexLabel1_Sigma->SetTextFont(42); // helvetica
   latexLabel1_Sigma->SetTextColor(kRed);
   
   TLatex *latexLabel2_Min = new TLatex();      
   latexLabel2_Min->SetTextSize(0.04);
   latexLabel2_Min->SetNDC();
   latexLabel2_Min->SetTextFont(42); // helvetica
   latexLabel2_Min->SetTextColor(kRed);

   TLatex *latexLabel2_Amp = new TLatex();      
   latexLabel2_Amp->SetTextSize(0.04);
   latexLabel2_Amp->SetNDC();
   latexLabel2_Amp->SetTextFont(42); // helvetica
   latexLabel2_Amp->SetTextColor(kRed);

   TLatex *latexLabel2_Eff = new TLatex();      
   latexLabel2_Eff->SetTextSize(0.04);
   latexLabel2_Eff->SetNDC();
   latexLabel2_Eff->SetTextFont(42); // helvetica
   latexLabel2_Eff->SetTextColor(kRed);

   TLatex *latexLabel2_Sigma = new TLatex();      
   latexLabel2_Sigma->SetTextSize(0.04);
   latexLabel2_Sigma->SetNDC();
   latexLabel2_Sigma->SetTextFont(42); // helvetica
   latexLabel2_Sigma->SetTextColor(kRed);
        
    if(iMCP != "MiB2"){
        g_Res_vs_eff_wrtMiB2->GetXaxis()->SetTitle("#epsilon");
        g_Res_vs_eff_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
        g_Res_vs_eff_wrtMiB2->SetMarkerStyle(20);
        g_Res_vs_eff_wrtMiB2->SetMarkerSize(0.7);
        g_Res_vs_eff_wrtMiB2->SetMarkerColor(kBlack);
        g_Res_vs_eff_wrtMiB2->SetLineColor(kBlack); 
        g_Res_vs_eff_wrtMiB2->GetXaxis()->SetRangeUser(efficiency_min_wrtMiB2-0.1,1.1);
        g_Res_vs_eff_wrtMiB2->GetYaxis()->SetRangeUser(points_ROC_wrtMiB2.at(0)-2.,points_ROC_wrtMiB2.at(points_ROC_wrtMiB2.size()-1)+2.);

        TCanvas* c9 = new TCanvas();
        c9->cd();
        c9->SetGrid();
        g_Res_vs_eff_wrtMiB2->Draw("AP");
        latexLabel1_Min->DrawLatex(0.35, 0.85,std::string(latex1_Min).c_str());
        latexLabel1_Amp->DrawLatex(0.35, 0.80,std::string(latex1_Amp).c_str());
        latexLabel1_Eff->DrawLatex(0.35, 0.75,std::string(latex1_Eff).c_str());
        latexLabel1_Sigma->DrawLatex(0.35, 0.70,std::string(latex1_Sigma).c_str());
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
    }

    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
        g_Res_vs_eff_wrtRm2->GetXaxis()->SetTitle("#epsilon");
        g_Res_vs_eff_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
        g_Res_vs_eff_wrtRm2->SetMarkerStyle(20);
        g_Res_vs_eff_wrtRm2->SetMarkerSize(0.7);
        g_Res_vs_eff_wrtRm2->SetMarkerColor(kBlack);
        g_Res_vs_eff_wrtRm2->SetLineColor(kBlack);
        g_Res_vs_eff_wrtRm2->GetXaxis()->SetRangeUser(efficiency_min_wrtRm2-0.1,1.1);
        g_Res_vs_eff_wrtRm2->GetYaxis()->SetRangeUser(points_ROC_wrtRm2.at(0)-2.,points_ROC_wrtRm2.at(points_ROC_wrtRm2.size()-1)+2.);

        TCanvas* c10 = new TCanvas();
        c10->cd();
        c10->SetGrid();
        g_Res_vs_eff_wrtRm2->Draw("AP");
        latexLabel2_Min->DrawLatex(0.35, 0.85,std::string(latex2_Min).c_str());
        latexLabel2_Amp->DrawLatex(0.35, 0.80,std::string(latex2_Amp).c_str());
        latexLabel2_Eff->DrawLatex(0.35, 0.75,std::string(latex2_Eff).c_str());
        latexLabel2_Sigma->DrawLatex(0.35, 0.70,std::string(latex2_Sigma).c_str());
        c10 -> Print(std::string("TimeResolution_vs_eff_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
        c10 -> Print(std::string("TimeResolution_vs_eff_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }    
   } 
}

void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan_Corr, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,ampMin,4000,timeMin);
    TH2F* timingCorrection_wrtRm2 = new TH2F("timingCorrection_wrtRm2","",nBins,ampMin,4000,timeMin);
    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
   
    std::string Selection1;   
    char Selection2 [1000];

    //if(iMCP == "M25") sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*amp_max[")+iMCP+std::string("]))<0.25")).c_str(),Params->at(0),Params->at(1));
    sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,"0.",false);
      Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection1,"0.",false);
      //Selection1 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection1,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
         //Selection1 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection1,"0.",false);
      }else{
         Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,"0.",false);
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
         //Selection1 = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection1,"0.",false);
      }   
    }
    std::string Selection3 = Selection1+" && "+Selection2+HodoSelection; 
    
    if(iMCP != "MiB2")
       h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtMiB2")).c_str(),Selection3.c_str());
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false)
       h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtRm2")).c_str(),Selection3.c_str());

    if(iMCP != "MiB2"){
       std::vector<float> points_wrtMiB2;
       timingCorrection_wrtMiB2->FitSlicesY();
       TH1F* timingCorrection_wrtMiB2_1 = (TH1F*)inputFile->Get("timingCorrection_wrtMiB2_1");
       timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time-time[MiB2]");

       for(int ii = 1; ii <= timingCorrection_wrtMiB2_1->GetNbinsX(); ii++)
           if(timingCorrection_wrtMiB2_1->GetBinContent(ii)!= 0) points_wrtMiB2.push_back(timingCorrection_wrtMiB2_1->GetBinContent(ii));
       std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());

       timingCorrection_wrtMiB2_1->SetAxisRange(points_wrtMiB2.at(0)-0.2,points_wrtMiB2.at(points_wrtMiB2.size()-1)+0.2, "Y");
       timingCorrection_wrtMiB2_1->SetMarkerStyle(20);
       timingCorrection_wrtMiB2_1->SetMarkerSize(0.9);
       timingCorrection_wrtMiB2_1->SetMarkerColor(kBlack);
       timingCorrection_wrtMiB2_1->SetLineColor(kBlack);
       
       TF1* fit_corr1;
       fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,3000.);
       timingCorrection_wrtMiB2_1->Fit("fit_corr1");
       if(doScan_Corr == false){
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
       }
    
       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       TCanvas* c1 = new TCanvas();
       c1->cd();
       timingCorrection_wrtMiB2_1->Draw();
       //timingCorrection_wrtMiB2->Draw();
       fit_corr1->Draw("same");
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
    }
 
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       std::vector<float> points_wrtRm2;
       timingCorrection_wrtRm2->FitSlicesY();
       TH1F* timingCorrection_wrtRm2_1 = (TH1F*)inputFile->Get("timingCorrection_wrtRm2_1");
       timingCorrection_wrtRm2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       timingCorrection_wrtRm2_1->GetYaxis()->SetTitle("time-time[Rm2]");

       for(int ii = 1; ii <= timingCorrection_wrtRm2_1->GetNbinsX(); ii++)
           if(timingCorrection_wrtRm2_1->GetBinContent(ii)!= 0) points_wrtRm2.push_back(timingCorrection_wrtRm2_1->GetBinContent(ii));
       std::sort(points_wrtRm2.begin(),points_wrtRm2.end());

       timingCorrection_wrtRm2_1->SetAxisRange(points_wrtRm2.at(0)-0.2,points_wrtRm2.at(points_wrtRm2.size()-1)+0.2, "Y");
       timingCorrection_wrtRm2_1->SetMarkerStyle(20);
       timingCorrection_wrtRm2_1->SetMarkerSize(0.9);
       timingCorrection_wrtRm2_1->SetMarkerColor(kBlack);
       timingCorrection_wrtRm2_1->SetLineColor(kBlack);
    
       TF1* fit_corr1;
       fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,3000.);
       timingCorrection_wrtRm2_1->Fit("fit_corr1");
       if(doScan_Corr == false){
         Params_wrtRm2->push_back(fit_corr1->GetParameter(0));
         Params_wrtRm2->push_back(fit_corr1->GetParameter(1));
         Params_wrtRm2->push_back(fit_corr1->GetParameter(2));
       }
    
       TCanvas* c2 = new TCanvas();
       c2->cd();
       timingCorrection_wrtRm2_1->Draw();
       //timingCorrection_wrtRm2->Draw();
       fit_corr1->Draw("same");
       c2 -> Print(std::string("timingCorrection_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c2 -> Print(std::string("timingCorrection_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }
      
    std::vector<TH1F*> resHist_wrtMiB2;
    resHist_wrtMiB2.resize(nBins);
    std::vector<TF1*> resFit_wrtMiB2;
    resFit_wrtMiB2.resize(nBins);
    std::vector<TF1*> resFitAlt_wrtMiB2;
    resFitAlt_wrtMiB2.resize(nBins);

    std::vector<TH1F*> resHist_wrtRm2;
    resHist_wrtRm2.resize(nBins);
    std::vector<TF1*> resFit_wrtRm2;
    resFit_wrtRm2.resize(nBins);
    std::vector<TF1*> resFitAlt_wrtRm2;
    resFitAlt_wrtRm2.resize(nBins);

    TGraphAsymmErrors* g_Res_vs_Amp_wrtMiB2 = new TGraphAsymmErrors(); 
    TGraphAsymmErrors* g_Res_vs_Amp_wrtRm2 = new TGraphAsymmErrors(); 

    std::vector<float> points_wrtMiB2;
    std::vector<float> points_wrtRm2;

    TCanvas* c3;
    TCanvas* c4;

    if(doScan_Corr == true){
      for(int ii = 0; ii < nBins; ii++)
      {
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_timingCorr_1_%d_wrtMiB2",ii);
        resHist_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,4000,-20.,20.);

        char Name_wrtRm2 [50];
        sprintf (Name_wrtRm2,"h_timingCorr_1_%d_wrtRm2",ii);
        resHist_wrtRm2[ii] = new TH1F(Name_wrtRm2,Name_wrtRm2,4000,-20.,20.);

        char cutMin [10];
        char cutMax [10];
        sprintf (cutMin,"%f",ampMin[ii]);
        sprintf (cutMax,"%f",ampMin[ii+1]);
        
        std::string Selection4 = Selection3+std::string(" && amp_max[")+iMCP+std::string("]>")+std::string(cutMin)+std::string(" && amp_max[")+iMCP+std::string("]<")+std::string(cutMax);
        std::string Selection5 = "time["+iMCP+iTiming+"]-time[MiB2] >> "+Name_wrtMiB2; 
        std::string Selection6 = "time["+iMCP+iTiming+"]-time[Rm2] >> "+Name_wrtRm2; 
        h4->Draw(Selection5.c_str(),Selection4.c_str());
        h4->Draw(Selection6.c_str(),Selection4.c_str());
      
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];

        std::string wrtMCP = "";
        if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

        sprintf (NameOutput_wrtMiB2_png,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(),wrtMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(),wrtMCP.c_str());
        sprintf (NameOutput_wrtRm2_png,"timingCorrection_%s_wrtRm2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtRm2_pdf,"timingCorrection_%s_wrtRm2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        
        if(iMCP != "MiB2"){
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"gaus",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.3,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.3);
           resFitAlt_wrtMiB2[ii]->SetParameters(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtMiB2[ii]->SetParameters(2,0.,2.);
           resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,2.);
           resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
   
           g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2,resFitAlt_wrtMiB2[ii]->GetParameter(1));
           g_Res_vs_Amp_wrtMiB2->SetPointError(ii,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,resFitAlt_wrtMiB2[ii]->GetParError(1),resFitAlt_wrtMiB2[ii]->GetParError(1));

           points_wrtMiB2.push_back(resFitAlt_wrtMiB2[ii]->GetParameter(1));
        
           c3 = new TCanvas();
           c3->cd();
           resHist_wrtMiB2[ii]->Draw("hist");
           resFitAlt_wrtMiB2[ii]->Draw("same");
           c3 -> Print(NameOutput_wrtMiB2_png,"png");
           c3 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c3;
        }

        if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.3,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.3);
           resFitAlt_wrtRm2[ii]->SetParameters(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtRm2[ii]->SetParameters(2,0.,2.);
           resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,2.);
           resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
        
           g_Res_vs_Amp_wrtRm2->SetPoint(ii,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2,resFitAlt_wrtRm2[ii]->GetParameter(1));
           g_Res_vs_Amp_wrtRm2->SetPointError(ii,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,resFitAlt_wrtRm2[ii]->GetParError(1),resFitAlt_wrtRm2[ii]->GetParError(1));

           points_wrtRm2.push_back(resFitAlt_wrtRm2[ii]->GetParameter(1));
        
           c4 = new TCanvas();
           c4->cd();
           resHist_wrtRm2[ii]->Draw("hist");
           resFitAlt_wrtRm2[ii]->Draw("same");
           c4 -> Print(NameOutput_wrtRm2_png,"png");
           c4 -> Print(NameOutput_wrtRm2_pdf,"pdf");
           delete c4;
        }
      }

      std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
      std::sort(points_wrtRm2.begin(),points_wrtRm2.end());
    
      if(iMCP != "MiB2"){
         g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
         g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#mu (ns)");
         g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
         g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
         g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
         g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
         g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-0.5,points_wrtMiB2.at(points_wrtMiB2.size()-1)+0.5);

         TF1* fit_corr3;
         if(iMCP == "M10") fit_corr3 = new TF1("fit_corr3","[0]+[1]*1/(x+[2])",0.,2500.);
         else fit_corr3 = new TF1("fit_corr3","pol1",0.,2500.);
         g_Res_vs_Amp_wrtMiB2->Fit("fit_corr3");
         Params_wrtMiB2->push_back(fit_corr3->GetParameter(0));
         Params_wrtMiB2->push_back(fit_corr3->GetParameter(1));
         if(iMCP == "M10")  Params_wrtMiB2->push_back(fit_corr3->GetParameter(2));

         std::string wrtMCP = "";
         if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

         TCanvas* c5 = new TCanvas();
         c5->cd();
         g_Res_vs_Amp_wrtMiB2->Draw("AP");
         fit_corr3->Draw("same");
         c5 -> Print(std::string("timingCorrection2_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
         c5 -> Print(std::string("timingCorrection2_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
      }
      
      if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
         g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
         g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#mu(ns)");
         g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
         g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
         g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
         g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
         g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-0.5,points_wrtRm2.at(points_wrtRm2.size()-1)+0.5);

         TF1* fit_corr4;
         if(iMCP == "M10") fit_corr4 = new TF1("fit_corr4","[0]+[1]*1/(x+[2])",0.,2500.);
         else fit_corr4 = new TF1("fit_corr4","pol1",0.,2500.);
         g_Res_vs_Amp_wrtRm2->Fit("fit_corr4");
         Params_wrtRm2->push_back(fit_corr4->GetParameter(0));
         Params_wrtRm2->push_back(fit_corr4->GetParameter(1));
         if(iMCP == "M10")  Params_wrtRm2->push_back(fit_corr4->GetParameter(2));

         TCanvas* c6 = new TCanvas();
         c6->cd();
         g_Res_vs_Amp_wrtRm2->Draw("AP");
         fit_corr4->Draw("same");
         c6 -> Print(std::string("timingCorrection2_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
         c6 -> Print(std::string("timingCorrection2_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
      }
    }

}

void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    TH2D* h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    TH2D* h2_time_maximum_vs_amp = new TH2D("h2_time_maximum_vs_amp","",300,0.,3000.,500,0.,5.);

    TF1* pol1_max = new TF1("pol1_max","pol1",0.,3000.);
    TF1* fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",0.,3000.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }   
    }

    Selection = Selection + HodoSelection;

    h4->Draw((std::string("time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]:amp_max[")+iMCP+std::string("] >> h2_time_max_vs_amp")).c_str(),Selection.c_str(),"goff");

    if(Timing == "CFD50"){
       fit_corr_max->SetParameters(2,0.,40.);
       fit_corr_max->SetParLimits(2,0.,40.);
    }else if(Timing == "LED50"){
       fit_corr_max->SetParameters(2,30.,70.);
       fit_corr_max->SetParLimits(2,30.,70.);
    }else if(Timing == "LED100"){
       fit_corr_max->SetParameters(2,80.,120.);
       fit_corr_max->SetParLimits(2,80.,120.);
    }else if(Timing == "LED150"){
       fit_corr_max->SetParameters(2,130.,170.);
       fit_corr_max->SetParLimits(2,130.,170.);
    }
    
    h2_time_max_vs_amp->GetYaxis()->SetTitle((std::string("time_max-time[")+iMCP+std::string("] (ns)")).c_str());
    h2_time_max_vs_amp->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("] (ns)")).c_str());
    //h2_time_max_vs_amp->SetAxisRange(0.,1750., "Z");
    //if(iMCP == "M25") h2_time_max_vs_amp->Fit("pol1_max");
    if(Timing == "CFD50") h2_time_max_vs_amp->Fit("fit_corr_max","B","",20.,3000.);
    else if(Timing == "LED50") h2_time_max_vs_amp->Fit("fit_corr_max","B","",50.,3000.);
    else if(Timing == "LED100") h2_time_max_vs_amp->Fit("fit_corr_max","B","",100.,3000.);
    else if(Timing == "LED150") h2_time_max_vs_amp->Fit("fit_corr_max","B","",150.,3000.);

    Params->push_back(fit_corr_max->GetParameter(0));
    Params->push_back(fit_corr_max->GetParameter(1));
    Params->push_back(fit_corr_max->GetParameter(2));
    
    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyMIB2";

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h2_time_max_vs_amp->Draw("COLZ");
    if(iMCP == "M25") fit_corr_max->Draw("same");
    else fit_corr_max->Draw("same");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
}

void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",300,-10,20,300,-1.,1.5,100.,3000.);
    TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time","",300,-10,20,300,-1.,1.5);

    std::string Selection;
    std::string iTiming = "";
  
    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }   
    }

    Selection = Selection+" && WF_ch == "+iMCP;//+HodoSelection;

    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time")).c_str(),Selection.c_str(),"goff");
    h4->Draw((std::string("WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> h2_amp_vs_time")).c_str(),Selection.c_str());

    TProfile* waveForm = h2_amp_vs_time->ProfileX(); 
    waveForm->SetName(std::string(nameiMCP+std::string("_waveform_prof")).c_str());
    
    p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    waveForm->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
    h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");

    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyMIB2";
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    h2_amp_vs_time->Draw("COLZ");
    c1 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+"_h2.png").c_str(),"png");
    c1 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+"_h2.pdf").c_str(),"pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    p2D_amp_vs_time->Draw("COLZ");
    c2 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c2 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");

    TCanvas* c3 = new TCanvas();
    c3->cd();
    waveForm->Draw("P");
    c3 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+"_profile.png").c_str(),"png");
    c3 -> Print(std::string("pulseShape_"+nameiMCP+wrtMCP+"_profile.pdf").c_str(),"pdf");

    TFile* output_Waveform = new TFile(std::string(nameiMCP+std::string("_Waveform.root")).c_str(),"RECREATE");
    output_Waveform->cd();
    waveForm->Write(std::string(nameiMCP+wrtMCP+"_waveform_prof").c_str());
    output_Waveform->Close();
}

void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2)
{
    TH2F* h2_hodoscope_Y_vs_X_noSelection = new TH2F("h2_hodoscope_Y_vs_X_noSelection","",30,0.,30.,30,0.,30.);
    TH2F* h2_hodoscope_Y_vs_X = new TH2F("h2_hodoscope_Y_vs_X","",30,0.,30.,30,0.,30.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }   
    }

    h4->Draw("Y:X >> h2_hodoscope_Y_vs_X_noSelection");
    h4->Draw("Y:X >> h2_hodoscope_Y_vs_X",Selection.c_str());

    h2_hodoscope_Y_vs_X_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X->GetYaxis()->SetTitle("Y");

    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMIB2";

    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetGrid();
    h2_hodoscope_Y_vs_X_noSelection->Draw("COLZ");
    c1 -> Print("Hodoscope_Y_vs_X_noSelection.png","png");
    c1 -> Print("Hodoscope_Y_vs_X_noSelection.pdf","pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetGrid();
    h2_hodoscope_Y_vs_X->Draw("COLZ");
    c2 -> Print(std::string("Hodoscope_Y_vs_X_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c2 -> Print(std::string("Hodoscope_Y_vs_X_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");
   
}

void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2)
{
    TH1F* h_timeChi2_noSelection = new TH1F("h_timeChi2_noSelection","",1500,1.,15000.);
    TH1F* h_timeChi2 = new TH1F("h_timeChi2","",1500,1.,15000.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }   
    } 

    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2_noSelection").c_str());
    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2").c_str(),Selection.c_str());
    
    h_timeChi2_noSelection->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2_noSelection->GetYaxis()->SetTitle("Events");
    h_timeChi2->GetYaxis()->SetTitle("Events");

    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMIB2";
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogx();
    c1->SetLogy();
    h_timeChi2_noSelection->Draw("H");
    c1 -> Print(std::string("timeChi2_noSelection_"+nameiMCP+".png").c_str(),"png");
    c1 -> Print(std::string("timeChi2_noSelection_"+nameiMCP+".pdf").c_str(),"pdf");
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogx();
    c2->SetLogy();
    h_timeChi2->Draw("H");
    c2 -> Print(std::string("timeChi2_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c2 -> Print(std::string("timeChi2_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");
}   

void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2)
{
    TH1F* h_ampMax_noSelection = new TH1F(std::string("h_ampMax_"+nameiMCP+"_noSelection").c_str(),"",450,0.,4500.);
    TH1F* h_ampMax = new TH1F(std::string("h_ampMax_"+nameiMCP).c_str(),"",450,0.,4500.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
      }   
    }

    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"_noSelection").c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP).c_str(),Selection.c_str());
    
    h_ampMax_noSelection->GetXaxis()->SetTitle("amp_max");
    h_ampMax->GetXaxis()->SetTitle("amp_max");
    h_ampMax_noSelection->GetYaxis()->SetTitle("Events");
    h_ampMax->GetYaxis()->SetTitle("Events");
 
    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMIB2";
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogy();
    h_ampMax_noSelection->Draw("H");
    c1 -> Print(std::string("ampMax_noSelection_"+nameiMCP+".png").c_str(),"png");
    c1 -> Print(std::string("ampMax_noSelection_"+nameiMCP+".pdf").c_str(),"pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    h_ampMax->Draw("H");
    c2 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c2 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");
    
}   

void CheckEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TH1F* num = new TH1F("num","",40,0.,4000.);
    TH1F* den = new TH1F("den","",40,0.,4000.);

    std::string HV = "";
    if(iMCP == "M25") HV = "HV25";
    else if(iMCP == "M10") HV = "HV10";
    else if(iMCP == "M8") HV = "HV8";
    else if(iMCP == "M5") HV = "HV5";
    else if(iMCP == "BINP1" && inputs.find("2378") != std::string::npos) HV = "HVBINP1NEG";
    else if(iMCP == "BINP1" && inputs.find("2378") == std::string::npos) HV = "HVBINP1";
    else if(iMCP == "BINP2" && inputs.find("2547") != std::string::npos) HV = "HVBINP2NEG";
    else if(iMCP == "BINP2" && inputs.find("2547") == std::string::npos) HV = "HVBINP2";
    else if(iMCP == "BINP3") HV = "HVBINP3";
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") != std::string::npos) HV = "HVBINP4NEG";
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") == std::string::npos) HV = "HVBINP4";


    h4->Draw(std::string(HV+" >> num").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[MiB2]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax).c_str());     

    TGraphAsymmErrors* eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    eff->Draw("AP");
    c1 -> Print(std::string("Efficiency_MiB2.png").c_str(),"png");
    c1 -> Print(std::string("Efficiency_MiB2.pdf").c_str(),"pdf");

    h4->Draw(std::string(HV+" >> num").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[Rm2]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax).c_str());     

    eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    eff->Draw("AP");
    c2 -> Print(std::string("Efficiency_Rm2.png").c_str(),"png");
    c2 -> Print(std::string("Efficiency_Rm2.pdf").c_str(),"pdf");

    h4->Draw(std::string(HV+" >> num").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[Rm2]>200. && amp_max[Rm2]<1200. && amp_max[MiB2]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[Rm2]>200. && amp_max[Rm2]<1200.").c_str());   

    eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c3 = new TCanvas();
    c3->cd();
    eff->Draw("AP");
    c3 -> Print(std::string("Efficiency_MiB2_wrtRm2.png").c_str(),"png");
    c3 -> Print(std::string("Efficiency_MiB2_wrtRm2.pdf").c_str(),"pdf");

    h4->Draw(std::string(HV+" >> num").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[MiB2]>200. && amp_max[Rm2]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[MiB2]>200.").c_str());   

    eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c4 = new TCanvas();
    c4->cd();
    eff->Draw("AP");
    c4 -> Print(std::string("Efficiency_Rm2_wrtMiB2.png").c_str(),"png");
    c4 -> Print(std::string("Efficiency_Rm2_wrtMiB2.pdf").c_str(),"pdf");


    h4->Draw(std::string(HV+" >> num").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[MiB2]>200. && amp_max[Rm2]>200. && amp_max[Rm2]<1200. && amp_max["+iMCP+"]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),std::string("adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && amp_max[MiB2]>200. && amp_max[Rm2]>200. && amp_max[Rm2]<1200.").c_str());   

    eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c5 = new TCanvas();
    c5->cd();
    eff->Draw("AP");
    c5 -> Print(std::string("Efficiency_"+iMCP+".png").c_str(),"png");
    c5 -> Print(std::string("Efficiency_"+iMCP+".pdf").c_str(),"pdf");
}   

void HodoValidation(TTree* h4, TTree* hodo, TTree* adc, int nParticle)
{
    std::vector<TH1F*> fiberEnergyX;
    fiberEnergyX.resize(32);
    std::vector<TH1F*> fiberEnergyY;
    fiberEnergyY.resize(32);

    std::vector<TH1F*> fiberEnergyX_g1000;
    fiberEnergyX_g1000.resize(32);
    std::vector<TH1F*> fiberEnergyY_g1000;
    fiberEnergyY_g1000.resize(32);

    std::string particleSel = "";
    if(nParticle == -1) particleSel = "adc_data[scint]<50."; 
    if(nParticle == 0) particleSel = "adc_data[scint]<150."; 
    if(nParticle == 1) particleSel = "adc_data[scint]>200. && adc_data[scint]<700."; 
    if(nParticle > 1) particleSel = "adc_data[scint]>200."; 
   
    char nPar[50];
    sprintf (nPar,"%d",nParticle);
   
    for(int iFib=0; iFib<32; iFib++)
    {
        char iFiber[50];
        sprintf (iFiber,"%d",iFib);

        if(nParticle==0) fiberEnergyX[iFib] = new TH1F(std::string("h_fiberEnergyX_"+std::string(iFiber)).c_str(),"",250,0.,2500.);
        else fiberEnergyX[iFib] = new TH1F(std::string("h_fiberEnergyX_"+std::string(iFiber)).c_str(),"",400,0.,4000.);
        h4->Draw(std::string("fiberEnergyX["+std::string(iFiber)+"] >> h_fiberEnergyX_"+std::string(iFiber)).c_str(),particleSel.c_str());
        if(nParticle==0) fiberEnergyX_g1000[iFib] = new TH1F(std::string("h_fiberEnergyX_"+std::string(iFiber)+"_g1000").c_str(),"",250,0.,2500.);
        else fiberEnergyX_g1000[iFib] = new TH1F(std::string("h_fiberEnergyX_"+std::string(iFiber)+"_g1000").c_str(),"",400,0.,4000.);
        h4->Draw(std::string("fiberEnergyX["+std::string(iFiber)+"] >> h_fiberEnergyX_"+std::string(iFiber)+"_g1000").c_str(),std::string(particleSel+" && fiberEnergyX["+std::string(iFiber)+"] >1150.").c_str());

        float iFrac = float(fiberEnergyX_g1000[iFib]->GetEntries())/float(fiberEnergyX[iFib]->GetEntries());
        char frac[100];
        TLatex *latexLabel1 = new TLatex();      
        sprintf (frac,"fraction > 1150 = %.1f",iFrac*100);
        latexLabel1->SetTextSize(0.05);
        latexLabel1->SetNDC();
        latexLabel1->SetTextFont(42); // helvetica
        
        TCanvas* c1 = new TCanvas();
        c1->cd();
        c1->SetLogy();
        fiberEnergyX[iFib]->Draw();
        latexLabel1->DrawLatex(0.56, 0.55,std::string(std::string(frac)+"%").c_str());
        c1 -> Print(std::string("fiberEnergyX_"+std::string(iFiber)+"_nPC"+std::string(nPar)+".png").c_str(),"png");
        c1 -> Print(std::string("fiberEnergyX_"+std::string(iFiber)+"_nPC_"+std::string(nPar)+".pdf").c_str(),"pdf");
         
        if(nParticle==0) fiberEnergyY[iFib] = new TH1F(std::string("h_fiberEnergyY_"+std::string(iFiber)).c_str(),"",250,0.,2500.);
        else fiberEnergyY[iFib] = new TH1F(std::string("h_fiberEnergyY_"+std::string(iFiber)).c_str(),"",400,0.,4000.);
        h4->Draw(std::string("fiberEnergyY["+std::string(iFiber)+"] >> h_fiberEnergyY_"+std::string(iFiber)).c_str(),particleSel.c_str());
        if(nParticle==0) fiberEnergyY_g1000[iFib] = new TH1F(std::string("h_fiberEnergyY_"+std::string(iFiber)+"_g1000").c_str(),"",250,0.,2500.);
        else fiberEnergyY_g1000[iFib] = new TH1F(std::string("h_fiberEnergyY_"+std::string(iFiber)+"_g1000").c_str(),"",400,0.,4000.);
        h4->Draw(std::string("fiberEnergyY["+std::string(iFiber)+"] >> h_fiberEnergyY_"+std::string(iFiber)+"_g1000").c_str(),std::string(particleSel+" && fiberEnergyY["+std::string(iFiber)+"] >1150.").c_str());

        iFrac = float(fiberEnergyY_g1000[iFib]->GetEntries())/float(fiberEnergyY[iFib]->GetEntries());
        TLatex *latexLabel2 = new TLatex();      
        sprintf (frac,"fraction > 1150 = %.1f",iFrac*100);
        latexLabel2->SetTextSize(0.05);
        latexLabel2->SetNDC();
        latexLabel2->SetTextFont(42); // helvetica
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        c2->SetLogy();
        fiberEnergyY[iFib]->Draw();
        latexLabel2->DrawLatex(0.56, 0.55,std::string(std::string(frac)+"%").c_str());
        c2 -> Print(std::string("fiberEnergyY_"+std::string(iFiber)+"_nPC_"+std::string(nPar)+".png").c_str(),"png");
        c2 -> Print(std::string("fiberEnergyY_"+std::string(iFiber)+"_nPC_"+std::string(nPar)+".pdf").c_str(),"pdf");

        delete c1;
        delete c2;
        delete latexLabel1;
        delete latexLabel2;
    }

    TH1F* fiberEnergyX_total = new TH1F("h_fiberEnergyX","",400,0.,4000.);
    TH1F* fiberEnergyY_total = new TH1F("h_fiberEnergyY","",400,0.,4000.);
    TH1F* fiberEnergyX_g1000_total = new TH1F("h_fiberEnergyX_g1000","",400,0.,4000.);
    TH1F* fiberEnergyY_g1000_total = new TH1F("h_fiberEnergyY_g1000","",400,0.,4000.);

    h4->Draw("fiberEnergyX >> h_fiberEnergyX",particleSel.c_str());
    h4->Draw("fiberEnergyX >> h_fiberEnergyX_g1000",std::string(particleSel+" && (fiberEnergyX[0] > 1150.  || fiberEnergyX[1] > 1150.  || fiberEnergyX[2] > 1150.  || fiberEnergyX[3] > 1150.  || fiberEnergyX[4] > 1150.  || fiberEnergyX[5] > 1150.  || fiberEnergyX[6] > 1150.  || fiberEnergyX[7] > 1150.  || fiberEnergyX[8] > 1150.  || fiberEnergyX[9] > 1150.  || fiberEnergyX[10] > 1150.  ||  fiberEnergyX[11] > 1150.  || fiberEnergyX[12] > 1150.  || fiberEnergyX[13] > 1150.  || fiberEnergyX[14] > 1150.  || fiberEnergyX[15] > 1150.  || fiberEnergyX[16] > 1150.  || fiberEnergyX[17] > 1150.  || fiberEnergyX[18] > 1150.  || fiberEnergyX[19] > 1150.  || fiberEnergyX[20] > 1150.  || fiberEnergyX[21] > 1150.  || fiberEnergyX[22] > 1150.  || fiberEnergyX[23] > 1150.  || fiberEnergyX[24] > 1150.  || fiberEnergyX[25] > 1150.  || fiberEnergyX[26] > 1150.  || fiberEnergyX[27] > 1150.  || fiberEnergyX[28] > 1150.  || fiberEnergyX[29] > 1150.  || fiberEnergyX[30] > 1150.  || fiberEnergyX[31] > 1150.)").c_str());

    float frac_total_X = float(fiberEnergyX_g1000_total->GetEntries())/float(fiberEnergyX_total->GetEntries());

    h4->Draw("fiberEnergyY >> h_fiberEnergyY",particleSel.c_str());
    h4->Draw("fiberEnergyY >> h_fiberEnergyY_g1000",std::string(particleSel+" && (fiberEnergyY[0] > 1150.  || fiberEnergyY[1] > 1150.  || fiberEnergyY[2] > 1150.  || fiberEnergyY[3] > 1150.  || fiberEnergyY[4] > 1150.  || fiberEnergyY[5] > 1150.  || fiberEnergyY[6] > 1150.  || fiberEnergyY[7] > 1150.  || fiberEnergyY[8] > 1150.  || fiberEnergyY[9] > 1150.  || fiberEnergyY[10] > 1150.  ||  fiberEnergyY[11] > 1150.  || fiberEnergyY[12] > 1150.  || fiberEnergyY[13] > 1150.  || fiberEnergyY[14] > 1150.  || fiberEnergyY[15] > 1150.  || fiberEnergyY[16] > 1150.  || fiberEnergyY[17] > 1150.  || fiberEnergyY[18] > 1150.  || fiberEnergyY[19] > 1150.  || fiberEnergyY[20] > 1150.  || fiberEnergyY[21] > 1150.  || fiberEnergyY[22] > 1150.  || fiberEnergyY[23] > 1150.  || fiberEnergyY[24] > 1150.  || fiberEnergyY[25] > 1150.  || fiberEnergyY[26] > 1150.  || fiberEnergyY[27] > 1150.  || fiberEnergyY[28] > 1150.  || fiberEnergyY[29] > 1150.  || fiberEnergyY[30] > 1150.  || fiberEnergyY[31] > 1150.)").c_str());

    float frac_total_Y = float(fiberEnergyY_g1000_total->GetEntries())/float(fiberEnergyY_total->GetEntries());

    std::cout << "Inclusive fraction X = " << frac_total_X*100. << "% " << std::endl;
    std::cout << "Inclusive fraction Y = " << frac_total_Y*100. << "% " << std::endl;

    TH1F* fiberEnergyMaxX = new TH1F("fiberEnergyMaxX","",400,0.,4000.);
    TH1F* fiberEnergyMaxY = new TH1F("fiberEnergyMaxY","",400,0.,4000.);
     
    //loop over the events of hodo tree
    ULong64_t     index;
    int           n_planes;
    int           n_hitsX;
    int           n_hitsY;
    float         X[1];   //[n_planes]
    float         Y[1];   //[n_planes]
    float         clusterX[1];   //[n_planes]
    float         clusterY[1];   //[n_planes]
    float         v_fiberEnergyX[32];
    float         v_fiberEnergyY[32];
    vector<float>   *adc_data;
    int           scint;

    TBranch        *b_index;   //!
    TBranch        *b_n_planes;   //!
    TBranch        *b_n_hitsX;   //!
    TBranch        *b_n_hitsY;   //!
    TBranch        *b_X;   //!
    TBranch        *b_Y;   //!
    TBranch        *b_clusterX;   //!
    TBranch        *b_clusterY;   //!
    TBranch        *b_fiberEnergyX;   //!
    TBranch        *b_fiberEnergyY;   //!
    TBranch        *b_adc_data;   //!
    TBranch        *b_scint;   //!

    hodo->SetBranchAddress("index", &index, &b_index);
    hodo->SetBranchAddress("n_planes", &n_planes, &b_n_planes);
    hodo->SetBranchAddress("n_hitsX", &n_hitsX, &b_n_hitsX);
    hodo->SetBranchAddress("n_hitsY", &n_hitsY, &b_n_hitsY);
    hodo->SetBranchAddress("X", X, &b_X);
    hodo->SetBranchAddress("Y", Y, &b_Y);
    hodo->SetBranchAddress("clusterX", clusterX, &b_clusterX);
    hodo->SetBranchAddress("clusterY", clusterY, &b_clusterY);
    hodo->SetBranchAddress("fiberEnergyX", v_fiberEnergyX, &b_fiberEnergyX);
    hodo->SetBranchAddress("fiberEnergyY", v_fiberEnergyY, &b_fiberEnergyY);
    adc->SetBranchAddress("adc_data", &adc_data, &b_adc_data);
    adc->SetBranchAddress("scint", &scint, &b_scint);

    for(int entry = 0; entry < hodo->GetEntries(); entry++){

        if(entry%1000==0) std::cout<<"--- Reading entry = "<<entry<<std::endl;
        hodo->GetEntry(entry);
        adc->GetEntry(entry);
        
        //std::cout << entry << " " << adc_data->at(scint) << std::endl;
        if(nParticle == -1 && adc_data->at(scint) > 50.) continue;
        if(nParticle == 0 && adc_data->at(scint) > 150.) continue;
        if(nParticle == 1 && (adc_data->at(scint) < 200. || adc_data->at(scint)>700.)) continue;
        if(nParticle > 1 && adc_data->at(scint) < 200.) continue;

        float energyMaxX = 0.;
        float energyMaxY = 0.;
        for(int iFib=0; iFib<32; iFib++)
        {
            if(v_fiberEnergyX[iFib]> energyMaxX) energyMaxX = v_fiberEnergyX[iFib];
            if(v_fiberEnergyY[iFib]> energyMaxY) energyMaxY = v_fiberEnergyY[iFib];
        }
        fiberEnergyMaxX->Fill(energyMaxX);
        fiberEnergyMaxY->Fill(energyMaxY);
    }

    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogy();
    fiberEnergyMaxX->Draw();
    c1 -> Print("fiberEnergyMaxX.png","png");
    c1 -> Print("fiberEnergyMaxX.pdf","pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    fiberEnergyMaxY->Draw();
    c2 -> Print("fiberEnergyMaxY.png","png");
    c2 -> Print("fiberEnergyMaxY.pdf","pdf");
}   

std::string AddSelection(TTree* h4, std::string Var, std::string Selection, std::string Cut = "0.", bool isCut = false)
{
    TH1F* h = new TH1F("h","h",4000,-20.,20.);
    TF1* g_fit = new TF1("g_fit","gaus",-20.,20.);
    h4->Draw((Var+std::string(" >> h")).c_str(),Selection.c_str());
    
    if(isCut == false){
      h->Fit("g_fit","","",h->GetMean()-3*h->GetRMS(),h->GetMean()+3*h->GetRMS());
      char Mean [100];
      sprintf(Mean,"%f",h->GetMean());
      std::string sMean = std::string(Mean);
      char Sigma [100];
      sprintf(Sigma,"%f",5*g_fit->GetParameter(2));
      std::string sSigma = std::string(Sigma);

      //std::cout << "Sigma = " << 5*g_fit->GetParameter(2) << std::endl;

      if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+std::string(" && fabs(")+Var+std::string("+")+std::string(sMean)+std::string(")<")+sSigma;
      }else{
       Selection = Selection+std::string(" && fabs(")+Var+std::string("-")+std::string(sMean)+std::string(")<")+sSigma;
      }
    }else{
      char Mean [100];
      sprintf(Mean,"%f",h->GetMean());
      std::string sMean = std::string(Mean);

      if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+std::string(" && fabs(")+Var+std::string("+")+std::string(sMean)+std::string(")<")+Cut;
      }else{
       Selection = Selection+std::string(" && fabs(")+Var+std::string("-")+std::string(sMean)+std::string(")<")+Cut;
      }
    }

    delete h;
    delete g_fit;
    return Selection;
}

std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel)
{
    TH1F* num = new TH1F("num","",40,0.,4000.);
    TH1F* den = new TH1F("den","",40,0.,4000.);

    std::string HV = "";
    if(iMCP == "M25") HV = "HV25";
    else if(iMCP == "M10") HV = "HV10";
    else if(iMCP == "M8") HV = "HV8";
    else if(iMCP == "M5") HV = "HV5";
    else if(iMCP == "BINP1" && inputs.find("2378") != std::string::npos) HV = "HVBINP1NEG";
    else if(iMCP == "BINP1" && inputs.find("2378") == std::string::npos) HV = "HVBINP1";
    else if(iMCP == "BINP2" && inputs.find("2547") != std::string::npos) HV = "HVBINP2NEG";
    else if(iMCP == "BINP2" && inputs.find("2547") == std::string::npos) HV = "HVBINP2";
    else if(iMCP == "BINP3") HV = "HVBINP3";
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") != std::string::npos) HV = "HVBINP4NEG";
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") == std::string::npos) HV = "HVBINP4";


    h4->Draw(std::string(HV+" >> num").c_str(),numSel.c_str());   
    h4->Draw(std::string(HV+" >> den").c_str(),denSel.c_str());     

    TGraphAsymmErrors* efficiency = new TGraphAsymmErrors(num,den);
    double hv;
    double eff;
    efficiency->GetPoint(0,hv,eff);
    float effErrorUp = efficiency->GetErrorYhigh(0);
    float effErrorDown = efficiency->GetErrorYlow(0);

   std::vector<float> finalEff;
   finalEff.push_back(float(eff)); 
   finalEff.push_back(effErrorDown); 
   finalEff.push_back(effErrorUp);

   delete num;
   delete den; 

   return finalEff;
}   


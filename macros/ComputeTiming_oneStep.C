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


int nBins = 34;
float ampMin[35] = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650, 700., 750., 800., 850., 900., 950., 1000., 1050., 1100., 1150., 1200., 1250, 1300., 1350., 1400., 1450, 1500., 1700., 2000.,2500.,3000.};
//float ampMin[20] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1700., 2000.,2500.,3000.};
float timeMin[4001];

//int nBins = 8;
//float stepSize = 250.; 
//float ampMin[9] = {0., 200., 400., 600., 800., 1000., 1500., 2000.,3000.};
//float timeMin[4001];

//MiB
//std::string amp_max_MiB2 = "200"
//std::string amp_max_Rm2 = "200"
//std::string b_rms = "4"

//BINP
std::string amp_max_MiB2 = "200.";
std::string time_max_MiB2 = "150.";
std::string time_max_Rm2 = "150.";
std::string amp_max_Rm2 = "200.";
std::string scintMin = "200.";
std::string scintMax = "700.";
std::string timeChi2 = "50";
bool doScan_Corr = false;


void FinalTiming(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, bool doScan, bool doDoubleGauss, std::string HodoSelection);
void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, bool doScan_Corr, std::string HodoSelection);
void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string HodoSelection);
void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string HodoSelection);
void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP);
void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP);
void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP);
std::string AddSelection(TTree*, std::string, std::string, std::string, bool);

void ComputeTiming_oneStep(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, bool doFirstStep = true, bool doPulseShapes = false, bool doScan = false, bool doDoubleGauss = false, std::string hodoXmin="0", std::string hodoXmax="30", std::string hodoYmin="0", std::string hodoYmax="30")
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    TTree* h4 = (TTree*)inputFile->Get("h4");
    
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
    
    if(doFirstStep == true){
       Hodoscope(h4, iMCP, nameiMCP, thresMCP);
       TimeChi2(h4, iMCP, nameiMCP, thresMCP);
       AmpMax(h4, iMCP, nameiMCP, thresMCP);
    }else if(doPulseShapes == true) PulseShapes(h4, iMCP, nameiMCP, thresMCP, HodoSelection);
    else{
       AmpVsTime_Selection(h4, iMCP, nameiMCP,Timing, Params, thresMCP,HodoSelection);
       TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP,doScan_Corr, HodoSelection);
       FinalTiming(h4, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, doScan, doDoubleGauss, HodoSelection);
    }
}

void FinalTiming(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, bool doScan, bool doDoubleGauss, std::string HodoSelection)
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
    if(iMCP != "Rm2" && iMCP != "BINP3"){
       sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
    }

    //no time correction
    std::string Selection6 = "time["+iMCP+iTiming+"]-time[MiB2] >>";
    std::string Selection7 = "time["+iMCP+iTiming+"]-time[Rm2] >>";
    std::string Selection8;

    if(iMCP == "MiB2"){
      Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection8 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection8,"0.",false);
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
      Selection8 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
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
       h4->Draw(Selection10.c_str(),Selection9.c_str()); 
       h4->Draw(Selection11.c_str(),Selection9.c_str()); 
       h4->Draw(Selection14.c_str(),Selection9.c_str()); 
    }
    if(iMCP != "Rm2" && iMCP != "BINP3"){
      h4->Draw(Selection12.c_str(),Selection9.c_str());
       h4->Draw(Selection13.c_str(),Selection9.c_str()); 
       h4->Draw(Selection15.c_str(),Selection9.c_str()); 
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
       if(doDoubleGauss == true){
          float f1 = g_res->GetParameter(0)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float f2 = g_res->GetParameter(3)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float sigma_eff = sqrt(f1*g_res->GetParameter(2)*g_res->GetParameter(2)+f2*g_res->GetParameter(5)*g_res->GetParameter(5) + f1*g_res->GetParameter(1)*g_res->GetParameter(1)+f2*g_res->GetParameter(4)*g_res->GetParameter(4) - (f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4))*(f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4)));
       
          sprintf (Sigma,"#sigma = %.0f ps",sigma_eff*1000);

          latexLabel->SetTextSize(0.05);
          latexLabel->SetNDC();
          latexLabel->SetTextFont(42); // helvetica
       }

       time_wrtMiB2->GetXaxis()->SetTitle("t-t_{ref} (ns)");
  
       TCanvas* c1 = new TCanvas();
       c1->cd();
       time_wrtMiB2->Draw("hist");
       if(doDoubleGauss == true) latexLabel->DrawLatex(0.72, 0.55,Sigma);
       g_res->Draw("same");
       c1 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
       delete g_res; 
    }

    if(iMCP != "Rm2" && iMCP != "BINP3"){ 
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
       if(doDoubleGauss == true){
          float f1 = g_res->GetParameter(0)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float f2 = g_res->GetParameter(3)/(g_res->GetParameter(0)+g_res->GetParameter(3));
          float sigma_eff = sqrt(f1*g_res->GetParameter(2)*g_res->GetParameter(2)+f2*g_res->GetParameter(5)*g_res->GetParameter(5) + f1*g_res->GetParameter(1)*g_res->GetParameter(1)+f2*g_res->GetParameter(4)*g_res->GetParameter(4) - (f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4))*(f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4)));

          sprintf (Sigma,"#sigma = %.0f ps",sigma_eff*1000);

          latexLabel->SetTextSize(0.05);
          latexLabel->SetNDC();
          latexLabel->SetTextFont(42); // helvetica
       }

       time_wrtRm2->GetXaxis()->SetTitle("t-t_{ref} (ns)");
        
       TCanvas* c2 = new TCanvas();
       c2->cd();
       time_wrtRm2->Draw("hist");
       if(doDoubleGauss == true) latexLabel->DrawLatex(0.72, 0.55,Sigma);
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

       TCanvas* c3 = new TCanvas();
       c3->cd();
       time_wrtMiB2_noCorrection->Draw("hist");
       g_res->Draw("same");
       c3 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
       delete g_res; 
    }
    
    if(iMCP != "Rm2" && iMCP != "BINP3"){
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
       c4 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_")+nameiMCP+std::string("_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
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

       TCanvas* c5 = new TCanvas();
       c5->cd();
       time_vs_amp_wrtMiB2_2->Draw();
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string("_auto.png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string("_auto.pdf")).c_str(),"pdf");
    }

    if(iMCP != "Rm2" && iMCP != "BINP3"){
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
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string("_auto.png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string("_auto.pdf")).c_str(),"pdf");
   }

   TGraphAsymmErrors* g_Res_vs_Amp_wrtMiB2 = new TGraphAsymmErrors(); 
   TGraphAsymmErrors* g_Res_vs_Amp_wrtRm2 = new TGraphAsymmErrors(); 

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

        h4->Draw(Selection17.c_str(),Selection16.c_str()); 
        h4->Draw(Selection18.c_str(),Selection16.c_str()); 
       
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];
        sprintf (NameOutput_wrtMiB2_png,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtRm2_png,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
         
        if(iMCP != "MiB2"){
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.1);
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
           if(doDoubleGauss == true){
              float f1 = resFitAlt_wrtMiB2[ii]->GetParameter(0)/(resFitAlt_wrtMiB2[ii]->GetParameter(0)+resFitAlt_wrtMiB2[ii]->GetParameter(3));
              float f2 = resFitAlt_wrtMiB2[ii]->GetParameter(3)/(resFitAlt_wrtMiB2[ii]->GetParameter(0)+resFitAlt_wrtMiB2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_wrtMiB2[ii]->GetParameter(2)*resFitAlt_wrtMiB2[ii]->GetParameter(2)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(5)*resFitAlt_wrtMiB2[ii]->GetParameter(5) + f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4)*resFitAlt_wrtMiB2[ii]->GetParameter(4) - (f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4))*(f1*resFitAlt_wrtMiB2[ii]->GetParameter(1)+f2*resFitAlt_wrtMiB2[ii]->GetParameter(4)));
           }else{
              sigma_eff = resFitAlt_wrtMiB2[ii]->GetParameter(2);
           }

           g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2,sigma_eff);
           //g_Res_vs_Amp_wrtMiB2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtMiB2[ii]->GetParError(2),resFitAlt_wrtMiB2[ii]->GetParError(2));
           g_Res_vs_Amp_wrtMiB2->SetPointError(ii,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,0.,0.);

           points_wrtMiB2.push_back(sigma_eff);
        
           TCanvas* c2 = new TCanvas();
           c2->cd();
           resHist_wrtMiB2[ii]->Draw("hist");
           resFitAlt_wrtMiB2[ii]->Draw("same");
           c2 -> Print(NameOutput_wrtMiB2_png,"png");
           c2 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c2;
        }

        if(iMCP != "Rm2" && iMCP != "BINP3"){
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.1);
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
           if(doDoubleGauss == true){
              float f1 = resFitAlt_wrtRm2[ii]->GetParameter(0)/(resFitAlt_wrtRm2[ii]->GetParameter(0)+resFitAlt_wrtRm2[ii]->GetParameter(3));
              float f2 = resFitAlt_wrtRm2[ii]->GetParameter(3)/(resFitAlt_wrtRm2[ii]->GetParameter(0)+resFitAlt_wrtRm2[ii]->GetParameter(3));
              sigma_eff = sqrt(f1*resFitAlt_wrtRm2[ii]->GetParameter(2)*resFitAlt_wrtRm2[ii]->GetParameter(2)+f2*resFitAlt_wrtRm2[ii]->GetParameter(5)*resFitAlt_wrtRm2[ii]->GetParameter(5) + f1*resFitAlt_wrtRm2[ii]->GetParameter(1)*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4)*resFitAlt_wrtRm2[ii]->GetParameter(4) - (f1*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4))*(f1*resFitAlt_wrtRm2[ii]->GetParameter(1)+f2*resFitAlt_wrtRm2[ii]->GetParameter(4)));
           }else{
              sigma_eff = resFitAlt_wrtRm2[ii]->GetParameter(2);    
           }

           g_Res_vs_Amp_wrtRm2->SetPoint(ii,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2.,sigma_eff);
           //g_Res_vs_Amp_wrtRm2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtRm2[ii]->GetParError(2),resFitAlt_wrtRm2[ii]->GetParError(2));
           g_Res_vs_Amp_wrtRm2->SetPointError(ii,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,0.,0.);

           points_wrtRm2.push_back(sigma_eff);
        
           TCanvas* c3 = new TCanvas();
           c3->cd();
           resHist_wrtRm2[ii]->Draw("hist");
           resFitAlt_wrtRm2[ii]->Draw("same");
           c3 -> Print(NameOutput_wrtRm2_png,"png");
           c3 -> Print(NameOutput_wrtRm2_pdf,"pdf");
           delete c3;
       }
    }

    std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    std::sort(points_wrtRm2.begin(),points_wrtRm2.end());
    
    if(iMCP != "MiB2"){
       g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
       g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
       g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
       g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-0.01,0.15);

       TCanvas* c7 = new TCanvas();
       c7->cd();
       g_Res_vs_Amp_wrtMiB2->Draw("AP");
       c7 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c7 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }
    
    if(iMCP != "Rm2" && iMCP != "BINP3"){
       g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
       g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
       g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
       g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
       g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
       g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-0.01,0.15);

       TCanvas* c8 = new TCanvas();
       c8->cd();
       g_Res_vs_Amp_wrtRm2->Draw("AP");
       c8 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c8 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }         
   } 

   TFile* output_ampMax = new TFile(std::string("TimeResolution_vs_amp_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".root").c_str(),"RECREATE");
   output_ampMax->cd();
   time_vs_amp_wrtMiB2_2->Write(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".root").c_str());
   time_vs_amp_wrtRm2_2->Write(std::string("TimeResolution_vs_amp_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".root").c_str());
   output_ampMax->Close(); 
}

void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, bool doScan_Corr, std::string HodoSelection)
{
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,ampMin,4000,timeMin);
    TH2F* timingCorrection_wrtRm2 = new TH2F("timingCorrection_wrtRm2","",nBins,ampMin,4000,timeMin);
    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
   
    std::string Selection1;   
    char Selection2 [1000];

    //if(iMCP == "M25") sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*amp_max[")+iMCP+std::string("]))<0.25")).c_str(),Params->at(0),Params->at(1));
    sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    if(iMCP == "MiB2"){
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection1 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,"0.",false);
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      Selection1 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }
    std::string Selection3 = Selection1+" && "+Selection2+HodoSelection; 
    
    if(iMCP != "MiB2")
       h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtMiB2")).c_str(),Selection3.c_str());
    if(iMCP != "Rm2" && iMCP != "BINP3")
       h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtRm2")).c_str(),Selection3.c_str());

    if(iMCP != "MiB2"){
       timingCorrection_wrtMiB2->FitSlicesY();
       TH2F* timingCorrection_wrtMiB2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtMiB2_1");
       timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time-time[MiB2]");
       timingCorrection_wrtMiB2_1->SetAxisRange(timingCorrection_wrtMiB2_1->GetMinimum()-0.1,timingCorrection_wrtMiB2_1->GetMaximum()+0.1, "Y");
       timingCorrection_wrtMiB2_1->SetMarkerStyle(20);
       timingCorrection_wrtMiB2_1->SetMarkerSize(0.9);
       timingCorrection_wrtMiB2_1->SetMarkerColor(kBlack);
       timingCorrection_wrtMiB2_1->SetLineColor(kBlack);
       
       //std::cout << "Extreme = " << timingCorrection_wrtMiB2_1->GetMinimum()-0.1 << " " << timingCorrection_wrtMiB2_1->GetMaximum()+0.1 << std::endl;
    
       TF1* fit_corr1;
       fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,3000.);
       timingCorrection_wrtMiB2_1->Fit("fit_corr1");
       if(doScan_Corr == false){
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
         Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
       }
    
       TCanvas* c1 = new TCanvas();
       c1->cd();
       timingCorrection_wrtMiB2_1->Draw();
       //timingCorrection_wrtMiB2->Draw();
       fit_corr1->Draw("same");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }
 
    if(iMCP != "Rm2" && iMCP != "BINP3"){
       timingCorrection_wrtRm2->FitSlicesY();
       TH2F* timingCorrection_wrtRm2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtRm2_1");
       timingCorrection_wrtRm2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       timingCorrection_wrtRm2_1->GetYaxis()->SetTitle("time-time[Rm2]");
       timingCorrection_wrtRm2_1->SetAxisRange(timingCorrection_wrtRm2_1->GetMinimum()-0.1,timingCorrection_wrtRm2_1->GetMaximum()+0.1, "Y");
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
       c2 -> Print((std::string("timingCorrection_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
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

        sprintf (NameOutput_wrtMiB2_png,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
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

        if(iMCP != "Rm2" && iMCP != "BINP3"){
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

         TCanvas* c5 = new TCanvas();
         c5->cd();
         g_Res_vs_Amp_wrtMiB2->Draw("AP");
         fit_corr3->Draw("same");
         c5 -> Print((std::string("timingCorrection2_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
         c5 -> Print((std::string("timingCorrection2_wrtMiB2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
      }
      
      if(iMCP!="Rm2" && iMCP != "BINP3"){
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
         c6 -> Print((std::string("timingCorrection2_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
         c6 -> Print((std::string("timingCorrection2_wrtRm2_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
      }
    }

}

void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string HodoSelection)
{
    TH2D* h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    TH2D* h2_time_maximum_vs_amp = new TH2D("h2_time_maximum_vs_amp","",300,0.,3000.,500,0.,5.);

    TF1* pol1_max = new TF1("pol1_max","pol1",0.,3000.);
    TF1* fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",0.,3000.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      if(iMCP == "M25") Selection = Selection + " && time_max["+iMCP+"]-time["+iMCP+iTiming+"] < 3.";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
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
    
    TCanvas* c6 = new TCanvas();
    c6->cd();
    h2_time_max_vs_amp->Draw("COLZ");
    if(iMCP == "M25") fit_corr_max->Draw("same");
    else fit_corr_max->Draw("same");
    c6 -> Print((std::string("deltaT_max_vs_amp_max_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
    c6 -> Print((std::string("deltaT_max_vs_amp_max_")+nameiMCP+std::string("_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
}

void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string HodoSelection)
{
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",300,-10,20,300,-1.,1.5,1000.,1500.);
    TProfile2D* p2D_amp_vs_time_noSelection = new TProfile2D("p2D_amp_vs_time_noSelection","",300,-10,20,300,-1.,1.5,1000.,1500.);
    TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time","",300,-10,20,300,-1.,1.5);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }

    Selection = Selection+" && WF_ch == "+iMCP;//+HodoSelection;

    //std::cout << "Selection: " << (std::string("WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> h2_amp_vs_time")).c_str() << "," << Selection << std::endl;

    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time_noSelection")).c_str(),std::string(std::string("WF_ch == ")+iMCP).c_str(),"goff");
    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time")).c_str(),Selection.c_str(),"goff");
    h4->Draw((std::string("WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> h2_amp_vs_time")).c_str(),Selection.c_str());

    TProfile* waveForm = h2_amp_vs_time->ProfileX(); 
    waveForm->SetName(std::string(nameiMCP+std::string("_waveform_prof")).c_str());
    
    p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time_noSelection->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time_noSelection->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
    p2D_amp_vs_time_noSelection->GetZaxis()->SetTitle("amp_max");

    TCanvas* c1 = new TCanvas();
    c1->cd();
    p2D_amp_vs_time_noSelection->Draw("COLZ");
    c1 -> Print(std::string(std::string("pulseShape_")+nameiMCP+std::string("_noSelection.png")).c_str(),"png");
    c1 -> Print(std::string(std::string("pulseShape_")+nameiMCP+std::string("_noSelection.pdf")).c_str(),"pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    p2D_amp_vs_time->Draw("COLZ");
    c2 -> Print(std::string(std::string("pulseShape_")+nameiMCP+std::string(".png")).c_str(),"png");
    c2 -> Print(std::string(std::string("pulseShape_")+nameiMCP+std::string(".pdf")).c_str(),"pdf");

    TFile* output_Waveform = new TFile(std::string(nameiMCP+std::string("_Waveform.root")).c_str(),"RECREATE");
    output_Waveform->cd();
    waveForm->Write(std::string(nameiMCP+std::string("_waveform_prof")).c_str());
    output_Waveform->Close();
}

void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP)
{
    TH2F* h2_hodoscope_Y_vs_X_noSelection = new TH2F("h2_hodoscope_Y_vs_X_noSelection","",30,0.,30.,30,0.,30.);
    TH2F* h2_hodoscope_Y_vs_X = new TH2F("h2_hodoscope_Y_vs_X","",30,0.,30.,30,0.,30.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && time_chi2["+iMCP+"]<"+timeChi2;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }

    std::cout << "Selection: " << Selection << std::endl;

    h4->Draw("Y:X >> h2_hodoscope_Y_vs_X_noSelection");
    h4->Draw("Y:X >> h2_hodoscope_Y_vs_X",Selection.c_str());

    h2_hodoscope_Y_vs_X_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X->GetYaxis()->SetTitle("Y");
    
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
    c2 -> Print(std::string(std::string("Hodoscope_Y_vs_X_")+nameiMCP+std::string(".png")).c_str(),"png");
    c2 -> Print(std::string(std::string("Hodoscope_Y_vs_X_")+nameiMCP+std::string(".pdf")).c_str(),"pdf");
}

void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP)
{
    TH1F* h_timeChi2_noSelection = new TH1F("h_timeChi2_noSelection","",1500,1.,15000.);
    TH1F* h_timeChi2 = new TH1F("h_timeChi2","",1500,1.,15000.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }   

    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2_noSelection").c_str());
    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2").c_str(),Selection.c_str());
    
    h_timeChi2_noSelection->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2_noSelection->GetYaxis()->SetTitle("Events");
    h_timeChi2->GetYaxis()->SetTitle("Events");
    
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
    c2 -> Print(std::string("timeChi2_"+nameiMCP+".png").c_str(),"png");
    c2 -> Print(std::string("timeChi2_"+nameiMCP+".pdf").c_str(),"pdf");
    
}   

void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP)
{
    TH1F* h_ampMax_noSelection = new TH1F(std::string("h_ampMax_"+nameiMCP+"_noSelection").c_str(),"",450,0.,4500.);
    TH1F* h_ampMax = new TH1F(std::string("h_ampMax_"+nameiMCP).c_str(),"",450,0.,4500.);

    std::string Selection;
    std::string iTiming = "";

    if(iMCP == "MiB2"){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else if(iMCP == "Rm2" || iMCP == "BINP3" ){
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }else{
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max["+iMCP+"]-time[")+iMCP+iTiming+std::string("]"),Selection,"0.",false);
    }   

    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"noSelection").c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP).c_str(),Selection.c_str());
    
    h_ampMax_noSelection->GetXaxis()->SetTitle("amp_max");
    h_ampMax->GetXaxis()->SetTitle("amp_max");
    h_ampMax_noSelection->GetYaxis()->SetTitle("Events");
    h_ampMax->GetYaxis()->SetTitle("Events");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogy();
    h_ampMax_noSelection->Draw("H");
    c1 -> Print(std::string("ampMax_noSelection_"+nameiMCP+".png").c_str(),"png");
    c1 -> Print(std::string("ampMax_noSelection_"+nameiMCP+".pdf").c_str(),"pdf");

    TFile* output_ampMax_noSelection = new TFile(std::string("ampMax_noSelection_"+nameiMCP+".root").c_str(),"RECREATE");
    output_ampMax_noSelection->cd();
    h_ampMax_noSelection->Write();
    output_ampMax_noSelection->Close();

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetLogy();
    h_ampMax->Draw("H");
    c2 -> Print(std::string("ampMax_"+nameiMCP+".png").c_str(),"png");
    c2 -> Print(std::string("ampMax_"+nameiMCP+".pdf").c_str(),"pdf");
    
    TFile* output_ampMax = new TFile(std::string("ampMax_"+nameiMCP+".root").c_str(),"RECREATE");
    output_ampMax->cd();
    h_ampMax->Write();
    output_ampMax->Close();

}   

std::string AddSelection(TTree* h4, std::string Var, std::string Selection, std::string Cut = "0.", bool isCut = false)
{
    TH1F* h = new TH1F("h","h",4000,-20.,20.);
    TF1* g_fit = new TF1("g_fit","gaus",-20.,20.);
    h4->Draw((Var+std::string(" >> h")).c_str(),Selection.c_str());
    
    if(isCut == false){
      h->Fit("g_fit","","",h->GetMean()-h->GetRMS()/2.,h->GetMean()+h->GetRMS()/2.);
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


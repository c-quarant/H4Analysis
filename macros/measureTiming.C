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

//int nBins = 10;
//float ampMin[11] = {0.};
int nBins = 8;
float ampMin[9] = {0.};
float *ptr_minVec = ampMin;
int nBinsFrac = 10;
float ampMaxFrac[10] = {0.};
float *ptr_minVecFrac = ampMaxFrac;
float timeMin[4001];

//BINP
std::string amp_max_MiB2 = "200.";
std::string time_max_MiB2 = "150.";
std::string time_max_Rm2 = "150.";
std::string amp_max_Rm2 = "200.";
std::string scintMin = "200.";
std::string scintMax = "700.";
std::string timeChi2 = "99999999.";
bool isEfficiencyRun = false;

void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan, bool doScanEff, bool doDoubleGauss, std::string HodoSelection, bool doOnlyWrtMiB2);
void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2);
void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2);
void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2);
void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2);
void CheckEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void HodoValidation(TTree* h4, TTree* hodo, TTree* adc, int nParticle);
std::string AddSelection(TTree*, std::string, std::string, std::string, bool);
std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel);
void CheckSelectionEfficiency(TTree* h4, std::string iMCP, std::string Selection);
std::vector<std::string> split(const std::string text, std::string sep);
void setBins(TTree* h4, std::string var, std::string Selection, std::string maxMCP, float step, float* ptr_minVec, int nbins);
void setBinsFraction(TTree* h4, std::string var, std::string Selection, std::string maxMCP, float step, float* ptr_minVec, int nbins);
std::pair<float,float> FWHMSigma(TH1F*,TF1*);
std::pair<float,float> EffectiveSigma(TF1*,int nSteps = 30000,float fraction = 0.34);

void measureTiming(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, std::string maxMCP, bool doOnlyWrtMiB2 = false, bool doFirstStep = true, bool doPulseShapes = false, bool doScan = false, bool doScanEff = false, bool doDoubleGauss = false, std::string hodoXmin="0", std::string hodoXmax="30", std::string hodoYmin="0", std::string hodoYmax="30")
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
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
       std::cout << "TIME CORRECTION!" << std::endl;
       TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, maxMCP, HodoSelection, doOnlyWrtMiB2);
       std::cout << "FINAL TIMING!" << std::endl;
       FinalTiming(h4, inputs, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, maxMCP, doScan, doScanEff, doDoubleGauss, HodoSelection, doOnlyWrtMiB2);
    }
}

void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, bool doScan, bool doScanEff, bool doDoubleGauss, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    //TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",200,-0.5,0.5);
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",100,-0.5,0.5);
    TH1F* time_wrtRm2 = new TH1F("time_wrtRm2","",200,-0.5,0.5);
    TH1F* time_wrtMiB2_noCorrection = new TH1F("time_wrtMiB2_noCorrection","",1600,-10.,10.);
    TH1F* time_wrtRm2_noCorrection = new TH1F("time_wrtRm2_noCorrection","",1600,-10.,10.);
    TH2F* time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins,ampMin,4000,timeMin);
    TH2F* time_vs_amp_wrtRm2 = new TH2F("time_vs_amp_wrtRm2","",nBins,ampMin,4000,timeMin);
    TH1F* time_wrtMiB2_FWHM = new TH1F("time_wrtMiB2_FWHM","",2000,-1.,1.);
    TH1F* time_wrtRm2_FWHM = new TH1F("time_wrtRm2_FWHM","",2000,-1.,1.);

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    //TF1 *g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-5.,5.);
    TF1 *g_res;
    TF1 *g_res_single;
    TF1 *g_res_double;

    char Selection1 [1000];  
    char Selection2 [1000];
    char Selection3 [1000];
    char Selection4 [1000];
    char Selection5 [1000];
    
    //(time-time_max) vs amp selection
    sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    //time correction
    if(iMCP != "MiB2"){
       //sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       //sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       //sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       //sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*exp((%f)*amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*exp((%f)*amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       //sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*amp_max[")+iMCP+std::string("]) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       //sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*amp_max[")+iMCP+std::string("]):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
    }
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       //sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       //sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       //sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       //sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       //sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-sqrt((%f)+(%f)*amp_max[")+iMCP+std::string("]+(%f)*amp_max[")+iMCP+std::string("]*amp_max[")+iMCP+std::string("]) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       //sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-sqrt((%f)+(%f)*amp_max[")+iMCP+std::string("]+(%f)*amp_max[")+iMCP+std::string("]*amp_max[")+iMCP+std::string("]):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       sprintf (Selection4,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-((%f)+(%f)*amp_max[")+iMCP+std::string("]) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
       sprintf (Selection5,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-((%f)+(%f)*amp_max[")+iMCP+std::string("]):amp_max[")+iMCP+std::string("] >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
    }

    //no time correction
    std::string Selection6 = "time["+iMCP+iTiming+"]-time[MiB2] >>";
    std::string Selection7 = "time["+iMCP+iTiming+"]-time[Rm2] >>";
    std::string Selection8;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection8 = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection8,"0.",false);
      Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection8,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
      }else{
         Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection8,"0.",false);
         Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
         Selection8 = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
      }   
    } 
    std::string Selection9 = Selection8+std::string(" && ")+std::string(Selection1)+HodoSelection; 

    std::cout << "Selection = " << Selection9 << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection9);

    std::string Selection10 = Selection2+std::string(" time_wrtMiB2");
    std::string Selection11 = Selection3+std::string(" time_vs_amp_wrtMiB2");
    std::string Selection12 = Selection4+std::string(" time_wrtRm2");
    std::string Selection13 = Selection5+std::string(" time_vs_amp_wrtRm2"); 
    std::string Selection14 = Selection6+std::string(" time_wrtMiB2_noCorrection");
    std::string Selection15 = Selection7+std::string(" time_wrtRm2_noCorrection");

    if(iMCP != "MiB2"){
       h4->Draw(Selection10.c_str(),Selection9.c_str()); 
       h4->Draw(Selection11.c_str(),Selection9.c_str()); 
       h4->Draw(Selection14.c_str(),Selection9.c_str()); 
    }
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       h4->Draw(Selection12.c_str(),Selection9.c_str());
       h4->Draw(Selection13.c_str(),Selection9.c_str()); 
       h4->Draw(Selection15.c_str(),Selection9.c_str()); 
    }

    if(time_wrtMiB2->GetEntries() < 2000) time_wrtMiB2->Rebin(4);
    if(time_wrtRm2->GetEntries() < 2000) time_wrtRm2->Rebin(4);

    float sub_wrtMiB2 = 0.;
    float sub_wrtMiB2_error = 0.;
    float sub_wrtRm2 = 0.;
    float sub_wrtRm2_error = 0.;

    if(isEfficiencyRun == false){
       sub_wrtMiB2 = 24.E-3;
       sub_wrtMiB2_error = 2.E-3;
       sub_wrtRm2 = 18.E-3;
       sub_wrtRm2_error = 3.E-3;
    }else{
       sub_wrtMiB2 = 26.E-3;
       sub_wrtMiB2_error = 3.E-3;
       sub_wrtRm2 = 24.E-3;
       sub_wrtRm2_error = 3.E-3;
    }

    TF1* g_res_final;

    if(iMCP != "MiB2"){

       float fwhm;
       float fwhm_error;

       if(doDoubleGauss == false) g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-3*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+3*time_wrtMiB2->GetRMS());  
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-3*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+3*time_wrtMiB2->GetRMS());  
       if(doDoubleGauss == true) g_res_single = new TF1("g_res_single","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-3*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+3*time_wrtMiB2->GetRMS());
       if(doDoubleGauss == true) g_res_double = new TF1("g_res_double","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-3*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+3*time_wrtMiB2->GetRMS());
       g_res->SetParameters(0,time_wrtMiB2->GetEntries()/2.);
       g_res->SetParameters(1,0.);
       g_res->SetParameters(2,0.04);
       g_res->SetParameters(3,0.);
       if(doDoubleGauss == true) g_res->SetParameters(3,time_wrtMiB2->GetEntries()/2.);
       if(doDoubleGauss == true) g_res->SetParameters(4,0.);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.1);
       //f(doDoubleGauss == true) g_res->SetParameters(5,0.04);
       if(doDoubleGauss == true) g_res->SetParameters(6,0.);
       g_res->SetParLimits(0,0.,time_wrtMiB2->GetEntries());
       g_res->SetParLimits(1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       //g_res->SetParLimits(2,0.,0.2);
       g_res->SetParLimits(2,0.,5.);
       g_res->SetParLimits(3,0.,10.);
       g_res->SetParLimits(3,0.,time_wrtMiB2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       //if(doDoubleGauss == true) g_res->SetParLimits(5,0.04,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(6,0.,10.);
       time_wrtMiB2->Fit("g_res","B");
       g_res->SetParameters(0,g_res->GetParameter(0));
       g_res->SetParameters(1,g_res->GetParameter(1));
       g_res->SetParameters(2,g_res->GetParameter(2));
       g_res->SetParameters(3,g_res->GetParameter(3));
       if(doDoubleGauss == true) g_res->SetParameters(4,g_res->GetParameter(4));
       if(doDoubleGauss == true) g_res->SetParameters(5,g_res->GetParameter(5));
       if(doDoubleGauss == true) g_res->SetParameters(6,g_res->GetParameter(6));
       g_res->SetParLimits(0,0.,time_wrtMiB2->GetEntries());
       g_res->SetParLimits(1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       //g_res->SetParLimits(2,0.,0.2);
       g_res->SetParLimits(2,0.,5.);
       g_res->SetParLimits(3,0.,10.);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtMiB2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
       //if(doDoubleGauss == true) g_res->SetParLimits(5,0.04,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(6,0.,10.);
       g_res->SetNpx(10000);
       time_wrtMiB2->Fit("g_res","B");
     
       if(doDoubleGauss == true){
          g_res_single->FixParameter(0,g_res->GetParameter(0));
          g_res_single->FixParameter(1,g_res->GetParameter(1));
          g_res_single->FixParameter(2,g_res->GetParameter(2));
          g_res_single->FixParameter(3,g_res->GetParameter(6));
          g_res_single->SetLineColor(kBlue+1);

          g_res_double->FixParameter(0,g_res->GetParameter(3));
          g_res_double->FixParameter(1,g_res->GetParameter(4));
          g_res_double->FixParameter(2,g_res->GetParameter(5));
          g_res_double->FixParameter(3,g_res->GetParameter(6));
          g_res_double->SetLineColor(kGreen+1);

          fwhm = FWHMSigma(time_wrtMiB2,g_res).first;
          fwhm_error = FWHMSigma(time_wrtMiB2,g_res).second;
       }

       char Sigma[100];
       TLatex *latexLabel = new TLatex();
       float sigma_eff;
       float s_sigma;
       float sigma_eff_sub;
       float s_sigma_sub;
       //float sub_wrtMiB2 = 24.E-3;
       //float sub_wrtMiB2_error = 2.E-3;

       sigma_eff = g_res->GetParameter(2);
       s_sigma = g_res->GetParError(2);
       if(doDoubleGauss == true){
          sigma_eff = fwhm/(2.*sqrt(2*log(2)));
          s_sigma = fwhm_error/(2*sqrt(2*log(2)));
       }

       if(iMCP != "MiB2" && iMCP != "Rm2"){
          sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtMiB2*sub_wrtMiB2);
          s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtMiB2*sub_wrtMiB2*sub_wrtMiB2_error*sub_wrtMiB2_error)/sigma_eff_sub;
       }else{
          sigma_eff_sub = sigma_eff; 
          s_sigma_sub = s_sigma;
       }
       sigma_eff = sigma_eff_sub;
       s_sigma = s_sigma_sub;

       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       /*TFile* output_ampMaxMiB2 = new TFile(std::string("Data_TimeResolution_vs_amp_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+"_Inclusive.root").c_str(),"RECREATE");
       output_ampMaxMiB2->cd();
       time_wrtMiB2->Write();
       g_res->Write();
       output_ampMaxMiB2->Close();*/
     
      std::cout << "sigma = " << sigma_eff*1000. << "+/-" << s_sigma*1000. << std::endl; 
      
       sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

       latexLabel->SetTextSize(0.05);
       latexLabel->SetNDC();
       latexLabel->SetTextFont(42); // helvetica

       time_wrtMiB2->GetXaxis()->SetTitle("t-t_{ref} (ns)");
  
       TCanvas* c1 = new TCanvas();
       c1->cd();
       time_wrtMiB2->Draw("hist");
       latexLabel->DrawLatex(0.72, 0.45,Sigma);
       if(doDoubleGauss == true) g_res_single->Draw("same");
       if(doDoubleGauss == true) g_res_double->Draw("same");
       g_res->Draw("same");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
       g_res_final = (TF1*)g_res->Clone("g_res_final");

       delete g_res; 
    }

    

    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){ 
       
       float fwhm;
       float fwhm_error;

       if(doDoubleGauss == false) g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-3*time_wrtRm2->GetRMS(),time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+3*time_wrtRm2->GetRMS());  
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-3*time_wrtRm2->GetRMS(),time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+3*time_wrtRm2->GetRMS());  
       if(doDoubleGauss == true) g_res_single = new TF1("g_res_single","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-3*time_wrtRm2->GetRMS(),time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+3*time_wrtRm2->GetRMS());
       if(doDoubleGauss == true) g_res_double = new TF1("g_res_double","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-3*time_wrtRm2->GetRMS(),time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+3*time_wrtRm2->GetRMS());
       g_res->SetParameters(0,time_wrtRm2->GetEntries()/2.);
       g_res->SetParameters(1,0.);
       g_res->SetParameters(2,0.04);
       g_res->SetParameters(3,0.);
       if(doDoubleGauss == true) g_res->SetParameters(3,time_wrtRm2->GetEntries()/2.);
       if(doDoubleGauss == true) g_res->SetParameters(4,0.);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.1);
       //f(doDoubleGauss == true) g_res->SetParameters(5,0.04);
       if(doDoubleGauss == true) g_res->SetParameters(6,0.);
       g_res->SetParLimits(0,0.,time_wrtRm2->GetEntries());
       g_res->SetParLimits(1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       //g_res->SetParLimits(2,0.,0.2);
       g_res->SetParLimits(2,0.,5.);
       g_res->SetParLimits(3,0.,10.);
       g_res->SetParLimits(3,0.,time_wrtRm2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       //if(doDoubleGauss == true) g_res->SetParLimits(5,0.04,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(6,0.,10.);
       time_wrtRm2->Fit("g_res","B");
       g_res->SetParameters(0,g_res->GetParameter(0));
       g_res->SetParameters(1,g_res->GetParameter(1));
       g_res->SetParameters(2,g_res->GetParameter(2));
       g_res->SetParameters(3,g_res->GetParameter(3));
       if(doDoubleGauss == true) g_res->SetParameters(4,g_res->GetParameter(4));
       if(doDoubleGauss == true) g_res->SetParameters(5,g_res->GetParameter(5));
       if(doDoubleGauss == true) g_res->SetParameters(6,g_res->GetParameter(6));
       g_res->SetParLimits(0,0.,time_wrtRm2->GetEntries());
       g_res->SetParLimits(1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       //g_res->SetParLimits(2,0.,0.2);
       g_res->SetParLimits(2,0.,5.);
       g_res->SetParLimits(3,0.,10.);
       if(doDoubleGauss == true) g_res->SetParLimits(3,0.,time_wrtRm2->GetEntries());
       if(doDoubleGauss == true) g_res->SetParLimits(4,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
       //if(doDoubleGauss == true) g_res->SetParLimits(5,0.04,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(5,0.,5.);
       if(doDoubleGauss == true) g_res->SetParLimits(6,0.,10.);
       time_wrtRm2->Fit("g_res","B");
     
       if(doDoubleGauss == true){
          g_res_single->FixParameter(0,g_res->GetParameter(0));
          g_res_single->FixParameter(1,g_res->GetParameter(1));
          g_res_single->FixParameter(2,g_res->GetParameter(2));
          g_res_single->FixParameter(3,g_res->GetParameter(6));
          g_res_single->SetLineColor(kBlue+1);

          g_res_double->FixParameter(0,g_res->GetParameter(3));
          g_res_double->FixParameter(1,g_res->GetParameter(4));
          g_res_double->FixParameter(2,g_res->GetParameter(5));
          g_res_double->FixParameter(3,g_res->GetParameter(6));
          g_res_double->SetLineColor(kGreen+1);

          fwhm = FWHMSigma(time_wrtRm2,g_res).first;
          fwhm_error = FWHMSigma(time_wrtRm2,g_res).second;
       }

       char Sigma[100];
       TLatex *latexLabel = new TLatex();
       float sigma_eff;
       float s_sigma;
       float sigma_eff_sub;
       float s_sigma_sub;
       //float sub_wrtRm2 = 24.E-3;
       //float sub_wrtRm2_error = 2.E-3;

       sigma_eff = g_res->GetParameter(2);
       s_sigma = g_res->GetParError(2);
       if(doDoubleGauss == true){
          sigma_eff = fwhm/(2.*sqrt(2*log(2)));
          s_sigma = fwhm_error/(2*sqrt(2*log(2)));
       }

       if(iMCP != "Rm2" && iMCP != "Rm2"){
          sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtRm2*sub_wrtRm2);
          s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtRm2*sub_wrtRm2*sub_wrtRm2_error*sub_wrtRm2_error)/sigma_eff_sub;
       }else{
          sigma_eff_sub = sigma_eff; 
          s_sigma_sub = s_sigma;
       }
       sigma_eff = sigma_eff_sub;
       s_sigma = s_sigma_sub;

       sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

       latexLabel->SetTextSize(0.05);
       latexLabel->SetNDC();
       latexLabel->SetTextFont(42); // helvetica

       time_wrtRm2->GetXaxis()->SetTitle("t-t_{ref} (ns)");

       TCanvas* c1 = new TCanvas();
       c1->cd();
       time_wrtRm2->Draw("hist");
       latexLabel->DrawLatex(0.72, 0.45,Sigma);
       if(doDoubleGauss == true) g_res_single->Draw("same");
       if(doDoubleGauss == true) g_res_double->Draw("same");
       g_res->Draw("same");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtRm2_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtRm2_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
       g_res_final = (TF1*)g_res->Clone("g_res_final");

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

   TGraphAsymmErrors* g_Res_vs_ROC_wrtMiB2 = new TGraphAsymmErrors(); 
   TGraphAsymmErrors* g_Res_vs_ROC_wrtRm2 = new TGraphAsymmErrors(); 

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

    std::string wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";
    TFile* file_timingCorrection_wrtMiB2 = new TFile(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".root").c_str(),"RECREATE");
    file_timingCorrection_wrtMiB2->cd();

    for(int ii = 0; ii < nBins; ii++)
    {
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_Res_1_%d_wrtMiB2",ii);
        resHist_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,100,-0.5,0.5);

        char Name_wrtRm2 [50];
        sprintf (Name_wrtRm2,"h_Res_1_%d_wrtRm2",ii);
        resHist_wrtRm2[ii] = new TH1F(Name_wrtRm2,Name_wrtRm2,100,-0.5,0.5);

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
         
        doDoubleGauss = false; 
        
        if(iMCP != "MiB2"){
           h4->Draw(Selection17.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           if(doDoubleGauss == false) resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-3*resHist_wrtMiB2[ii]->GetRMS(),resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+3*resHist_wrtMiB2[ii]->GetRMS()); 
           else resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-5*resHist_wrtMiB2[ii]->GetRMS(),resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+5*resHist_wrtMiB2[ii]->GetRMS());  
           resFitAlt_wrtMiB2[ii]->SetParameters(0,resHist_wrtMiB2[ii]->GetEntries()/2.);
           resFitAlt_wrtMiB2[ii]->SetParameters(1,0.);
           if(doDoubleGauss == false) resFitAlt_wrtMiB2[ii]->SetParameters(2,0.04); 
           else resFitAlt_wrtMiB2[ii]->SetParameters(2,0.04);
           //else resFitAlt_wrtMiB2[ii]->SetParameters(2,0.1);
           resFitAlt_wrtMiB2[ii]->SetParameters(3,0.); 
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(3,resHist_wrtMiB2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(4,0.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(5,0.1); 
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(6,0.); 
           resFitAlt_wrtMiB2[ii]->SetParLimits(0,0.,resHist_wrtMiB2[ii]->GetEntries());
           resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == false) resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,0.3); 
           else resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,0.5);
           resFitAlt_wrtMiB2[ii]->SetParLimits(3,0.,10.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(3,0.,resHist_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(4,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(5,0.06,5.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(6,0.,10.);
           resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
           resFitAlt_wrtMiB2[ii]->SetParameters(0,resFitAlt_wrtMiB2[ii]->GetParameter(0));
           resFitAlt_wrtMiB2[ii]->SetParameters(1,resFitAlt_wrtMiB2[ii]->GetParameter(1));
           resFitAlt_wrtMiB2[ii]->SetParameters(2,resFitAlt_wrtMiB2[ii]->GetParameter(2));
           resFitAlt_wrtMiB2[ii]->SetParameters(3,resFitAlt_wrtMiB2[ii]->GetParameter(3));
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(4,resFitAlt_wrtMiB2[ii]->GetParameter(4));
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(5,resFitAlt_wrtMiB2[ii]->GetParameter(5)); 
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParameters(6,resFitAlt_wrtMiB2[ii]->GetParameter(6)); 
           resFitAlt_wrtMiB2[ii]->SetParLimits(0,0.,resHist_wrtMiB2[ii]->GetEntries());
           resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == false) resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,1.);
           else resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,0.3);
           //else resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,0.5);
           resFitAlt_wrtMiB2[ii]->SetParLimits(3,0.,10.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(3,0.,resHist_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(4,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(5,0.06,5.);
           if(doDoubleGauss == true) resFitAlt_wrtMiB2[ii]->SetParLimits(6,0.,10.);
           resHist_wrtMiB2[ii]->Write(resHist_wrtMiB2[ii]->GetName());
           resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");

           float sigma_eff = resFitAlt_wrtMiB2[ii]->GetParameter(2);
           float s_sigma = resFitAlt_wrtMiB2[ii]->GetParError(2);
           float sigma_eff_sub;
           float s_sigma_sub;
    
           if(iMCP != "MiB2"){
              sigma_eff_sub = sqrt(sigma_eff*sigma_eff - sub_wrtMiB2*sub_wrtMiB2);
              s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma + sub_wrtMiB2*sub_wrtMiB2*sub_wrtMiB2_error*sub_wrtMiB2_error)/sigma_eff_sub;
              if(sigma_eff < sub_wrtMiB2) sigma_eff_sub = 0.;
           }else{
              sigma_eff_sub = sigma_eff;
              s_sigma_sub = s_sigma;
           }
           sigma_eff = sigma_eff_sub;
           s_sigma = s_sigma_sub;
          
           if(sigma_eff_sub*1000 < 10) continue;
           g_Res_vs_Amp_wrtMiB2->SetPoint(iPoint_MiB2,ampMin[ii]+(ampMin[ii+1]-ampMin[ii])/2,sigma_eff*1000.);
           g_Res_vs_Amp_wrtMiB2->SetPointError(iPoint_MiB2,(ampMin[ii+1]-ampMin[ii])/2,(ampMin[ii+1]-ampMin[ii])/2,s_sigma*1000.,s_sigma*1000.);

           iPoint_MiB2++;

           points_wrtMiB2.push_back(sigma_eff*1000.);
        
           char Sigma[100];
           sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

           std::cout << "sigma = " << sigma_eff*1000. << "+/-" << s_sigma*1000. << std::endl; 

           TLatex *latexLabel = new TLatex();
           latexLabel->SetTextSize(0.05);
           latexLabel->SetNDC();
           latexLabel->SetTextFont(42); // helvetica

           resHist_wrtMiB2[ii]->GetXaxis()->SetTitle("t-t_{ref} (ns)");

           TCanvas* c2 = new TCanvas();
           c2->cd();
           resHist_wrtMiB2[ii]->Draw("hist");
           latexLabel->DrawLatex(0.72, 0.55,Sigma);
           resFitAlt_wrtMiB2[ii]->Draw("same");
           c2 -> Print(NameOutput_wrtMiB2_png,"png");
           c2 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c2;
        }

        if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
           h4->Draw(Selection18.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           if(doDoubleGauss == false) resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-3*resHist_wrtRm2[ii]->GetRMS(),resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+3*resHist_wrtRm2[ii]->GetRMS()); 
           else resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-3*resHist_wrtRm2[ii]->GetRMS(),resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+3*resHist_wrtRm2[ii]->GetRMS());  
           resFitAlt_wrtRm2[ii]->SetParameters(0,resHist_wrtRm2[ii]->GetEntries()/2.);
           resFitAlt_wrtRm2[ii]->SetParameters(1,0.);
           resFitAlt_wrtRm2[ii]->SetParameters(2,0.04);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(3,resHist_wrtRm2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(4,0.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(5,0.1); 
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(6,0.); 
           resFitAlt_wrtRm2[ii]->SetParLimits(0,0.,resHist_wrtRm2[ii]->GetEntries());
           resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,0.5);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(3,0.,resHist_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(4,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(5,0.1,5.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(6,0.,10.);
           resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
           resFitAlt_wrtRm2[ii]->SetParameters(0,resFitAlt_wrtRm2[ii]->GetParameter(0));
           resFitAlt_wrtRm2[ii]->SetParameters(1,resFitAlt_wrtRm2[ii]->GetParameter(1));
           resFitAlt_wrtRm2[ii]->SetParameters(2,resFitAlt_wrtRm2[ii]->GetParameter(2));
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(3,resFitAlt_wrtRm2[ii]->GetParameter(3));
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(4,resFitAlt_wrtRm2[ii]->GetParameter(4));
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(5,resFitAlt_wrtRm2[ii]->GetParameter(5)); 
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParameters(6,resFitAlt_wrtRm2[ii]->GetParameter(6)); 
           resFitAlt_wrtRm2[ii]->SetParLimits(0,0.,resHist_wrtRm2[ii]->GetEntries());
           resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,0.5);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(3,0.,resHist_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(4,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(5,0.1,5.);
           if(doDoubleGauss == true) resFitAlt_wrtRm2[ii]->SetParLimits(6,0.,10.);
           resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");

           float sigma_eff = resFitAlt_wrtRm2[ii]->GetParameter(2);
           float s_sigma = resFitAlt_wrtRm2[ii]->GetParError(2);
           float sigma_eff_sub;
           float s_sigma_sub;
    
           if(iMCP != "Rm2"){
              sigma_eff_sub = sqrt(sigma_eff*sigma_eff - sub_wrtRm2*sub_wrtRm2);
              s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma + sub_wrtRm2*sub_wrtRm2*sub_wrtRm2_error*sub_wrtRm2_error)/sigma_eff_sub;
           }else{
              sigma_eff_sub = sigma_eff;
              s_sigma_sub = s_sigma;
           }
           sigma_eff = sigma_eff_sub;
           s_sigma = s_sigma_sub;
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
    file_timingCorrection_wrtMiB2->Close();

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

   int nBinsROC = nBinsFrac; ;

   if(doScanEff == true){
    
    std::vector<float> points_wrtMiB2;
    std::vector<float> points_wrtRm2;
    
    std::vector<TH1F*> resHist_ROC_wrtMiB2;
    resHist_ROC_wrtMiB2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_wrtMiB2;
    resFitAlt_ROC_wrtMiB2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_single_wrtMiB2;
    resFitAlt_ROC_single_wrtMiB2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_double_wrtMiB2;
    resFitAlt_ROC_double_wrtMiB2.resize(nBinsROC);


    std::vector<TH1F*> resHist_ROC_wrtRm2;
    resHist_ROC_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFit_ROC_wrtRm2;
    resFit_ROC_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_wrtRm2;
    resFitAlt_ROC_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_single_wrtRm2;
    resFitAlt_ROC_single_wrtRm2.resize(nBinsROC);
    std::vector<TF1*> resFitAlt_ROC_double_wrtRm2;
    resFitAlt_ROC_double_wrtRm2.resize(nBinsROC);

    int iPoint_MiB2 = 0;
    int iPoint_Rm2 = 0;

    wrtMCP = "";
    if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

    TFile* output_ampMax = new TFile(std::string("Data_TimeResolution_vs_eff_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".root").c_str(),"RECREATE");
    output_ampMax->cd();
    for(int ii = 0; ii < nBinsROC; ii++)
    {

        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];

        sprintf (NameOutput_wrtMiB2_png,"TimeResolution_vs_fraction_%s_wrtMiB2_%d_%s_thres%s%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(), wrtMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_vs_fraction_%s_wrtMiB2_%d_%s_thres%s%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str(), wrtMCP.c_str());
        sprintf (NameOutput_wrtRm2_png,"TimeResolution_vs_fraction_%s_wrtRm2_%d_%s_thres%s.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_vs_fraction_%s_wrtRm2_%d_%s_thres%s.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
         
        doDoubleGauss = true; 

        if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";
        else wrtMCP = "";
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_Res_%d_ROC_wrtMiB2%s",ii,wrtMCP.c_str());
        resHist_ROC_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,100,-0.5,0.5);

        char Name_wrtRm2 [50];
        sprintf (Name_wrtRm2,"h_Res_%d_ROC_wrtRm2",ii);
        resHist_ROC_wrtRm2[ii] = new TH1F(Name_wrtRm2,Name_wrtRm2,100,-0.5,0.5);

        char cutMin [10];
        char cutMax [10];
        sprintf (cutMin,"%f",ampMaxFrac[ii]);

        std::cout << "CutMin = " << cutMin << std::endl;
        
        std::string Selection16 = Selection9+" && amp_max["+iMCP+"]>"+cutMin;
        std::string Selection17 = Selection2+std::string(Name_wrtMiB2);
        std::string Selection18 = Selection4+std::string(Name_wrtRm2);

        float step = 100./nBinsFrac;

        if(iMCP != "MiB2"){
           float fwhm = 0.;
           float fwhm_error = 0.;

           h4->Draw(Selection17.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtMiB2 [100];
           sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
           if(doDoubleGauss == false) resFitAlt_ROC_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-5*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+5*resHist_ROC_wrtMiB2[ii]->GetRMS());  
           else resFitAlt_ROC_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-5*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+5*resHist_ROC_wrtMiB2[ii]->GetRMS());  
           if(doDoubleGauss == true) resFitAlt_ROC_single_wrtMiB2[ii] = new TF1((std::string(NameFitAlt_wrtMiB2)+std::string("_single")).c_str(),"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-5*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+5*resHist_ROC_wrtMiB2[ii]->GetRMS());
           if(doDoubleGauss == true) resFitAlt_ROC_double_wrtMiB2[ii] = new TF1((std::string(NameFitAlt_wrtMiB2)+std::string("_double")).c_str(),"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-5*resHist_ROC_wrtMiB2[ii]->GetRMS(),resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+5*resHist_ROC_wrtMiB2[ii]->GetRMS());
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(0,resHist_ROC_wrtMiB2[ii]->GetEntries()/2.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(1,0.);
           //resFitAlt_ROC_wrtMiB2[ii]->SetParameters(2,0.04);
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(2,0.1);
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(3,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(3,resHist_ROC_wrtMiB2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(4,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(5,0.1);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(5,0.04);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(6,0.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(0,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(1,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           //resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(2,0.,0.2);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(2,0.,5.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(3,0.,10.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(3,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(4,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(5,0.04,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(5,0.,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(6,0.,10.);
           resHist_ROC_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(0,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(0));
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(1,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1));
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(2,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2));
           resFitAlt_ROC_wrtMiB2[ii]->SetParameters(3,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(3));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(4,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(5,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(5));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParameters(6,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(6));
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(0,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(1,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           //resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(2,0.,0.2);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(2,0.,5.);
           resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(3,0.,10.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(3,0.,resHist_ROC_wrtMiB2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(4,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtMiB2[ii]->GetBinCenter(resHist_ROC_wrtMiB2[ii]->GetMaximumBin())+0.2);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(5,0.04,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(5,0.,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtMiB2[ii]->SetParLimits(6,0.,10.);
           resHist_ROC_wrtMiB2[ii]->Write();
           resHist_ROC_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
           resFitAlt_ROC_wrtMiB2[ii]->SetNpx(10000);
           
           if(doDoubleGauss == true){
              resFitAlt_ROC_single_wrtMiB2[ii]->FixParameter(0,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(0));
              resFitAlt_ROC_single_wrtMiB2[ii]->FixParameter(1,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(1));
              resFitAlt_ROC_single_wrtMiB2[ii]->FixParameter(2,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2));
              resFitAlt_ROC_single_wrtMiB2[ii]->FixParameter(3,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(6));
              resFitAlt_ROC_single_wrtMiB2[ii]->SetLineColor(kBlue+1);

              resFitAlt_ROC_double_wrtMiB2[ii]->FixParameter(0,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(3));
              resFitAlt_ROC_double_wrtMiB2[ii]->FixParameter(1,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(4));
              resFitAlt_ROC_double_wrtMiB2[ii]->FixParameter(2,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(5));
              resFitAlt_ROC_double_wrtMiB2[ii]->FixParameter(3,resFitAlt_ROC_wrtMiB2[ii]->GetParameter(6));
              resFitAlt_ROC_double_wrtMiB2[ii]->SetLineColor(kGreen+1);

             fwhm = FWHMSigma(resHist_ROC_wrtMiB2[ii],resFitAlt_ROC_wrtMiB2[ii]).first;
             fwhm_error = FWHMSigma(resHist_ROC_wrtMiB2[ii],resFitAlt_ROC_wrtMiB2[ii]).second;
           }

           char Sigma[100];
           float sigma_eff;
           float s_sigma;
           float sigma_eff_sub;
           float s_sigma_sub;
           //float sub_wrtRm2 = 24.E-3;
           //float sub_wrtRm2_error = 2.E-3;

           sigma_eff = resFitAlt_ROC_wrtMiB2[ii]->GetParameter(2);
           s_sigma = resFitAlt_ROC_wrtMiB2[ii]->GetParError(2);
           if(doDoubleGauss == true){
             sigma_eff = fwhm/(2.*sqrt(2*log(2)));
             s_sigma = fwhm_error/(2*sqrt(2*log(2)));
           }

           if(iMCP != "Rm2" && iMCP != "Rm2"){
             sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtMiB2*sub_wrtMiB2);
             s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtRm2*sub_wrtRm2*sub_wrtRm2_error*sub_wrtRm2_error)/sigma_eff_sub;
           }else{
             sigma_eff_sub = sigma_eff; 
             s_sigma_sub = s_sigma;
           }
           sigma_eff = sigma_eff_sub;
           s_sigma = s_sigma_sub;

           std::cout << ii << " - sigma = " << sigma_eff*1000. << "+/-" << s_sigma*1000. << std::endl; 

           if(sigma_eff_sub*1000 < 10) continue;
           g_Res_vs_ROC_wrtMiB2->SetPoint(iPoint_MiB2,step*(ii+1),sigma_eff*1000.);
           g_Res_vs_ROC_wrtMiB2->SetPointError(iPoint_MiB2,step/2.,step/2.,s_sigma*1000.,s_sigma*1000.);

           iPoint_MiB2++;

           points_wrtMiB2.push_back(sigma_eff*1000.);
        
           sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

           TLatex *latexLabel = new TLatex();
           latexLabel->SetTextSize(0.05);
           latexLabel->SetNDC();
           latexLabel->SetTextFont(42); // helvetica

           resHist_ROC_wrtMiB2[ii]->GetXaxis()->SetTitle("t-t_{ref} (ns)");

           TCanvas* c2 = new TCanvas();
           c2->cd();
           resHist_ROC_wrtMiB2[ii]->Draw("hist");
           latexLabel->DrawLatex(0.72, 0.55,Sigma);
           resFitAlt_ROC_wrtMiB2[ii]->Draw("same");
           resFitAlt_ROC_single_wrtMiB2[ii]->Draw("same");
           resFitAlt_ROC_double_wrtMiB2[ii]->Draw("same");
           c2 -> Print(NameOutput_wrtMiB2_png,"png");
           c2 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c2;
        }

        if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
           float fwhm = 0.;
           float fwhm_error = 0.;

           h4->Draw(Selection17.c_str(),Selection16.c_str()); 
           char NameFitAlt_wrtRm2 [100];
           sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
           if(doDoubleGauss == false) resFitAlt_ROC_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-5*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+5*resHist_ROC_wrtRm2[ii]->GetRMS());  
           else resFitAlt_ROC_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-5*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+5*resHist_ROC_wrtRm2[ii]->GetRMS());  
           if(doDoubleGauss == true) resFitAlt_ROC_single_wrtRm2[ii] = new TF1((std::string(NameFitAlt_wrtRm2)+std::string("_single")).c_str(),"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-5*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+5*resHist_ROC_wrtRm2[ii]->GetRMS());
           if(doDoubleGauss == true) resFitAlt_ROC_double_wrtRm2[ii] = new TF1((std::string(NameFitAlt_wrtRm2)+std::string("_double")).c_str(),"[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-5*resHist_ROC_wrtRm2[ii]->GetRMS(),resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+5*resHist_ROC_wrtRm2[ii]->GetRMS());
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(0,resHist_ROC_wrtRm2[ii]->GetEntries()/2.);
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(1,0.);
           //resFitAlt_ROC_wrtRm2[ii]->SetParameters(2,0.04);
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(2,0.1);
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(3,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(3,resHist_ROC_wrtRm2[ii]->GetEntries()/2.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(4,0.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(5,0.1);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(5,0.04);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(6,0.);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(0,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(1,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           //resFitAlt_ROC_wrtRm2[ii]->SetParLimits(2,0.,0.2);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(2,0.,5.);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(3,0.,10.);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(3,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(4,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(5,0.04,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(5,0.,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(6,0.,10.);
           resHist_ROC_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(0,resFitAlt_ROC_wrtRm2[ii]->GetParameter(0));
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(1,resFitAlt_ROC_wrtRm2[ii]->GetParameter(1));
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(2,resFitAlt_ROC_wrtRm2[ii]->GetParameter(2));
           resFitAlt_ROC_wrtRm2[ii]->SetParameters(3,resFitAlt_ROC_wrtRm2[ii]->GetParameter(3));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(4,resFitAlt_ROC_wrtRm2[ii]->GetParameter(4));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(5,resFitAlt_ROC_wrtRm2[ii]->GetParameter(5));
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParameters(6,resFitAlt_ROC_wrtRm2[ii]->GetParameter(6));
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(0,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(1,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           //resFitAlt_ROC_wrtRm2[ii]->SetParLimits(2,0.,0.2);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(2,0.,5.);
           resFitAlt_ROC_wrtRm2[ii]->SetParLimits(3,0.,10.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(3,0.,resHist_ROC_wrtRm2[ii]->GetEntries());
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(4,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_ROC_wrtRm2[ii]->GetBinCenter(resHist_ROC_wrtRm2[ii]->GetMaximumBin())+0.2);
           //if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(5,0.04,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(5,0.,5.);
           if(doDoubleGauss == true) resFitAlt_ROC_wrtRm2[ii]->SetParLimits(6,0.,10.);
           resHist_ROC_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
           resFitAlt_ROC_wrtRm2[ii]->SetNpx(10000);
           
           if(doDoubleGauss == true){
              resFitAlt_ROC_single_wrtRm2[ii]->FixParameter(0,resFitAlt_ROC_wrtRm2[ii]->GetParameter(0));
              resFitAlt_ROC_single_wrtRm2[ii]->FixParameter(1,resFitAlt_ROC_wrtRm2[ii]->GetParameter(1));
              resFitAlt_ROC_single_wrtRm2[ii]->FixParameter(2,resFitAlt_ROC_wrtRm2[ii]->GetParameter(2));
              resFitAlt_ROC_single_wrtRm2[ii]->FixParameter(3,resFitAlt_ROC_wrtRm2[ii]->GetParameter(6));
              resFitAlt_ROC_single_wrtRm2[ii]->SetLineColor(kBlue+1);

              resFitAlt_ROC_double_wrtRm2[ii]->FixParameter(0,resFitAlt_ROC_wrtRm2[ii]->GetParameter(3));
              resFitAlt_ROC_double_wrtRm2[ii]->FixParameter(1,resFitAlt_ROC_wrtRm2[ii]->GetParameter(4));
              resFitAlt_ROC_double_wrtRm2[ii]->FixParameter(2,resFitAlt_ROC_wrtRm2[ii]->GetParameter(5));
              resFitAlt_ROC_double_wrtRm2[ii]->FixParameter(3,resFitAlt_ROC_wrtRm2[ii]->GetParameter(6));
              resFitAlt_ROC_double_wrtRm2[ii]->SetLineColor(kGreen+1);

             fwhm = FWHMSigma(resHist_ROC_wrtRm2[ii],resFitAlt_ROC_wrtRm2[ii]).first;
             fwhm_error = FWHMSigma(resHist_ROC_wrtRm2[ii],resFitAlt_ROC_wrtRm2[ii]).second;
           }

           char Sigma[100];
           float sigma_eff;
           float s_sigma;
           float sigma_eff_sub;
           float s_sigma_sub;
           //float sub_wrtRm2 = 24.E-3;
           //float sub_wrtRm2_error = 2.E-3;

           sigma_eff = resFitAlt_ROC_wrtRm2[ii]->GetParameter(2);
           s_sigma = resFitAlt_ROC_wrtRm2[ii]->GetParError(2);
           if(doDoubleGauss == true){
             sigma_eff = fwhm/(2.*sqrt(2*log(2)));
             s_sigma = fwhm_error/(2*sqrt(2*log(2)));
           }

           if(iMCP != "Rm2" && iMCP != "Rm2"){
             sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtRm2*sub_wrtRm2);
             s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtRm2*sub_wrtRm2*sub_wrtRm2_error*sub_wrtRm2_error)/sigma_eff_sub;
           }else{
             sigma_eff_sub = sigma_eff; 
             s_sigma_sub = s_sigma;
           }
           sigma_eff = sigma_eff_sub;
           s_sigma = s_sigma_sub;

           if(sigma_eff_sub*1000 < 10) continue;
           g_Res_vs_ROC_wrtRm2->SetPoint(iPoint_Rm2,step*(ii+1),sigma_eff*1000.);
           g_Res_vs_ROC_wrtRm2->SetPointError(iPoint_Rm2,step/2.,step/2.,s_sigma*1000.,s_sigma*1000.);

           iPoint_Rm2++;

           points_wrtRm2.push_back(sigma_eff*1000.);
        
           sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

           TLatex *latexLabel = new TLatex();
           latexLabel->SetTextSize(0.05);
           latexLabel->SetNDC();
           latexLabel->SetTextFont(42); // helvetica

           resHist_ROC_wrtRm2[ii]->GetXaxis()->SetTitle("t-t_{ref} (ns)");

           TCanvas* c2 = new TCanvas();
           c2->cd();
           resHist_ROC_wrtRm2[ii]->Draw("hist");
           latexLabel->DrawLatex(0.72, 0.55,Sigma);
           resFitAlt_ROC_wrtRm2[ii]->Draw("same");
           resFitAlt_ROC_single_wrtRm2[ii]->Draw("same");
           resFitAlt_ROC_double_wrtRm2[ii]->Draw("same");
           c2 -> Print(NameOutput_wrtRm2_png,"png");
           c2 -> Print(NameOutput_wrtRm2_pdf,"pdf");
           delete c2;
        }
    }
    output_ampMax->Close();

    if(iMCP != "MiB2") std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false) std::sort(points_wrtRm2.begin(),points_wrtRm2.end());

        
    if(iMCP != "MiB2"){
        g_Res_vs_ROC_wrtMiB2->GetXaxis()->SetTitle("#epsilon");
        g_Res_vs_ROC_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
        g_Res_vs_ROC_wrtMiB2->SetMarkerStyle(20);
        g_Res_vs_ROC_wrtMiB2->SetMarkerSize(0.7);
        g_Res_vs_ROC_wrtMiB2->SetMarkerColor(kBlack);
        g_Res_vs_ROC_wrtMiB2->SetLineColor(kBlack); 
        g_Res_vs_ROC_wrtMiB2->GetXaxis()->SetRangeUser(0.,100.);
        g_Res_vs_ROC_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-2.,points_wrtMiB2.at(points_wrtMiB2.size()-1)+2.);

        TCanvas* c9 = new TCanvas();
        c9->cd();
        c9->SetGrid();
        g_Res_vs_ROC_wrtMiB2->Draw("AP");
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
    }   

    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
        g_Res_vs_ROC_wrtRm2->GetXaxis()->SetTitle("#epsilon");
        g_Res_vs_ROC_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
        g_Res_vs_ROC_wrtRm2->SetMarkerStyle(20);
        g_Res_vs_ROC_wrtRm2->SetMarkerSize(0.7);
        g_Res_vs_ROC_wrtRm2->SetMarkerColor(kBlack);
        g_Res_vs_ROC_wrtRm2->SetLineColor(kBlack); 
        g_Res_vs_ROC_wrtRm2->GetXaxis()->SetRangeUser(0.,100.);
        g_Res_vs_ROC_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-2.,points_wrtRm2.at(points_wrtRm2.size()-1)+2.);

        TCanvas* c9 = new TCanvas();
        c9->cd();
        c9->SetGrid();
        g_Res_vs_ROC_wrtRm2->Draw("AP");
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
        c9 -> Print(std::string("TimeResolution_vs_eff_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
    }   
     
    
    /*g_Res_vs_ROC_wrtMiB2->Write(std::string("TimeResolution_vs_eff_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP).c_str());
    time_wrtMiB2->Write();
    if(doOnlyWrtMiB2 == false) g_Res_vs_ROC_wrtRm2->Write(std::string("TimeResolution_vs_eff_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP).c_str());
    if(doOnlyWrtMiB2 == false) time_wrtRm2->Write(); 
    output_ampMax->Close();*/

   } 
}

void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
   
    std::string Selection1;   
    char Selection2 [1000];

    sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection1 = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,"0.",false);
      Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection1,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
      }else{
         Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,"0.",false);
         Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
         Selection1 = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
      }   
    } 
    std::string Selection3 = Selection1+" && "+Selection2+HodoSelection; 

    std::cout << "Selection = " << Selection3 << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection3);
 
    setBins(h4, std::string("amp_max["+iMCP+"]"), Selection3, maxMCP, 10., ptr_minVec, nBins);
    setBinsFraction(h4, std::string("amp_max["+iMCP+"]"), Selection3, maxMCP, 1., ptr_minVecFrac, nBinsFrac);
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,ampMin,4000,timeMin);
    TH2F* timingCorrection_wrtRm2 = new TH2F("timingCorrection_wrtRm2","",nBins,ampMin,4000,timeMin);

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

       TF1* fit_corr1;
       //fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,atof(maxMCP.c_str()));
       //fit_corr1->SetParLimits(2,0.,999999999.);
       //fit_corr1 = new TF1("fit_corr1","[0]+[1]*log(x+[2])",0.,atof(maxMCP.c_str()));
       fit_corr1 = new TF1("fit_corr1","[0]+[1]*exp([2]*x)",0.,atof(maxMCP.c_str()));
       //fit_corr1 = new TF1("fit_corr1","sqrt([0]+[1]*x+[2]*x*x)",0.,atof(maxMCP.c_str()));
       //fit_corr1 = new TF1("fit_corr1","pol1",0.,atof(maxMCP.c_str()));
       timingCorrection_wrtMiB2_1->Fit("fit_corr1","B");
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));

       std::cout << "Time correction = " << Params_wrtMiB2->at(0) << " " << Params_wrtMiB2->at(1) << " " << Params_wrtMiB2->at(2) << std::endl;

       for(int ii = 1; ii <= timingCorrection_wrtMiB2_1->GetNbinsX(); ii++)
       {
           if(timingCorrection_wrtMiB2_1->GetBinContent(ii)!= 0) points_wrtMiB2.push_back(timingCorrection_wrtMiB2_1->GetBinContent(ii)-fit_corr1->GetParameter(0));
           timingCorrection_wrtMiB2_1->SetBinContent(ii,timingCorrection_wrtMiB2_1->GetBinContent(ii)-fit_corr1->GetParameter(0));
       }
       std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());

       timingCorrection_wrtMiB2_1->SetAxisRange(points_wrtMiB2.at(0)-0.2,points_wrtMiB2.at(points_wrtMiB2.size()-1)+0.2, "Y");
       timingCorrection_wrtMiB2_1->SetMarkerStyle(20);
       timingCorrection_wrtMiB2_1->SetMarkerSize(0.9);
       timingCorrection_wrtMiB2_1->SetMarkerColor(kBlack);
       timingCorrection_wrtMiB2_1->SetLineColor(kBlack);

       fit_corr1->FixParameter(0,0.);
    
       std::string wrtMCP = "";
       if(doOnlyWrtMiB2) wrtMCP = "_onlyWrtMiB2";

       TCanvas* c1 = new TCanvas();
       c1->cd();
       timingCorrection_wrtMiB2_1->Draw();
       //timingCorrection_wrtMiB2->Draw();
       fit_corr1->Draw("same");
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");

       /*TFile* file_timingCorrection_wrtMiB2 = new TFile(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".root").c_str(),"RECREATE");
       file_timingCorrection_wrtMiB2->cd();
       timingCorrection_wrtMiB2_1->Write("timingCorrection_wrtMiB2");
       fit_corr1->Write("timingCorrection_wrtMiB2_f");
       file_timingCorrection_wrtMiB2->Close();*/
    }
 
    if(iMCP != "Rm2" && doOnlyWrtMiB2 == false){
       std::vector<float> points_wrtRm2;
       timingCorrection_wrtRm2->FitSlicesY();
       TH1F* timingCorrection_wrtRm2_1 = (TH1F*)inputFile->Get("timingCorrection_wrtRm2_1");
       timingCorrection_wrtRm2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
       timingCorrection_wrtRm2_1->GetYaxis()->SetTitle("time-time[Rm2]");

       TF1* fit_corr1;
       //fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,atof(maxMCP.c_str()));
       //fit_corr1->SetParLimits(2,0.,999999999.);
       //fit_corr1 = new TF1("fit_corr1","[0]+[1]*log(x+[2])",0.,atof(maxMCP.c_str()));
       //fit_corr1 = new TF1("fit_corr1","sqrt([0]+[1]*x+[2]x*x)",0.,atof(maxMCP.c_str()));
       fit_corr1 = new TF1("fit_corr1","pol1",0.,atof(maxMCP.c_str()));
       timingCorrection_wrtRm2_1->Fit("fit_corr1");
       Params_wrtRm2->push_back(fit_corr1->GetParameter(0));
       Params_wrtRm2->push_back(fit_corr1->GetParameter(1));
       Params_wrtRm2->push_back(fit_corr1->GetParameter(2));

       for(int ii = 1; ii <= timingCorrection_wrtRm2_1->GetNbinsX(); ii++)
       {
           if(timingCorrection_wrtRm2_1->GetBinContent(ii)!= 0) points_wrtRm2.push_back(timingCorrection_wrtRm2_1->GetBinContent(ii)-fit_corr1->GetParameter(0));
           timingCorrection_wrtRm2_1->SetBinContent(ii,timingCorrection_wrtRm2_1->GetBinContent(ii)-fit_corr1->GetParameter(0));
       }
       std::sort(points_wrtRm2.begin(),points_wrtRm2.end());

       timingCorrection_wrtRm2_1->SetAxisRange(points_wrtRm2.at(0)-0.2,points_wrtRm2.at(points_wrtRm2.size()-1)+0.2, "Y");
       timingCorrection_wrtRm2_1->SetMarkerStyle(20);
       timingCorrection_wrtRm2_1->SetMarkerSize(0.9);
       timingCorrection_wrtRm2_1->SetMarkerColor(kBlack);
       timingCorrection_wrtRm2_1->SetLineColor(kBlack);

       //fit_corr1->FixParameter(0,0.);
    
       TCanvas* c2 = new TCanvas();
       c2->cd();
       timingCorrection_wrtRm2_1->Draw();
       //timingCorrection_wrtRm2->Draw();
       fit_corr1->Draw("same");
       c2 -> Print(std::string("timingCorrection_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c2 -> Print(std::string("timingCorrection_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");

       /*TFile* file_timingCorrection_wrtRm2 = new TFile(std::string("timingCorrection_wrtRm2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".root").c_str(),"RECREATE");
       file_timingCorrection_wrtRm2->cd();
       timingCorrection_wrtRm2_1->Write("timingCorrection_wrtRm2");
       fit_corr1->Write("timingCorrection_wrtRm2_f");
       file_timingCorrection_wrtRm2->Close();*/
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
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }   
    } 

    Selection = Selection + HodoSelection;

    std::cout << "Selection = " << Selection << std::endl;

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
    fit_corr_max->Draw("same");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
}

void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP, std::string HodoSelection, bool doOnlyWrtMiB2)
{
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",300,-10,20,300,-1.,1.5,100.,3000.);
    TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time","",300,-10,20,300,-1.,1.5);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;
  
    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }   
    } 

    Selection = Selection+" && WF_ch == "+iMCP;//+HodoSelection;

    //std::cout << "Selection = " << Selection << std::endl;

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
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
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
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
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
    TH1F* h_ampMax_2part = new TH1F(std::string("h_ampMax_"+nameiMCP+"_2part").c_str(),"",450,0.,4500.);
    TH1F* h_ampMax_3part = new TH1F(std::string("h_ampMax_"+nameiMCP+"_3part").c_str(),"",450,0.,4500.);

    std::string Selection;
    std::string Selection_2part;
    std::string Selection_3part;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
      Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }else{
         Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>"+scintMin+" && adc_data[scint]<"+scintMax+" && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,"0.",false);
         Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
         Selection = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
      }   
    } 
    
    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection_2part = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>700 && adc_data[scint]<1500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection_2part = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection_2part,"0.",false);
      Selection_2part = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection_2part,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection_2part = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>700 && adc_data[scint]<1500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection_2part = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection_2part,"1",true);
         Selection_2part = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection_2part,"1",true);
      }else{
         Selection_2part = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>700 && adc_data[scint]<1500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection_2part = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection_2part,"0.",false);
         Selection_2part = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection_2part,"1",true);
         Selection_2part = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection_2part,"1",true);
         Selection_2part = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection_2part,"1",true);
         Selection_2part = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection_2part,"1",true);
      }   
    } 

    if(iMCP == "MiB2" || iMCP == "Rm2"){
      Selection_3part = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>1500 && adc_data[scint]<2500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
      Selection_3part = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection_3part,"0.",false);
      Selection_3part = AddSelection(h4,std::string("time_max[MiB2]-time_max[Rm2]"),Selection_3part,"0.",false);
    }else{ 
      if(doOnlyWrtMiB2){
         Selection_3part = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && adc_data[scint]>1500 && adc_data[scint]<2500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection_3part = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection_3part,"1",true);
         Selection_3part = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection_3part,"1",true);
      }else{
         Selection_3part = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && fabs(time_max[Rm2])<"+time_max_Rm2+" && amp_max[Rm2]< 1200 && adc_data[scint]>1500 && adc_data[scint]<2500 && n_hitsX> 0 && n_hitsX<3 && n_hitsY> 0 && n_hitsY<3";
         Selection_3part = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection_3part,"0.",false);
         Selection_3part = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection_3part,"1",true);
         Selection_3part = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection_3part,"1",true);
         Selection_3part = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection_3part,"1",true);
         Selection_3part = AddSelection(h4,std::string("time_max[Rm2]-time_max[")+iMCP+std::string("]"),Selection_3part,"1",true);
      }   
    } 

    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"_noSelection").c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP).c_str(),Selection.c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"_2part").c_str(),Selection_2part.c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"_3part").c_str(),Selection_3part.c_str());
    
    h_ampMax_noSelection->GetXaxis()->SetTitle("amp_max");
    h_ampMax->GetXaxis()->SetTitle("amp_max");
    h_ampMax_2part->GetXaxis()->SetTitle("amp_max");
    h_ampMax_3part->GetXaxis()->SetTitle("amp_max");
    h_ampMax_noSelection->GetYaxis()->SetTitle("Events");
    h_ampMax->GetYaxis()->SetTitle("Events");
    h_ampMax_2part->GetYaxis()->SetTitle("Events");
    h_ampMax_3part->GetYaxis()->SetTitle("Events");
 
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

    TCanvas* c3 = new TCanvas();
    c3->cd();
    c3->SetLogy();
    h_ampMax_2part->Draw("H");
    c3 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+"_2part.png").c_str(),"png");
    c3 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+"_2part.pdf").c_str(),"pdf");

    TCanvas* c4 = new TCanvas();
    c4->cd();
    c4->SetLogy();
    h_ampMax_3part->Draw("H");
    c4 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+"_3part.png").c_str(),"png");
    c4 -> Print(std::string("ampMax_"+nameiMCP+wrtMCP+"_3part.pdf").c_str(),"pdf");
    
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
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") == std::string::npos) HV = "HVBINP2NEG";


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

void CheckSelectionEfficiency(TTree* h4, std::string iMCP, std::string Selection)
{
    //std::cout << "String = " << Selection << std::endl;
    TH1F* h = new TH1F("h","",10000,-99999999.,99999999);
    std::string splittedSelection = "";
    h4->Draw(std::string("time["+iMCP+"] >> h").c_str(),splittedSelection.c_str());
    int Ntot = h->GetEntries();
    for(unsigned int ii = 0; ii < split(Selection, std::string("&&")).size(); ii++){ 
        if(ii == 0) splittedSelection = split(Selection, std::string("&&")).at(ii);
        else splittedSelection = splittedSelection + "&" + split(Selection, std::string("&&")).at(ii);
        h4->Draw(std::string("time["+iMCP+"] >> h").c_str(),splittedSelection.c_str());
        std::cout << "Selection = " << splittedSelection << " - " << float(h->GetEntries())/float(Ntot)*100. << "%" << std::endl;
    }  
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
    else if((iMCP == "BINP4" || iMCP == "MiB2" || iMCP == "Rm2") && inputs.find("2378") == std::string::npos) HV = "HVBINP2NEG";


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

void setBins(TTree* h4, std::string var, std::string Selection, std::string maxMCP, float step, float* ptr_minVec, int nbins)
{
    TH1F* h_var = new TH1F("h_var","",9999999,-99999.,99999);
    
    h4->Draw(std::string(var+">>h_var").c_str(),Selection.c_str());   
    float Total_Events = h_var->Integral(h_var->FindBin(0.),h_var->FindBin(atof(maxMCP.c_str())));
    float nEvents = (float)Total_Events/(float)nbins; 
    
    float binMax = 0.;
    float binMin = 0.;
    for(int ii = 0; ii < nbins; ii++){
      do{
        binMax+= step; 
        char min[100];
        sprintf(min,"%f",binMin);
        char max[100];
        sprintf(max,"%f",binMax);
        //h4->Draw(std::string(var+">>h_var").c_str(),std::string(Selection+ " && "+var+">"+min+" && "+var+"<"+max).c_str());  
        //std::cout << ii << " - binMin " << min << " binMax = " << max << " " << atof(maxMCP.c_str()) << " " << h_var->Integral(h_var->FindBin(binMin),h_var->FindBin(binMax)) << " " << nEvents << std::endl;   
        if(binMax >= atof(maxMCP.c_str())) break; 
      }while(h_var->Integral(h_var->FindBin(binMin),h_var->FindBin(binMax))<nEvents);
      *(ptr_minVec+ii+1) = binMax;
      binMin = binMax;
    }
    for(int ii = 0; ii <= nbins; ii++)
      std::cout << ii << " " << *(ptr_minVec+ii) << std::endl;
}

void setBinsFraction(TTree* h4, std::string var, std::string Selection, std::string maxMCP, float step, float* ptr_minVec, int nbins)
{
    TH1F* h_var = new TH1F("h_var","",9999999,-99999.,99999);
    
    h4->Draw(std::string(var+">>h_var").c_str(),Selection.c_str());   
    float Total_Events = h_var->Integral(h_var->FindBin(0.),h_var->FindBin(atof(maxMCP.c_str())));
    float nEvents = (float)Total_Events/(float)nbins; 
    
    float binMax = 0.;
    float binMin = 0.;
    for(int ii = 1; ii < nbins; ii++){
      do{
        binMax+= step; 
        char min[100];
        sprintf(min,"%f",binMin);
        char max[100];
        sprintf(max,"%f",binMax);
        if(binMax >= atof(maxMCP.c_str())) break; 
      }while(h_var->Integral(h_var->FindBin(binMin),h_var->FindBin(binMax))<nEvents*ii);
      if(ii != nbins) *(ptr_minVec+ii) = binMax;
      else *(ptr_minVec+ii) = atof(maxMCP.c_str());
      binMin = 0.;
    }
    for(int ii = 0; ii < nbins; ii++)
      std::cout << ii << " " << *(ptr_minVec+ii) << std::endl;
}

std::vector<std::string> split(const std::string text, std::string sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep.c_str(), start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

std::pair<float,float> FWHMSigma(TH1F* h,TF1* f)
{
   float fwhm = 0.;
   float fwhm_up = 0.;
   float fwhm_down = 0.;
   
   TH1F* time_FWHM = new TH1F("time_FWHM","",40000,-10.,10.);
   time_FWHM->FillRandom(f->GetName(),10000000);
   
   int bin1 = time_FWHM->FindFirstBinAbove(time_FWHM->GetMaximum()/2);
   int bin2 = time_FWHM->FindLastBinAbove(time_FWHM->GetMaximum()/2);
   fwhm = time_FWHM->GetBinCenter(bin2) - time_FWHM->GetBinCenter(bin1);
   
   /*int bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
   int bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
   fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);*/
   
   TF1* f_new = new TF1("f_new","([0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6])",-10,10);   
   f_new->FixParameter(0,f->GetParameter(0));  
   f_new->FixParameter(1,f->GetParameter(1));  
   f_new->FixParameter(2,f->GetParameter(2)+f->GetParError(2));  
   f_new->FixParameter(3,f->GetParameter(3));  
   f_new->FixParameter(4,f->GetParameter(4));
   f_new->FixParameter(5,f->GetParameter(5)+f->GetParError(5));  
   f_new->FixParameter(6,f->GetParameter(6));   
   
   delete time_FWHM;
    
   time_FWHM = new TH1F("time_FWHM","",40000,-10.,10.);
   time_FWHM->FillRandom("f_new",10000000);
   bin1 = time_FWHM->FindFirstBinAbove(time_FWHM->GetMaximum()/2);
   bin2 = time_FWHM->FindLastBinAbove(time_FWHM->GetMaximum()/2);
   fwhm_up = time_FWHM->GetBinCenter(bin2) - time_FWHM->GetBinCenter(bin1);

   f_new->FixParameter(0,f->GetParameter(0));  
   f_new->FixParameter(1,f->GetParameter(1));  
   f_new->FixParameter(2,f->GetParameter(2)-f->GetParError(2));  
   f_new->FixParameter(3,f->GetParameter(3));  
   f_new->FixParameter(4,f->GetParameter(4));
   f_new->FixParameter(5,f->GetParameter(5)-f->GetParError(5));  
   f_new->FixParameter(6,f->GetParameter(6));   
 
   delete time_FWHM;
    
   time_FWHM = new TH1F("time_FWHM","",40000,-10.,10.);
   time_FWHM->FillRandom("f_new",10000000);
   bin1 = time_FWHM->FindFirstBinAbove(time_FWHM->GetMaximum()/2);
   bin2 = time_FWHM->FindLastBinAbove(time_FWHM->GetMaximum()/2);
   fwhm_down = time_FWHM->GetBinCenter(bin2) - time_FWHM->GetBinCenter(bin1);

   delete time_FWHM;
   
   std::pair<float,float> outpair = std::make_pair(fwhm,fabs(fwhm_up-fwhm_down)/2.);
   return outpair;
}

std::pair<float,float> EffectiveSigma(TF1* f, int nSteps, float fraction)
{
   float integral = f->Integral(-10.,10.);
   float mean = f->GetMaximumX(-10.,10.);
   
   float sigma_low = 0.;
   for(int iStep = 0; iStep < nSteps; iStep++){
       sigma_low = mean-0.0001*iStep; 
       if(f->Integral(mean-0.0001*iStep,mean)/integral >= fraction) break;
   }
   
   float sigma_high = 0.;
   for(int iStep = 0; iStep < nSteps; iStep++){
       sigma_high = mean+0.0001*iStep; 
       if(f->Integral(mean,mean+0.0001*iStep)/integral >= fraction) break;
   }
   
   std::pair<float,float> outpair = std::make_pair (fabs(sigma_low-sigma_high)/2,fabs(sigma_low-sigma_high)/2);
   return outpair ;
}

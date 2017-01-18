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
#include "TPaveStats.h"

#include<iostream>
#include<string>
#include<fstream>

int nBins = 12;
float minVec[13] = {0.};
float *ptr_minVec = minVec;
//ptr_minVec = minVec; 
float timeMin[4001];
float ampMin[3001];

//BINP
std::string minMiB2 = "200.";
std::string maxMiB2 = "3000.";
std::string time_max_MiB2 = "150.";
std::string time_max_Rm2 = "150.";
std::string amp_max_Rm2 = "200.";
std::string scintMin = "200.";
std::string scintMax = "700.";
std::string timeChi2 = "99999999.";
std::string hodoX = "-999";
std::string hodoY = "-1000";
bool isSaturated = true;
float sShift=0.;

void FinalTiming(TTree* h4,TTree* digi,TTree* hodo, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, bool doDoubleGauss, std::string* ptr_timeShifted);
void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, std::string* ptr_timeShifted);
void SaturationCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP);
void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP);
void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void CheckEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
std::string AddSelection(TTree*, std::string, std::string, std::string, bool);
std::string shiftVar(TTree* h4, std::string Var, std::string Selection);
std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel);
void setBins(TTree* h4, std::string var, std::string Selection, std::string maxMCP,float step, float* ptr_minVec);
void CheckSelectionEfficiency(TTree* h4, std::string iMCP, std::string Selection);
int getHV(TTree* h4, std::string HV);
std::vector<std::string> split(const std::string text, std::string sep);


void measureTiming_H42014(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, std::string maxMCP, bool doFirstStep = true, bool doPulseShapes = false, bool doDoubleGauss = false)
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open(inputs.c_str());
    TTree* h4 = (TTree*)inputFile->Get("h4");
    TTree* digi = (TTree*)inputFile->Get("digi");
    TTree* hodo = (TTree*)inputFile->Get("hodo");
    
    std::vector<float>* Params = new std::vector<float>;
    std::vector<float>* Params_wrtMiB2 = new std::vector<float>;
    
    for(int ii = 0; ii < 4001; ii++)
        timeMin[ii] = 0.01*ii-20.;
    for(int ii = 0; ii < 3001; ii++)
        ampMin[ii] = ii;

    std::string nameiMCP = iMCP;
    std::string timeShifted;
    std::string* ptr_timeShifted = &timeShifted;
    
    if(doFirstStep == true){
       TimeChi2(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       AmpMax(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       Hodoscope(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       CheckEfficiency(h4, inputs, iMCP, nameiMCP, thresMCP, maxMCP);
    }else if(doPulseShapes == true) PulseShapes(h4, iMCP, nameiMCP, thresMCP, maxMCP);
    else{
       //SaturationCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP);
       if(!isSaturated) AmpVsTime_Selection(h4, iMCP, nameiMCP,Timing, Params, thresMCP, maxMCP);
       TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, ptr_timeShifted);
       FinalTiming(h4,digi,hodo, inputs, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, doDoubleGauss, ptr_timeShifted);
    }
}

void FinalTiming(TTree* h4,TTree* digi,TTree* hodo, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, bool doDoubleGauss, std::string* ptr_timeShifted)
{
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",400,-1.,1.);
    TH2F* time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins,minVec,4000,timeMin);

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    TF1 *g_res;

    // Declaration of leaf types
    int           SEE;
    int           ZS1;
    int           ZS2;
    int           MiB3;
    int           MiB2;
    //int           void;
    int           CFD;
    int           LED50;
    int           LED100;
    int           LED150;
    int           LED200;
    int           LED300;
    int           LED400;
    int           LED500;
    int           LED600;
    ULong64_t       index;
    unsigned int          n_channels;
    unsigned int          n_timetypes;
    int           n_planes;
    int           n_hitsX;
    int           n_hitsY;
    float         X[2];   //[n_planes]
    float         Y[2];   //[n_planes]
    float         b_charge[6];   //[n_channels]
    float         b_slope[6];   //[n_channels]
    float         b_rms[6];   //[n_channels]
    float         time[54];   //[n_timetypes]
    float         time_mirror[54];   //[n_timetypes]
    float         time_chi2[54];   //[n_timetypes]
    float         maximum[6];   //[n_channels]
    float         time_maximum[6];   //[n_channels]
    float         amp_max[6];   //[n_channels]
    float         time_max[6];   //[n_channels]
    float         chi2_max[6];   //[n_channels]
    float         charge_tot[6];   //[n_channels]
    float         charge_sig[6];   //[n_channels]
    float         fit_ampl[6];   //[n_channels]
    float         fit_time[6];   //[n_channels]
    float         fit_chi2[6];   //[n_channels]
    float         calibration[6];   //[n_channels]

    // List of branches
    TBranch        *b_SEE;   //!
    TBranch        *b_ZS1;   //!
    TBranch        *b_ZS2;   //!
    TBranch        *b_MiB3;   //!
    TBranch        *b_MiB2;   //!
    TBranch        *b_void;   //!
    TBranch        *b_CFD;   //!
    TBranch        *b_LED50;   //!
    TBranch        *b_LED100;   //!
    TBranch        *b_LED150;   //!
    TBranch        *b_LED200;   //!
    TBranch        *b_LED300;   //!
    TBranch        *b_LED400;   //!
    TBranch        *b_LED500;   //!
    TBranch        *b_LED600;   //!
    TBranch        *b_index;   //!
    TBranch        *b_n_channels;   //!
    TBranch        *b_n_timetypes;   //!
    TBranch        *b_b_charge;   //!
    TBranch        *b_b_slope;   //!
    TBranch        *b_b_rms;   //!
    TBranch        *b_time;   //!
    TBranch        *b_time_mirror;   //!
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
    TBranch        *b_n_planes;   //!
    TBranch        *b_n_hitsX;   //!
    TBranch        *b_n_hitsY;   //!
    TBranch        *b_X;   //!
    TBranch        *b_Y;   //!

    digi->SetBranchAddress("SEE", &SEE, &b_SEE);
    digi->SetBranchAddress("ZS1", &ZS1, &b_ZS1);
    digi->SetBranchAddress("ZS2", &ZS2, &b_ZS2);
    digi->SetBranchAddress("MiB3", &MiB3, &b_MiB3);
    digi->SetBranchAddress("MiB2", &MiB2, &b_MiB2);
    //digi->SetBranchAddress("void", &void, &b_void);
    digi->SetBranchAddress("CFD", &CFD, &b_CFD);
    digi->SetBranchAddress("LED50", &LED50, &b_LED50);
    digi->SetBranchAddress("LED100", &LED100, &b_LED100);
    digi->SetBranchAddress("LED150", &LED150, &b_LED150);
    digi->SetBranchAddress("LED200", &LED200, &b_LED200);
    digi->SetBranchAddress("LED300", &LED300, &b_LED300);
    digi->SetBranchAddress("LED400", &LED400, &b_LED400);
    digi->SetBranchAddress("LED500", &LED500, &b_LED500);
    digi->SetBranchAddress("LED600", &LED600, &b_LED600);
    digi->SetBranchAddress("index", &index, &b_index);
    digi->SetBranchAddress("n_channels", &n_channels, &b_n_channels);
    digi->SetBranchAddress("n_timetypes", &n_timetypes, &b_n_timetypes);
    digi->SetBranchAddress("b_charge", b_charge, &b_b_charge);
    digi->SetBranchAddress("b_slope", b_slope, &b_b_slope);
    digi->SetBranchAddress("b_rms", b_rms, &b_b_rms);
    digi->SetBranchAddress("time", time, &b_time);
    digi->SetBranchAddress("time_mirror", time_mirror, &b_time_mirror);
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
    hodo->SetBranchAddress("index", &index, &b_index);
    hodo->SetBranchAddress("n_planes", &n_planes, &b_n_planes);
    hodo->SetBranchAddress("n_hitsX", &n_hitsX, &b_n_hitsX);
    hodo->SetBranchAddress("n_hitsY", &n_hitsY, &b_n_hitsY);
    hodo->SetBranchAddress("X", X, &b_X);
    hodo->SetBranchAddress("Y", Y, &b_Y);

    std::string Selection;
    Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
    Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
    
    for(int entry = 0; entry < h4->GetEntries(); entry++){

        if(entry%1000==0) std::cout<<"--- Reading entry = "<<entry<<std::endl;
        h4->GetEntry(entry);
        digi->GetEntry(entry);
        hodo->GetEntry(entry);

        int imcp = 0;
        if(iMCP == "SEE") imcp = SEE;
        else if(iMCP == "ZS1") imcp = ZS1;
        else if(iMCP == "ZS2") imcp = ZS2;
        else if(iMCP == "MiB3") imcp = MiB3;
        else if(iMCP == "MiB2") imcp = MiB2;

        int timing = 0;
        if(Timing == "CFD50") timing = CFD; 
        else if(Timing == "LED50") timing = LED50; 
        else if(Timing == "LED100") timing = LED100; 
        else if(Timing == "LED150") timing = LED150; 
        else if(Timing == "LED200") timing = LED200; 
        else if(Timing == "LED300") timing = LED300; 
        else if(Timing == "LED400") timing = LED400; 
        else if(Timing == "LED500") timing = LED500; 
        else if(Timing == "LED600") timing = LED600; 
        
        //if(amp_max[imcp]<atof(thresMCP.c_str()) || amp_max[imcp]> atof(maxMCP.c_str()) || amp_max[MiB2]<atof(minMiB2.c_str()) || amp_max[MiB2]> atof(maxMiB2.c_str()) || X[0]<atoi(hodoX.c_str()) || Y[0]<atoi(hodoY.c_str())) continue;

        if(amp_max[imcp]<atof(thresMCP.c_str())) continue;
        if(amp_max[imcp]> atof(maxMCP.c_str())) continue;
        if(amp_max[MiB2]<atof(minMiB2.c_str())) continue;
        if(amp_max[MiB2]> atof(maxMiB2.c_str())) continue;
        if(X[0]<atoi(hodoX.c_str()) || Y[0]<atoi(hodoY.c_str())) continue;
        if(fabs(time[MiB2]-time[imcp+timing]-(sShift))>1) continue;
        if(!isSaturated && fabs(time_max[imcp]-time[imcp+timing]-(Params->at(0)+(Params->at(1))*log(Params->at(2)+amp_max[imcp])))>0.25) continue;   //(time-time_max) vs amp selection
        
    
        //time correction
        float time_corr = 0.;
        float amp_corr = amp_max[imcp];
        float deltaCross = time_mirror[imcp+timing]-time[imcp+timing];
        if(iMCP == "ZS1" && isSaturated) amp_corr = 488.8 + 53.87*deltaCross + 2.998*deltaCross*deltaCross + 8.715*deltaCross*deltaCross*deltaCross - 0.0308*deltaCross*deltaCross*deltaCross*deltaCross;
        if(iMCP == "ZS2" && isSaturated) amp_corr = 493.7 + 32.73*deltaCross + 61.56*deltaCross*deltaCross - 27.54*deltaCross*deltaCross*deltaCross + 7.141*deltaCross*deltaCross*deltaCross*deltaCross;
        if(isSaturated) time_corr = Params_wrtMiB2->at(0) + Params_wrtMiB2->at(1)*amp_corr + Params_wrtMiB2->at(2)*amp_corr*amp_corr;
        else time_corr = Params_wrtMiB2->at(0) + Params_wrtMiB2->at(1)/(amp_corr + Params_wrtMiB2->at(2)); 

        time_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr);
        time_vs_amp_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr,amp_corr);
    }
    
    if(time_wrtMiB2->GetEntries() < 2000) time_wrtMiB2->Rebin(8);
    
    float sub_wrtMiB2 = 24.E-3;
    float sub_wrtMiB2_error = 2.E-3;
    if(doDoubleGauss == false) g_res = new TF1("g_res","gaus",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-1*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+1*time_wrtMiB2->GetRMS());
    else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-1*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+1*time_wrtMiB2->GetRMS());  
    //else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)",-1.,1.);  
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
    float sigma_eff_sub;
    float s_sigma_sub;
    if(doDoubleGauss == true){
       float f1 = g_res->GetParameter(0)/(g_res->GetParameter(0)+g_res->GetParameter(3));
       float f2 = g_res->GetParameter(3)/(g_res->GetParameter(0)+g_res->GetParameter(3));
       sigma_eff = sqrt(f1*g_res->GetParameter(2)*g_res->GetParameter(2)+f2*g_res->GetParameter(5)*g_res->GetParameter(5) + f1*g_res->GetParameter(1)*g_res->GetParameter(1)+f2*g_res->GetParameter(4)*g_res->GetParameter(4) - (f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4))*(f1*g_res->GetParameter(1)+f2*g_res->GetParameter(4)));
       s_sigma = sqrt((f1*f1*g_res->GetParameter(2)*g_res->GetParameter(2)*g_res->GetParError(2)*g_res->GetParError(2)+f2*f2*g_res->GetParameter(4)*g_res->GetParameter(4)*g_res->GetParError(4)*g_res->GetParError(4)/(sigma_eff*sigma_eff))); 
    }else{
       sigma_eff = g_res->GetParameter(2); 
       s_sigma = g_res->GetParError(2);
    }

    if(iMCP != "MiB2" && !isSaturated){
       sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtMiB2*sub_wrtMiB2);
       s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtMiB2*sub_wrtMiB2*sub_wrtMiB2_error*sub_wrtMiB2_error)/sigma_eff_sub;
    }else{
       sigma_eff_sub = g_res->GetParameter(2); 
       s_sigma_sub = g_res->GetParError(2);
    }
    sigma_eff = sigma_eff_sub;
    s_sigma = s_sigma_sub;

    sprintf (Sigma,"#sigma = %.0f+/-%.0f ps",sigma_eff*1000.,s_sigma*1000.);

    latexLabel->SetTextSize(0.05);
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica

    time_wrtMiB2->GetXaxis()->SetTitle("t-t_{ref} (ns)");

    TCanvas* c1 = new TCanvas();
    c1->cd();
    time_wrtMiB2->Draw("hist");
    latexLabel->DrawLatex(0.72, 0.55,Sigma);
    g_res->Draw("same");
    if(iMCP != "MiB2"){
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }else{
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }
    delete g_res; 

    int HV = getHV(h4, std::string("HV"+iMCP));

    std::ofstream res_txt;
    res_txt.open("res_vs_HV.dat", ios::out | ios::app);    
    if (res_txt.is_open())
    {
        res_txt << HV << " " << sigma_eff*1000 << " " << s_sigma*1000 << "\n";
        res_txt.close();
    }else{
        ofstream res_txt_new("res_vs_HV.dat");
        res_txt_new << HV << " " << sigma_eff*1000 << " " << s_sigma*1000 << "\n";
        res_txt_new.close();
    }
}

void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, std::string* ptr_timeShifted)
{    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
   
    std::string Selection1;   
    char Selection2 [1000];

    if(!isSaturated) sprintf (Selection2,std::string("fabs(time_max["+iMCP+"]-time["+iMCP+iTiming+"]-(%f+(%f)*log(%f+amp_max["+iMCP+"])))<0.25").c_str(),Params->at(0),Params->at(1),Params->at(2));
   
    if(iMCP == "MiB2"){
      Selection1 = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
    } 
    std::string Selection3;
    if(!isSaturated) Selection3 = Selection1+" && "+Selection2;
    else Selection3 = Selection1;
    std::cout << "Selection = " << Selection3 << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection3);
  
    std::string varAmpMax = "0.";
    if(iMCP == "ZS1" && isSaturated) varAmpMax = "(488.8+53.87*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+2.998*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+8.715*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])-0.0308*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]))";
    else if(iMCP == "ZS2" && isSaturated) varAmpMax = "(493.7+32.73*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+61.56*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])-27.54*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+7.141*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]))";
    //std::cout << "new AmpMax = " << varAmpMax << std::endl;
    
    if(isSaturated) setBins(h4, varAmpMax, Selection3,maxMCP,10,ptr_minVec);
    else setBins(h4, std::string("amp_max["+iMCP+"]"), Selection3, maxMCP,10.,ptr_minVec);
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,minVec,4000,timeMin);
   
    if(isSaturated) h4->Draw(std::string("time["+iMCP+iTiming+"]-time[MiB2]:"+varAmpMax+" >> timingCorrection_wrtMiB2").c_str(),Selection3.c_str());
    else h4->Draw(std::string("time["+iMCP+iTiming+"]-time[MiB2]:amp_max["+iMCP+"] >> timingCorrection_wrtMiB2").c_str(),Selection3.c_str());
   
    std::vector<float> points_wrtMiB2;
    timingCorrection_wrtMiB2->FitSlicesY();
    TH1F* timingCorrection_wrtMiB2_1 = (TH1F*)inputFile->Get("timingCorrection_wrtMiB2_1");
    timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle(std::string("amp_max["+iMCP+"]").c_str());
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
    if(isSaturated){
        fit_corr1 = new TF1("fit_corr1","pol2",0.,7000.);
    }else{
        fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,3000.);
        fit_corr1->SetParLimits(2,0.,999999999.);
    }
    timingCorrection_wrtMiB2_1->Fit("fit_corr1","B");
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    timingCorrection_wrtMiB2_1->Draw();
    fit_corr1->Draw("same");
    if(iMCP != "MiB2"){
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("timingCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }else{
       c1 -> Print(std::string("timingCorrection_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("timingCorrection_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }
}

void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP)
{
    TH2D* h2_time_max_vs_amp;
    if(isSaturated) h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",700,0.,7000.,500,0.,5.);
    else h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    
    TF1* fit_corr_max;
    if(isSaturated) fit_corr_max = new TF1("fit_corr_max","[0]+[1]/(x+[2])",0.,7000.);
    else fit_corr_max = new TF1("fit_corr_max","[0]+[1]/(x+[2])",0.,3000.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      if(isSaturated) Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      else Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);
    
    if(iMCP == "ZS1" && isSaturated) h4->Draw((std::string("time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]:(488.8+53.87*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+2.998*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+8.715*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])-0.0308*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])) >> h2_time_max_vs_amp")).c_str(),Selection.c_str(),"goff");
    else if(iMCP == "ZS2" && isSaturated) h4->Draw((std::string("time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]:(493.7+32.73*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+61.56*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])-27.54*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])+7.141*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])) >> h2_time_max_vs_amp")).c_str(),Selection.c_str(),"goff");
    else h4->Draw((std::string("time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]:amp_max[")+iMCP+std::string("] >> h2_time_max_vs_amp")).c_str(),Selection.c_str(),"goff");

    if(Timing == "CFD50"){
       fit_corr_max->SetParameters(2,0.,40.);
       fit_corr_max->SetParLimits(2,0.,40.);
    }else if(Timing == "LED50"){
       fit_corr_max->SetParameters(2,45.,55.);
       fit_corr_max->SetParLimits(2,45.,55.);
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
    else h2_time_max_vs_amp->Fit("fit_corr_max","B","",0.,7000.);

    Params->push_back(fit_corr_max->GetParameter(0));
    Params->push_back(fit_corr_max->GetParameter(1));
    Params->push_back(fit_corr_max->GetParameter(2));
    
    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h2_time_max_vs_amp->Draw("COLZ");
    fit_corr_max->Draw("same");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
}

void SaturationCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP)
{    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
     
    std::string min = "20.";
    if(Timing == "LED50") min = "50.";
    else if(Timing == "LED100") min = "100.";
    else if(Timing == "LED150") min = "150.";
    else if(Timing == "LED500") min = "500.";

    std::string Selection = "amp_max["+iMCP+"]>"+min+" && amp_max["+iMCP+"]<3000.";
    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);

    for(int ii = 0; ii<nBins+1;ii++)
        minVec[ii] = ii*0.714;
    //setBins(h4, std::string("time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]"), Selection,std::string("6."),0.01,ptr_minVec);
    //TH2F* saturationCorrection_wrtMiB2 = new TH2F("saturationCorrection_wrtMiB2","",nBins,minVec,3000,ampMin);
    TH2F* saturationCorrection_wrtMiB2 = new TH2F("saturationCorrection_wrtMiB2","",nBins,minVec,3000,ampMin);
    
    h4->Draw(std::string("amp_max["+iMCP+"]:(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) >> saturationCorrection_wrtMiB2").c_str(),Selection.c_str());

    TCanvas* c2 = new TCanvas();
    c2->cd();
    saturationCorrection_wrtMiB2->Draw("COLZ");
    c2 -> Print(std::string("h2_saturationCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    c2 -> Print(std::string("h2_saturationCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
   
    TProfile *saturationCorrection_wrtMiB2_profileX = saturationCorrection_wrtMiB2->ProfileX();
    saturationCorrection_wrtMiB2_profileX->SetMarkerStyle(20);
    saturationCorrection_wrtMiB2_profileX->SetMarkerSize(0.9);
    saturationCorrection_wrtMiB2_profileX->SetMarkerColor(kBlack);
    saturationCorrection_wrtMiB2_profileX->SetLineColor(kBlack);
    saturationCorrection_wrtMiB2_profileX->GetXaxis()->SetTitle("time_mirror-time (ns)");
    saturationCorrection_wrtMiB2_profileX->GetYaxis()->SetTitle("amp_max (ADC)");
    
    TF1* fit_corr1;
    fit_corr1 = new TF1("fit_corr1","pol4",0.,10.);
    saturationCorrection_wrtMiB2_profileX->Fit("fit_corr1","B");
   
    TCanvas* c1 = new TCanvas();
    c1->cd();
    saturationCorrection_wrtMiB2_profileX->Draw();
    fit_corr1->Draw("same");
    c1->Update();
    TPaveStats *stats = (TPaveStats*)c1->GetPrimitive("stats");
    //stats->SetName("stats");
    stats->SetX1NDC(.15);
    stats->SetX2NDC(.55);
    stats->SetY1NDC(.65);
    stats->SetY2NDC(.85);
    if(iMCP != "MiB2"){
       c1 -> Print(std::string("saturationCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("saturationCorrection_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }else{
       c1 -> Print(std::string("saturationCorrection_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
       c1 -> Print(std::string("saturationCorrection_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
    }
}


void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",300,-10,20,300,-1.,1.5,100.,3000.);
    TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time","",300,-10,20,300,-1.,1.5);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;
  
    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

    Selection = Selection+" && WF_ch == "+iMCP;

    h4->Draw(std::string("amp_max["+iMCP+"]:WF_val/amp_max["+iMCP+"]:WF_time-time["+iMCP+"] >> p2D_amp_vs_time").c_str(),Selection.c_str(),"goff");
    h4->Draw(std::string("WF_val/amp_max["+iMCP+"]:WF_time-time["+iMCP+"] >> h2_amp_vs_time").c_str(),Selection.c_str());

    TProfile* waveForm = h2_amp_vs_time->ProfileX(); 
    waveForm->SetName(std::string(nameiMCP+"_waveform_prof").c_str());
    
    p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    waveForm->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
    h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");

    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

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

void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TH2F* h2_hodoscope_Y_vs_X_plane1_noSelection = new TH2F("h2_hodoscope_Y_vs_X_plane1_noSelection","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane1 = new TH2F("h2_hodoscope_Y_vs_X_plane1","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane2_noSelection = new TH2F("h2_hodoscope_Y_vs_X_plane2_noSelection","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane2 = new TH2F("h2_hodoscope_Y_vs_X_plane2","",32,-16,16,32,-16,16);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

    h4->Draw("Y[0]:X[0] >> h2_hodoscope_Y_vs_X_plane1_noSelection");
    h4->Draw("Y[0]:X[0] >> h2_hodoscope_Y_vs_X_plane1",Selection.c_str());
    h4->Draw("Y[1]:X[1] >> h2_hodoscope_Y_vs_X_plane2_noSelection");
    h4->Draw("Y[1]:X[1] >> h2_hodoscope_Y_vs_X_plane2",Selection.c_str());

    h2_hodoscope_Y_vs_X_plane1_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane1->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane1_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane1->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane2_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane2->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane2_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane2->GetYaxis()->SetTitle("Y");

    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetGrid();
    h2_hodoscope_Y_vs_X_plane1_noSelection->Draw("COLZ");
    c1 -> Print("Hodoscope_Y_vs_X_plane1_noSelection.png","png");
    c1 -> Print("Hodoscope_Y_vs_X_plane1_noSelection.pdf","pdf");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    c2->SetGrid();
    h2_hodoscope_Y_vs_X_plane1->Draw("COLZ");
    c2 -> Print(std::string("Hodoscope_Y_vs_X_plane1_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c2 -> Print(std::string("Hodoscope_Y_vs_X_plane1_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");

    TCanvas* c3 = new TCanvas();
    c3->cd();
    c3->SetGrid();
    h2_hodoscope_Y_vs_X_plane2_noSelection->Draw("COLZ");
    c3 -> Print("Hodoscope_Y_vs_X_plane2_noSelection.png","png");
    c3 -> Print("Hodoscope_Y_vs_X_plane2_noSelection.pdf","pdf");

    TCanvas* c4 = new TCanvas();
    c4->cd();
    c4->SetGrid();
    h2_hodoscope_Y_vs_X_plane2->Draw("COLZ");
    c4 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c4 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");
   
}

void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TH1F* h_timeChi2_noSelection = new TH1F("h_timeChi2_noSelection","",1500,1.,15000.);
    TH1F* h_timeChi2 = new TH1F("h_timeChi2","",1500,1.,15000.);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2_noSelection").c_str());
    h4->Draw(std::string("time_chi2["+iMCP+"] >> h_timeChi2").c_str(),Selection.c_str());
    
    h_timeChi2_noSelection->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2->GetXaxis()->SetTitle("time_chi2");
    h_timeChi2_noSelection->GetYaxis()->SetTitle("Events");
    h_timeChi2->GetYaxis()->SetTitle("Events");

    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

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

void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TH1F* h_ampMax_noSelection = new TH1F(std::string("h_ampMax_"+nameiMCP+"_noSelection").c_str(),"",450,0.,4500.);
    TH1F* h_ampMax = new TH1F(std::string("h_ampMax_"+nameiMCP).c_str(),"",450,0.,4500.);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 
    
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP+"_noSelection").c_str());
    h4->Draw(std::string("amp_max["+iMCP+"] >> h_ampMax_"+nameiMCP).c_str(),Selection.c_str());
    
    h_ampMax_noSelection->GetXaxis()->SetTitle("amp_max");
    h_ampMax->GetXaxis()->SetTitle("amp_max");
    h_ampMax_noSelection->GetYaxis()->SetTitle("Events");
    h_ampMax->GetYaxis()->SetTitle("Events");

    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";
 
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
    if(iMCP == "ZS1") HV = "HVZS1";
    else if(iMCP == "ZS2") HV = "HVZS2";
    else if(iMCP == "SEE") HV = "HVSEE";
    else if(iMCP == "MiB3") HV = "HVMiB3";
    else if(iMCP == "MiB2") HV = "HVZS2";


    h4->Draw(std::string(HV+" >> num").c_str(),std::string("amp_max[MiB2]>20.").c_str());   
    h4->Draw(std::string(HV+" >> den").c_str());     

    TGraphAsymmErrors* eff = new TGraphAsymmErrors(num,den);
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    eff->Draw("AP");
    c1 -> Print(std::string("Efficiency_MiB2.png").c_str(),"png");
    c1 -> Print(std::string("Efficiency_MiB2.pdf").c_str(),"pdf");

    if(iMCP != "MiB2")
    {
       h4->Draw(std::string(HV+" >> num").c_str(),std::string("amp_max[MiB2]>200. && amp_max["+iMCP+"]>20.").c_str());   
       h4->Draw(std::string(HV+" >> den").c_str(),std::string("amp_max[MiB2]>200.").c_str());   

       eff = new TGraphAsymmErrors(num,den);
    
       TCanvas* c2 = new TCanvas();
       c2->cd();
       eff->Draw("AP");
       c2 -> Print(std::string("Efficiency_"+iMCP+".png").c_str(),"png");
       c2 -> Print(std::string("Efficiency_"+iMCP+".pdf").c_str(),"pdf");
    }
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
    std::string sMean;

    if(isCut == false){
      h->Fit("g_fit","","",h->GetMean()-3*h->GetRMS(),h->GetMean()+3*h->GetRMS());
      char Mean [100];
      sprintf(Mean,"%f",h->GetMean());
      sMean = std::string(Mean);
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
      sMean = std::string(Mean);

      if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+std::string(" && fabs(")+Var+std::string("+")+std::string(sMean)+std::string(")<")+Cut;
      }else{
       Selection = Selection+std::string(" && fabs(")+Var+std::string("-")+std::string(sMean)+std::string(")<")+Cut;
      }
    }

    sShift = atof(sMean.c_str());

    delete h;
    delete g_fit;
    return Selection;
}

std::string shiftVar(TTree* h4, std::string Var, std::string Selection)
{
    TH1F* h = new TH1F("h","h",2000000,-1000.,1000.);
    h4->Draw((Var+std::string(" >> h")).c_str(),Selection.c_str());
    
    char Mean [100];
    sprintf(Mean,"%f",h->GetMean());
    std::string sMean = std::string(Mean);

    delete h;
    
    return std::string(Var+"-("+Mean+")");
}

std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel)
{
    TH1F* num = new TH1F("num","",40,0.,4000.);
    TH1F* den = new TH1F("den","",40,0.,4000.);

    std::string HV = "";
    if(iMCP == "ZS1") HV = "HVZS1";
    else if(iMCP == "ZS2") HV = "HVZS2";
    else if(iMCP == "SEE") HV = "HVSEE";
    else if(iMCP == "MiB3") HV = "HVMiB3";
    else if(iMCP == "MiB2") HV = "HVMiB2";


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

int getHV(TTree* h4, std::string HV)
{
    TH1F* h_HV = new TH1F("h_HV","",38,0,3800.);
    
    h4->Draw(std::string(HV+" >> h_HV").c_str());   
    int hv = 0;
    for(int ii = 1; ii < h_HV->GetNbinsX(); ii++)
        if(h_HV->GetBinContent(ii) > 0) hv = (int)h_HV->GetBinCenter(ii)+50;
    
   return hv;
}   

void setBins(TTree* h4, std::string var, std::string Selection, std::string maxMCP, float step, float* ptr_minVec)
{
    TH1F* h_var = new TH1F("h_var","",99999,-99999.,99999);
    
    h4->Draw(std::string(var+">>h_var").c_str(),std::string(Selection+" && "+var+"<"+maxMCP).c_str());   
    int Total_Events = h_var->GetEntries();
    int nEvents = (float)Total_Events/nBins;
    
    float binMax = 0.;
    float binMin = 0.;
    for(int ii = 0; ii < nBins; ii++){
      do{
        binMax+= step; 
        char min[100];
        sprintf(min,"%f",binMin);
        char max[100];
        sprintf(max,"%f",binMax);
        h4->Draw(std::string(var+">>h_var").c_str(),std::string(Selection+ " && "+var+">"+min+" && "+var+"<"+max).c_str());  
        //std::cout << ii << " - binMin " << min << " binMax = " << max << " " << atof(maxMCP.c_str()) << " " << h_var->GetEntries() << " " << nEvents << std::endl;   
        if(binMax >= atof(maxMCP.c_str())) break; 
      }while(h_var->GetEntries()<nEvents);
      *(ptr_minVec+ii+1) = binMax;
      binMin = binMax;
    }
    for(int ii = 0; ii <= nBins; ii++)
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

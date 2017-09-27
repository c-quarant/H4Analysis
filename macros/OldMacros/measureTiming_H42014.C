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

int nBins = 7;
float minVec[8] = {0.};
//int nBins = 10;
//float minVec[11] = {0.};
int nBins2 = 10;
float minVec2[11] = {0.};
float *ptr_minVec = minVec;
float *ptr_minVec2 = minVec2;
float timeMin[4001];
float ampMin[3001];

//BINP
std::string minMiB2 = "200.";
std::string maxMiB2 = "3400.";
std::string time_max_MiB2 = "150.";
std::string time_max_Rm2 = "150.";
std::string amp_max_Rm2 = "200.";
std::string deltaCrossMin = "0.";
//std::string deltaCrossMin = "0.8";
std::string deltaCrossMax = "9999.";
//std::string deltaCrossMax = "5.6";
std::string scintMin = "200.";
std::string scintMax = "700.";
std::string timeChi2 = "99999999.";
std::string hodoX = "-800";
std::string hodoY = "-800";
std::string hodoXmin = "-9";
std::string hodoXmax = "1";
std::string hodoYmin = "-6";
std::string hodoYmax = "5";
std::string thresSaturation = "3000.";

bool isSaturated = true;
bool doDoubleGauss = true;
bool doDoubleSetBins = false;
float sShift=0.;

void FinalTiming(TTree* h4,TTree* digi,TTree* hodo, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, std::string* ptr_timeShifted);
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
void setBins(TTree*, std::string, std::string, std::string, std::string, float, float*, int);
void CheckSelectionEfficiency(TTree* h4, std::string iMCP, std::string Selection);
int getHV(TTree* h4, std::string HV);
std::vector<std::string> split(const std::string text, std::string sep);
std::pair<float,float> FWHMSigma(TH1F* h,TF1* f);

void measureTiming_H42014(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, std::string maxMCP, bool doFirstStep = true, bool doPulseShapes = false, bool saturationCorr = false)
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
       Hodoscope(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       TimeChi2(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       AmpMax(h4, iMCP, nameiMCP, thresMCP, maxMCP);
       CheckEfficiency(h4, inputs, iMCP, nameiMCP, thresMCP, maxMCP);
    }else if(doPulseShapes == true) PulseShapes(h4, iMCP, nameiMCP, thresMCP, maxMCP);
    else{
       if(saturationCorr) SaturationCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP);
       else{
          if(!isSaturated) AmpVsTime_Selection(h4, iMCP, nameiMCP,Timing, Params, thresMCP, maxMCP);
          TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, ptr_timeShifted);
          FinalTiming(h4,digi,hodo, inputs, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, ptr_timeShifted);
       }
    }
}

void FinalTiming(TTree* h4,TTree* digi,TTree* hodo, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, std::string* ptr_timeShifted)
{
    //TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",200,-0.5,0.5);
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",100,-0.5,0.5);
    TH2F* time_vs_amp_wrtMiB2;
    TH1F* time_wrtMiB2_FWHM = new TH1F("time_wrtMiB2_FWHM","",2000,-1.,1.);

    TGraphAsymmErrors* g_Res_vs_Amp_wrtMiB2 = new TGraphAsymmErrors(); 
    std::vector<float> points_wrtMiB2;
    std::vector<TH1F*> resHist_wrtMiB2;
    resHist_wrtMiB2.resize(nBins);
    std::vector<TF1*> resFitAlt_wrtMiB2;
    resFitAlt_wrtMiB2.resize(nBins);
    int iPoint_MiB2 = 0;
    int iPoint_Rm2 = 0;

    for(int ii = 0; ii < nBins; ii++)
    {
        char Name_wrtMiB2 [50];
        sprintf (Name_wrtMiB2,"h_Res_1_%d_wrtMiB2",ii);
        //resHist_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,100,-0.5,0.5); 
        resHist_wrtMiB2[ii] = new TH1F(Name_wrtMiB2,Name_wrtMiB2,50,-0.5,0.5);       
    }

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    TF1 *g_res;
    TF1 *g_res_single;
    TF1 *g_res_double;

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
    if(isSaturated) Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax+" && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])>"+deltaCrossMin+" && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])<"+deltaCrossMax;
    else Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
    Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);

    if(doDoubleSetBins){
       setBins(h4, std::string("amp_max["+iMCP+"]"), Selection,thresMCP,maxMCP,1.,ptr_minVec2, nBins2);
       time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins2,minVec2,4000,timeMin);
    }else{
       time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins,minVec,4000,timeMin);
    }

    float thres, par0, par1, par2, par3, par4, par5, par6, par7;
    ifstream infile;  
    if(iMCP == "ZS1" && isSaturated) infile.open("saturation_correction_ZS1_LED100.dat"); 
    else if(iMCP == "ZS2" && isSaturated) infile.open("saturation_correction_ZS2_LED100.dat");     
    if(!infile.fail()  && isSaturated)
    {
        while(!infile.eof())
        {
             if(infile.eof())
             break;
             infile >> thres >> par0 >> par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> par7;
             std::cout << thres << " " << par0 << " " << par1 << " " << par2 << " " << par3 << " " << par4 << " " << par5 << " " << par6 << " " << par7 << std::endl;
        }    
    }
    if(isSaturated) infile.close(); 

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
        
        if(!isSaturated){
           if(amp_max[imcp]<atof(thresMCP.c_str())) continue;
        }else if(isSaturated){
           if(amp_max[imcp]<20.) continue;
        }

        if(amp_max[imcp]> atof(maxMCP.c_str())) continue;
        if(amp_max[MiB2]<atof(minMiB2.c_str())) continue;
        if(amp_max[MiB2]> atof(maxMiB2.c_str())) continue;
        if(fabs(time[MiB2]-time[imcp+timing]-(sShift))>1.) continue;
        if(X[0]<atoi(hodoXmin.c_str()) || X[0]>atoi(hodoXmax.c_str()) || Y[0]<atoi(hodoYmin.c_str()) || Y[0]>atoi(hodoYmax.c_str())) continue;
        if(isSaturated && (time_mirror[imcp+timing]-time[imcp+timing])<0.) continue;    
        if(isSaturated && (time_mirror[imcp+timing]-time[imcp+timing])>atof(deltaCrossMax.c_str())) continue;
        if(!isSaturated && fabs(time_max[imcp]-time[imcp+timing]-(Params->at(0)+(Params->at(1))*log(Params->at(2)+amp_max[imcp])))>0.25) continue;   //(time-time_max) vs amp selection
        
    
        //time correction
        float time_corr = 0.;
        float amp_corr = amp_max[imcp];
        float deltaCross = time_mirror[imcp+timing]-time[imcp+timing];
        
        if(isSaturated){
           if(amp_max[imcp]>atof(thresSaturation.c_str())) amp_corr = par0 + par1*deltaCross + par2*deltaCross*deltaCross + par3*deltaCross*deltaCross*deltaCross + par4*deltaCross*deltaCross*deltaCross*deltaCross + par5*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross + par6*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross + par7*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross*deltaCross;
           else if(amp_max[imcp]<atof(thresSaturation.c_str())) amp_corr = amp_max[imcp];
           if(amp_max[imcp]>atof(thresSaturation.c_str())) time_corr = Params_wrtMiB2->at(0) + Params_wrtMiB2->at(1)*1/(amp_corr+Params_wrtMiB2->at(2)); 
           else if(amp_max[imcp]<atof(thresSaturation.c_str())) time_corr = 0.4149 + (-84.0554)*1/(amp_corr+189.473); //1X0
           //else if(amp_max[imcp]<atof(thresSaturation.c_str())) time_corr = 0.3984 + (-129.806)*1/(amp_corr+496.510); //2X0
        }else time_corr = Params_wrtMiB2->at(0) + Params_wrtMiB2->at(1)*1/(amp_corr + Params_wrtMiB2->at(2)); 
       
        if(isSaturated){
           if(amp_max[imcp]>atof(thresSaturation.c_str())){
              time_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr);
              time_vs_amp_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr,amp_corr);
           }else{
              time_wrtMiB2->Fill(time[imcp+CFD]-time[MiB2]-time_corr);
              time_vs_amp_wrtMiB2->Fill(time[imcp+CFD]-time[MiB2]-time_corr,amp_corr); 
           }
        }else{
           time_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr);
           time_vs_amp_wrtMiB2->Fill(time[imcp+timing]-time[MiB2]-time_corr,amp_corr);
        }

        for(int ii = 0; ii < nBins; ii++){
            if(isSaturated){
               if(amp_max[imcp]>atof(thresSaturation.c_str())){
                  if(doDoubleSetBins && (amp_corr>minVec2[ii] && amp_corr<minVec2[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+timing]-time[MiB2]-time_corr);
                  else if(!doDoubleSetBins && (amp_corr>minVec[ii] && amp_corr<minVec[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+timing]-time[MiB2]-time_corr);
               }else{
                  if(doDoubleSetBins && (amp_corr>minVec2[ii] && amp_corr<minVec2[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+CFD]-time[MiB2]-time_corr);
                  else if(!doDoubleSetBins && (amp_corr>minVec[ii] && amp_corr<minVec[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+CFD]-time[MiB2]-time_corr);
               } 
            }else{
               if(doDoubleSetBins && (amp_corr>minVec2[ii] && amp_corr<minVec2[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+timing]-time[MiB2]-time_corr);
               else if(!doDoubleSetBins && (amp_corr>minVec[ii] && amp_corr<minVec[ii+1])) resHist_wrtMiB2[ii]->Fill(time[imcp+timing]-time[MiB2]-time_corr);
            } 
        }
    }
    
    if(time_wrtMiB2->GetEntries() < 2000) time_wrtMiB2->Rebin(2);
    
    float sigma_eff;
    float s_sigma;

    if(iMCP != "MiB2"){

       float fwhm;
       float fwhm_error;

       if(doDoubleGauss == false) g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-7*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+7*time_wrtMiB2->GetRMS());  
       else g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-7*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+7*time_wrtMiB2->GetRMS());  
       if(doDoubleGauss == true) g_res_single = new TF1("g_res_single","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-7*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+7*time_wrtMiB2->GetRMS());
       if(doDoubleGauss == true) g_res_double = new TF1("g_res_double","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-7*time_wrtMiB2->GetRMS(),time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+7*time_wrtMiB2->GetRMS());
       g_res->SetParameters(0,time_wrtMiB2->GetEntries()/2.);
       g_res->SetParameters(1,0.);
       //g_res->SetParameters(2,0.04);
       g_res->SetParameters(2,0.1);
       g_res->SetParameters(3,0.);
       if(doDoubleGauss == true) g_res->SetParameters(3,time_wrtMiB2->GetEntries()/2.);
       if(doDoubleGauss == true) g_res->SetParameters(4,0.);
       if(doDoubleGauss == true) g_res->SetParameters(5,0.1);
       //if(doDoubleGauss == true) g_res->SetParameters(5,0.04);
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
       float sigma_eff_sub;
       float s_sigma_sub;
       float sub_wrtMiB2 = 24.E-3;
       float sub_wrtMiB2_error = 2.E-3;

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

       TFile* output_time_wrtMiB2 = new TFile("output_time_wrtMiB2.root","RECREATE");
       output_time_wrtMiB2->cd();
       time_wrtMiB2->Write();
       output_time_wrtMiB2->Close();

       std::cout << "sigma_eff = " << sigma_eff*1000. << "+/-" << s_sigma*1000. << std::endl;

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
       if(iMCP != "MiB2"){
          c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
          c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_wrtMiB2_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
       }else{
          c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".png").c_str(),"png");
          c1 -> Print(std::string("TimeResolution_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".pdf").c_str(),"pdf");
      }
      delete g_res;
    }

    std::cout << "sigma_eff = " << sigma_eff << "+/-" << s_sigma << std::endl;

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

    TH2F* time_vs_amp_wrtMiB2_2;
    
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

       TCanvas* c2 = new TCanvas();
       c2->cd();
       time_vs_amp_wrtMiB2_2->Draw();
       c2 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_wrtMiB2_auto.png").c_str(),"png");
       c2 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_wrtMiB2_auto.pdf").c_str(),"pdf");
    }

    doDoubleGauss = true;

    int nBinsTotal = 0;
    if(doDoubleSetBins) nBinsTotal = nBins2;
    else nBinsTotal = nBins;
 
    for(int ii = 0; ii<nBinsTotal; ii++)
    {
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        
        sprintf (NameOutput_wrtMiB2_png,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s_wrtMiB2.png",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s_wrtMiB2.pdf",nameiMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());

        if(iMCP != "MiB2"){
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
           resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");

           float sigma_eff = resFitAlt_wrtMiB2[ii]->GetParameter(2);
           float s_sigma = resFitAlt_wrtMiB2[ii]->GetParError(2);
           float sigma_eff_sub;
           float s_sigma_sub;
           float sub_wrtMiB2 = 26.E-3;
           float sub_wrtMiB2_error = 3.E-3;
    
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
           
           if(!doDoubleSetBins)
           {
              g_Res_vs_Amp_wrtMiB2->SetPoint(iPoint_MiB2,minVec[ii]+(minVec[ii+1]-minVec[ii])/2,sigma_eff*1000.);
              g_Res_vs_Amp_wrtMiB2->SetPointError(iPoint_MiB2,(minVec[ii+1]-minVec[ii])/2,(minVec[ii+1]-minVec[ii])/2,s_sigma*1000.,s_sigma*1000.);
           }else{
              g_Res_vs_Amp_wrtMiB2->SetPoint(iPoint_MiB2,minVec2[ii]+(minVec2[ii+1]-minVec2[ii])/2,sigma_eff*1000.);
              g_Res_vs_Amp_wrtMiB2->SetPointError(iPoint_MiB2,(minVec2[ii+1]-minVec2[ii])/2,(minVec2[ii+1]-minVec2[ii])/2,s_sigma*1000.,s_sigma*1000.);
           }
           iPoint_MiB2++;

           std::cout << ii << " - sigma_eff = " << sigma_eff*1000. << "+/-" << s_sigma*1000. << std::endl; 
           points_wrtMiB2.push_back(sigma_eff*1000.);
        
           TCanvas* c3 = new TCanvas();
           c3->cd();
           resHist_wrtMiB2[ii]->Draw("hist");
           resFitAlt_wrtMiB2[ii]->Draw("same");
           c3 -> Print(NameOutput_wrtMiB2_png,"png");
           c3 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
           delete c3;
        }
    }
    if(points_wrtMiB2.size() != 0) std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    
    if(iMCP != "MiB2"){
       g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ps)");
       g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
       g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
       g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
       g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-15,points_wrtMiB2.at(points_wrtMiB2.size()-1)+15);
       //g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-15,120);

       TCanvas* c4 = new TCanvas();
       c4->cd();
       g_Res_vs_Amp_wrtMiB2->Draw("AP");
       c4 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_wrtMiB2.png").c_str(),"png");
       c4 -> Print(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+"_wrtMiB2.pdf").c_str(),"pdf");
    }

    /*TFile* output_ampMax = new TFile(std::string("Data_TimeResolution_vs_amp_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+".root").c_str(),"RECREATE");
    output_ampMax->cd();
    g_Res_vs_Amp_wrtMiB2->Write(std::string("TimeResolution_vs_amp_wrtMiB2_"+nameiMCP+"_"+Timing+"_thres"+thresMCP).c_str());
    output_ampMax->Close();*/
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
      if(!isSaturated) Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
      else Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax+" && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])>"+deltaCrossMin+" && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])<"+deltaCrossMax;
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
    } 
    std::string Selection3;
    if(!isSaturated) Selection3 = Selection1+" && "+Selection2;
    else Selection3 = Selection1;
    
    std::cout << "Selection = " << Selection3 << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection3);

    std::string varAmpMax = "0.";
    ifstream infile;  
    if(iMCP == "ZS1" && isSaturated) infile.open("saturation_correction_ZS1_LED100.dat"); 
    else if(iMCP == "ZS2" && isSaturated) infile.open("saturation_correction_ZS2_LED100.dat");     
    if(!infile.fail()  && isSaturated)
    {
       while(!infile.eof())
       {
          if(infile.eof())
             break;
          char thres[50], par0[50], par1[50], par2[50], par3[50], par4[50], par5[50], par6[50], par7[50];
          infile >> thres >> par0 >> par1 >> par2 >> par3 >> par4 >> par5 >> par6 >> par7;
          std::cout << thres << " " << par0 << " " << par1 << " " << par2 << " " << par3 << " " << par4 << " " << par5 << " " << par6 << " " << par7 << std::endl;
          if(isSaturated) varAmpMax = "(("+std::string(par0)+") + ("+std::string(par1)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par2)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par3)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par4)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par5)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par6)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]) + ("+std::string(par7)+")*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])*(time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"]))";
       }    
    }
    if(isSaturated) infile.close(); 

    std::cout << "new AmpMax = " << varAmpMax << std::endl;
    
    if(isSaturated) setBins(h4, varAmpMax, Selection3,thresSaturation,maxMCP,1.,ptr_minVec, nBins);
    else setBins(h4, std::string("amp_max["+iMCP+"]"), Selection3,thresMCP,maxMCP,1.,ptr_minVec, nBins);
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,minVec,4000,timeMin);

    //for(int ii = 0; ii < 16; ii++)
    //    minVec[ii] = ii*466; 
    //TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,minVec,4000,timeMin);
   
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
        //fit_corr1 = new TF1("fit_corr1","pol2",0.,7000.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]*log([2]+x)",0.,7000.);
        fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,8000.);
        //fit_corr1->SetParLimits(2,0.,999999999.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]*exp(x*[2])",0.,13000.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]/(([2]+x)*([3]+x)*([4]+x)*([5]+x))",0.,13000.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]/(1+[2]*exp([3]*x))",0.,13000.);
        //fit_corr1->SetParLimits(3,-999999999.,0.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]/(([2]+x)*([3]+x)*([4]+x)*([5]+x))",0.,7000.);
        //fit_corr1 = new TF1("fit_corr1","[0]+sqrt([1]/(([2]+x)*([3]+x)*([4]+x))+[5]*x)",0.,7000.);
    }else{
        fit_corr1 = new TF1("fit_corr1","[0]+[1]/(x+[2])",0.,3400.);
        //fit_corr1->SetParLimits(2,0.,999999999.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]*exp([2]*x)",0.,3400.);
        //fit_corr1 = new TF1("fit_corr1","[0]+[1]/(([2]+x)*([3]+x)*([4]+x)*([5]+x))",0.,3400.);
    }
    timingCorrection_wrtMiB2_1->Fit("fit_corr1","B");
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(3));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(4));
    Params_wrtMiB2->push_back(fit_corr1->GetParameter(5));
    
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
    fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",0.,3400.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
    }else{ 
      if(isSaturated) Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      else Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);
    
    if(!isSaturated) h4->Draw((std::string("time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]:amp_max[")+iMCP+std::string("] >> h2_time_max_vs_amp")).c_str(),Selection.c_str(),"goff");

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
    if(Timing == "CFD50") h2_time_max_vs_amp->Fit("fit_corr_max","B","",20.,3400.);
    else if(Timing == "LED50") h2_time_max_vs_amp->Fit("fit_corr_max","B","",50.,3400.);
    else if(Timing == "LED100") h2_time_max_vs_amp->Fit("fit_corr_max","B","",100.,3400.);
    else if(Timing == "LED150") h2_time_max_vs_amp->Fit("fit_corr_max","B","",150.,3400.);
    
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
    else if(Timing == "LED200") min = "200.";
    else if(Timing == "LED300") min = "300.";
    else if(Timing == "LED400") min = "400.";
    else if(Timing == "LED500") min = "500.";
    else if(Timing == "LED600") min = "600.";
    else if(Timing == "LED700") min = "700.";
    else if(Timing == "LED800") min = "800.";
    else if(Timing == "LED900") min = "900.";
    else if(Timing == "LED1000") min = "1000.";
   

    std::string Selection = "amp_max["+iMCP+"]>"+min+" && amp_max["+iMCP+"]<"+maxMCP+" && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])>"+deltaCrossMin+"  && (time_mirror["+iMCP+iTiming+"]-time["+iMCP+iTiming+"])<"+deltaCrossMax;
    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);

    TH2F* saturationCorrection_wrtMiB2 = new TH2F("saturationCorrection_wrtMiB2","",70,0.,7.,340,0.,3400.);
    
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
    saturationCorrection_wrtMiB2_profileX->GetXaxis()->SetRangeUser(0.,6.);

    TF1* fit_corr1;
    fit_corr1 = new TF1("fit_corr1","pol7",0.,7.);
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
    
    int thres = 0;
    if(Timing == "LED50") thres = 50;
    else if(Timing == "LED100") thres = 100;
    else if(Timing == "LED150") thres = 150;
    else if(Timing == "LED200") thres = 200;
    else if(Timing == "LED300") thres = 300;
    else if(Timing == "LED400") thres = 400;
    else if(Timing == "LED500") thres = 500;
    else if(Timing == "LED600") thres = 600;
    else if(Timing == "LED700") thres = 700;
    else if(Timing == "LED800") thres = 800;
    else if(Timing == "LED900") thres = 900;
    else if(Timing == "LED1000") thres = 1000;

    std::ofstream sat_txt;
    sat_txt.open(std::string("saturation_correction_"+iMCP+".dat").c_str(), ios::out | ios::app);    
    if (sat_txt.is_open())
    {
        sat_txt << thres << " " << fit_corr1->GetParameter(0) << " " << fit_corr1->GetParameter(1) << " " << fit_corr1->GetParameter(2) << " " << fit_corr1->GetParameter(3) << " " << fit_corr1->GetParameter(4) << " " << fit_corr1->GetParameter(5) << " " << fit_corr1->GetParameter(6) << " " << fit_corr1->GetParameter(7) << "\n";
        sat_txt.close();
    }else{
        ofstream res_txt_new(std::string("saturation_correction_"+iMCP+".dat").c_str());
        sat_txt << iTiming.c_str() << " " << fit_corr1->GetParameter(0) << " " << fit_corr1->GetParameter(1) << " " << fit_corr1->GetParameter(2) << " " << fit_corr1->GetParameter(3) << " " << fit_corr1->GetParameter(4) << " " << fit_corr1->GetParameter(5) << " " << fit_corr1->GetParameter(6) << " " << fit_corr1->GetParameter(7) << "\n";
        sat_txt.close();
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
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 
    
    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);

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
    TH2F* h2_hodoscope_Y_vs_X_plane1_efficiency = new TH2F("h2_hodoscope_Y_vs_X_plane1_effiency","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane2_noSelection = new TH2F("h2_hodoscope_Y_vs_X_plane2_noSelection","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane2 = new TH2F("h2_hodoscope_Y_vs_X_plane2","",32,-16,16,32,-16,16);
    TH2F* h2_hodoscope_Y_vs_X_plane2_efficiency = new TH2F("h2_hodoscope_Y_vs_X_plane2_effiency","",32,-16,16,32,-16,16);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+"&& X[0]>"+hodoX+" && Y[0]>"+hodoY;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoX+" && Y[0]>"+hodoY;
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 
    
    std::cout << "Selection = " << Selection << std::endl;
    CheckSelectionEfficiency(h4,iMCP,Selection);

    h4->Draw("Y[0]:X[0] >> h2_hodoscope_Y_vs_X_plane1_noSelection");
    h4->Draw("Y[0]:X[0] >> h2_hodoscope_Y_vs_X_plane1",Selection.c_str());
    h4->Draw("Y[1]:X[1] >> h2_hodoscope_Y_vs_X_plane2_noSelection");
    h4->Draw("Y[1]:X[1] >> h2_hodoscope_Y_vs_X_plane2",Selection.c_str());

    h2_hodoscope_Y_vs_X_plane1_efficiency->Divide(h2_hodoscope_Y_vs_X_plane1,h2_hodoscope_Y_vs_X_plane1_noSelection,1,1,"B");
    h2_hodoscope_Y_vs_X_plane2_efficiency->Divide(h2_hodoscope_Y_vs_X_plane2,h2_hodoscope_Y_vs_X_plane2_noSelection,1,1,"B");

    h2_hodoscope_Y_vs_X_plane1_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane1->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane1_efficiency->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane1_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane1->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane1_efficiency->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane2_noSelection->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane2->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane2_efficiency->GetXaxis()->SetTitle("X");
    h2_hodoscope_Y_vs_X_plane2_noSelection->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane2->GetYaxis()->SetTitle("Y");
    h2_hodoscope_Y_vs_X_plane2_efficiency->GetYaxis()->SetTitle("Y");

    std::string wrtMCP = "";
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

    /*TCanvas* c1 = new TCanvas();
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
    c2 -> Print(std::string("Hodoscope_Y_vs_X_plane1_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");*/

    TCanvas* c3 = new TCanvas();
    c3->cd();
    c3->SetGrid();
    h2_hodoscope_Y_vs_X_plane1_efficiency->Draw("COLZ");
    c3 -> Print(std::string("Hodoscope_Y_vs_X_plane1_"+nameiMCP+wrtMCP+"_efficiency.png").c_str(),"png");
    c3 -> Print(std::string("Hodoscope_Y_vs_X_plane1_"+nameiMCP+wrtMCP+"_efficiency.pdf").c_str(),"pdf");

    /*TCanvas* c4 = new TCanvas();
    c4->cd();
    c4->SetGrid();
    h2_hodoscope_Y_vs_X_plane2_noSelection->Draw("COLZ");
    c4 -> Print("Hodoscope_Y_vs_X_plane2_noSelection.png","png");
    c4 -> Print("Hodoscope_Y_vs_X_plane2_noSelection.pdf","pdf");

    TCanvas* c5 = new TCanvas();
    c5->cd();
    c5->SetGrid();
    h2_hodoscope_Y_vs_X_plane2->Draw("COLZ");
    c5 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+".png").c_str(),"png");
    c5 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+".pdf").c_str(),"pdf");*/

    TCanvas* c6 = new TCanvas();
    c6->cd();
    c6->SetGrid();
    h2_hodoscope_Y_vs_X_plane2_efficiency->Draw("COLZ");
    c6 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+"_efficiency.png").c_str(),"png");
    c6 -> Print(std::string("Hodoscope_Y_vs_X_plane2_"+nameiMCP+wrtMCP+"_efficiency.pdf").c_str(),"pdf");
   
}

void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TH1F* h_timeChi2_noSelection = new TH1F("h_timeChi2_noSelection","",1500,1.,15000.);
    TH1F* h_timeChi2 = new TH1F("h_timeChi2","",1500,1.,15000.);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
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
      Selection = "amp_max[MiB2]>"+minMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+minMiB2+" && amp_max[MiB2]<"+maxMiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>"+hodoXmin+" && X[0]<"+hodoXmax+" && Y[0]>"+hodoYmin+" && Y[0]<"+hodoYmax;
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

void setBins(TTree* h4, std::string var, std::string Selection, std::string thresMCP, std::string maxMCP, float step, float* ptr_minVec, int nbins)
{
    TH1F* h_var = new TH1F("h_var","",9999999,-99999.,99999);
    
    h4->Draw(std::string(var+">>h_var").c_str(),Selection.c_str());   
    float Total_Events = h_var->Integral(h_var->FindBin(atof(thresMCP.c_str())),h_var->FindBin(atof(maxMCP.c_str())));
    float nEvents = (float)Total_Events/(float)nbins; 
    
    float binMax = atof(thresMCP.c_str());
    float binMin = atof(thresMCP.c_str());
    for(int ii = 0; ii < nbins; ii++){
      do{
        binMax+= step; 
        char min[100];
        sprintf(min,"%f",binMin);
        char max[100];
        sprintf(max,"%f",binMax);
        //h4->Draw(std::string(var+">>h_var").c_str(),std::string(Selection+ " && "+var+">"+min+" && "+var+"<"+max).c_str());  
        std::cout << ii << " - binMin " << min << " binMax = " << max << " " << atof(maxMCP.c_str()) << " " << h_var->Integral(h_var->FindBin(binMin),h_var->FindBin(binMax)) << " " << nEvents << std::endl;   
        if(binMax >= atof(maxMCP.c_str())) break; 
      }while(h_var->Integral(h_var->FindBin(binMin),h_var->FindBin(binMax))<nEvents);
      *(ptr_minVec+ii+1) = binMax;
      binMin = binMax;
    }
    *(ptr_minVec+0) = atof(thresMCP.c_str());
    for(int ii = 0; ii <= nbins; ii++)
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
   
   std::pair<float,float> outpair = std::make_pair(fwhm,fabs(fwhm_up-fwhm_down)/2.);
   return outpair;
}


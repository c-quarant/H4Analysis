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

#include<iostream>
#include<string>
#include<fstream>

//int nBins = 26;
//float ampMin[27] = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650, 700., 750., 800., 850., 900., 950., 1000., 1250, 1500., 1750., 2000.,2500.,3000.};
int nBins = 11;
float ampMin[12] = {0.,100.,200.,300.,400.,500.,800.,1200.,1700.,2000.,2500.,3000.};
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
bool isEfficiencyRun = false;

void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, bool doDoubleGauss, std::string* ptr_timeShifted);
void TimeCorrection(TTree* h4, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, std::string* ptr_timeShifted);
void AmpVsTime_Selection(TTree* h4, std::string iMCP, std::string nameiMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP, std::string maxMCP);
void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void Hodoscope(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void TimeChi2(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void AmpMax(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
void CheckEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP);
std::string AddSelection(TTree*, std::string, std::string, std::string, bool);
std::string shiftVar(TTree* h4, std::string Var, std::string Selection);
std::vector<float> ComputeEfficiency(TTree* h4, std::string inputs, std::string iMCP, std::string numSel, std::string denSel);
void CheckSelectionEfficiency(TTree* h4, std::string iMCP, std::string Selection);
int getHV(TTree* h4, std::string HV);
std::vector<std::string> split(const std::string text, std::string sep);

void ComputeTiming_oneStep_H42014(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, std::string maxMCP, bool doFirstStep = true, bool doPulseShapes = false, bool doDoubleGauss = false)
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
    
    for(int ii = 0; ii < 4001; ii++)
        timeMin[ii] = 0.01*ii-20.;

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
       AmpVsTime_Selection(h4, iMCP, nameiMCP,Timing, Params, thresMCP, maxMCP);
       TimeCorrection(h4, iMCP, nameiMCP,inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, ptr_timeShifted);
       FinalTiming(h4, inputs, iMCP, nameiMCP, inputFile, Timing, Params, Params_wrtMiB2, thresMCP, maxMCP, doDoubleGauss, ptr_timeShifted);
    }
}

void FinalTiming(TTree* h4, std::string inputs, std::string iMCP, std::string nameiMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::string thresMCP, std::string maxMCP, bool doDoubleGauss, std::string* ptr_timeShifted)
{
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",400,-1.,1.);
    TH2F* time_vs_amp_wrtMiB2 = new TH2F("time_vs_amp_wrtMiB2","",nBins,ampMin,4000,timeMin);

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    TF1 *g_res;

    char Selection1 [1000];  
    char Selection2 [1000];
    char Selection3 [1000];
    char Selection4 [1000];
    char Selection5 [1000];

    std::string newTime = *ptr_timeShifted;
    
    //(time-time_max) vs amp selection
    sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    //time correction
    if(iMCP != "MiB2"){
       sprintf (Selection2,std::string("time["+iMCP+iTiming+"]-time[MiB2]-(%f+(%f)*1/(%f+amp_max["+iMCP+"])) >>").c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       sprintf (Selection3,std::string("time["+iMCP+iTiming+"]-time[MiB2]-(%f+(%f)*1/(%f+amp_max["+iMCP+"])):amp_max["+iMCP+"] >>").c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
    }else{
       sprintf (Selection2,std::string(newTime+"-(%f+(%f)*1/(%f+amp_max["+iMCP+"])) >>").c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
       sprintf (Selection3,std::string(newTime+"-(%f+(%f)*1/(%f+amp_max["+iMCP+"])):amp_max["+iMCP+"] >>").c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
    }

    std::string Selection8;

    if(iMCP == "MiB2"){
      Selection8 = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection8 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
      Selection8 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection8,"1",true);
      //Selection8 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection8,"1",true);
    } 
    std::string Selection9 = Selection8+std::string(" && ")+std::string(Selection1); 

    //std::cout << "Selection = " << Selection9 << std::endl;
    //CheckSelectionEfficiency(h4,iMCP,Selection9);

    std::string Selection10 = Selection2+std::string(" time_wrtMiB2");
    std::string Selection11 = Selection3+std::string(" time_vs_amp_wrtMiB2");
    
    //std::cout << "Selection10 = " << Selection10 << std::endl;
    h4->Draw(Selection10.c_str(),Selection9.c_str()); 
    h4->Draw(Selection11.c_str(),Selection9.c_str()); 
    
    if(time_wrtMiB2->GetEntries() < 2000) time_wrtMiB2->Rebin(8);
    
    float sub_wrtMiB2 = 24.E-3;
    float sub_wrtMiB2_error = 2.E-3;

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

    if(iMCP != "MiB2"){
       sigma_eff_sub = sqrt(sigma_eff*sigma_eff-sub_wrtMiB2*sub_wrtMiB2);
       s_sigma_sub = sqrt(sigma_eff*sigma_eff*s_sigma*s_sigma+sub_wrtMiB2*sub_wrtMiB2*sub_wrtMiB2_error*sub_wrtMiB2_error)/sigma_eff_sub;
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
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",nBins,ampMin,4000,timeMin);
    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
   
    std::string Selection1;   
    char Selection2 [1000];

    sprintf (Selection2,std::string("fabs(time_max["+iMCP+"]-time["+iMCP+iTiming+"]-(%f+(%f)*log(%f+amp_max["+iMCP+"])))<0.25").c_str(),Params->at(0),Params->at(1),Params->at(2));
   
    if(iMCP == "MiB2"){
      Selection1 = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection1 = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
      Selection1 = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection1,"1",true);
      //Selection1 = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection1,"1",true);
    }  
    std::string Selection3 = Selection1+" && "+Selection2;
   
    std::string newTime;
    if(iMCP != "MiB2")
       h4->Draw(std::string("time["+iMCP+iTiming+"]-time[MiB2]:amp_max["+iMCP+"] >> timingCorrection_wrtMiB2").c_str(),Selection3.c_str());
    else{
       newTime = shiftVar(h4,std::string("time[MiB2"+iTiming+"]"),Selection3);
       *ptr_timeShifted = newTime;
       h4->Draw(std::string(newTime+":amp_max["+iMCP+"] >> timingCorrection_wrtMiB2").c_str(),Selection3.c_str());
    }
    std::cout << "timeShifted = " << newTime << std::endl;
    std::vector<float> points_wrtMiB2;
    timingCorrection_wrtMiB2->FitSlicesY();
    TH1F* timingCorrection_wrtMiB2_1 = (TH1F*)inputFile->Get("timingCorrection_wrtMiB2_1");
    timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle(std::string("amp_max["+iMCP+"]").c_str());
    if(iMCP != "MiB2") timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time-time[MiB2]");
    else timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time[MiB2]");

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
    fit_corr1->SetParLimits(2,0.,999999999.);
    //fit_corr1 = new TF1("fit_corr1","[0]+[1]*log([2]+x)",0.,3000.);
    timingCorrection_wrtMiB2_1->Fit("fit_corr1","B");
    if(doScan_Corr == false){
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
       Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
    }
    
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
    TH2D* h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    TH2D* h2_time_maximum_vs_amp = new TH2D("h2_time_maximum_vs_amp","",300,0.,3000.,500,0.,5.);

    TF1* pol1_max = new TF1("pol1_max","pol1",0.,3000.);
    TF1* fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",0.,3000.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
      Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,"1",true);
      //Selection = AddSelection(h4,std::string("time_max[MiB2]-time_max[")+iMCP+std::string("]"),Selection,"1",true);
    } 

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
    if(iMCP != "MiB2") wrtMCP = "_wrtMIB2";

    TCanvas* c1 = new TCanvas();
    c1->cd();
    h2_time_max_vs_amp->Draw("COLZ");
    fit_corr_max->Draw("same");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".png").c_str(),"png");
    c1 -> Print(std::string("deltaT_max_vs_amp_max_"+nameiMCP+"_"+Timing+"_thres"+thresMCP+wrtMCP+".pdf").c_str(),"pdf");
}

void PulseShapes(TTree* h4, std::string iMCP, std::string nameiMCP, std::string thresMCP, std::string maxMCP)
{
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",300,-10,20,300,-1.,1.5,100.,3000.);
    TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time","",300,-10,20,300,-1.,1.5);

    std::string Selection;
    std::string iTiming = "";
    //if(Timing != "CFD50") iTiming = "+"+Timing;
  
    if(iMCP == "MiB2"){
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
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
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
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
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
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
      Selection = "amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2;
    }else{ 
      Selection = "amp_max["+iMCP+"]>"+thresMCP+" && amp_max["+iMCP+"]<"+maxMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && fabs(time_max[MiB2])<"+time_max_MiB2+" && X[0]>-999 && Y[0]>-999";
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
    TH1F* h_HV = new TH1F("h_HV","",80,-25,3975.);
    
    h4->Draw(std::string(HV+" >> h_HV").c_str());   
    int hv = 0;
    for(int ii = 1; ii < h_HV->GetNbinsX(); ii++)
        if(h_HV->GetBinContent(ii) > 0)hv = (int)h_HV->GetBinCenter(ii);
    
   return hv;
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


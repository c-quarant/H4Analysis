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

int nBins = 5;
float stepSize = 680.; 

//MiB
//std::string amp_max_MiB2 = "200"
//std::string amp_max_Rm2 = "200"
//std::string b_rms = "4"

//BINP
std::string amp_max_MiB2 = "200";
std::string amp_max_Rm2 = "200";
std::string b_rms = "3.5";

void FinalTiming(TTree* h4, std::string iMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP);
void TimeCorrection(TTree* h4, std::string iMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP);
void PulseShapes(TTree* h4, std::string iMCP);
void TimeMax_vs_Amp(TTree* h4, std::string iMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP);
std::string AddSelection(TTree* h4, std::string Var, std::string Selection, std::string Cut);

void ComputeTiming_oneStep(std::string inputs, std::string iMCP, std::string Timing, std::string thresMCP, bool doScan = false, bool doPulseShapes = false)
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

    if(doPulseShapes == true) PulseShapes(h4, iMCP);
    TimeMax_vs_Amp(h4, iMCP, Timing, Params, thresMCP);
    TimeCorrection(h4, iMCP, inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP);
    FinalTiming(h4, iMCP, inputFile, Timing, Params, Params_wrtMiB2, Params_wrtRm2, thresMCP, doScan);

}

void FinalTiming(TTree* h4, std::string iMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP, bool doScan)
{
    TH1F* time_wrtMiB2 = new TH1F("time_wrtMiB2","",400,-1.,1.);
    TH1F* time_wrtRm2 = new TH1F("time_wrtRm2","",400,-1.,1.);
    TH1F* time_wrtMiB2_noCorrection = new TH1F("time_wrtMiB2_noCorrection","",400,-5.,5.);
    TH1F* time_wrtRm2_noCorrection = new TH1F("time_wrtRm2_noCorrection","",400,-5.,5.);

    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;

    //TF1 *g_res = new TF1("g_res","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",-5.,5.);
    TF1 *g_res;

    char Selection1 [1000];  
    char Selection2 [1000];
    char Selection3 [1000];
    std::string Selection8;
    std::string Selection9;
    if(iMCP == "M25") sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*amp_max[")+iMCP+std::string("]))<0.25")).c_str(),Params->at(0),Params->at(1));
    else sprintf (Selection1,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));
    if(iMCP == "M25" || iMCP == "M10") sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1),Params_wrtMiB2->at(2));
    else sprintf (Selection2,(std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]-(%f+(%f)*amp_max[")+iMCP+std::string("]) >>")).c_str(),Params_wrtMiB2->at(0),Params_wrtMiB2->at(1));  
    if(iMCP == "M25" || iMCP == "M10") sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*1/(%f+amp_max[")+iMCP+std::string("])) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1),Params_wrtRm2->at(2));
    else sprintf (Selection3,(std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]-(%f+(%f)*amp_max[")+iMCP+std::string("]) >>")).c_str(),Params_wrtRm2->at(0),Params_wrtRm2->at(1));  
    Selection8 = "time["+iMCP+iTiming+"]-time[MiB2] >>";
    Selection9 = "time["+iMCP+iTiming+"]-time[Rm2] >>";

    std::string Selection4;
    std::string Selection5;
    std::string Selection6 = Selection2+std::string(" time_wrtMiB2");
    std::string Selection7 = Selection3+std::string(" time_wrtRm2");
    std::string Selection10 = Selection8+std::string(" time_wrtMiB2_noCorrection");
    std::string Selection11 = Selection9+std::string(" time_wrtRm2_noCorrection");

    if(iMCP == "M25") Selection4 = "amp_max[M25]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M25] <= "+b_rms;
    else if(iMCP == "M10") Selection4 = "amp_max[M10]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M10] <= "+b_rms;
    else if(iMCP == "MiB2") Selection4 = "amp_max[MiB2]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_MiB2+" && b_rms[MiB2] <= "+b_rms;
    else if(iMCP == "BINP1") Selection4 = "time_max[BINP1]-time[BINP1"+iTiming+"]<3. && amp_max[BINP1]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP1] <= "+b_rms;
    else if(iMCP == "BINP2") Selection4 = "time_max[BINP2]-time[BINP2"+iTiming+"]<3. && amp_max[BINP2]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP2] <= "+b_rms;
    else if(iMCP == "BINP3") Selection4 = "time_max[BINP3]-time[BINP3"+iTiming+"]<3. && amp_max[BINP3]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP3] <= "+b_rms;
    else if(iMCP == "BINP4") Selection4 = "time_max[BINP4]-time[BINP4"+iTiming+"]<3. && amp_max[BINP4]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP4] <= "+b_rms;

    Selection4 = AddSelection(h4,(std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]")).c_str(),Selection4,std::string("1."));
    Selection4 = AddSelection(h4,(std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]")).c_str(),Selection4,std::string("1."));
    Selection4 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection4,std::string("0.2"));
    Selection5 = Selection4+std::string(" && ")+std::string(Selection1); 

    //std::cout << Selection6 << "," << Selection5 << std::endl;

    h4->Draw(Selection6.c_str(),Selection5.c_str()); 
    h4->Draw(Selection7.c_str(),Selection5.c_str()); 
    h4->Draw(Selection10.c_str(),Selection5.c_str()); 
    h4->Draw(Selection11.c_str(),Selection5.c_str()); 

    g_res = new TF1("g_res","gaus",time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.1);
    g_res->SetParameters(0,0.,time_wrtMiB2->GetEntries());
    g_res->SetParameters(1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
    g_res->SetParameters(2,0.,0.5);
    //g_res->SetParameters(3,0.,time_wrtMiB2->GetEntries());
    //g_res->SetParameters(4,0.,2.);
    g_res->SetParLimits(0,0.,time_wrtMiB2->GetEntries());
    g_res->SetParLimits(1,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())-0.2,time_wrtMiB2->GetBinCenter(time_wrtMiB2->GetMaximumBin())+0.2);
    g_res->SetParLimits(2,0.,0.5);
    //g_res->SetParLimits(3,0.,time_wrtMiB2->GetEntries());
    //g_res->SetParLimits(4,0.,2.);
    time_wrtMiB2->Fit("g_res","B"); 

    TCanvas* c1 = new TCanvas();
    c1->cd();
    time_wrtMiB2->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c1 -> Print((std::string("TimeResolution_MiB_25mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_MiB_25mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c1 -> Print((std::string("TimeResolution_MiB_10mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_MiB_10mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c1 -> Print((std::string("TimeResolution_BINP1_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_BINP1_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c1 -> Print((std::string("TimeResolution_BINP2_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_BINP2_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c1 -> Print((std::string("TimeResolution_BINP3_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_BINP3_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c1 -> Print((std::string("TimeResolution_BINP4_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("TimeResolution_BINP4_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

    delete g_res; 
    g_res = new TF1("g_res","gaus",time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.1);
    g_res->SetParameters(0,0.,time_wrtRm2->GetEntries());
    g_res->SetParameters(1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
    g_res->SetParameters(2,0.,0.5);
    //g_res->SetParameters(3,0.,time_wrtRm2->GetEntries());
    //g_res->SetParameters(4,0.,2.);
    g_res->SetParLimits(0,0.,time_wrtRm2->GetEntries());
    g_res->SetParLimits(1,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())-0.2,time_wrtRm2->GetBinCenter(time_wrtRm2->GetMaximumBin())+0.2);
    g_res->SetParLimits(2,0.,0.5);
    //g_res->SetParLimits(3,0.,time_wrtRm2->GetEntries());
    //g_res->SetParLimits(4,0.,2.);
    time_wrtRm2->Fit("g_res","B");

    TCanvas* c2 = new TCanvas();
    c2->cd();
    time_wrtMiB2->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c2 -> Print((std::string("TimeResolution_MiB_25mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_MiB_25mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c2 -> Print((std::string("TimeResolution_MiB_10mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_MiB_10mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c2 -> Print((std::string("TimeResolution_BINP1_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_BINP1_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c2 -> Print((std::string("TimeResolution_BINP2_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_BINP2_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c2 -> Print((std::string("TimeResolution_BINP3_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_BINP3_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c2 -> Print((std::string("TimeResolution_BINP4_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("TimeResolution_BINP4_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

    delete g_res; 
    g_res = new TF1("g_res","gaus",time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.1);
    g_res->SetParameters(0,0.,time_wrtMiB2_noCorrection->GetEntries());
    g_res->SetParameters(1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
    g_res->SetParameters(2,0.,0.5);
    //g_res->SetParameters(3,0.,time_wrtMiB2_noCorrection->GetEntries());
    //g_res->SetParameters(4,0.,2.);
    g_res->SetParLimits(0,0.,time_wrtMiB2_noCorrection->GetEntries());
    g_res->SetParLimits(1,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())-0.2,time_wrtMiB2_noCorrection->GetBinCenter(time_wrtMiB2_noCorrection->GetMaximumBin())+0.2);
    g_res->SetParLimits(2,0.,0.5);
    //g_res->SetParLimits(3,0.,time_wrtMiB2_noCorrection->GetEntries());
    //g_res->SetParLimits(4,0.,2.);
    time_wrtMiB2_noCorrection->Fit("g_res","B"); 

    TCanvas* c3 = new TCanvas();
    c3->cd();
    time_wrtMiB2_noCorrection->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c3 -> Print((std::string("TimeResolution_MiB_25mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_MiB_25mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c3 -> Print((std::string("TimeResolution_MiB_10mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_MiB_10mu_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c3 -> Print((std::string("TimeResolution_BINP1_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_BINP1_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c3 -> Print((std::string("TimeResolution_BINP2_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_BINP2_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c3 -> Print((std::string("TimeResolution_BINP3_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_BINP3_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c3 -> Print((std::string("TimeResolution_BINP4_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c3 -> Print((std::string("TimeResolution_BINP4_wrtMiB2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }

    delete g_res; 
    g_res = new TF1("g_res","gaus",time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.1);
    g_res->SetParameters(0,0.,time_wrtRm2_noCorrection->GetEntries());
    g_res->SetParameters(1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
    g_res->SetParameters(2,0.,0.5);
    //g_res->SetParameters(3,0.,time_wrtRm2_noCorrection->GetEntries());
    //g_res->SetParameters(4,0.,2.);
    g_res->SetParLimits(0,0.,time_wrtRm2_noCorrection->GetEntries());
    g_res->SetParLimits(1,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())-0.2,time_wrtRm2_noCorrection->GetBinCenter(time_wrtRm2_noCorrection->GetMaximumBin())+0.2);
    g_res->SetParLimits(2,0.,0.5);
    //g_res->SetParLimits(3,0.,time_wrtRm2_noCorrection->GetEntries());
    //g_res->SetParLimits(4,0.,2.);
    time_wrtRm2_noCorrection->Fit("g_res","B");

    TCanvas* c4 = new TCanvas();
    c4->cd();
    time_wrtRm2_noCorrection->Draw("hist");
    g_res->Draw("same");
    if(iMCP == "M25"){
       c4 -> Print((std::string("TimeResolution_MiB_25mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_MiB_25mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c4 -> Print((std::string("TimeResolution_MiB_10mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_MiB_10mu_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c4 -> Print((std::string("TimeResolution_BINP1_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_BINP1_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c4 -> Print((std::string("TimeResolution_BINP2_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_BINP2_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c4 -> Print((std::string("TimeResolution_BINP3_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_BINP3_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c4 -> Print((std::string("TimeResolution_BINP4_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.png")).c_str(),"png");
       c4 -> Print((std::string("TimeResolution_BINP4_wrtRm2_")+Timing+std::string("_thres")+thresMCP+std::string("_noCorrection.pdf")).c_str(),"pdf");
    }

   if(doScan == true){
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
        sprintf (cutMin,"%f",ii*stepSize);
        sprintf (cutMax,"%f",(ii+1)*stepSize);
        
        std::string Selection8;
        Selection8 = Selection5+std::string(" && amp_max[")+iMCP+std::string("]>")+std::string(cutMin)+std::string(" && amp_max[")+iMCP+std::string("]<")+std::string(cutMax);
        Selection6 = Selection2+std::string(Name_wrtMiB2);
        Selection7 = Selection3+std::string(Name_wrtRm2);

        h4->Draw(Selection6.c_str(),Selection8.c_str()); 
        h4->Draw(Selection7.c_str(),Selection8.c_str()); 
       
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];
        if(iMCP == "M25"){
           sprintf (NameOutput_wrtMiB2_png,"TimeResolution_MiB_25mu_wrtMiB2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_MiB_25mu_wrtMiB2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_png,"TimeResolution_MiB_25mu_wrtRm2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_MiB_25mu_wrtRm2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
        }else if(iMCP == "M10"){
           sprintf (NameOutput_wrtMiB2_png,"TimeResolution_MiB_10mu_wrtMiB2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_MiB_10mu_wrtMiB2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_png,"TimeResolution_MiB_10mu_wrtRm2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_MiB_10mu_wrtRm2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
        }else{
           sprintf (NameOutput_wrtMiB2_png,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s.png",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"TimeResolution_%s_wrtMiB2_%d_%s_thres%s.pdf",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());  
           sprintf (NameOutput_wrtRm2_png,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.png",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"TimeResolution_%s_wrtRm2_%d_%s_thres%s.pdf",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        }

        char NameFitAlt_wrtMiB2 [100];
        sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
        resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"gaus",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.1);
        resFitAlt_wrtMiB2[ii]->SetParameters(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtMiB2[ii]->SetParameters(2,0.,2.);
        resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,2.);
        resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
   
        g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt_wrtMiB2[ii]->GetParameter(2));
        g_Res_vs_Amp_wrtMiB2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtMiB2[ii]->GetParError(2),resFitAlt_wrtMiB2[ii]->GetParError(2));

        points_wrtMiB2.push_back(resFitAlt_wrtMiB2[ii]->GetParameter(2));
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        resFitAlt_wrtMiB2[ii]->Draw("same");
        c2 -> Print(NameOutput_wrtMiB2_png,"png");
        c2 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
        delete c2;

        char NameFitAlt_wrtRm2 [100];
        sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
        resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.1);
        resFitAlt_wrtRm2[ii]->SetParameters(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtRm2[ii]->SetParameters(2,0.,2.);
        resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,2.);
        resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
        
        g_Res_vs_Amp_wrtRm2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt_wrtRm2[ii]->GetParameter(2));
        g_Res_vs_Amp_wrtRm2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtRm2[ii]->GetParError(2),resFitAlt_wrtRm2[ii]->GetParError(2));

        points_wrtRm2.push_back(resFitAlt_wrtRm2[ii]->GetParameter(2));
        
        TCanvas* c3 = new TCanvas();
        c3->cd();
        resHist_wrtRm2[ii]->Draw("hist");
        resFitAlt_wrtRm2[ii]->Draw("same");
        c3 -> Print(NameOutput_wrtRm2_png,"png");
        c3 -> Print(NameOutput_wrtRm2_pdf,"pdf");
        delete c3;
    }

    std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    std::sort(points_wrtRm2.begin(),points_wrtRm2.end());
    
    g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-0.1,points_wrtMiB2.at(points_wrtMiB2.size()-1)+0.1);

    TCanvas* c5 = new TCanvas();
    c5->cd();
    g_Res_vs_Amp_wrtMiB2->Draw("AP");
    if(iMCP == "M25"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("TimeResolution_vs_amp_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

    g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#sigma_{t}(ns)");
    g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-0.1,points_wrtRm2.at(points_wrtRm2.size()-1)+0.1);

    TCanvas* c6 = new TCanvas();
    c6->cd();
    g_Res_vs_Amp_wrtRm2->Draw("AP");
    if(iMCP == "M25"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("TimeResolution_vs_amp_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }
   } 
}

void TimeCorrection(TTree* h4, std::string iMCP, TFile* inputFile, std::string Timing, std::vector<float>* Params, std::vector<float>* Params_wrtMiB2, std::vector<float>* Params_wrtRm2, std::string thresMCP)
{
    TH2F* timingCorrection_wrtMiB2 = new TH2F("timingCorrection_wrtMiB2","",30,0.,3000.,4000,-20.,20.);
    TH2F* timingCorrection_wrtRm2 = new TH2F("timingCorrection_wrtRm2","",30,0.,3000.,4000,-20.,20.);
    
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
    /*std::string thresMCP = "20";
    if(Timing == "CFD50" && iMCP == "M10") thresMCP = "50";
    else if(Timing == "LED50") thresMCP = "50";
    else if(Timing == "LED100") thresMCP = "100";
    else if(Timing == "LED150") thresMCP = "150";*/

    std::string Selection1;   
    char Selection2 [1000];
    if(iMCP == "M25") sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*amp_max[")+iMCP+std::string("]))<0.25")).c_str(),Params->at(0),Params->at(1));
    else sprintf (Selection2,(std::string("fabs(time_max[")+iMCP+std::string("]-time[")+iMCP+iTiming+std::string("]-(%f+(%f)*log(%f+amp_max[")+iMCP+std::string("])))<0.25")).c_str(),Params->at(0),Params->at(1),Params->at(2));

    if(iMCP == "M25") Selection1 = "amp_max[M25]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M25] <= "+b_rms;
    else if(iMCP == "M10") Selection1 = "amp_max[M10]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M10] <= "+b_rms;
    else if(iMCP == "MiB2") Selection1 = "amp_max[MiB2]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_MiB2+" && b_rms[MiB2] <= "+b_rms;
    else if(iMCP == "BINP1") Selection1 = "time_max[BINP1]-time[BINP1"+iTiming+"]<3. && amp_max[BINP1]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP1] <= "+b_rms;
    else if(iMCP == "BINP2") Selection1 = "time_max[BINP2]-time[BINP2"+iTiming+"]<3. && amp_max[BINP2]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP2] <= "+b_rms;
    else if(iMCP == "BINP3") Selection1 = "time_max[BINP3]-time[BINP3"+iTiming+"]<3. && amp_max[BINP3]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP3] <= "+b_rms;
    else if(iMCP == "BINP4") Selection1 = "time_max[BINP4]-time[BINP4"+iTiming+"]<3. && amp_max[BINP4]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP4] <= "+b_rms;

    Selection1 = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection1,std::string("0.2"));
    std::string Selection3 = Selection1+" && "+Selection2; 
    h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtMiB2")).c_str(),Selection3.c_str());
    h4->Draw((std::string("time[")+iMCP+iTiming+std::string("]-time[Rm2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtRm2")).c_str(),Selection3.c_str());

    //std::cout << (std::string("time[")+iMCP+iTiming+std::string("]-time[MiB2]:amp_max[")+iMCP+std::string("] >> timingCorrection_wrtMiB2")).c_str() << "," << Selection3.c_str() << std::endl;
    
    timingCorrection_wrtMiB2->FitSlicesY();
    TH2F* timingCorrection_wrtMiB2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtMiB2_1");
    timingCorrection_wrtMiB2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
    timingCorrection_wrtMiB2_1->GetYaxis()->SetTitle("time-time[MiB2]");
    timingCorrection_wrtMiB2_1->SetAxisRange(-15.,15., "Y");
    timingCorrection_wrtMiB2_1->SetMarkerStyle(20);
    timingCorrection_wrtMiB2_1->SetMarkerSize(0.9);
    timingCorrection_wrtMiB2_1->SetMarkerColor(kBlack);
    timingCorrection_wrtMiB2_1->SetLineColor(kBlack);
    
    TF1* fit_corr1;
    if(iMCP == "M25" || iMCP == "M10") fit_corr1 = new TF1("fit_corr1","[0]+[1]*1/(x+[2])",0.,3000.);
    else fit_corr1 = new TF1("fit_corr1","pol1",0.,3000.);
    timingCorrection_wrtMiB2_1->Fit("fit_corr1");
    //Params_wrtMiB2->push_back(fit_corr1->GetParameter(0));
    //Params_wrtMiB2->push_back(fit_corr1->GetParameter(1));
    //if(iMCP == "M25" || iMCP == "M10")  Params_wrtMiB2->push_back(fit_corr1->GetParameter(2));
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    //timingCorrection_wrtMiB2_1->Draw();
    timingCorrection_wrtMiB2->Draw();
    fit_corr1->Draw("same");
    if(iMCP == "M25"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c1 -> Print((std::string("timingCorrection_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

    timingCorrection_wrtRm2->FitSlicesY();
    TH2F* timingCorrection_wrtRm2_1 = (TH2F*)inputFile->Get("timingCorrection_wrtRm2_1");
    timingCorrection_wrtRm2_1->GetXaxis()->SetTitle((std::string("amp_max[")+iMCP+std::string("]")).c_str());
    timingCorrection_wrtRm2_1->GetYaxis()->SetTitle("time-time[Rm2]");
    timingCorrection_wrtRm2_1->SetAxisRange(-15.,15., "Y");
    timingCorrection_wrtRm2_1->SetMarkerStyle(20);
    timingCorrection_wrtRm2_1->SetMarkerSize(0.9);
    timingCorrection_wrtRm2_1->SetMarkerColor(kBlack);
    timingCorrection_wrtRm2_1->SetLineColor(kBlack);

    TF1* fit_corr2;
    if(iMCP == "M25" || iMCP == "M10") fit_corr2 = new TF1("fit_corr2","[0]+[1]*1/(x+[2])",0.,3000.);
    else fit_corr2 = new TF1("fit_corr2","pol1",0.,3000.);
    timingCorrection_wrtRm2_1->Fit("fit_corr2");
    //Params_wrtRm2->push_back(fit_corr2->GetParameter(0));
    //Params_wrtRm2->push_back(fit_corr2->GetParameter(1));
    //if(iMCP == "M25" || iMCP == "M10") Params_wrtRm2->push_back(fit_corr2->GetParameter(2));

    TCanvas* c2 = new TCanvas();
    c2->cd();
    timingCorrection_wrtRm2_1->Draw();
    fit_corr2->Draw("same");
    if(iMCP == "M25"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");   
    }else if(iMCP == "BINP4"){
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c2 -> Print((std::string("timingCorrection_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
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

    TCanvas* c3;
    TCanvas* c4;

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
        sprintf (cutMin,"%f",ii*stepSize);
        sprintf (cutMax,"%f",(ii+1)*stepSize);
        
        std::string Selection4 = Selection3+std::string(" && amp_max[")+iMCP+std::string("]>")+std::string(cutMin)+std::string(" && amp_max[")+iMCP+std::string("]<")+std::string(cutMax);
        std::string Selection5 = "time["+iMCP+iTiming+"]-time[MiB2] >> "+Name_wrtMiB2; 
        std::string Selection6 = "time["+iMCP+iTiming+"]-time[Rm2] >> "+Name_wrtRm2; 
        h4->Draw(Selection5.c_str(),Selection4.c_str());
        h4->Draw(Selection6.c_str(),Selection4.c_str());
      
        char NameOutput_wrtMiB2_png [500];
        char NameOutput_wrtMiB2_pdf [500];
        char NameOutput_wrtRm2_png [500];
        char NameOutput_wrtRm2_pdf [500];
        if(iMCP == "M25"){
           sprintf (NameOutput_wrtMiB2_png,"timingCorrection_MiB_25mu_wrtMiB2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"timingCorrection_MiB_25mu_wrtMiB2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_png,"timingCorrection_MiB_25mu_wrtRm2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"timingCorrection_MiB_25mu_wrtRm2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
        }else if(iMCP == "M10"){
           sprintf (NameOutput_wrtMiB2_png,"timingCorrection_MiB_10mu_wrtMiB2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"timingCorrection_MiB_10mu_wrtMiB2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_png,"timingCorrection_MiB_10mu_wrtRm2_%d_%s_thres%s.png",ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"timingCorrection_MiB_10mu_wrtRm2_%d_%s_thres%s.pdf",ii+1,Timing.c_str(),thresMCP.c_str());
        }else{
           sprintf (NameOutput_wrtMiB2_png,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s.png",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtMiB2_pdf,"timingCorrection_%s_wrtMiB2_%d_%s_thres%s.pdf",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());  
           sprintf (NameOutput_wrtRm2_png,"timingCorrection_%s_wrtRm2_%d_%s_thres%s.png",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
           sprintf (NameOutput_wrtRm2_pdf,"timingCorrection_%s_wrtRm2_%d_%s_thres%s.pdf",iMCP.c_str(),ii+1,Timing.c_str(),thresMCP.c_str());
        }

        char NameFitAlt_wrtMiB2 [100];
        sprintf (NameFitAlt_wrtMiB2,"f_ResAlt_2_%d_wrtMiB2",ii);
        resFitAlt_wrtMiB2[ii] = new TF1(NameFitAlt_wrtMiB2,"gaus",resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.3,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.3);
        resFitAlt_wrtMiB2[ii]->SetParameters(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtMiB2[ii]->SetParameters(2,0.,2.);
        resFitAlt_wrtMiB2[ii]->SetParLimits(1,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())-0.2,resHist_wrtMiB2[ii]->GetBinCenter(resHist_wrtMiB2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtMiB2[ii]->SetParLimits(2,0.,2.);
        resHist_wrtMiB2[ii]->Fit(NameFitAlt_wrtMiB2,"B");
   
        g_Res_vs_Amp_wrtMiB2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt_wrtMiB2[ii]->GetParameter(1));
        g_Res_vs_Amp_wrtMiB2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtMiB2[ii]->GetParError(1),resFitAlt_wrtMiB2[ii]->GetParError(1));

        points_wrtMiB2.push_back(resFitAlt_wrtMiB2[ii]->GetParameter(1));
        
        c3 = new TCanvas();
        c3->cd();
        resHist_wrtMiB2[ii]->Draw("hist");
        resFitAlt_wrtMiB2[ii]->Draw("same");
        c3 -> Print(NameOutput_wrtMiB2_png,"png");
        c3 -> Print(NameOutput_wrtMiB2_pdf,"pdf");
        delete c3;

        char NameFitAlt_wrtRm2 [100];
        sprintf (NameFitAlt_wrtRm2,"f_ResAlt_2_%d_wrtRm2",ii);
        resFitAlt_wrtRm2[ii] = new TF1(NameFitAlt_wrtRm2,"gaus",resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.3,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.3);
        resFitAlt_wrtRm2[ii]->SetParameters(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtRm2[ii]->SetParameters(2,0.,2.);
        resFitAlt_wrtRm2[ii]->SetParLimits(1,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())-0.2,resHist_wrtRm2[ii]->GetBinCenter(resHist_wrtRm2[ii]->GetMaximumBin())+0.2);
        resFitAlt_wrtRm2[ii]->SetParLimits(2,0.,2.);
        resHist_wrtRm2[ii]->Fit(NameFitAlt_wrtRm2,"B");
        
        g_Res_vs_Amp_wrtRm2->SetPoint(ii,ii*stepSize+stepSize/2.,resFitAlt_wrtRm2[ii]->GetParameter(1));
        g_Res_vs_Amp_wrtRm2->SetPointError(ii,stepSize/2.,stepSize/2.,resFitAlt_wrtRm2[ii]->GetParError(1),resFitAlt_wrtRm2[ii]->GetParError(1));

        points_wrtRm2.push_back(resFitAlt_wrtRm2[ii]->GetParameter(2));
        
        c4 = new TCanvas();
        c4->cd();
        resHist_wrtRm2[ii]->Draw("hist");
        resFitAlt_wrtRm2[ii]->Draw("same");
        c4 -> Print(NameOutput_wrtRm2_png,"png");
        c4 -> Print(NameOutput_wrtRm2_pdf,"pdf");
        delete c4;
    }

    std::sort(points_wrtMiB2.begin(),points_wrtMiB2.end());
    std::sort(points_wrtRm2.begin(),points_wrtRm2.end());
    
    g_Res_vs_Amp_wrtMiB2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetTitle("#mu (ns)");
    g_Res_vs_Amp_wrtMiB2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtMiB2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtMiB2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtMiB2->GetYaxis()->SetRangeUser(points_wrtMiB2.at(0)-0.1,points_wrtMiB2.at(points_wrtMiB2.size()-1)+0.1);

    TF1* fit_corr3;
    if(iMCP == "M25" || iMCP == "M10") fit_corr3 = new TF1("fit_corr3","[0]+[1]*1/(x+[2])",0.,3000.);
    else fit_corr3 = new TF1("fit_corr3","pol1",0.,3000.);
    g_Res_vs_Amp_wrtMiB2->Fit("fit_corr3");
    Params_wrtMiB2->push_back(fit_corr3->GetParameter(0));
    Params_wrtMiB2->push_back(fit_corr3->GetParameter(1));
    if(iMCP == "M25" || iMCP == "M10")  Params_wrtMiB2->push_back(fit_corr3->GetParameter(2));

    TCanvas* c5 = new TCanvas();
    c5->cd();
    g_Res_vs_Amp_wrtMiB2->Draw("AP");
    fit_corr3->Draw("same");
    if(iMCP == "M25"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c5 -> Print((std::string("timingCorrection2_wrtMiB2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

    g_Res_vs_Amp_wrtRm2->GetXaxis()->SetTitle("amp_max");
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetTitle("#mu(ns)");
    g_Res_vs_Amp_wrtRm2->SetMarkerStyle(20);
    g_Res_vs_Amp_wrtRm2->SetMarkerSize(0.7);
    g_Res_vs_Amp_wrtRm2->SetMarkerColor(kBlack);
    g_Res_vs_Amp_wrtRm2->SetLineColor(kBlack);
    g_Res_vs_Amp_wrtRm2->GetYaxis()->SetRangeUser(points_wrtRm2.at(0)-0.1,points_wrtRm2.at(points_wrtRm2.size()-1)+0.1);

    TF1* fit_corr4;
    if(iMCP == "M25" || iMCP == "M10") fit_corr4 = new TF1("fit_corr4","[0]+[1]*1/(x+[2])",0.,3000.);
    else fit_corr4 = new TF1("fit_corr4","pol1",0.,3000.);
    g_Res_vs_Amp_wrtRm2->Fit("fit_corr4");
    Params_wrtRm2->push_back(fit_corr4->GetParameter(0));
    Params_wrtRm2->push_back(fit_corr4->GetParameter(1));
    if(iMCP == "M25" || iMCP == "M10")  Params_wrtRm2->push_back(fit_corr4->GetParameter(2));

    TCanvas* c6 = new TCanvas();
    c6->cd();
    g_Res_vs_Amp_wrtRm2->Draw("AP");
    fit_corr4->Draw("same");
    if(iMCP == "M25"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("timingCorrection2_wrtRm2_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }

}

void PulseShapes(TTree* h4, std::string iMCP)
{
  
    TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time","",2250,-50,175,300,-1.,1.5,1000.,1500.);
    TProfile2D* p2D_amp_vs_time_zoom1 = new TProfile2D("p2D_amp_vs_time_zoom1","",300,-10,20,300,-1.,1.5,1000.,1500.);
    TProfile2D* p2D_amp_vs_time_zoom2 = new TProfile2D("p2D_amp_vs_time_zoom2","",300,-10,20,160,0.5,1.1,1000.,1500.);

    std::string Selection;
    std::string thresMCP = "20";
    std::string iTiming = "";
    
    if(iMCP == "M25") Selection = "amp_max[M25]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M25] <= "+b_rms;
    else if(iMCP == "M10") Selection = "amp_max[M10]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M10] <= "+b_rms;
    else if(iMCP == "MiB2") Selection = "amp_max[MiB2]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_MiB2+" && b_rms[MiB2] <= "+b_rms;
    else if(iMCP == "BINP1") Selection = "time_max[BINP1]-time[BINP1"+iTiming+"]<3. && amp_max[BINP1]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP1] <= "+b_rms;
    else if(iMCP == "BINP2") Selection = "time_max[BINP2]-time[BINP2"+iTiming+"]<3. && amp_max[BINP2]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP2] <= "+b_rms;
    else if(iMCP == "BINP3") Selection = "time_max[BINP3]-time[BINP3"+iTiming+"]<3. && amp_max[BINP3]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP3] <= "+b_rms;
    else if(iMCP == "BINP4") Selection = "time_max[BINP4]-time[BINP4"+iTiming+"]<3. && amp_max[BINP4]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP4] <= "+b_rms;
    
    Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,std::string("1."));
    Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,std::string("1."));
    Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,std::string("0.2"));

    Selection = Selection+" && WF_ch == "+iMCP;

    //std::cout << Selection << std::endl; 

    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time")).c_str(),Selection.c_str(),"goff");
    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time_zoom1")).c_str(),Selection.c_str(),"goff");
    h4->Draw((std::string("amp_max[")+iMCP+std::string("]:WF_val/amp_max[")+iMCP+std::string("]:WF_time-time[")+iMCP+std::string("] >> p2D_amp_vs_time_zoom2")).c_str(),Selection.c_str(),"goff");
    
    p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
    p2D_amp_vs_time_zoom1->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time_zoom1->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time_zoom1->GetZaxis()->SetTitle("amp_max");
    p2D_amp_vs_time_zoom2->GetXaxis()->SetTitle((std::string("WF_time-time[")+iMCP+std::string("] (ns)")).c_str());
    p2D_amp_vs_time_zoom2->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+iMCP+std::string("]")).c_str());
    p2D_amp_vs_time_zoom2->GetZaxis()->SetTitle("amp_max");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    p2D_amp_vs_time->Draw("COLZ");
    if(iMCP == "M25"){
       c1 -> Print("pulseShape_MiB_25mu.png","png");
       c1 -> Print("pulseShape_MiB_25mu.pdf","pdf");
    }else if(iMCP == "M10"){
       c1 -> Print("pulseShape_MiB_10mu.png","png");
       c1 -> Print("pulseShape_MiB_10mu.pdf","pdf");
    }else if(iMCP == "BINP1"){
       c1 -> Print("pulseShape_BINP1.png","png");
       c1 -> Print("pulseShape_BINP1.pdf","pdf");
    }else if(iMCP == "BINP2"){
       c1 -> Print("pulseShape_BINP2.png","png");
       c1 -> Print("pulseShape_BINP2.pdf","pdf");
    }else if(iMCP == "BINP3"){
       c1 -> Print("pulseShape_BINP3.png","png");
       c1 -> Print("pulseShape_BINP3.pdf","pdf");
    }else if(iMCP == "BINP4"){
       c1 -> Print("pulseShape_BINP4.png","png");
       c1 -> Print("pulseShape_BINP4.pdf","pdf");
    }

    TCanvas* c2 = new TCanvas();
    c2->cd();
    p2D_amp_vs_time_zoom1->Draw("COLZ");
    if(iMCP == "M25"){
       c2 -> Print("pulseShape_MiB_25mu_zoom1.png","png");
       c2 -> Print("pulseShape_MiB_25mu_zoom1.pdf","pdf");
    }else if(iMCP == "M10"){
       c2 -> Print("pulseShape_MiB_10mu_zoom1.png","png");
       c2 -> Print("pulseShape_MiB_10mu_zoom1.pdf","pdf");
    }else if(iMCP == "BINP1"){
       c2 -> Print("pulseShape_BINP1_zoom1.png","png");
       c2 -> Print("pulseShape_BINP1_zoom1.pdf","pdf");
    }else if(iMCP == "BINP2"){
       c2 -> Print("pulseShape_BINP2_zoom1.png","png");
       c2 -> Print("pulseShape_BINP2_zoom1.pdf","pdf");
    }else if(iMCP == "BINP3"){
       c2 -> Print("pulseShape_BINP3_zoom1.png","png");
       c2 -> Print("pulseShape_BINP3_zoom1.pdf","pdf");
    }else if(iMCP == "BINP4"){
       c2 -> Print("pulseShape_BINP4_zoom1.png","png");
       c2 -> Print("pulseShape_BINP4_zoom1.pdf","pdf");
    }

    TCanvas* c3 = new TCanvas();
    c3->cd();
    p2D_amp_vs_time_zoom2->Draw("COLZ");
    if(iMCP == "M25"){
       c3 -> Print("pulseShape_MiB_25mu_zoom2.png","png");
       c3 -> Print("pulseShape_MiB_25mu_zoom2.pdf","pdf");
    }else if(iMCP == "M10"){
       c3 -> Print("pulseShape_MiB_10mu_zoom2.png","png");
       c3 -> Print("pulseShape_MiB_10mu_zoom2.pdf","pdf");
    }else if(iMCP == "BINP1"){
       c3 -> Print("pulseShape_BINP1_zoom2.png","png");
       c3 -> Print("pulseShape_BINP1_zoom2.pdf","pdf");
    }else if(iMCP == "BINP2"){
       c3 -> Print("pulseShape_BINP2_zoom2.png","png");
       c3 -> Print("pulseShape_BINP2_zoom2.pdf","pdf");
    }else if(iMCP == "BINP3"){
       c3 -> Print("pulseShape_BINP3_zoom2.png","png");
       c3 -> Print("pulseShape_BINP3_zoom2.pdf","pdf");
    }else if(iMCP == "BINP4"){
       c3 -> Print("pulseShape_BINP4_zoom2.png","png");
       c3 -> Print("pulseShape_BINP4_zoom2.pdf","pdf");
    }

}

void TimeMax_vs_Amp(TTree* h4, std::string iMCP, std::string Timing, std::vector<float>* Params, std::string thresMCP)
{
    TH2D* h2_time_max_vs_amp = new TH2D("h2_time_max_vs_amp","",300,0.,3000.,500,0.,5.);
    TH2D* h2_time_maximum_vs_amp = new TH2D("h2_time_maximum_vs_amp","",300,0.,3000.,500,0.,5.);

    TF1* pol1_max = new TF1("pol1_max","pol1",0.,3000.);
    TF1* fit_corr_max = new TF1("fit_corr_max","[0]+[1]*log(x+[2])",0.,3000.);

    std::string Selection;
    std::string iTiming = "";
    if(Timing != "CFD50") iTiming = "+"+Timing;
    /*std::string thresMCP = "20";
    if(Timing == "CFD50" && iMCP == "M10") thresMCP = "50";
    else if(Timing == "LED50") thresMCP = "50";
    else if(Timing == "LED100") thresMCP = "100";
    else if(Timing == "LED150") thresMCP = "150";*/

    if(iMCP == "M25") Selection = "amp_max[M25]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M25] <= "+b_rms;
    else if(iMCP == "M10") Selection = "amp_max[M10]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[M10] <= "+b_rms;
    else if(iMCP == "MiB2") Selection = "amp_max[MiB2]>"+thresMCP+" && amp_max[Rm2]>"+amp_max_MiB2+" && b_rms[MiB2] <= "+b_rms;
    else if(iMCP == "BINP1") Selection = "time_max[BINP1]-time[BINP1"+iTiming+"]<3. && amp_max[BINP1]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP1] <= "+b_rms;
    else if(iMCP == "BINP2") Selection = "time_max[BINP2]-time[BINP2"+iTiming+"]<3. && amp_max[BINP2]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP2] <= "+b_rms;
    else if(iMCP == "BINP3") Selection = "time_max[BINP3]-time[BINP3"+iTiming+"]<3. && amp_max[BINP3]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP3] <= "+b_rms;
    else if(iMCP == "BINP4") Selection = "time_max[BINP4]-time[BINP4"+iTiming+"]<3. && amp_max[BINP4]>"+thresMCP+" && amp_max[MiB2]>"+amp_max_MiB2+" && amp_max[Rm2]>"+amp_max_Rm2+" && b_rms[BINP4] <= "+b_rms;
    
    Selection = AddSelection(h4,std::string("time[MiB2]-time[")+iMCP+iTiming+std::string("]"),Selection,std::string("1."));
    Selection = AddSelection(h4,std::string("time[Rm2]-time[")+iMCP+iTiming+std::string("]"),Selection,std::string("1."));
    Selection = AddSelection(h4,std::string("time[MiB2]-time[Rm2]"),Selection,std::string("0.2"));
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
    if(iMCP == "M25") h2_time_max_vs_amp->Fit("pol1_max");
    else if(Timing == "CFD50") h2_time_max_vs_amp->Fit("fit_corr_max","B","",20.,3000.);
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
    
    if(iMCP == "M25"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB_25mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "M10"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB_10mu_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "MiB2"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_MiB2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP1"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP1_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP2"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP2_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP3"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP3_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }else if(iMCP == "BINP4"){
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".png")).c_str(),"png");
       c6 -> Print((std::string("deltaT_max_vs_amp_max_BINP4_")+Timing+std::string("_thres")+thresMCP+std::string(".pdf")).c_str(),"pdf");
    }
}

std::string AddSelection(TTree* h4, std::string Var, std::string Selection, std::string Cut)
{
    TH1F* h = new TH1F("h","h",4000,-20.,20.);
    h4->Draw((Var+std::string(" >> h")).c_str(),Selection.c_str());
    //std::cout << "Draw(" << (Var+std::string(" >> h")).c_str() << "," << Selection.c_str() << ")" << std::endl;
    char Mean [100];
    sprintf(Mean,"%f",h->GetMean());
    std::string sMean = std::string(Mean);

    if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+std::string(" && fabs(")+Var+std::string("+")+std::string(sMean)+std::string(")<")+Cut;
    }else{
       Selection = Selection+std::string(" && fabs(")+Var+std::string("-")+std::string(sMean)+std::string(")<")+Cut;
    }

    //std::cout << "Selection = " << Selection << std::endl;

    delete h;
    return Selection;
}
	

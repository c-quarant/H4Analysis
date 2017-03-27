#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

#include<iostream>
#include<string>
#include<fstream>

void draw_PulseShape()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    gStyle->SetOptFit(1); 
    gStyle->SetErrorX(0);

    TFile* inputFile = TFile::Open("root://eoscms.cern.ch//store/group/dpg_ecal/comm_ecal/upgrade/testbeam/TimingTB_BTF_Jun2016/ntuples/v4/ntuples_v4/btf2016_RU4_2363.root");
    TTree* h4 = (TTree*)inputFile->Get("h4");
    TTree* digi = (TTree*)inputFile->Get("digi");
    TTree* wf = (TTree*)inputFile->Get("wf");

    // Declaration of leaf types
    ULong64_t       index;
    ULong64_t       start_time;
    ULong64_t       time_stamp;
    UInt_t          run;
    UInt_t          spill;
    UInt_t          event;
    int           BINP3;
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
    int           WF_samples;
    int           WF_ch[6144];   //[WF_samples]
    float         WF_time[6144];   //[WF_samples]
    float         WF_val[6144];   //[WF_samples]

    // List of branches
    TBranch        *b_start_time;   //!
    TBranch        *b_time_stamp;   //!
    TBranch        *b_run;   //!
    TBranch        *b_spill;   //!
    TBranch        *b_event;   //!
    TBranch        *b_BINP3;   //!
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
    TBranch        *b_WF_samples;   //!
    TBranch        *b_WF_ch;   //!
    TBranch        *b_WF_time;   //!
    TBranch        *b_WF_val;   //!

    h4->SetBranchAddress("start_time", &start_time, &b_start_time);
    h4->SetBranchAddress("time_stamp", &time_stamp, &b_time_stamp);
    h4->SetBranchAddress("run", &run, &b_run);
    h4->SetBranchAddress("spill", &spill, &b_spill);
    h4->SetBranchAddress("event", &event, &b_event);
    digi->SetBranchAddress("BINP3", &BINP3, &b_BINP3);
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
    wf->SetBranchAddress("WF_samples", &WF_samples, &b_WF_samples);
    wf->SetBranchAddress("WF_ch", WF_ch, &b_WF_ch);
    wf->SetBranchAddress("WF_time", WF_time, &b_WF_time);
    wf->SetBranchAddress("WF_val", WF_val, &b_WF_val);

    TH1F* pulseShape = new TH1F("pulseShape","",150,-10,20);
    TGraph* g_pulseShape = new TGraph();

    for(int entry = 0; entry < h4->GetEntries(); entry++){

        if(entry%1000==0) std::cout<<"--- Reading entry = "<<entry<<std::endl;
        h4->GetEntry(entry);
        digi->GetEntry(entry);
        wf->GetEntry(entry);
            
        //if(event != 722 || spill != 4) continue;
        if(event != 58 || spill != 4) continue;
        for(int ii = 0; ii < WF_samples; ii++)
        {
            if(WF_ch[ii] != BINP3) continue;
            pulseShape->SetBinContent(pulseShape->FindBin(WF_time[ii]-time[BINP3]),WF_val[ii]);

            //std::cout << "amp_max" << amp_max[BINP3] << std::endl;
        }
    }

    for(int ii = 1; ii < pulseShape->GetNbinsX(); ii++)
        g_pulseShape->SetPoint(ii-1,pulseShape->GetBinCenter(ii),pulseShape->GetBinContent(ii));

    g_pulseShape->SetLineWidth(2);

    setStyle();

    //TH2F* H2 = new TH2F("H2","",100,-5,15,1800,-300,1500);
    TH2F* H2 = new TH2F("H2","",100,-5,15,1800,-300.,2200.);
    H2->GetXaxis()->SetTitle("t-t_{CFD} (ns)");
    H2->GetYaxis()->SetTitle("ADC counts");

    //TLine *vLine1 = new TLine(0,-280,0,420);
    TLine *vLine1 = new TLine(0,-280,0,1000.);
    vLine1->SetLineColor(kGreen+2);
    vLine1->SetLineWidth(2);

    //TLine *hLine1 = new TLine(-5,420,0,420);
    TLine *hLine1 = new TLine(-5,1000.,0,1000.);
    hLine1->SetLineColor(kGreen+2);
    hLine1->SetLineWidth(2);

    //TLine *vLine2 = new TLine(0.75,-300,0.75,1200);
    TLine *vLine2 = new TLine(0.73,-300,0.73,2000.);
    vLine2->SetLineColor(kRed+1);
    vLine2->SetLineWidth(2);

    //TLine *hLine2 = new TLine(-5,1200,0.75,1200);
    TLine *hLine2 = new TLine(-5,2000,0.75,2000);
    hLine2->SetLineColor(kRed+1);
    hLine2->SetLineWidth(2);

    TGraphAsymmErrors* g_tmp1 = new TGraphAsymmErrors();
    g_tmp1->SetLineColor(kGreen+2);
    g_tmp1->SetMarkerColor(kGreen+2);
    g_tmp1->SetLineWidth(2);

    TGraphAsymmErrors* g_tmp2 = new TGraphAsymmErrors();
    g_tmp2->SetLineColor(kRed+1);
    g_tmp2->SetMarkerColor(kRed+1);
    g_tmp2->SetLineWidth(2);

    TLegend* legend = new TLegend(0.45, 0.67, 0.79, 0.79);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(g_tmp1,"Time at half amplitude (CFD)","L");
    legend -> AddEntry(g_tmp2,"Time at maximum amplitude","L");

    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    H2->Draw();
    g_pulseShape->Draw("C");
    legend->Draw("same");
    vLine1->Draw("same");
    hLine1->Draw("same");
    vLine2->Draw("same");
    hLine2->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw();
    c1 -> Print(std::string("pulseShape.png").c_str(),"png");
    c1 -> Print(std::string("pulseShape.pdf").c_str(),"pdf");
}

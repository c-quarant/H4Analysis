#include "TLeaf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "MyLib.h"
#include <string>
#include <fstream>
#include <math.h>


void WFPulseShape(std::string FileIn, std::string detector, Float_t bound, std::string runNum, std::string MCP)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, hodoCfg;
	Float_t X[2], Y[2], Ampl[6];
	Float_t X1true, Y1true, Xshift, Yshift, Xavg, Yavg;
	Float_t Gain_val, Energy_val;
	float XMax, YMax;

	float fitRange=5;
	//setStyle();
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();

	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = "E"+Energy+"_G"+Gain;

	//Calculating shift of hodo plane 1 respect to hodo plane 0 faro X and Y axis
	Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");

	std::string AmplMCPSel = "amp_max[MCP1]>200 && amp_max[MCP1]<2000";
	AmplitudeProfilesFit(h4, detector, AmplMCPSel, pathToOutput, RunStats, bound, &XMax, &YMax, Xshift, Yshift, "amp_max");
	AmplitudeMaps(h4, detector, AmplMCPSel, std::to_string(Xshift), std::to_string(Yshift), XMax, YMax, bound, RunStats, "amp_max");
	EfficiencyMaps(h4, detector, AmplMCPSel, std::to_string(Xshift), std::to_string(Yshift), XMax, YMax, bound, RunStats);

	if(fabs(XMax)<100 && fabs(YMax)<100)
	{
		//PulseShapes(h4, detector, 0, XMax, YMax, bound, pathToOutput, RunStats, runNum, MCP);
		WFPulseShapes(h4, detector, 0, XMax, YMax, bound, pathToOutput, RunStats, runNum, MCP);
	}
	else cout<< "SBALLATO!" << endl;
	
}





//TProfile2D* p2D_amp_vs_time_no_cut = new TProfile2D("p2D_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	//TH2F* h2_amp_vs_time_no_cut = new TH2F("h2_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5);

/*
	//Drawing WITHOUT CUTS on HODOSCOPE
	Selection = "amp_max[MCP1]>200 && fabs(time[MCP1]-("+TimeMCP+"))<8 && X[0]>-800 && Y[0]>-800";

	h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+")>>p2D_amp_vs_time_no_cut").c_str(),Selection.c_str());
        //std::cout << (std::string("amp_max[")+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+") >> p2D_amp_vs_time") << std::endl;	
	//cout << Selection << endl;
	cout << "Draw 1" << endl;

	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time[MCP1]"+TimeShift+" >> h2_amp_vs_time_no_cut").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

    	TProfile* waveForm_no_cut = h2_amp_vs_time_no_cut->ProfileX(); 
    	waveForm_no_cut->SetName(std::string(detector+std::string("_waveform_prof_no_cut")).c_str());

	p2D_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	h2_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	waveForm_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");	
	*/

/*
	TCanvas* c4 = new TCanvas();
    	c5->cd();
    	h2_amp_vs_time_no_cut->Draw("COLZ");
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.png").c_str());
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.pdf").c_str());

    	TCanvas* c5 = new TCanvas();
    	c6->cd();
    	p2D_amp_vs_time_no_cut->Draw("COLZ");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.png").c_str(),"png");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.pdf").c_str(),"pdf");

    	TCanvas* c6 = new TCanvas();
    	c7->cd();
    	waveForm_no_cut->Draw("P");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.png").c_str(),"png");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.pdf").c_str(),"pdf");
	*/




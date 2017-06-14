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
#include "AmplitudeHist.C"
#include "AmplitudeProfilesFit.C"
#include <string>
#include <fstream>
#include <math.h>

void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats);
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats);
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats);

void AveragePulseShape(std::string FileIn, std::string detector, Float_t bound)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, runNum, hodoCfg;
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
	runNum = h4->GetLeaf("run")->GetValue(0);

	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain+"_"+std::to_string(runNum);

	AmplitudeProfilesFit(FileIn, detector, bound, &XMax, &YMax);
	
	if(fabs(XMax)<3 && fabs(YMax)<3)
	{
		PulseShapes(h4, detector, 0, XMax, YMax, 2, pathToOutput, RunStats);
	}
	else cout<< "SBALLATO!" << endl;
	
}


void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats)
{
	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpThresh;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max[MCP1]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats));
	Selection = Selection + " && fabs(time[MCP1]-("+TimeMCP+"))<7";

	AmpThresh = std::to_string(AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats));	
	Selection = Selection + " && amp_max["+detector+"]>"+AmpThresh;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats));
	Selection = "WF_ch == " + detector + " && " + Selection;	

	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

    	TProfile* waveForm = h2_amp_vs_time->ProfileX(); 
    	waveForm->SetName("Waveform");
	
	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	TCanvas* c0 = new TCanvas();
    	c0->cd();
    	aSlices[1]->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"TestFitSlices/PS_"+detector+"_"+RunStats+"_h2_fit_slicesY.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"TestFitSlices/PS_"+detector+"_"+RunStats+"_h2_fit_slicesY.pdf").c_str());

	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	waveForm->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");	
    	    	
	gStyle->SetOptStat(0);
	
    	TCanvas* c1 = new TCanvas();
    	c1->cd();
    	h2_amp_vs_time->Draw("COLZ");
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_Amp.pdf").c_str(),"pdf");
	
    	TCanvas* c3 = new TCanvas();
    	c3->cd();
    	waveForm->Draw("P");
    	c3 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_profile.png").c_str(),"png");
    	c3 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_profile.pdf").c_str(),"pdf");
	
    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms");
	h2_amp_vs_time->SetName("H2Waveforms");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close(); 
}

float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats)
{
	TH1F* raw_time_dist = new TH1F("raw_time_dist", "", 200, -3, 3);

	h4->Draw((std::string("time["+detector+"]-time[MCP1]>>raw_time_dist")).c_str(), Selection.c_str());

	raw_time_dist->GetXaxis()->SetTitle((std::string("time["+detector+"]-time[MCP1] (ns)")).c_str());
	raw_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* c0 = new TCanvas();
    	c0->cd();
	raw_time_dist->Fit("gaus", "", "", -3, 3);
    	raw_time_dist->Draw();
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist_"+detector+"_"+RunStats+".png").c_str());
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist_"+detector+"_"+RunStats+".pdf").c_str());

	cout << raw_time_dist->GetMean() << endl;
	cout << raw_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return raw_time_dist->GetFunction("gaus")->GetMaximumX();
}

float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats)
{
	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);

	h4->Draw((std::string("time[MCP1]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle("time[MCP1] (ns)");
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas();
    	ca->cd();
	MCP_time_dist->Fit("gaus", "", "", 0, 50);
    	MCP_time_dist->Draw();
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/MCP_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/MCP_TimeDist_"+RunStats+".pdf").c_str());

	cout << MCP_time_dist->GetMean() << endl;
	cout << MCP_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return MCP_time_dist->GetFunction("gaus")->GetMaximumX();
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




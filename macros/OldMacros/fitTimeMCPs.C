#include "TFile.h"
#include "TStyle.h"
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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "MyLib.h"
#include <string>
#include <fstream>

GaussPar fitTimeMCPs(std::string FileIn)
{
	GaussPar TimePar;
	std::string TimeMCP1, TimeMCP2, TimeShift; 

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");

	gStyle->SetOptFit();
	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain;

	//plot detector time distribution
	std::string Selection = "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && amp_max[MCP2]>200 && amp_max[MCP2]<2000";
	
	TimeMCP1 = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, "MCP1"));
	TimeMCP2 = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, "MCP2"));
	Selection = Selection + " && fabs(time[MCP1]-("+TimeMCP1+"))<7 && fabs(time[MCP2]-("+TimeMCP2+"))<7";

	TH1F* tD = new TH1F("tD", "", 2500, -20, 20);
	h4->Draw("(time[MCP1]-time[MCP2])>>tD", Selection.c_str());

	
	float Xfirst = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())-0.5;
	float Xlast = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())+0.5;

	TCanvas* c0 = new TCanvas("c0", "c0");
	tD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD->GetXaxis()->SetTitle("time[MCP1]-time[MCP2] (ns)");
	tD->GetYaxis()->SetTitle("events");
	tD->Fit("gaus", "", "", Xfirst, Xlast);
	tD->Draw();

	TimePar.Mean = tD->GetFunction("gaus")->GetParameter(1);
	TimePar.MeanErr = tD->GetFunction("gaus")->GetParError(1);
	TimePar.Sigma = tD->GetFunction("gaus")->GetParameter(2);
	TimePar.SigmaErr = tD->GetFunction("gaus")->GetParError(2); 
	
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_MCP1-MCP2_"+RunStats+".png").c_str());
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_MCP1-MCP2_"+RunStats+".pdf").c_str());

	return TimePar;
}
	

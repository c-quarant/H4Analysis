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

void PlotSelectedEventsSep(std::string NtupleFile, std::string detector, int bound, std::string MCP)
{

	float XMax, YMax, Xshift, Yshift, AmpMean, AmpSigma;
	std::string TimeMCP, TimeShift, AmpMean_str, AmpSigma_str; 
	TFile *f = TFile::Open(NtupleFile.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain;

	//plot detector time distribution
	Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");
	
	AmplitudeProfilesFit(h4, detector, pathToOutput, RunStats, bound, &XMax, &YMax);

	AmplitudeHist(h4, detector, "", pathToOutput + "SelectedEvents/", "NoCut_"+RunStats, &AmpMean, &AmpSigma, "A");
	cout << "ciao" << endl;
	std::string Selection = "(fabs(X[0]-("+std::to_string(XMax)+"))<"+std::to_string(bound)+" || fabs(X[1]-("+std::to_string(XMax)+")-("+std::to_string(Xshift)+"))<"+std::to_string(bound)+") && (fabs(Y[0]-("+std::to_string(YMax)+"))<"+std::to_string(bound)+" || fabs(Y[1]-("+std::to_string(YMax)+")-("+std::to_string(Yshift)+"))<"+std::to_string(bound)+") && amp_max["+MCP+"]>100";

	AmplitudeHist(h4, detector, Selection, pathToOutput + "SelectedEvents/", "PosCut_"+RunStats, &AmpMean, &AmpSigma, "B");
 
	Selection = "amp_max["+MCP+"]>200"; 

	AmplitudeHist(h4, detector, Selection, pathToOutput + "SelectedEvents/", "AmpMCPCut_"+RunStats, &AmpMean, &AmpSigma, "C");
	
	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, MCP));
	Selection = "fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput + "SelectedEvents/", "TimeMCPcut"+RunStats, &AmpMean, &AmpSigma, "D");

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma, "GOOD");
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = "fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;	
	//Selection = "WF_ch == " + detector + " && " + Selection;

	AmplitudeHist(h4, detector, Selection, pathToOutput + "SelectedEvents/", "AmpDetCut_"+RunStats, &AmpMean, &AmpSigma, "E");	
	
	cout << Selection << endl;

	TH1F* tD = new TH1F("tD", "", 2000, -20, 20);
	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]>>tD").c_str(), Selection.c_str());
	
	float Xfirst = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())-1;
	float Xlast = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())+1;

	TCanvas* c0 = new TCanvas("c0", "c0");
	tD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD->GetXaxis()->SetTitle(("time["+detector+"]-time["+MCP+"] (ns)").c_str());
	tD->GetYaxis()->SetTitle("events");
	tD->Fit("gaus", "", "", Xfirst, Xlast);
	tD->Draw();
}

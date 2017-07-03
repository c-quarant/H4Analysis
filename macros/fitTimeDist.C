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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "MyLib.h"
#include <string>
#include <fstream>

void fitTimeDist(std::string FileIn, std::string detector, Float_t bound, std::string MCP, std::string fWinMin, std::string fWinMax)
{
	float XMax, YMax, Xshift, Yshift, AmpMean, AmpSigma;
	std::string TimeMCP, TimeShift, AmpMean_str, AmpSigma_str; 
	TFile *f = TFile::Open(FileIn.c_str());
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

	//std::string Selection = "(fabs(X[0]-("+std::to_string(XMax)+"))<"+std::to_string(bound)+" || fabs(X[1]-("+std::to_string(XMax)+")-("+std::to_string(Xshift)+"))<"+std::to_string(bound)+") && (fabs(Y[0]-("+std::to_string(YMax)+"))<"+std::to_string(bound)+" || fabs(Y[1]-("+std::to_string(YMax)+")-("+std::to_string(Yshift)+"))<"+std::to_string(bound)+") && amp_max["+MCP+"]>100";
	
	std::string Selection = "fabs(X[0]-("+std::to_string(XMax)+"))<"+std::to_string(bound)+" && fabs(Y[0]-("+std::to_string(YMax)+"))<"+std::to_string(bound)+  " && amp_max["+MCP+"]>100";

	//std::string Selection = "fabs(X[1]-("+std::to_string(XMax)+")-("+std::to_string(Xshift)+"))<"+std::to_string(bound)+" && fabs(Y[1]-("+std::to_string(YMax)+")-("+std::to_string(Yshift)+"))<"+std::to_string(bound)+" && amp_max["+MCP+"]>100";

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;	
	//Selection = "WF_ch == " + detector + " && " + Selection;	
	
	cout << Selection << endl;

	TH1F* tD = new TH1F("tD", "", 500, -5, 10);
	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]>>tD").c_str(), Selection.c_str());
	
	float Xfirst = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())-1;
	float Xlast = tD->GetXaxis()->GetBinCenter(tD->GetMaximumBin())+1;

	TCanvas* c0 = new TCanvas("c0", "c0");
	tD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD->Fit("gaus", "", "", Xfirst, Xlast);
	tD->Draw();
	c0->SaveAs((pathToOutput+"fitTimeDist/fitWindow/Time_"+detector+"-"+MCP+"_"+RunStats+"_fWmin_"+fWinMin+"_fWmax_"+fWimMax+".png").c_str());
	c0->SaveAs((pathToOutput+"fitTimeDist/fitWindow/Time_"+detector+"-"+MCP+"_"+RunStats+"_fWmin_"+fWinMin+"_fWmax_"+fWimMax+".pdf").c_str());
	/*
	TH1F* tD1 = new TH1F("tD1", "", 500, -5, 10);
	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]>>tD1").c_str());
	
	TCanvas* c1 = new TCanvas("c1", "c1");
	tD1->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD1->Fit("gaus", "", "", Xfirst, Xlast);
	tD1->Draw();
	c1->SaveAs((pathToOutput+"fitTimeDist/RawTimeDistribution/TimeNoCut_"+detector+"-"+MCP+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"fitTimeDist/RawTimeDistribution/TimeNoCut_"+detector+"-"+MCP+"_"+RunStats+".pdf").c_str());
	*/
	
}
	

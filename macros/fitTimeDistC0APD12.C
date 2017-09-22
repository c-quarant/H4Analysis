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

void fitTimeDistC0APD12(std::string FileIn, std::string detector1, std::string detector2, Float_t bound)
{
	float XMax, YMax, Xshift, Yshift, AmpMean, AmpSigma;
	std::string TimeMCP1, TimeMCP2, AmpMean_str, AmpSigma_str; 
	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain;

	//Shift between Hodo planes
	Xshift = HodoPlaneShift(h4, detector1, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector1, pathToOutput, RunStats, "Y");
	
	//Selection on distance of hitting position from the center of the detector 1
	AmplitudeProfilesFit(h4, detector1, pathToOutput, RunStats, bound, &XMax, &YMax);

	std::string XMax_str = std::to_string(XMax);
	std::string YMax_str = std::to_string(YMax);
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);


	std::string PosSel1 = "(fabs(X[0]-("+XMax_str+"))<"+bound_str+" || fabs(X[1]-("+XMax_str+")-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]-("+YMax_str+"))<"+bound_str+" || fabs(Y[1]-("+YMax_str+")-("+Yshift_str+"))<"+bound_str+") && (fabs(X[0])<5 || fabs(X[1]-"+Xshift_str+")<5) && (fabs(Y[0])<5 || fabs(Y[1]-"+Yshift_str+")<5)";
	
	//Selection on Amplitude detector 1
	AmplitudeHist(h4, detector1, PosSel1, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	std::string AmpSel1 = "fabs(fit_ampl["+detector1+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;	

	//Selection on distance of hitting position from the center of the detector 2
	AmplitudeProfilesFit(h4, detector2, pathToOutput, RunStats, bound, &XMax, &YMax);

	XMax_str = std::to_string(XMax);
	YMax_str = std::to_string(YMax);
	
	std::string PosSel2 = "(fabs(X[0]-("+XMax_str+"))<"+bound_str+" || fabs(X[1]-("+XMax_str+")-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]-("+YMax_str+"))<"+bound_str+" || fabs(Y[1]-("+YMax_str+")-("+Yshift_str+"))<"+bound_str+") && (fabs(X[0])<5 || fabs(X[1]-"+Xshift_str+")<5) && (fabs(Y[0])<5 || fabs(Y[1]-"+Yshift_str+")<5)";
	
	//Selection on Amplitude detector 2
	AmplitudeHist(h4, detector2, PosSel2, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	std::string AmpSel2 = "fabs(fit_ampl["+detector2+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;	
	
	//plot and fit APD1_APD2 time distribution
	std::string tD_APD1_APD2_Sel = PosSel1 + " && " + AmpSel1 + " && " + PosSel2 + " && " + AmpSel2;

	TH1F* tD_APD1_APD2 = new TH1F("tD_APD1_APD2", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD1_APD2->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD1_APD2->SetBins(1500, -40, 40);
	h4->Draw(("fit_time["+detector1+"]-time["+detector2+"]>>tD_APD1_APD2").c_str(), tD_APD1_APD2_Sel.c_str());
	
	float Xfirst = tD_APD1_APD2->GetXaxis()->GetBinCenter(tD_APD1_APD2->GetMaximumBin())-1;
	float Xlast = tD_APD1_APD2->GetXaxis()->GetBinCenter(tD_APD1_APD2->GetMaximumBin())+1;

	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD1_APD2->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD1_APD2->GetXaxis()->SetTitle(("time["+detector1+"]-time["+detector2+"] (ns)").c_str());
	tD_APD1_APD2->GetYaxis()->SetTitle("events");
	tD_APD1_APD2->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD1_APD2->Draw();
	
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector1+"-"+detector2+"_"+RunStats+".png").c_str());
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector1+"-"+detector2+"_"+RunStats+".pdf").c_str());

}



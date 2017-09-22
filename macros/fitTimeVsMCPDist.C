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

TimeParameters fitTimeVsMCPDist(std::string FileIn, std::string detector, Float_t bound)
{
	TimeParameters FitResults;
	float XCentre, YCentre, Xshift, Yshift, AmpMean, AmpSigma, ANoiseRatio, ANoiseRatioErr;
	std::string TimeMCP1, TimeMCP2, AmpMean_str, AmpSigma_str; 
	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4->GetEntry(0);

	FitResults.Run = (int)h4->GetLeaf("run")->GetValue(0);
	FitResults.Energy = (int)h4->GetLeaf("Energy")->GetValue(0);
	FitResults.Gain = (int)h4->GetLeaf("CHGain")->GetValue(0);

	std::string Gain = std::to_string(FitResults.Gain);
	std::string Energy = std::to_string(FitResults.Energy);
	std::string RunStats = Energy+"Gev_G"+Gain;
	
	//Shift between Hodo planes
	Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");

	//Selection on distance of hitting position from the center of the detector
	std::string AmplMCPSel = "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && amp_max[MCP2]>200 && amp_max[MCP2]<2000";
	AmplitudeProfilesFit(h4, detector, AmplMCPSel, pathToOutput, RunStats, bound, &XCentre, &YCentre, Xshift, Yshift);

	std::string XCentre_str = std::to_string(XCentre);
	std::string YCentre_str = std::to_string(YCentre);
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	//Drawing Amplitude Maps
	AmplitudeMaps(h4, detector, AmplMCPSel, Xshift_str, Yshift_str, XCentre, YCentre, bound, RunStats);

	std::string PosSel = "(fabs(X[0]-("+XCentre_str+"))<"+bound_str+" || fabs(X[1]-("+XCentre_str+")-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]-("+YCentre_str+"))<"+bound_str+" || fabs(Y[1]-("+YCentre_str+")-("+Yshift_str+"))<"+bound_str+") && (fabs(X[0])<5 || fabs(X[1]-"+Xshift_str+")<5) && (fabs(Y[0])<5 || fabs(Y[1]-"+Yshift_str+")<5)";


////////// MCP1 AS TIME REFERENCE/////////////

	//Selection on MCP1 time & amplitude
	TimeMCP1 = std::to_string(MeanTimeMCP(h4, PosSel, pathToOutput+"fitTimeDist/", RunStats, "MCP1"));
	std::string MCP1Sel = "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && fabs(time[MCP1]-("+TimeMCP1+"))<7";

	//Selection on Amplitude
	GaussPar SignalAmp;
	AmplitudeHistPar(h4, detector, PosSel+" && "+MCP1Sel, pathToOutput,  RunStats, &SignalAmp);
	AmpMean_str = std::to_string(SignalAmp.Mean);
	AmpSigma_str = std::to_string(SignalAmp.Sigma);
	std::string AmpSel = "fabs(amp_max["+detector+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;

	/*
	AmplitudeHist(h4, detector, PosSel+" && "+MCP1Sel, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);
	std::string AmpSel = "fabs(amp_max["+detector+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;	
	*/

	std::string tD_APD_MCP1_Sel = PosSel + " && " + MCP1Sel + " && " + AmpSel;
	//cout << tD_APD_MCP1_Sel << endl;

	//TemplateFit Chi2 distribution
	PlotTemplateChi2(h4, detector, "MCP1", tD_APD_MCP1_Sel, RunStats);

	//define APD_MCP1 time distribution
	TH1F* tD_APD_MCP1 = new TH1F("tD_APD_MCP1", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_MCP1->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_MCP1->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+detector+"]-time[MCP1]>>tD_APD_MCP1").c_str(), tD_APD_MCP1_Sel.c_str());
	
	float Xfirst = tD_APD_MCP1->GetXaxis()->GetBinCenter(tD_APD_MCP1->GetMaximumBin())-1;
	float Xlast = tD_APD_MCP1->GetXaxis()->GetBinCenter(tD_APD_MCP1->GetMaximumBin())+1;

	//plot and fit APD_MCP1 time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_MCP1->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP1->GetXaxis()->SetTitle(("time["+detector+"]-time[MCP1] (ns)").c_str());
	tD_APD_MCP1->GetYaxis()->SetTitle("events");
	tD_APD_MCP1->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_MCP1->Draw();

	FitResults.TimeMean[0] = tD_APD_MCP1->GetFunction("gaus")->GetParameter(1);
	FitResults.TimeMeanErr[0] = tD_APD_MCP1->GetFunction("gaus")->GetParError(1);
	FitResults.TimeSigma[0] = tD_APD_MCP1->GetFunction("gaus")->GetParameter(2);
	FitResults.TimeSigmaErr[0] = tD_APD_MCP1->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP1_"+RunStats+".png").c_str());
	c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP1_"+RunStats+".pdf").c_str());

	//Compute Ratio Signal_Amplitude/Noise_RMS
	GaussPar Useless;
	ComputeAvsNoise(h4, detector, bound, "MCP1", Xshift, Yshift, XCentre, YCentre, &Useless, &ANoiseRatio, &ANoiseRatioErr);

	FitResults.ANoiseR[0] = ANoiseRatio;
	FitResults.ANoiseRErr[0] = ANoiseRatioErr;

	//Drawing Time Maps MCP1
	TimeMaps(h4, detector, "MCP1", MCP1Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XCentre, YCentre, bound, RunStats, tD_APD_MCP1->GetFunction("gaus")->GetParameter(1), tD_APD_MCP1->GetFunction("gaus")->GetParameter(2));
	
	
	


////////// MCP2 AS TIME REFERENCE/////////////////////////

	//Selection on MCP2 time && amplitude
	TimeMCP2 = std::to_string(MeanTimeMCP(h4, PosSel, pathToOutput+"fitTimeDist/", RunStats, "MCP2"));
	std::string MCP2Sel = "amp_max[MCP2]>200 && amp_max[MCP2]<2000 && fabs(time[MCP2]-("+TimeMCP2+"))<7";

	std::string tD_APD_MCP2_Sel = PosSel + " && " + MCP2Sel + " && " + AmpSel;
	//cout << tD_APD_MCP2_Sel << endl;
	
	//Drawing Time Maps MCP1
	TimeMaps(h4, detector, "MCP2", MCP2Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XCentre, YCentre, bound, RunStats,0,0);

	//TemplateFit Chi2 distribution
	PlotTemplateChi2(h4, detector, "MCP1", tD_APD_MCP2_Sel, RunStats);
	
	
	//Define APD_MCP2 time distribution
	TH1F* tD_APD_MCP2 = new TH1F("tD_APD_MCP2", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_MCP2->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_MCP2->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+detector+"]-time[MCP2]>>tD_APD_MCP2").c_str(), tD_APD_MCP2_Sel.c_str());

	Xfirst = tD_APD_MCP2->GetXaxis()->GetBinCenter(tD_APD_MCP2->GetMaximumBin())-1;
	Xlast = tD_APD_MCP2->GetXaxis()->GetBinCenter(tD_APD_MCP2->GetMaximumBin())+1;
	
	//plot and fit APD_MCP2 time distribution
	gStyle->SetOptStat();
	TCanvas* c2 = new TCanvas("c2", "c2");
	tD_APD_MCP2->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP2->GetXaxis()->SetTitle(("time["+detector+"]-time[MCP2] (ns)").c_str());
	tD_APD_MCP2->GetYaxis()->SetTitle("events");
	tD_APD_MCP2->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_MCP2->Draw();

	FitResults.TimeMean[1] = tD_APD_MCP2->GetFunction("gaus")->GetParameter(1);
	FitResults.TimeMeanErr[1] = tD_APD_MCP2->GetFunction("gaus")->GetParError(1);
	FitResults.TimeSigma[1] = tD_APD_MCP2->GetFunction("gaus")->GetParameter(2);
	FitResults.TimeSigmaErr[1] = tD_APD_MCP2->GetFunction("gaus")->GetParError(2);

	c2->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP2_"+RunStats+".png").c_str());
	c2->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP2_"+RunStats+".pdf").c_str());

	//Compute Ratio Signal_Amplitude/Noise_RMS
	ComputeAvsNoise(h4, detector, bound, "MCP2", Xshift, Yshift, XCentre, YCentre, &Useless, &ANoiseRatio, &ANoiseRatioErr);

	FitResults.ANoiseR[1] = ANoiseRatio;
	FitResults.ANoiseRErr[1] = ANoiseRatioErr;	

	//Drawing Time Maps MCP2
	TimeMaps(h4, detector, "MCP2", MCP1Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XCentre, YCentre, bound, RunStats, tD_APD_MCP2->GetFunction("gaus")->GetParameter(1), tD_APD_MCP2->GetFunction("gaus")->GetParameter(2));





////////// MCP MEAN AS TIME REFERENCE ///////////////////

	//Selection on both MCP time & amplitude
	std::string tD_APD_MCP_Mean_Sel = PosSel + " && " + MCP1Sel + " && " + MCP2Sel + " && " + AmpSel;
	//cout << tD_APD_MCP_Mean_Sel
	
	//Define APD_MCP_Mean time distribution
	TH1F* tD_APD_MCP_Mean = new TH1F("tD_APD_MCP_Mean", "", 4000, -40, 40);
	if(Energy == "20" && Gain!="200") tD_APD_MCP_Mean->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_MCP_Mean->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+detector+"]-0.5*(time[MCP1]+time[MCP2])>>tD_APD_MCP_Mean").c_str(), tD_APD_MCP_Mean_Sel.c_str());
	
	Xfirst = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())-1;
	Xlast = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())+1;

	//plot and fit APD_MCP_Mean time distribution
	gStyle->SetOptStat();
	TCanvas* c3 = new TCanvas("c3", "c3");
	tD_APD_MCP_Mean->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP_Mean->GetXaxis()->SetTitle(("time["+detector+"]-timeMean (ns)").c_str());
	tD_APD_MCP_Mean->GetYaxis()->SetTitle("events");
	tD_APD_MCP_Mean->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_MCP_Mean->Draw();

	FitResults.TimeMean[2] = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(1);
	FitResults.TimeMeanErr[2] = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(1);
	FitResults.TimeSigma[2] = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(2);
	FitResults.TimeSigmaErr[2] = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(2);
	FitResults.ANoiseR[2] = ANoiseRatio;
	FitResults.ANoiseRErr[2] = ANoiseRatioErr;
	
	c3->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP_Mean_"+RunStats+".png").c_str());
	c3->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP_Mean_"+RunStats+".pdf").c_str());

	return FitResults;

}



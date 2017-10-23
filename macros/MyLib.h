#ifndef __MYLIB_H__
#define __MYLIB_H__

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
#include <iostream>
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "TGraphErrors.h"
#include <string>
#include <fstream>
#include <math.h>
#include <vector>

class GaussPar{
	public:
	float Mean, Sigma, MeanErr, SigmaErr;	
};

class TimeParameters{
	public:
	int Run;
	int Energy;
	int Gain;

	float TimeMean[3];
	float TimeMeanErr[3];
	float TimeSigma[3];
	float TimeSigmaErr[3];

	float ANoiseR[3];
	float ANoiseRErr[3];
};

//returns the mean time shift between detector and MCP selected
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOutput, std::string RunStats, std::string MCP)
{
	TH1F* raw_time_dist = new TH1F("raw_time_dist", "", 200, -3, 3);

	h4->Draw((std::string("time["+detector+"]-time["+MCP+"]>>raw_time_dist")).c_str(), Selection.c_str());

	raw_time_dist->GetXaxis()->SetTitle((std::string("time["+detector+"]-time["+MCP+"] (ns)")).c_str());
	raw_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* c0 = new TCanvas();
    	c0->cd();
	raw_time_dist->Fit("gaus", "", "", -3, 3);
    	raw_time_dist->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".pdf").c_str());

	cout << raw_time_dist->GetMean() << endl;
	cout << raw_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return raw_time_dist->GetFunction("gaus")->GetMaximumX();
}

//returns the mean time measured by MCP selected
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOutput, std::string RunStats, std::string MCP)
{
	float MCPTimeMean;

	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);

	h4->Draw((std::string("time["+MCP+"]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle(("time["+MCP+"] (ns)").c_str());
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas("ca", "ca");
    	ca->cd();
	MCP_time_dist->Fit("gaus", "Q", "", 0, 50);
    	MCP_time_dist->Draw();

	MCPTimeMean = MCP_time_dist->GetFunction("gaus")->GetMaximumX();

    	ca -> SaveAs(std::string(pathToOutput+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string(pathToOutput+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".pdf").c_str());

	ca->~TCanvas();
	MCP_time_dist->~TH1F();	

	return MCPTimeMean;
}

//Draw Amplitude histogram and fit it with gaus
void AmplitudeHist2(TTree *h4, std::string detector, std::string Selection, std::string pathToOutput,  std::string RunStats, float* AmpMean, float* AmpSigma, std::string HAmpName)
{
	TH1F* HAmp = new TH1F(HAmpName.c_str(), "", 2500, -50, 5000);

	h4->Draw(("amp_max["+detector+"]>>"+HAmpName+"").c_str(), Selection.c_str());
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.80), (int)(HAmp->GetMaximumBin()*1.10));


	TCanvas* c1 = new TCanvas("c1", "c1");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->GetXaxis()->SetTitle("Amplitude Max (ADC counts)");
	HAmp->GetYaxis()->SetTitle("events");
	HAmp->Draw();

	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());
	/*
	HAmp->Fit("gaus");
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());
	
	*AmpMean = HAmp->GetFunction("gaus")->GetParameter(1);
	*AmpSigma = HAmp->GetFunction("gaus")->GetParameter(2);
	*/
}

//Draw Amplitude histogram and fit it with gaus
void AmplitudeHist(TTree *h4, std::string detector, std::string Selection, std::string pathToOutput,  std::string RunStats, float* AmpMean, float* AmpSigma, std::string AmpType)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2000, +50, 5000);

	h4->Draw((AmpType+"["+detector+"]>>HAmp").c_str(), Selection.c_str());
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.80), (int)(HAmp->GetMaximumBin()*1.10));


	TCanvas* c0 = new TCanvas("c0", "c0");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus", "Q");	
	HAmp->Draw();
	
	c0->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c0->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());
	
	*AmpMean = HAmp->GetFunction("gaus")->GetParameter(1);
	*AmpSigma = HAmp->GetFunction("gaus")->GetParameter(2);

	//c0->~TCanvas();
	
}

//Draw Amplitude histogram and fit it with gaus
void AmplitudeHistPar(TTree *h4, std::string detector, std::string Selection, std::string pathToOutput,  std::string RunStats, GaussPar* gPar)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2000, -50, 5000);
	float tmpMean, tmpSigma;
	float Xfirst, Xlast, XMax;

	h4->Draw(("fit_ampl["+detector+"]>>HAmp").c_str(), Selection.c_str());
	HAmp->GetXaxis()->SetRangeUser(50, 4000);

	XMax = HAmp->GetBinCenter(HAmp->GetMaximumBin());
	Xfirst = XMax*0.6;
	Xlast = XMax*1.2;

	HAmp->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]*TMath::Exp(-(x-[4])*(x-[4])/(2*[5]*[5]))", Xfirst, Xlast);

	if(Selection!=""){
		fitFunc->SetParLimits(0, HAmp->GetMaximum()*0.8, HAmp->GetMaximum()*1.39);
		fitFunc->SetParLimits(1, XMax*0.95, XMax*1.05);
		fitFunc->SetParLimits(2, 0, 150);
		fitFunc->SetParLimits(3, 0, HAmp->GetMaximum()*0.7);
		fitFunc->SetParLimits(4, XMax*0.85, XMax*0.95);
		fitFunc->SetParLimits(5, 5, 500);

		fitFunc->SetParameter(0, HAmp->GetMaximum());
		fitFunc->SetParameter(1, XMax);
		fitFunc->SetParameter(2, 40);
		fitFunc->SetParameter(3, HAmp->GetMaximum()/10);
		fitFunc->SetParameter(4, XMax*0.9);
		fitFunc->SetParameter(5, 80);	

		fitFunc->SetParNames("A0", "Mean0", "Sigma0", "A1", "Mean1", "Sigma1");

		HAmp->Fit("fitFunc", "Q");	
	}
	else if(Selection == "")
	{
		Xfirst=-30;
		Xlast = 30;
		HAmp->GetXaxis()->SetRangeUser(Xfirst, Xlast);
		XMax = HAmp->GetBinCenter(HAmp->GetMaximumBin());
		
		fitFunc->SetParLimits(0, HAmp->GetMaximum()*0.75, HAmp->GetMaximum()*1.2);
		fitFunc->SetParLimits(1, XMax*0.95-1, XMax*1.05+1);
		fitFunc->SetParLimits(2, 0, 100);
		fitFunc->SetParLimits(3, 0, HAmp->GetMaximum()*0.45);
		fitFunc->SetParLimits(4, XMax*0.8-5, XMax);
		fitFunc->SetParLimits(5, 0, 200);

		fitFunc->SetParameter(0, HAmp->GetMaximum());
		fitFunc->SetParameter(1, XMax);
		fitFunc->SetParameter(2, 4);
		fitFunc->FixParameter(3, 0);
		fitFunc->FixParameter(4, 0);
		fitFunc->FixParameter(5, 0);

		fitFunc->SetParNames("A0", "Mean0", "Sigma0", "A1", "Mean1", "Sigma1");
		HAmp->Fit("fitFunc", "RQ", "", XMax-20, XMax+20);		
	}

	TCanvas* c1 = new TCanvas("c1", "c1");

	tmpMean = HAmp->GetFunction("fitFunc")->GetParameter(1);
	tmpSigma = HAmp->GetFunction("fitFunc")->GetParameter(2);
	
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->GetXaxis()->SetTitle("Amplitude (ADC counts)");
	HAmp->GetYaxis()->SetTitle("events");
	HAmp->Draw();
	
	//draw components
	TF1* gaus0 = new TF1("gaus0", "gaus", Xfirst+0.6, Xlast-0.6);
	TF1* gaus1 = new TF1("gaus1", "gaus", Xfirst+0.6, Xlast-0.6);

	gaus0->SetParameter(0, fitFunc->GetParameter("A0"));
	gaus0->SetParameter(1, fitFunc->GetParameter("Mean0"));
	gaus0->SetParameter(2, fitFunc->GetParameter("Sigma0"));

	gaus1->SetParameter(0, fitFunc->GetParameter("A1"));
	gaus1->SetParameter(1, fitFunc->GetParameter("Mean1"));
	gaus1->SetParameter(2, fitFunc->GetParameter("Sigma1"));

	gaus0->SetLineColor(kBlue);
	gaus0->Draw("SAME");

	gaus1->SetLineColor(kGreen);
	gaus1->Draw("SAME");

	gPar->Mean = HAmp->GetFunction("fitFunc")->GetParameter(1);
	gPar->MeanErr = HAmp->GetFunction("fitFunc")->GetParError(1);
	gPar->Sigma = HAmp->GetFunction("fitFunc")->GetParameter(2);
	gPar->SigmaErr = HAmp->GetFunction("fitFunc")->GetParError(2);

	c1->SaveAs((pathToOutput+"Amp_plot/FitAmp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"Amp_plot/FitAmp_"+detector+"_"+RunStats+".pdf").c_str());

	c1->~TCanvas();
	
}

//Fit amplitude histogram with gaus and returns mean value
float AmplitudeMean(TTree *h4, std::string detector, std::string Selection)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -200, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus");	
	
	return HAmp->GetFunction("gaus")->GetParameter(1);
}

//Fit amplitude histogram with gaus and returns mean value's error
float AmplitudeMeanErr(TTree *h4, std::string detector, std::string Selection)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -200, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus");	
	
	return HAmp->GetFunction("gaus")->GetParError(1);
}

//Fit amplitude histogram with gaus and returns sigma
float AmplitudeSigma(TTree *h4, std::string detector, std::string Selection)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -200, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus");	
	
	return HAmp->GetFunction("gaus")->GetParameter(2);
}

//Draw Amplitude histogram and fit
void DrawAmplitudeHist(TTree *h4, std::string detector, std::string Selection, std::string pathToOutput,  std::string RunStats)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -50, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	TCanvas* c1 = new TCanvas("c1", "c1");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus");	
	HAmp->Draw();
	
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());

	c1->~TCanvas();	
}

//Fit amplitude histogram with gaus and returns sigma's error
float AmplitudeSigmaErr(TTree *h4, std::string detector, std::string Selection)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -200, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus");	
	
	return HAmp->GetFunction("gaus")->GetParError(2);
}


//returns the shift of plane 1 from plane 0 of hodoscope along selected axis
float HodoPlaneShift(TTree* h4, std::string detector, std::string pathToOutput, std::string RunStats, std::string axis)
{
	auto *DXvsX = new TProfile(("D"+axis+"vs"+axis+"").c_str(), "", 128, -16, 16, -10, 10);

	//Filling DeltaX histogram
	h4->Draw(("("+axis+"[1]-"+axis+"[0]):"+axis+"[0]>>D"+axis+"vs"+axis+"").c_str(), (axis+"[0]>-800 && "+axis+"[1]>-800").c_str());

	//Fitting and Drawing DeltaX
	TCanvas *cDeltaX = new TCanvas(("cDelta"+axis+"").c_str(), ("cDelta"+axis+"").c_str());
	TH2F* Hset = new TH2F(("Hset"+axis).c_str(),"", 128, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	Hset->GetYaxis()->SetTitle((axis+"[1]-"+axis+"[0]").c_str());
  
	Hset->Draw();	
	DXvsX->Fit("pol1", "Q", "", -4, 4);
	DXvsX->Draw("SAME");
	
	std::string fileOutpdf = pathToOutput + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = pathToOutput + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".png";
    
  	cDeltaX->SaveAs(fileOutpdf.c_str(), "Q");
  	cDeltaX->SaveAs(fileOutpng.c_str(), "Q");

	cDeltaX->~TCanvas();

	return DXvsX->GetFunction("pol1")->GetParameter(0);
}

//Draw amplitude profiles in X and Y calculate maximum X
void AmplitudeProfilesFit(TTree *h4, std::string detector, std::string AmpMCPSel, std::string pathToOutput, std::string RunStats, Float_t bound, float* XMax, float* YMax, float Xshift, float Yshift, std::string AmplType)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, hodoCfg;
	Float_t X1true, Y1true, Xavg, Yavg;
	Float_t Gain_val, Energy_val;
	float fitRange=5;
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string runNum = std::to_string((int)h4->GetLeaf("run")->GetValue(0));

	gStyle->SetOptStat();
	gStyle->SetOptFit();

	auto *AmpXavg = new TProfile("AmpXavg", "", 128, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 128, -16, 16, 0, 10000);
	
	//Filling Amplitude profile histograms Draw method
	std::string varexp, selection;	
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	cout << "\n\nXshift = " << Xshift << "    Yshift = " << Yshift << endl << endl;

	varexp = AmplType+"[" + detector + "]:0.5*(X[0]+X[1]-(" + Xshift_str + "))>>AmpXavg";
	selection = AmpMCPSel + " && 0.5*(Y[0]+Y[1]-(" + Yshift_str + "))>-" + bound_str + " && 0.5*(Y[0]+Y[1]-(" + Yshift_str + "))<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = AmplType+"[" + detector + "]:0.5*(Y[0]+Y[1]-(" + Yshift_str + "))>>AmpYavg";
	selection = AmpMCPSel + " && 0.5*(X[0]+X[1]-(" + Xshift_str + "))>-" + bound_str + " && 0.5*(X[0]+X[1]-(" + Xshift_str + "))<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	TH2F* H1 = new TH2F("H1","", 128, -16, 16, 50, 0, (AmpXavg->GetMaximum())*1.25);
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpXavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	std::string fileOutpdf = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpXAVG_" + runNum + "_" + detector + "_G" + Gain + ".pdf";
  	std::string fileOutpng = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpXAVG_" + runNum + "_" + detector + "_G" + Gain + ".png";
  
	cout << "\n\nX Centre Position = " << fitResX->GetMaximumX() << "\t";  

  	c5->SaveAs(fileOutpdf.c_str(), "Q");
  	c5->SaveAs(fileOutpng.c_str(), "Q");
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpYAVG_" + runNum + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpYAVG_" + runNum + "_" + detector + "_G" + Gain + ".png";

	cout << "Y Centre Position = " << fitResY->GetMaximumX() << "\n" << endl;      

  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());
	/*
	if(fitResX->GetMaximumX()-bound<-5) *XMax = bound-5;
	else if(fitResX->GetMaximumX()+bound>5) *XMax = 5-bound;
	else *XMax = fitResX->GetMaximumX();
	
	if(fitResY->GetMaximumX()-bound<-5) *YMax = bound-5;
	else if(fitResY->GetMaximumX()+bound>5) *YMax = 5-bound;
	else *YMax = fitResY->GetMaximumX();
	*/
	
	*XMax = fitResX->GetMaximumX();
	*YMax = fitResY->GetMaximumX();
	
	H1->~TH2F();
	c5->~TCanvas();
	c6->~TCanvas();
}

void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP)
{
	int i;
	float AmpMean, AmpSigma;

	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpMean_str;
	std::string AmpSigma_str;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max["+MCP+"]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma, "amp_max");
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = "WF_ch == " + detector + " && " + Selection;	
	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *waveForm = new TProfile("waveForm", "", 300, -10, 40);
	waveForm = (TProfile*)aSlices[1];

	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
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
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+detector+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.pdf").c_str(),"pdf");
	
	TCanvas* c0 = new TCanvas();
    	c0->cd();
	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName("Waveform_");
    	waveForm->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.pdf").c_str());    	

    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms_");
	h2_amp_vs_time->SetName("H2Waveforms_");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close(); 
}

void WFPulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP)
{
	int i;
	float AmpMean, AmpSigma;

	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpMean_str;
	std::string AmpSigma_str;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 2048, -204.8, 204.8 , 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 2048, -204.8, 204.8, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max["+MCP+"]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma, "amp_max");
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = "WF_ch == " + detector + " && " + Selection;	
	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *waveForm = new TProfile(("XTAL_"+detector+"_"+RunStats+"_prof_").c_str(), "", 2048, -204.8, 204.8);
	waveForm = (TProfile*)aSlices[1];
	
	int NBins=2048, C=0;
	for(i=1; i<NBins; i++)
	{
		if((waveForm->GetBinCenter(i)>0 || waveForm->GetBinCenter(i)<50) && waveForm->GetBinError(i)>0.3)
		{
			waveForm->SetBinError(i, 0.01);
			C++;
		}
		if((waveForm->GetBinCenter(i)<-25 || waveForm->GetBinCenter(i)>165) && (waveForm->GetBinContent(i)>0.08 || waveForm->GetBinError(i)>0.1))
		{
			waveForm->SetBinContent(i,0);
			waveForm->SetBinError(i, 0.01);
			C++;
		}
	}
	cout << "CORRECTIONS:   " << C << endl;

	p2D_amp_vs_time->GetXaxis()->SetTitle("WF_time (ns)");
    	h2_amp_vs_time->GetXaxis()->SetTitle("WF_time (ns)");
    	waveForm->GetXaxis()->SetTitle("WF_time");
    	p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");	
    	    	
	gStyle->SetOptStat(0);
	
    	TCanvas* c1 = new TCanvas();
    	c1->cd();
    	h2_amp_vs_time->Draw("COLZ");
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+runNum+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+detector+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.pdf").c_str(),"pdf");
	
	TCanvas* c0 = new TCanvas();
    	c0->cd();
	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName(("XTAL_"+detector+"_"+RunStats+"_prof_").c_str());
    	waveForm->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+runNum+"_profile.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/WFPulseShapes/WFPS_"+detector+"_"+RunStats+"_"+runNum+"_profile.pdf").c_str());    	

    	TFile* output_Waveform = new TFile(std::string("/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/templateProf/WF_"+detector+"_"+RunStats+".root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms_");
	h2_amp_vs_time->SetName("H2Waveforms_");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close();
	cout << waveForm->GetName() << endl; 
}

void PulseShapesMCP(TTree* h4, std::string detector, int plane, std::string MCP2, float XMax, float YMax, float bound, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP)
{
	int i;
	float AmpMean, AmpSigma;

	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpMean_str;
	std::string AmpSigma_str;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(bound) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(bound) + " && amp_max["+MCP+"]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma, "amp_max");
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = "WF_ch == " + detector + " && " + Selection;	
	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *waveForm = new TProfile("waveForm", "", 300, -10, 40);
	waveForm = (TProfile*)aSlices[1];

	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
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
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+detector+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.pdf").c_str(),"pdf");
	
	TCanvas* c0 = new TCanvas();
    	c0->cd();
	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName("Waveform_");
    	waveForm->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.pdf").c_str());    	

    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms_");
	h2_amp_vs_time->SetName("H2Waveforms_");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close(); 
}

float NumericalRootFinder(TF1* f, float Xfirst, float Xlast)
{
	// Create the function and wrap it
   	ROOT::Math::WrappedTF1 wf1(*f);
 
   	// Create the Integrator
 	ROOT::Math::Roots::Bisection brf;
 
	// Set parameters of the method
   	brf.SetFunction( wf1, Xfirst, Xlast );
   	brf.Solve();
 
	return brf.Root();
}


void TimeMaps(TTree* h4, std::string detector, std::string MCP, std::string Selection, std::string Xshift_str, std::string Yshift_str, float XCenter, float YCenter, float bound, std::string RunStats, float TimeMean, float TimeSigma)
{	
	gStyle->SetOptStat(0);
  	int i, Nentries, CH, runNum; 
  	float AmpTemp;
  	
  	h4->GetEntry(0);
 	CH=h4->GetLeaf(detector.c_str())->GetValue(0);

	
	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

		
  	//2DHist definition 
  	auto *TimeXY0 = new TProfile2D("TimeXY0","", 32, -16, 16, 32, -16, 16, 4, 15); 
  	auto *TimeXY1 = new TProfile2D("TimeXY1","", 32, -16, 16, 32, -16, 16, 4, 15); 
  	auto *TimeXYM = new TProfile2D("TimeXYM","", 32, -16, 16, 32, -16, 16, 0, 20); 
  	
  	Nentries = h4->GetEntries();
  	

   	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]:Y[0]:X[0]>>TimeXY0").c_str(), (Selection + " && X[0]>-16 && Y[0]>-16").c_str());
  	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>TimeXY1").c_str(), (Selection + " && X[1]>-16 && Y[1]>-16").c_str());
  	h4->Draw(("fit_time["+detector+"]-time["+MCP+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>TimeXYM").c_str(), (Selection + " && X[0]>-16 && Y[0]>-16 && X[1]>-16 && Y[1]>-16").c_str());

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	//FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	TimeXY0->GetZaxis()->SetTitle(("fit_time["+detector+"] (ADC counts)").c_str());
  	
  	H1->Draw();
	TimeXY0->GetZaxis()->SetRangeUser(TimeMean-3*TimeSigma, TimeMean+3*TimeSigma);	
  	TimeXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(XCenter-bound, YCenter-bound, XCenter-bound, YCenter+bound);
	line_left->SetLineColor(kRed);
	line_left->Draw();

	TLine *line_right = new TLine(XCenter+bound, YCenter-bound, XCenter+bound, YCenter+bound);
	line_right->SetLineColor(kRed);
	line_right->Draw();

	TLine *line_up = new TLine(XCenter-bound, YCenter+bound, XCenter+bound, YCenter+bound);
	line_up->SetLineColor(kRed);
	line_up->Draw();

	TLine *line_down = new TLine(XCenter-bound, YCenter-bound, XCenter+bound, YCenter-bound);
	line_down->SetLineColor(kRed);
	line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane0/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane0/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".png";
  	
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  
  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();
	TimeXY1->GetZaxis()->SetRangeUser(TimeMean-3*TimeSigma, TimeMean+3*TimeSigma);	
  	TimeXY1->Draw("COLZ SAME");

	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane1/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane1/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".png";
  	
	
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	
  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");

	c3->cd();
  	
  	H1->Draw();
	TimeXYM->GetZaxis()->SetRangeUser(TimeMean-3*TimeSigma, TimeMean+3*TimeSigma);	
  	TimeXYM->Draw("COLZ SAME");
  	
	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	
	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/PlaneAVG/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/PlaneAVG/TimeXY_" + detector + "-" + MCP + "_" + RunStats + ".png";
  	
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
	TimeXY0->~TProfile2D();
	TimeXY1->~TProfile2D();
	TimeXYM->~TProfile2D();
	H1->~TH2F(); 
}

void AmplitudeMaps(TTree* h4, std::string detector, std::string Selection, std::string Xshift_str, std::string Yshift_str, float XCenter, float YCenter, float bound, std::string RunStats, std::string AmpType)
{	
	gStyle->SetOptStat(0);
  	int i, Nentries, CH, runNum; 
  	float AmpTemp;
    	
  	h4->GetEntry(0);
 	CH=h4->GetLeaf(detector.c_str())->GetValue(0);

	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	
  	Nentries = h4->GetEntries();
  	

   	h4->Draw((AmpType+"["+detector+"]:Y[0]:X[0]>>AmpXY0").c_str(), (Selection + " && X[0]>-800 && Y[0]>-800").c_str());
  	h4->Draw((AmpType+"["+detector+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>AmpXY1").c_str(), (Selection + " && X[1]>-800 && Y[1]>-800").c_str());
  	h4->Draw((AmpType+"["+detector+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>AmpXYM").c_str(), (Selection + " && X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800").c_str());

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	//FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+detector+"] (ADC counts)").c_str());
  	
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(XCenter-bound, YCenter-bound, XCenter-bound, YCenter+bound);
	line_left->SetLineColor(kRed);
	line_left->Draw();

	TLine *line_right = new TLine(XCenter+bound, YCenter-bound, XCenter+bound, YCenter+bound);
	line_right->SetLineColor(kRed);
	line_right->Draw();

	TLine *line_up = new TLine(XCenter-bound, YCenter+bound, XCenter+bound, YCenter+bound);
	line_up->SetLineColor(kRed);
	line_up->Draw();

	TLine *line_down = new TLine(XCenter-bound, YCenter-bound, XCenter+bound, YCenter-bound);
	line_down->SetLineColor(kRed);
	line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + detector + "_" + RunStats + ".png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  

  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");

	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + detector + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + detector + "_" + RunStats + ".png";

  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	

  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + detector + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + detector + "_" + RunStats + ".png";  
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	H1->~TH2F();
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
}

void PlotTemplateChi2(TTree* h4, std::string detector, std::string MCP, std::string Selection, std::string RunStats)
{
	TH1F* tD_Chi2 = new TH1F("tD_Chi2", "", 50, 0, 8);
	h4->Draw(("fit_chi2["+detector+"]>>tD_Chi2").c_str(), Selection.c_str());
	
	TCanvas *c_Chi2 = new TCanvas("c_Chi2", "c_Chi2");
	tD_Chi2->GetXaxis()->SetTitle("#chi^2");
	tD_Chi2->GetYaxis()->SetTitle("events");
	tD_Chi2->Draw();
	c_Chi2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TemplateChi2/Chi2_"+detector+"-"+MCP+"_"+RunStats+".png").c_str());
	c_Chi2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TemplateChi2/Chi2_"+detector+"-"+MCP+"_"+RunStats+".pdf").c_str());

	tD_Chi2->~TH1F();
	c_Chi2->~TCanvas();	
}

int ComputeAvsNoise(TTree* h4, std::string detector, int bound, std::string MCP, float Xshift, float Yshift, float XCentre, float YCentre, GaussPar* AeNoise, float* ANoiseR, float* ANoiseRErr)
{
	GaussPar AmpPar, NoisePar;
	float XMax, YMax;
	std::string TimeMCP, TimeShift, NoiseNtuple; 
		
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain;

	if(detector == "C3")
	{
		if(Gain == "50") NoiseNtuple = "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/ntuples/H42016_Pedestal_5896.root";
		else if(Gain == "100") NoiseNtuple = "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/ntuples/H42016_Pedestal_5894.root";
	}
	else
	{ 
		if(Gain == "50") NoiseNtuple = "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/ntuples/H42016_C0_Pedestal_5902.root";
		else if(Gain == "100") NoiseNtuple = "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/ntuples/H42016_C0_Pedestal_5900.root";
		else if(Gain == "200") NoiseNtuple = "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/ntuples/H42016_C0_Pedestal_5898.root";
	}

	TFile *fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree *h4Noise = (TTree*)fNoise->Get("h4");

	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4Noise->GetEntry(0);
	std::string NoiseGain = std::to_string((int)h4Noise->GetLeaf("CHGain")->GetValue(0));

	if(NoiseGain!=Gain)
	{
		cout << endl << "!!!!!!!!!!!!!!!!!!!! Noise Gain =/= Run Gain !!!!!!!!!!!!!!!!!!!!!!11111111" << endl << endl;
		return -1;
	}

	std::string Selection = "(fabs(X[0]-("+std::to_string(XCentre)+"))<"+std::to_string(bound)+" || fabs(X[1]-("+std::to_string(XCentre)+")-("+std::to_string(Xshift)+"))<"+std::to_string(bound)+") && (fabs(Y[0]-("+std::to_string(YCentre)+"))<"+std::to_string(bound)+" || fabs(Y[1]-("+std::to_string(YCentre)+")-("+std::to_string(Yshift)+"))<"+std::to_string(bound)+") && amp_max["+MCP+"]>100";
	
	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHistPar(h4, detector, Selection, pathToOutput, RunStats, &AmpPar);
	Selection = "";
	AmplitudeHistPar(h4Noise, detector, Selection, pathToOutput, "Noise_"+Gain, &NoisePar);
	/*
	cout << "Signal Amplitude = " << AmpPar.Mean << " +- " << AmpPar.MeanErr << endl;
	cout << "Noise RMS = " << NoisePar.Sigma << " +- " << NoisePar.SigmaErr << endl;
	cout << "A/Noise = " << AmpPar.Mean/NoisePar.Sigma << " +- " << TMath::Sqrt(pow(AmpPar.MeanErr/AmpPar.Mean, 2) + pow(NoisePar.SigmaErr/NoisePar.Sigma, 2))*AmpPar.Mean/NoisePar.Sigma << endl;
	*/
	AeNoise->Mean = AmpPar.Mean;
	AeNoise->MeanErr = AmpPar.MeanErr;
	AeNoise->Sigma = NoisePar.Sigma;
	AeNoise->SigmaErr = NoisePar.SigmaErr;

	*ANoiseR = AmpPar.Mean/NoisePar.Sigma;
	*ANoiseRErr = TMath::Sqrt(pow(AmpPar.MeanErr/AmpPar.Mean, 2) + pow(NoisePar.SigmaErr/NoisePar.Sigma, 2))*AmpPar.Mean/NoisePar.Sigma;

	return 0;
}

void FitRes(TTree* TimeRes, int SelDetector, int SelMCP){
	gStyle->SetOptFit();	

	int i=0, Nentries=TimeRes->GetEntries();
	int EntryGain, EntryDetector;
	float ANoiseRatio[3], ANoiseRatioErr[3], Res[3], ResErr[3];
	std::vector<Float_t> X0, X0Err, Y0, Y0Err, X1, X1Err, Y1, Y1Err, X2, X2Err, Y2, Y2Err;

	TimeRes->SetBranchAddress("Detector", &EntryDetector);
	TimeRes->SetBranchAddress("Gain", &EntryGain);
	
	if(SelDetector==3)
	{
		TimeRes->SetBranchAddress("C0APDs_time_sigma", Res);
		TimeRes->SetBranchAddress("C0APDs_time_sigma_error", ResErr);
	}
	else
	{
		TimeRes->SetBranchAddress("time_sigma", Res);
		TimeRes->SetBranchAddress("time_sigma_error", ResErr);
	}
		
	TimeRes->SetBranchAddress("Amplitude_Noise_Ratio", ANoiseRatio);
	TimeRes->SetBranchAddress("Amplitude_Noise_Ratio_error", ANoiseRatioErr);

	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + [1]*[1] )", 15, 1000);

	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, 10, 100);
	//fitFunc->SetParLimits(2, 0, 2000);  

	fitFunc->SetParameter(0, 5000);
	fitFunc->SetParameter(1, 33);
	//fitFunc->SetParameter(2, 0);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");
	//fitFunc->SetParName(2, "~stochastic");
	


	for(i=0; i<Nentries; i++)
	{
	TimeRes->GetEntry(i);
		if(EntryDetector == SelDetector && EntryGain == 50)
		{
			X0.push_back(ANoiseRatio[SelMCP]);
			X0Err.push_back(ANoiseRatioErr[SelMCP]);
			Y0.push_back(Res[SelMCP]*1000);
			Y0Err.push_back(ResErr[SelMCP]*1000);
		}
		if(EntryDetector == SelDetector && EntryGain == 100)
		{
			X1.push_back(ANoiseRatio[SelMCP]);
			X1Err.push_back(ANoiseRatioErr[SelMCP]);
			Y1.push_back(Res[SelMCP]*1000);
			Y1Err.push_back(ResErr[SelMCP]*1000);
		}
		if(EntryDetector == SelDetector && EntryGain == 200)
		{
			X2.push_back(ANoiseRatio[SelMCP]);
			X2Err.push_back(ANoiseRatioErr[SelMCP]);
			Y2.push_back(Res[SelMCP]*1000);
			Y2Err.push_back(ResErr[SelMCP]*1000);
		}
	}

	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->cd();

	TGraphErrors* g50 = new TGraphErrors(X0.size(), &X0[0], &Y0[0], &X0Err[0], &Y0Err[0]);
	g50->SetTitle("");

	g50->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g50->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	g50->SetMarkerStyle(kFullCircle);
	g50->SetMarkerSize(1);
	fitFunc->SetLineColor(kBlack);	
	g50->Fit("fitFunc");
	g50->Draw("AP");
	
	std::string MCP_str;
	if(SelMCP==0) MCP_str="MCP1";
	else if(SelMCP==1) MCP_str="MCP2";
	else MCP_str="MCP_Mean";

	std::string Detector_str;
	if(SelDetector==0) Detector_str="C3";
	else if(SelDetector==1) Detector_str="C0APD1";
	else if(SelDetector==2) Detector_str="C0APD2";
	else{ Detector_str="C0APD1";  MCP_str="C0APD2"; }

	std::string Info = Detector_str + "-" + MCP_str + "_G50";

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());

	


	TCanvas* c1 = new TCanvas("c1", "c1");
	c1->cd();

	TGraphErrors* g100 = new TGraphErrors(X1.size(), &X1[0], &Y1[0], &X1Err[0], &Y1Err[0]);
	g100->SetTitle("");
		
	g100->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g100->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	g100->SetLineColor(kBlue);
	g100->SetMarkerStyle(kFullSquare);
	g100->SetMarkerSize(1);
	g100->SetMarkerColor(kBlue);
	fitFunc->SetLineColor(kBlue);
	g100->Fit("fitFunc");
       	g100->Draw("AP");

	Info = Detector_str + "-" + MCP_str + "_G" + std::to_string(100);

	c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
	c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());




	TCanvas* c2 = new TCanvas("c2", "c2");
	c2->cd();

	TGraphErrors* g200 = new TGraphErrors(X2.size(), &X2[0], &Y2[0], &X2Err[0], &Y2Err[0]);
	g200->SetTitle("");	

	g200->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g200->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	g200->SetLineColor(kViolet);
	g200->SetMarkerStyle(kFullTriangleDown);
	g200->SetMarkerSize(1);
	g200->SetMarkerColor(6);
	fitFunc->SetLineColor(kViolet);	
	g200->Fit("fitFunc");

        g200->Draw("AP");

	Info = Detector_str + "-" + MCP_str + "_G" + std::to_string(200);
	if(SelDetector!=0)
	{
		c2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
		c2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());
	}

	

	TCanvas* c3 = new TCanvas("c3", "c3");
	if(SelDetector!=0)
	{
		g200->GetYaxis()->SetRangeUser(0, (*max_element(Y0.begin(),Y0.end()))*1.1);
		g200->Draw("AP");
		g100->Draw("P");
	}
	else
	{
		g100->GetYaxis()->SetRangeUser(0, (*max_element(Y0.begin(),Y0.end()))*1.1);
		g100->Draw("AP");	
	}	
	g50->Draw("P");

	c3->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Detector_str+"-"+MCP_str+"AllGain.png").c_str());
	c3->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Detector_str+"-"+MCP_str+"AllGain.pdf").c_str());

			
}

void EfficiencyMaps(TTree* h4, std::string detector, std::string Selection, std::string Xshift_str, std::string Yshift_str, float XCenter, float YCenter, float bound, std::string RunStats)
{	
	gStyle->SetOptStat(0);
  	int i, Nentries, CH, runNum; 
  	float AmpTemp;
    	
  	h4->GetEntry(0);
 	CH=h4->GetLeaf(detector.c_str())->GetValue(0);

	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

	//2DHist definition 
  	auto *AmpXY0 = new TH2D("AmpXY0","", 32, -16, 16, 32, -16, 16); 
  	auto *AmpXY1 = new TH2D("AmpXY1","", 32, -16, 16, 32, -16, 16); 
  	auto *AmpXYM = new TH2D("AmpXYM","", 32, -16, 16, 32, -16, 16); 
  	
  	Nentries = h4->GetEntries();
  	

   	h4->Draw("Y[0]:X[0]>>AmpXY0", (Selection + " && X[0]>-800 && Y[0]>-800").c_str());
  	h4->Draw(("Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>AmpXY1").c_str(), (Selection + " && X[1]>-800 && Y[1]>-800").c_str());
  	h4->Draw(("(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>AmpXYM").c_str(), (Selection + " && X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800").c_str());

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	//FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+detector+"] (ADC counts)").c_str());
  	
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(XCenter-bound, YCenter-bound, XCenter-bound, YCenter+bound);
	line_left->SetLineColor(kRed);
	line_left->Draw();

	TLine *line_right = new TLine(XCenter+bound, YCenter-bound, XCenter+bound, YCenter+bound);
	line_right->SetLineColor(kRed);
	line_right->Draw();

	TLine *line_up = new TLine(XCenter-bound, YCenter+bound, XCenter+bound, YCenter+bound);
	line_up->SetLineColor(kRed);
	line_up->Draw();

	TLine *line_down = new TLine(XCenter-bound, YCenter-bound, XCenter+bound, YCenter-bound);
	line_down->SetLineColor(kRed);
	line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/EventsXY_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/EventsXY_" + detector + "_" + RunStats + ".png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  

  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");

	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/EventsXY_" + detector + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/EventsXY_" + detector + "_" + RunStats + ".png";

  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	

  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/EventsXY_" + detector + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/EventsXY_" + detector + "_" + RunStats + ".png";  
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	H1->~TH2F();
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
}
#endif	


#include "TFile.h"
#include "TMath.h"
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

int ComputeAvsNoise(std::string RunNtuple, std::string NoiseNtuple, std::string detector, int bound, std::string MCP)
{
	GaussPar AmpPar1, AmpPar2, NoisePar;
	float XMax, YMax, Xshift, Yshift;
	std::string TimeMCP, TimeShift; 
	TFile *f = TFile::Open(RunNtuple.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	TFile *fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree *h4Noise = (TTree*)fNoise->Get("h4");


	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain;

	h4Noise->GetEntry(0);
	std::string NoiseGain = std::to_string((int)h4Noise->GetLeaf("CHGain")->GetValue(0));

	if(NoiseGain!=Gain)
	{
		cout << "Noise Gain =/= Run Gain" << endl;
		return -1;
	}

	//plot detector time distribution
	Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");
	
	AmplitudeProfilesFit(h4, detector, pathToOutput, RunStats, bound, &XMax, &YMax);

	std::string Selection = "(fabs(X[0]-("+std::to_string(XMax)+"))<"+std::to_string(bound)+" || fabs(X[1]-("+std::to_string(XMax)+")-("+std::to_string(Xshift)+"))<"+std::to_string(bound)+") && (fabs(Y[0]-("+std::to_string(YMax)+"))<"+std::to_string(bound)+" || fabs(Y[1]-("+std::to_string(YMax)+")-("+std::to_string(Yshift)+"))<"+std::to_string(bound)+") && amp_max["+MCP+"]>100";
	
	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"fitTimeDist/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHistPar(h4, detector, Selection, pathToOutput, RunStats, &AmpPar);
	Selection = "";
	AmplitudeHistPar(h4Noise, detector, Selection, pathToOutput, "Noise_"+Gain, &NoisePar);

	cout << AmpPar.Mean << " +- " << AmpPar.MeanErr << endl;
	cout << NoisePar.Sigma << " +- " << NoisePar.SigmaErr << endl;
	cout << AmpPar.Mean/NoisePar.Sigma << " +- " << TMath::Sqrt(pow(AmpPar.MeanErr/AmpPar.Mean, 2) + pow(NoisePar.SigmaErr/NoisePar.Sigma, 2))*AmpPar.Mean/NoisePar.Sigma << endl;
	return 0;
}
	




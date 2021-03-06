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
#include "Mylib.h"
#include <string>
#include <fstream>
#include <math.h>

void AveragePulseShapeMCP(std::string FileIn, std::string MCP1, std::string MCP2, std::string runNum)
{
	Int_t Nentries, i, CH;
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

	if(fabs(XMax)<100 && fabs(YMax)<100)
	{
		PulseShapesMCP(h4, detector, 0, XMax, YMax, 2, pathToOutput, RunStats, runNum, MCP);
	}
	else cout<< "SBALLATO!" << endl;
	
}










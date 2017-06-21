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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>
#include <math.h>

void ComparePulseShapes(std::string detector, std::string Energy, std::string run50, std::string run100, std::string run200)
{
	gStyle->SetOptStat(0);

	std::string WF_50 = "WaveForms/"+detector+"_"+Energy+"Gev_G50_"+run50+"_Waveform.root";
	std::string WF_100 = "WaveForms/"+detector+"_"+Energy+"Gev_G100_"+run100+"_Waveform.root";
	std::string WF_200 = "WaveForms/"+detector+"_"+Energy+"Gev_G200_"+run200+"_Waveform.root";

	TFile *f1 = TFile::Open(WF_50.c_str());
	TFile *f2 = TFile::Open(WF_100.c_str());
	TFile *f3 = TFile::Open(WF_200.c_str());

	TF1 *w1 = (TF1*)f1->Get("Waveform_");
	w1->SetTitle("");
	TH1D *w2 = (TH1D*)f2->Get("Waveform_");
	w2->SetTitle("");
	TH1D *w3 = (TH1D*)f3->Get("Waveform_");
	w3->SetTitle("");

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();

	w1->Draw("L");
	w2->SetLineColor(2);
	w2->SetMarkerColor(2);
	w2->Draw("SAME");
	w3->SetLineColor(3);
	w3->Draw("SAME L");

	TLegend* legend = new TLegend(0.60, 0.7, 0.88, 0.88);
	legend->AddEntry(w1, "Gain 50");
	legend->AddEntry(w2, "Gain 100");
	legend->AddEntry(w3, "Gain 200");
	legend->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_"+detector+"_"+Energy+"Gev.pdf").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_"+detector+"_"+Energy+"Gev.png").c_str());


}

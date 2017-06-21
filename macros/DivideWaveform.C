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
//#include "FPCanvasStyle.C"
//#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>
#include <math.h>

void DivideWaveform(std::string detector, std::string Gain, std::string WaveRef, std::string WaveToDiv1, std::string WaveToDiv2, std::string WaveToDiv3)
{
	gStyle->SetOptStat(0);
	
	TFile *f0 = TFile::Open(WaveRef.c_str());
	TFile *f1 = TFile::Open(WaveToDiv1.c_str());
	TFile *f2 = TFile::Open(WaveToDiv2.c_str());
	TFile *f3 = TFile::Open(WaveToDiv3.c_str());
	
	TH1D *w0 = (TH1D*)f0->Get("Waveform_");
	w0->SetTitle("");	
	TH1D *w1 = (TH1D*)f1->Get("Waveform_");
	w1->SetTitle("");
	TH1D *w2 = (TH1D*)f2->Get("Waveform_");
	w2->SetTitle("");	
	TH1D *w3 = (TH1D*)f3->Get("Waveform_");
	w3->SetTitle("");


	w1->Divide(w0);
	w2->Divide(w0);
	w3->Divide(w0);
	
	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();

	w1->GetXaxis()->SetRangeUser(-5, 20);
	w1->SetLineColor(2);
	w1->Draw();

	w2->GetXaxis()->SetRangeUser(-5, 20);
	w2->SetLineColor(3);
	w2->Draw("SAME");

	w3->GetXaxis()->SetRangeUser(-5, 20);
	w3->SetLineColor(4);
	//w3->Draw("SAME");

	TLine *line = new TLine(-5, 1, 20, 1);
	line->SetLineColor(1);
  	line->Draw("SAME");
	
	TLegend* legend = new TLegend(0.5, 0.65, 0.88, 0.88);
	legend->SetHeader(("Detector "+detector+" Gain "+Gain).c_str());
	legend->AddEntry(w1, "20Gev/100Gev");
	legend->AddEntry(w2, "50Gev/100Gev");	
	//legend->AddEntry(w3, "200Gev/100Gev");		
	legend->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".pdf").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".png").c_str());
	
}

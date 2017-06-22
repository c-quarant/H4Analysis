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

void DivideWaveform(std::string detector, std::string Gain)
{
	gStyle->SetOptStat(0);

	std::string WF_20 = "WaveForms/"+detector+"_20Gev_G"+Gain+"_Waveform.root";
	std::string WF_50 = "WaveForms/"+detector+"_50Gev_G"+Gain+"_Waveform.root";
	std::string WF_REF = "WaveForms/"+detector+"_100Gev_G"+Gain+"_Waveform.root";
	std::string WF_150 = "WaveForms/"+detector+"_150Gev_G"+Gain+"_Waveform.root";
	std::string WF_200 = "WaveForms/"+detector+"_200Gev_G"+Gain+"_Waveform.root";
	
	TFile *f0 = TFile::Open(WF_REF.c_str());
	TFile *f1 = TFile::Open(WF_20.c_str());
	TFile *f2 = TFile::Open(WF_50.c_str());
	TFile *f3 = TFile::Open(WF_150.c_str());
	TFile *f4 = TFile::Open(WF_200.c_str());

	TH1D *w0 = (TH1D*)f0->Get("Waveform_");
	w0->SetTitle("");	
	TH1D *w1 = (TH1D*)f1->Get("Waveform_");
	w1->SetTitle("");
	TH1D *w2 = (TH1D*)f2->Get("Waveform_");
	w2->SetTitle("");	
	TH1D *w3 = (TH1D*)f3->Get("Waveform_");
	w3->SetTitle("");
	TH1D *w4 = (TH1D*)f4->Get("Waveform_");
	w4->SetTitle("");

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();
	
	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);

	cUp->SetBottomMargin(0.01); 
	cDown->SetTopMargin(0.01); 
    
    	cUp->Draw();
	cDown->Draw();

	cUp->cd();
	w0->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w0->SetLineColor(kBlack);
	w0->Draw();

	cUp->cd();
	w1->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w1->SetLineColor(kRed);
	w1->Draw();

	cUp->cd();
	w2->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w2->SetLineColor(kBlue);
	w2->Draw("SAME");

	cUp->cd();
	w3->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w3->SetLineColor(kGreen);
	w3->Draw("SAME");

	cUp->cd();
	w4->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w4->SetLineColor(kViolet);
	w4->Draw("SAME");

	w1->Divide(w0);
	w2->Divide(w0);
	w3->Divide(w0);
	w4->Divide(w0);

	cDown->cd();

	w1->GetYaxis()->SetRangeUser(0.9, 1.1);
	w1->SetLineColor(kRed);
	w1->Draw();

	w2->GetYaxis()->SetRangeUser(0.9, 1.1);
	w2->SetLineColor(kBlue);
	w2->Draw("SAME");

	w3->GetYaxis()->SetRangeUser(0.9, 1.1);
	w3->SetLineColor(4);
	w3->Draw("SAME");

	w4->GetYaxis()->SetRangeUser(0.9, 1.1);
	w4->SetLineColor(kRed);
	w4->Draw();

	TLine *line = new TLine(-5, 1, 20, 1);
	line->SetLineColor(1);
  	line->Draw("SAME");
	
	TLegend* legend = new TLegend(0.5, 0.65, 0.88, 0.88);
	legend->SetHeader(("Detector "+detector+" Gain "+Gain).c_str());
	legend->AddEntry(w1, "20Gev/100Gev");
	legend->AddEntry(w2, "50Gev/100Gev");	
	legend->AddEntry(w3, "150Gev/100Gev");	
	legend->AddEntry(w4, "200Gev/100Gev");			
	legend->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".pdf").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".png").c_str());
	
}

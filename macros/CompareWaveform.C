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

void CompareWaveform(std::string detector, std::string Gain, float xmin, float xmax)
{
	int i;
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
	
	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.5,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.5);

	cUp->SetBottomMargin(0.01); 
	cDown->SetTopMargin(0.01); 
    
    	cUp->Draw();
	cDown->Draw();

	cUp->cd();	
	w0->GetYaxis()->SetTitleSize(0.05);	
	w0->GetXaxis()->SetRangeUser(xmin, xmax);	
	w0->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w0->SetLineColor(kBlack);
	w0->Draw();

	w1->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w1->SetLineColor(kRed);
	w1->Draw("SAME");

	w2->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w2->SetLineColor(kBlue);
	w2->Draw("SAME");

	w3->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w3->SetLineColor(kGreen);
	w3->Draw("SAME");

	w4->GetYaxis()->SetRangeUser(-0.05, 1.05);
	w4->SetLineColor(kViolet);
	w4->Draw("SAME");

	TLegend* legend = new TLegend(0.65, 0.4, 0.89, 0.88);
	legend->SetHeader(("Detector "+detector+" Gain "+Gain).c_str());
	legend->AddEntry(w1, "Energy 20 Gev");
	legend->AddEntry(w2, "Energy 50 Gev");	
	legend->AddEntry(w0, "Energy 100 Gev");	
	legend->AddEntry(w3, "Energy 150 Gev");	
	legend->AddEntry(w4, "Energy 200 Gev");			
	legend->Draw("SAME");
	
	TH1D *R1 = new TH1D("R1", "", 300, -10, 40); 
	for(i=0; i<300; i++)
	{
		R1->SetBinContent(i, w1->GetBinContent(i));
		R1->SetBinError(i, w1->GetBinError(i));
	}
	R1->SetTitle("");
	TH1D *R2 = new TH1D("R2", "", 300, -10, 40); 
	for(i=0; i<300; i++)
	{
		R2->SetBinContent(i, w2->GetBinContent(i));
		R2->SetBinError(i, w2->GetBinError(i));
	}
	R2->SetTitle("");	
	TH1D *R3 = new TH1D("R3", "", 300, -10, 40); 
	for(i=0; i<300; i++)
	{
		R3->SetBinContent(i, w3->GetBinContent(i));
		R3->SetBinError(i, w3->GetBinError(i));
	}
	R3->SetTitle("");
	TH1D *R4 = new TH1D("R4", "", 300, -10, 40); 
	for(i=0; i<300; i++)
	{
		R4->SetBinContent(i, w4->GetBinContent(i));
		R4->SetBinError(i, w4->GetBinError(i));
	}
	R4->SetTitle("");

	R1->Divide(w0);
	R2->Divide(w0);
	R3->Divide(w0);
	R4->Divide(w0);

	cDown->cd();

	R1->GetXaxis()->SetTitle("Time_MCP1 (ns)");
	R1->GetXaxis()->SetTitleSize(0.06);
	R1->GetXaxis()->SetRangeUser(xmin, xmax);	
	R1->GetYaxis()->SetRangeUser(0.9, 1.1);
	R1->SetLineColor(kRed);
	R1->Draw();

	R2->GetYaxis()->SetRangeUser(0.9, 1.1);
	R2->SetLineColor(kBlue);
	R2->Draw("SAME");

	R3->GetYaxis()->SetRangeUser(0.9, 1.1);
	R3->SetLineColor(kGreen);
	R3->Draw("SAME");

	R4->GetYaxis()->SetRangeUser(0.9, 1.1);
	R4->SetLineColor(kViolet);
	R4->Draw("SAME");

	TLine *line = new TLine(xmin, 1, xmax, 1);
	line->SetLineColor(1);
  	line->Draw("SAME");
	
	
	TLegend *Rlegend = new TLegend(0.5, 0.88, 0.88, 0.96);
	Rlegend->SetHeader("Ratio plot whit respect to 100 Gev");
	Rlegend->SetMargin(0.05);		
	Rlegend->Draw("SAME");
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".pdf").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Gain"+Gain+".png").c_str());
	
}

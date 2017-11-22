#include "TFile.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include <iostream>
#include "TGaxis.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include <string>
#include <fstream>
#include <math.h>

void Compare2PulseShapes(std::string Waveform1, std::string Waveform2)
{
	gStyle->SetOptStat(0);

	TFile *f1 = TFile::Open(Waveform1.c_str());
	TFile *f2 = TFile::Open(Waveform2.c_str());
	
	TH1D* w1 = (TH1D*)f1->Get("XTAL_D3_E100_G50_prof_");
	w1->SetTitle("");
	w1->GetXaxis()->SetTitle("time (ns)");
	w1->GetXaxis()->SetRangeUser(-10, 40);
	w1->GetYaxis()->SetTitle("ADC counts");
	w1->GetYaxis()->SetTitleSize(0.055);
	w1->GetYaxis()->SetTitleOffset(0.75);

	TH1D *w2 = (TH1D*)f2->Get("XTAL_C3_E100_G50_prof_");
	w2->SetTitle("");
	w2->SetLineColor(kRed);
	w2->SetMarkerColor(kRed);

	TH1D* Ratio = new TH1D();
	Ratio = (TH1D*)w1->Clone();	
	Ratio->GetXaxis()->SetTitle("t_{MCP1} (ns)");
	Ratio->GetXaxis()->SetTitleSize(0.055);
	Ratio->GetXaxis()->SetTitleOffset(0.75);
	Ratio->GetYaxis()->SetRangeUser(0.8, 1.1);
	Ratio->GetYaxis()->SetTitle("A_{D3}/A_{C3}");
	Ratio->GetYaxis()->SetTitleSize(0.055);
	Ratio->GetYaxis()->SetTitleOffset(0.75);
	Ratio->Divide(w2);
	Ratio->SetLineColor(kRed);
	
	TLine *line = new TLine(-10, 1, 40, 1);
	line->SetLineColor(kBlack);

	TLegend* legend = new TLegend(0.55, 0.7, 0.88, 0.88);
	legend->AddEntry(w1, w1->GetName());
	legend->AddEntry(w2, w2->GetName());

	TLegend* Rlegend = new TLegend(0.55, 0.76, 0.88, 0.94);
	Rlegend->AddEntry(Ratio, "D3 over C3 pulse shape");

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();
	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.5,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.5);

	cUp->SetBottomMargin(0.01); 	
    	cUp->Draw();
	cDown->SetTopMargin(0.01);
	cDown->Draw();

	cUp->cd();

	w1->Draw();
	w2->Draw("SAME");
	legend->Draw("SAME");	

	cDown->cd();	

	Ratio->Draw();
	line->Draw("SAME");
	Rlegend->Draw("SAME");

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_D3_C3_100Gev_G50.pdf"));
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_D3_C3_100Gev_G50.png"));
}
